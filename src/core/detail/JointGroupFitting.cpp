#include "GroupPotentialFitting.hpp"
#include "data/detail/PostJointLocalValidation.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <limits>

namespace rhbm_gem::core::detail {
namespace {
bool PositiveDefinite(const Eigen::MatrixXd & matrix)
{
    if (!matrix.allFinite()) return false;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(matrix);
    return spectrum.info() == Eigen::Success && spectrum.eigenvalues()(0) >
        std::numeric_limits<double>::epsilon() * static_cast<double>(matrix.rows()) * spectrum.eigenvalues().maxCoeff();
}

bool HasCurrentPostJointLocal(
    const AtomLocalPotentialView & local, const LocalStageEstimate & joint, int atom_id)
{
    const auto & result = local.GetPostJointLocalResult();
    if (!result) return false;
    const auto & peeling = local.GetPostFitPeeling();
    const auto assessment = data_internal::AssessPostJointLocalResult(
        *result, atom_id, joint,
        peeling ? &*peeling : nullptr, local.GetRawSamplingEntriesView());
    if (!assessment.valid) throw std::invalid_argument(assessment.reason);
    return assessment.group_eligible;
}

GaussianModel3DWithUncertainty Decode(
    const Eigen::VectorXd & z, const Eigen::MatrixXd & covariance)
{
    const auto decoded = linearization_service::DecodeParameterVector(z, covariance);
    GaussianModel3D::RequireFinitePositiveWidthModel(decoded.GetModel());
    const auto & sd = decoded.GetStandardDeviationModel();
    if (!std::isfinite(sd.GetAmplitude()) || !std::isfinite(sd.GetWidth()))
        throw std::runtime_error("Nonfinite posterior uncertainty.");
    return {decoded.GetModel().WithOffset(0.0),
        GaussianModel3DUncertainty(sd.GetAmplitude(), sd.GetWidth(), std::numeric_limits<double>::quiet_NaN())};
}
}

void RunJointGroupPotentialFitting(ModelObject & model, const FitOptions & options)
{
    const auto view = model.GetAnalysisView();
    auto editor = model.EditAnalysis();
    const auto keys = view.CollectAtomGroupKeys();
    std::vector<std::vector<GaussianModel3D>> training;
    std::string training_run_id;
    for (const auto key : keys)
    {
        std::vector<GaussianModel3D> members;
        const auto & atoms = view.GetAtomObjectList(key);
        for (const auto * atom : atoms)
        {
            const auto local = AtomLocalPotentialView::For(*atom);
            const auto & joint = local.GetStageEstimate(FittingStage::Second);
            if (HasCurrentPostJointLocal(local, joint, atom->GetSerialID()))
            {
                if (training_run_id.empty()) training_run_id = joint.source.run_id;
                else if (training_run_id != joint.source.run_id)
                    throw std::invalid_argument("Post-Joint alpha training mixes Joint runs.");
                members.push_back(local.GetPostJointLocalResult()->fit.mdpde.GetModel());
            }
        }
        if (members.size() >= 10 && atoms.front()->IsMainChainAtom()) training.push_back(std::move(members));
    }
    const double alpha = TrainAlphaG(training, options);
    const RHBMExecutionOptions execution{.quiet_mode = options.quiet_mode, .thread_size = options.thread_size};
    for (const auto key : keys)
    {
        editor.SetAtomGroupAlphaG(key, alpha);
        GroupParameterSummary summary;
        summary.status = EvidenceStatus::Unavailable;
        std::vector<RHBMParameterVector> beta;
        std::vector<RHBMInformation> information;
        std::vector<JointFixedChargeProfile> fixed_charges;
        Eigen::Vector2d sums = Eigen::Vector2d::Zero();
        std::size_t target_count{};
        for (const auto * atom : view.GetAtomObjectList(key))
        {
            const auto local = AtomLocalPotentialView::For(*atom);
            const auto & joint = local.GetStageEstimate(FittingStage::Second);
            if (joint.source.method != EstimateMethod::JointComponents || joint.source.role != FittingRole::Target) continue;
            ++target_count;
            if (summary.source_id.empty()) summary.source_id = joint.source.run_id;
            else if (summary.source_id != joint.source.run_id)
                throw std::invalid_argument("Post-Joint group fitting mixes Joint runs.");
            if (!HasCurrentPostJointLocal(local, joint, atom->GetSerialID())) continue;
            const auto & post_joint = *local.GetPostJointLocalResult();
            const auto & point = post_joint.fit.mdpde.GetModel();
            const auto encoded = linearization_service::EncodeGaussianToParameterVector(point);
            if (!encoded.allFinite()) throw std::invalid_argument("Invalid post-Joint local parameter encoding.");
            beta.push_back(encoded);
            information.push_back(*post_joint.information);
            fixed_charges.push_back(post_joint.fixed_charge);
            sums += Eigen::Vector2d(point.GetAmplitude(), point.GetWidth());
            summary.member_ids.push_back(atom->GetSerialID());
        }
        summary.point_count = beta.size();
        summary.eligible_count = beta.size();
        summary.excluded_count = target_count - beta.size();
        if (!beta.empty())
            summary.descriptive_mean = GaussianModel3D(
                sums(0) / static_cast<double>(beta.size()), sums(1) / static_cast<double>(beta.size()),
                0.0);
        if (!fixed_charges.empty())
        {
            std::vector<double> coefficients;
            coefficients.reserve(fixed_charges.size());
            for (const auto & profile : fixed_charges) coefficients.push_back(profile.coefficient);
            summary.descriptive_charge_median = array_helper::ComputeMedian(coefficients);
        }
        if (beta.size() < 2)
            summary.reason = beta.empty() ? "no-eligible-post-joint-local-results" : "single-member";
        else
        {
            const auto beta_matrix = rhbm_helper::BuildBetaMatrix(beta);
            const auto mu = rhbm_helper::EstimateMuMDPDE(alpha, beta_matrix, execution);
            if (mu.status != RHBMEstimationStatus::SUCCESS) summary.reason = "group-estimation-failed";
            else if (!mu.covariance_available || !PositiveDefinite(mu.capital_lambda))
                summary.reason = "singular-group-covariance";
            else
            {
                try
                {
                    const auto web = rhbm_helper::EstimateWEBFromInformation(
                        information, mu.mu_mdpde, mu.member_capital_lambda_list, execution);
                    if (web.status != RHBMEstimationStatus::SUCCESS) throw std::runtime_error("WEB failed.");
                    GroupGaussianResult group;
                    group.alpha_g = alpha;
                    group.mean = linearization_service::DecodeParameterVector(mu.mu_mean).WithOffset(0.0);
                    group.mdpde = linearization_service::DecodeParameterVector(mu.mu_mdpde).WithOffset(0.0);
                    group.prior = Decode(web.mu_prior, mu.capital_lambda);
                    const auto distances = rhbm_helper::CalculateMemberStatisticalDistance(
                        web.mu_prior, mu.capital_lambda, web.beta_posterior_matrix);
                    const auto flags = rhbm_helper::CalculateOutlierMemberFlag(2, distances);
                    for (std::size_t i = 0; i < beta.size(); ++i)
                    {
                        const auto column = static_cast<Eigen::Index>(i);
                        GroupGaussianMemberResult member;
                        member.posterior = Decode(web.beta_posterior_matrix.col(column),
                            web.capital_sigma_posterior_list[i]);
                        member.is_outlier = flags(column);
                        member.statistical_distance = distances(column);
                        member.evidence_source_id = summary.source_id;
                        member.charge_inferred = false;
                        member.parameter_covariance = web.capital_sigma_posterior_list[i];
                        member.fixed_charge = fixed_charges[i];
                        group.member_results.push_back(std::move(member));
                    }
                    summary.inference = std::move(group);
                    summary.status = EvidenceStatus::Available;
                }
                catch (const std::exception &)
                { summary.reason = "posterior-estimation-failed"; }
            }
        }
        editor.ApplyAtomGroupParameterSummary(key, std::move(summary));
    }
}
}
