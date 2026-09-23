#include "JointUncertainty.hpp"
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
GaussianModel3DWithUncertainty Decode(const Eigen::VectorXd & z, const Eigen::MatrixXd & covariance, double charge)
{
    const auto decoded = linearization_service::DecodeParameterVector(z, covariance);
    GaussianModel3D::RequireFinitePositiveWidthModel(decoded.GetModel());
    const auto & sd = decoded.GetStandardDeviationModel();
    if (!std::isfinite(sd.GetAmplitude()) || !std::isfinite(sd.GetWidth()))
        throw std::runtime_error("Nonfinite posterior uncertainty.");
    return {decoded.GetModel().WithOffset(charge),
        GaussianModel3DUncertainty(sd.GetAmplitude(), sd.GetWidth(), std::numeric_limits<double>::quiet_NaN())};
}
}
void RunJointGroupPotentialFitting(ModelObject & model, const FitOptions & options)
{
    const auto view = model.GetAnalysisView();
    auto editor = model.EditAnalysis();
    const auto keys = view.CollectAtomGroupKeys();
    std::vector<std::vector<GaussianModel3D>> training;
    for (const auto key : keys)
    {
        std::vector<GaussianModel3D> members;
        const auto & atoms = view.GetAtomObjectList(key);
        for (const auto * atom : atoms)
        {
            const auto local = AtomLocalPotentialView::For(*atom);
            const auto & evidence = local.GetGroupEvidence();
            if (evidence && evidence->status == EvidenceStatus::Available)
                members.push_back(local.GetFinalModel(FittingStage::Second));
        }
        if (members.size() >= 10 && atoms.front()->IsMainChainAtom()) training.push_back(std::move(members));
    }
    const double alpha = TrainAlphaG(training, options);
    const RHBMExecutionOptions execution{.quiet_mode = options.quiet_mode, .thread_size = options.thread_size};
    for (const auto key : keys)
    {
        GroupParameterSummary summary;
        summary.status = EvidenceStatus::Unavailable;
        std::vector<RHBMParameterVector> z;
        std::vector<RHBMInformation> information;
        std::vector<double> charges, all_charges;
        Eigen::Vector2d sums = Eigen::Vector2d::Zero();
        for (const auto * atom : view.GetAtomObjectList(key))
        {
            const auto local = AtomLocalPotentialView::For(*atom);
            const auto & stage = local.GetStageEstimate(FittingStage::Second);
            if (stage.source.role != FittingRole::Target) continue;
            summary.source_id = stage.source.run_id;
            if (stage.point)
            {
                ++summary.point_count;
                sums += Eigen::Vector2d(stage.point->GetAmplitude(), stage.point->GetWidth());
                all_charges.push_back(stage.point->GetOffset());
            }
            const auto & evidence = local.GetGroupEvidence();
            if (!evidence || evidence->status != EvidenceStatus::Available || !evidence->estimate || !evidence->covariance)
            { ++summary.excluded_count; continue; }
            if (evidence->atom_id != stage.source.atom_id || evidence->source_id != stage.source.run_id)
                throw std::invalid_argument("Joint group evidence source mismatch.");
            const Eigen::Matrix2d precision = evidence->covariance->ldlt().solve(Eigen::Matrix2d::Identity());
            z.push_back(*evidence->estimate);
            information.push_back({precision, precision * *evidence->estimate});
            charges.push_back(stage.point->GetOffset());
            summary.member_ids.push_back(atom->GetSerialID());
        }
        summary.eligible_count = z.size();
        if (summary.point_count)
            summary.descriptive_mean = GaussianModel3D(sums(0) / static_cast<double>(summary.point_count), sums(1) / static_cast<double>(summary.point_count),
                array_helper::ComputeMedian(all_charges));
        if (z.size() < 2) summary.reason = z.empty() ? "no-eligible-evidence" : "single-member";
        else
        {
            const auto mu = rhbm_helper::EstimateMuMDPDE(alpha, rhbm_helper::BuildBetaMatrix(z), execution);
            if (mu.status != RHBMEstimationStatus::SUCCESS) summary.reason = "group-estimation-failed";
            else if (!mu.covariance_available || !PositiveDefinite(mu.capital_lambda)) summary.reason = "singular-group-covariance";
            else
            {
                try
                {
                    const auto web = rhbm_helper::EstimateWEBFromInformation(information, mu.mu_mdpde,
                        mu.member_capital_lambda_list, execution);
                    if (web.status != RHBMEstimationStatus::SUCCESS) throw std::runtime_error("WEB failed.");
                    const double charge = array_helper::ComputeMedian(charges);
                    GroupGaussianResult result;
                    result.alpha_g = alpha;
                    result.mean = linearization_service::DecodeParameterVector(mu.mu_mean).WithOffset(charge);
                    result.mdpde = linearization_service::DecodeParameterVector(mu.mu_mdpde).WithOffset(charge);
                    result.prior = Decode(web.mu_prior, mu.capital_lambda, charge);
                    const auto distances = rhbm_helper::CalculateMemberStatisticalDistance(web.mu_prior,
                        mu.capital_lambda, web.beta_posterior_matrix);
                    const auto flags = rhbm_helper::CalculateOutlierMemberFlag(2, distances);
                    for (std::size_t i = 0; i < z.size(); ++i)
                    {
                        const auto column = static_cast<Eigen::Index>(i);
                        GroupGaussianMemberResult member;
                        member.posterior = Decode(web.beta_posterior_matrix.col(column), web.capital_sigma_posterior_list[i], charges[i]);
                        member.is_outlier = flags(column);
                        member.statistical_distance = distances(column);
                        member.evidence_source_id = summary.source_id;
                        member.charge_inferred = false;
                        member.parameter_covariance = web.capital_sigma_posterior_list[i];
                        result.member_results.push_back(std::move(member));
                    }
                    summary.inference = std::move(result);
                    summary.status = EvidenceStatus::Available;
                }
                catch (const std::exception &) { summary.reason = "posterior-estimation-failed"; }
            }
        }
        editor.ApplyAtomGroupParameterSummary(key, std::move(summary));
    }
}
}
