#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include "detail/JointCliProgressReporter.hpp"
#include "detail/FirstStageInitialization.hpp"
#include "detail/FittingWorkset.hpp"
#include "detail/PotentialFittingWorkflow.hpp"
#include "detail/joint_component/Problem.hpp"
#include "detail/StageSummary.hpp"
#include "detail/PostFitPeeling.hpp"
#include "detail/LocalConditionalCovariance.hpp"
#include "detail/GroupPotentialFitting.hpp"
#include "detail/gaussian_fit/FittingRanges.hpp"
#include "detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "detail/second_stage/IterationProcess.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem::core::detail {
namespace {
bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role;
}

PostJointLocalResult FitPostJointLocalTarget(
    const MapObject & map, const ModelObject & model, const JointProblem & problem,
    const JointAnalysisResult & snapshot, int atom_id, const PostFitPeelingResult & peeling,
    const FitOptions & options)
{
    const auto atom = model.FindAtomPtr(atom_id);
    const auto local = AtomLocalPotentialView::For(*atom);
    const auto & joint = local.GetStageEstimate(FittingStage::Second);
    PostJointLocalResult result;
    result.source = joint.source;
    result.sample_count = local.GetRawSamplingEntriesView().size();
    result.status = EvidenceStatus::Ineligible;
    if (joint.source.method != EstimateMethod::JointComponents || joint.source.role != FittingRole::Target)
    { result.reason = "not-a-joint-target"; return result; }
    if (!joint.point) { result.reason = "missing-joint-endpoint"; return result; }
    result.fixed_charge = {joint.point->GetOffset(), joint.point->GetWidth()};
    if (joint.convergence != JointCheckStatus::Passed)
    { result.reason = "joint-endpoint-not-qualified"; return result; }
    if (!local.HasSampleGeometry())
    { result.status = EvidenceStatus::Unavailable; result.reason = "sample-geometry-unavailable"; return result; }
    if (!SameSource(joint.source, peeling.source))
        throw std::invalid_argument("Post-Joint local fitting source mismatch.");

    const auto prepared_samples = BuildPostJointLocalSamples(map, model, problem, atom_id, peeling, snapshot);
    result.fixed_charge = prepared_samples.fixed_charge;
    for (const auto & sample : prepared_samples.samples)
        if (!sample.reason.empty())
            result.sample_failures.push_back({sample.raw_sample_index, sample.distance, sample.reason});
    const auto & raw = local.GetRawSamplingEntriesView();
    auto gaussian_samples = prepared_samples.BuildGaussianSamples(raw);
    std::vector<double> responses;
    responses.reserve(gaussian_samples.size());
    for (const auto & sample : gaussian_samples) responses.push_back(sample.response);

    const GaussianModel3D no_charge{0.0, 1.0, 0.0};
    const PreparedLocalGaussianDesign design{gaussian_samples, 0.0, kSignalDistanceMax};
    const auto dataset = design.BuildDataset(responses, no_charge);
    result.effective_sample_count = static_cast<std::size_t>(dataset.y.size());
    if (dataset.y.size() <= 2)
    { result.status = EvidenceStatus::Unavailable; result.reason = "insufficient-positive-responses"; return result; }
    Eigen::ColPivHouseholderQR<RHBMDesignMatrix> design_qr(dataset.X);
    if (dataset.X.cols() != 2 || design_qr.rank() != 2)
    { result.status = EvidenceStatus::Unavailable; result.reason = "rank-deficient-local-design"; return result; }

    const double alpha = TrainAlphaR(gaussian_samples, options);
    if (!std::isfinite(alpha) || alpha < 0.0)
    { result.status = EvidenceStatus::Unavailable; result.reason = "local-alpha-training-failed"; return result; }
    result.fit = EstimateLocalGaussian(gaussian_samples, alpha, options, no_charge);
    result.fit.alpha_r = alpha;
    if (!result.fit.fit_result)
    { result.status = EvidenceStatus::Unavailable; result.reason = "local-mdpde-result-unavailable"; return result; }
    auto & fit = *result.fit.fit_result;
    result.fit.diagnostics = LocalFitDiagnostics{fit.status, fit.sigma_square, fit.diagnostics, fit.refinement};
    const auto reject = [&result](std::string reason) {
        result.status = EvidenceStatus::Unavailable;
        result.reason = std::move(reason);
        result.fit.fit_result.reset();
        return result;
    };
    const auto & ols = result.fit.ols.GetModel();
    const auto & mdpde = result.fit.mdpde.GetModel();
    if (fit.Qualification() == RHBMSolveQualification::Unqualified ||
        fit.status != RHBMEstimationStatus::SUCCESS)
        return reject("local-mdpde-unqualified");
    if (!fit.beta_ols.allFinite() || !fit.beta_mdpde.allFinite() ||
        !std::isfinite(ols.GetAmplitude()) || ols.GetAmplitude() <= 0.0 ||
        !std::isfinite(ols.GetWidth()) || ols.GetWidth() <= 0.0 ||
        !std::isfinite(mdpde.GetAmplitude()) || mdpde.GetAmplitude() <= 0.0 ||
        !std::isfinite(mdpde.GetWidth()) || mdpde.GetWidth() <= 0.0)
        return reject("invalid-local-parameters");
    if (!std::isfinite(fit.sigma_square) || fit.sigma_square <= 0.0 ||
        fit.sigma_square == std::numeric_limits<double>::max() ||
        fit.data_weight.diagonal().size() != dataset.y.size() ||
        fit.data_covariance.diagonal().size() != dataset.y.size() ||
        !fit.data_weight.diagonal().allFinite() || !fit.data_covariance.diagonal().allFinite() ||
        (fit.data_weight.diagonal().array() <= 0.0).any() ||
        (fit.data_covariance.diagonal().array() <= 0.0).any())
        return reject("invalid-local-working-covariance");

    const Eigen::VectorXd inverse_variance = fit.data_covariance.diagonal().cwiseInverse();
    RHBMInformation information;
    information.precision = dataset.X.transpose() * inverse_variance.asDiagonal() * dataset.X;
    information.moment = dataset.X.transpose() * (inverse_variance.array() * dataset.y.array()).matrix();
    const Eigen::Matrix2d precision = information.precision;
    Eigen::LLT<Eigen::Matrix2d> precision_factor(precision);
    if (!information.precision.allFinite() || !information.moment.allFinite() ||
        precision_factor.info() != Eigen::Success)
        return reject("invalid-local-measurement-information");

    result.information = std::move(information);
    result.status = EvidenceStatus::Available;
    result.reason.clear();
    result.fit.fit_result.reset();
    return result;
}

void LogJointInfo(const std::string & message)
{
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message);
}

}

void RunTwoStageFromPreparedSamples(
    ModelObject & model, const std::vector<AtomObject *> & atoms, const FitOptions & options)
{
    model.EditAnalysis().InitializeLocalFittingSeedModels();
    RunBatchFirstStageFromPreparedSamples(model, atoms, options);
    RunSecondStageIterations(model, options);

    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model));
    RunGroupAlphaTraining(model, options);
    RunGroupPotentialFitting(model, options);
}

void RunTwoStageFromPreparedSamples(ModelObject & model, const FitOptions & options)
{
    const auto workset = MakeTwoStageFittingWorkset(model);
    RunTwoStageFromPreparedSamples(model, workset.contributors, options);
}

void RunTwoStageWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    model.EditAnalysis().InitializeFromSelection();
    const auto workset = MakeTwoStageFittingWorkset(model);
    const auto first_stage_atoms = CollectFirstStageAtoms(workset);
    RunPotentialSamplingWorkflow(map, model, first_stage_atoms, options.sampling_method,
        options.thread_size, " Sampling", options.quiet_mode);
    RunTwoStageFromPreparedSamples(model, workset.contributors, options);
}

void RunJointComponentWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    if (!options.quiet_mode) LogJointInfo("[Joint] Building joint problem...");
    const auto construction_start = std::chrono::steady_clock::now();
    const auto problem = BuildJointProblem(map, model);
    const auto construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - construction_start).count();
    if (!options.quiet_mode)
    {
        const auto & input = problem.Input();
        const auto targets = input.selection_domain ? input.selection_domain->target_indices.size() : 0;
        std::ostringstream message;
        message << "[Joint] Problem ready: observations=" << input.observations.size()
            << ", contributors=" << input.atom_ids.size()
            << ", targets=" << targets
            << ", halo=" << input.atom_ids.size() - targets
            << ", elapsed=" << std::fixed << std::setprecision(2) << construction_seconds << 's';
        LogJointInfo(message.str());
        if(input.selection_domain && input.selection_domain->halo_eligibility)
        {
            const auto & eligibility=*input.selection_domain->halo_eligibility;
            std::size_t insufficient_rows{},insufficient_radial{};
            for(const auto & excluded:eligibility.excluded_halos)
            {
                insufficient_rows+=excluded.reason=="insufficient-observations";
                insufficient_radial+=excluded.reason=="insufficient-radial-information";
            }
            std::ostringstream halo_message;
            halo_message << "[Joint] Halo eligibility: candidates=" << eligibility.candidate_halos
                << ", retained=" << eligibility.retained_halos
                << ", excluded=" << eligibility.excluded_halos.size()
                << ", excluded_insufficient_rows=" << insufficient_rows
                << ", excluded_insufficient_radial_information=" << insufficient_radial;
            LogJointInfo(halo_message.str());
        }
    }
    const auto initialization_start = std::chrono::steady_clock::now();
    const auto initialization = [&] {
        const auto workset = MakeJointFittingWorkset(model, problem);
        model.EditAnalysis().InitializeFromSelection();
        return PrepareJointWidthInitialization(map, model, workset, options);
    }();
    const auto initialization_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - initialization_start).count();
    std::optional<JointCliProgressReporter> reporter;
    joint_component::JointProgressObserver observer;
    if (!options.quiet_mode)
    {
        reporter.emplace();
        observer = reporter->Observer();
    }
    auto snapshot = [&] {
        auto fit = joint_component::FitFixedNeighborComponents(problem, initialization.b, {}, observer);
        fit.costs.construction_seconds = construction_seconds;
        fit.costs.initialization_seconds = initialization_seconds;
        fit.initialization = initialization;
        return CaptureJointAnalysisResult(fit);
    }();
    data_internal::ApplyJointStageEstimates(model, snapshot, boost::uuids::to_string(boost::uuids::random_generator()()));
    if (!options.quiet_mode) LogJointInfo("[Joint] Computing post-fit outputs...");
    auto peelings = BuildPostFitPeelingSamples(
        map, model, problem, problem.Input().selection_domain->target_indices, &snapshot);
    std::vector<std::pair<int, PostFitPeelingResult>> peeling_updates;
    std::vector<std::pair<int, GroupParameterEvidence>> evidence_updates;
    peeling_updates.reserve(peelings.size());
    evidence_updates.reserve(peelings.size());
    for (auto & [id, peeling] : peelings)
    {
        auto & atom = *model.FindAtomPtr(id);
        const auto local = AtomLocalPotentialView::For(atom);
        const auto & stage = local.GetStageEstimate(FittingStage::Second);
        const auto & raw_samples = local.GetRawSamplingEntriesView();
        LocalConditionalCovariance covariance;
        covariance = ComputeLocalConditionalCovariance(stage, raw_samples, peeling, local.HasSampleGeometry());
        covariance.source = stage.source;
        auto evidence = BuildLocalGroupEvidence(stage, &peeling, covariance);
        peeling_updates.emplace_back(id, std::move(peeling));
        evidence_updates.emplace_back(id, std::move(evidence));
    }
    model.EditAnalysis().ApplyJointPostFitEvidenceBatch(std::move(peeling_updates), std::move(evidence_updates));
    model.EditAnalysis().SetJointResult(std::move(snapshot));
    const auto & joint_snapshot = *model.GetAnalysisView().GetJointResult();
    std::vector<std::pair<int, PostJointLocalResult>> local_updates;
    for (const auto * atom : model.GetSelectedAtoms())
    {
        const auto local = AtomLocalPotentialView::For(*atom);
        const auto & stage = local.GetStageEstimate(FittingStage::Second);
        if (stage.source.method != EstimateMethod::JointComponents || stage.source.role != FittingRole::Target) continue;
        const auto & peeling = local.GetPostFitPeeling();
        PostJointLocalResult result;
        if (peeling) result = FitPostJointLocalTarget(map, model, problem, joint_snapshot, atom->GetSerialID(), *peeling, options);
        else
        {
            result.source = stage.source;
            result.status = EvidenceStatus::Unavailable;
            result.reason = "peeling-unavailable";
        }
        local_updates.emplace_back(atom->GetSerialID(), std::move(result));
    }
    model.EditAnalysis().ApplyPostJointLocalResultsBatch(std::move(local_updates));
    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model));
    RunGroupPotentialFitting(model, options);
    if (!options.quiet_mode) LogJointInfo("[Joint] Completed.");
}

} // namespace rhbm_gem::core::detail

namespace rhbm_gem::core {

void RunPotentialFittingWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    switch (options.estimator)
    {
    case PotentialEstimator::TWO_STAGE:
        detail::RunTwoStageWorkflow(map, model, options);
        return;
    case PotentialEstimator::JOINT_COMPONENTS:
        detail::RunJointComponentWorkflow(map, model, options);
        return;
    }
    throw std::invalid_argument("Unsupported potential estimator.");
}

} // namespace rhbm_gem::core
