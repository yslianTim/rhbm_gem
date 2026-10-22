#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include "detail/joint_component/JointCliProgressReporter.hpp"
#include "detail/FirstStageInitialization.hpp"
#include "detail/FittingWorkset.hpp"
#include "detail/PotentialFittingWorkflow.hpp"
#include "detail/joint_component/HaloMedianReconstruction.hpp"
#include "detail/StageSummary.hpp"
#include "detail/joint_component/PostFitPeeling.hpp"
#include "detail/GroupPotentialFitting.hpp"
#include "detail/joint_component/ResourceWork.hpp"
#include "data/detail/PostJointLocalValidation.hpp"
#include "detail/gaussian_fit/FittingRanges.hpp"
#include "detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "detail/second_stage/IterationProcess.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <chrono>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem::core::detail {
namespace {
bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role;
}

PostJointLocalResult FitPostJointLocalTarget(
    const JointPostprocessingContext & postprocessing, const ModelObject & model,
    int atom_id, const PostFitPeelingResult & peeling,
    const FitOptions & options, JointPostprocessingMetrics & metrics)
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

    auto phase_start = std::chrono::steady_clock::now();
    const auto prepared_samples = postprocessing.BuildPostJointLocalSamples(atom_id, peeling);
    metrics.charge_correction_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start).count();
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

    phase_start = std::chrono::steady_clock::now();
    const double alpha = TrainAlphaR(gaussian_samples, options);
    metrics.local_alpha_training_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start).count();
    if (!std::isfinite(alpha) || alpha < 0.0)
    { result.status = EvidenceStatus::Unavailable; result.reason = "local-alpha-training-failed"; return result; }
    phase_start = std::chrono::steady_clock::now();
    result.fit = EstimateLocalGaussian(gaussian_samples, alpha, options, no_charge);
    metrics.local_mdpde_fitting_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start).count();
    result.fit.alpha_r = alpha;
    if (!result.fit.fit_result)
    { result.status = EvidenceStatus::Unavailable; result.reason = "local-mdpde-result-unavailable"; return result; }
    auto & fit = *result.fit.fit_result;
    result.fit.diagnostics = LocalFitDiagnostics{fit.status, fit.sigma_square, fit.diagnostics,
        fit.refinement, fit.failed_only_recovery};
    const auto reject = [&result](std::string reason) {
        result.status = EvidenceStatus::Unavailable;
        result.reason = std::move(reason);
        result.information.reset();
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
    const auto assessment = data_internal::AssessPostJointLocalResult(
        result, atom_id, joint, &peeling, raw);
    if (!assessment.valid) return reject(assessment.reason);
    return result;
}

void LogJointInfo(const std::string & message)
{
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message);
}

void LogJointDebug(const std::string & message)
{
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Debug, message);
}

void LogJointWarning(const std::string & message)
{
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Warning, message);
}

}

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
JointWorkflowFitHook & JointWorkflowFitHookForTesting()
{
    static thread_local JointWorkflowFitHook hook{};
    return hook;
}
#endif

PostJointLocalResult FitPostJointLocalTargetForResearch(
    const JointPostprocessingContext & postprocessing, const ModelObject & model, int atom_id,
    const PostFitPeelingResult & peeling, const FitOptions & options, JointPostprocessingMetrics & metrics)
{
    return FitPostJointLocalTarget(postprocessing, model, atom_id, peeling, options, metrics);
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
    RunPotentialSamplingWorkflow(map, model, first_stage_atoms, options.thread_size,
        " Sampling", options.quiet_mode);
    RunTwoStageFromPreparedSamples(model, workset.contributors, options);
}

void RunJointComponentWorkflow(
    MapObject & map, ModelObject & model, const FitOptions & options,
    JointPostprocessingMetrics * supplied_metrics)
{
    JointPostprocessingMetrics local_metrics;
    auto & postprocessing_metrics = supplied_metrics ? *supplied_metrics : local_metrics;
    postprocessing_metrics = {};
    const auto construction_start = std::chrono::steady_clock::now();
    const auto problem = BuildPureHaloSupportExclusionProblem(map, model);
    const auto construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - construction_start).count();
    if (!options.quiet_mode)
    {
        const auto & input = problem.Input();
        const auto targets = input.selection_domain ? input.selection_domain->target_indices.size() : 0;
        std::ostringstream message;
        message << "[Joint] Problem ready | observations=" << input.observations.size()
            << " | targets=" << targets
            << " | halo=" << input.atom_ids.size() - targets;
        LogJointInfo(message.str());
        std::ostringstream detail;
        detail << "[Joint] Problem build elapsed=" << construction_seconds << 's'
            << " | contributors=" << input.atom_ids.size();
        LogJointDebug(detail.str());
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
        auto fit = FitJointComponentsInitializedWithProgress(problem, initialization, observer);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if (const auto hook = JointWorkflowFitHookForTesting()) hook(fit);
#endif
        fit.costs.construction_seconds = construction_seconds;
        fit.costs.initialization_seconds = initialization_seconds;
        return CaptureJointAnalysisResult(fit);
    }();
    const auto postprocessing_start = std::chrono::steady_clock::now();
    postprocessing_metrics.rss_before_postprocessing_bytes = joint_component::CurrentRssBytes();
    ++postprocessing_metrics.joint_stage_mapping_builds;
    const auto joint_run_id=boost::uuids::to_string(boost::uuids::random_generator()());
    const auto finish_without_downstream = [&] {
        postprocessing_metrics.total_postprocessing_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - postprocessing_start).count();
    };
    if (!IsQualifiedJointAnalysisResult(snapshot))
    {
        LogJointWarning("[Joint] Qualification failed; preserving diagnostics and skipping downstream processing.");
        model.EditAnalysis().SetJointResult(std::move(snapshot));
        finish_without_downstream();
        return;
    }
    bool halo_reconstruction_unavailable{};
    std::optional<HaloMedianBuildResult> halo_median_build;
    const auto & original=problem.OriginalInput();
    const auto target_count=original.selection_domain ? original.selection_domain->target_indices.size() : 0;
    const bool original_has_halos=original.atom_ids.size()>target_count;
    try
    {
        halo_median_build=BuildHaloMedianReconstruction(model,problem,snapshot,joint_run_id);
        if (original_has_halos)
            snapshot.halo_median_reconstruction=halo_median_build->TakeSnapshot();
        postprocessing_metrics.halo_median_valid_target_count=halo_median_build->audit.valid_target_count;
        postprocessing_metrics.halo_median_geometry_local_halos=halo_median_build->audit.geometry_local_halos;
    }
    catch (const std::invalid_argument & error)
    {
        if (std::string_view(error.what()) != "halo-median-no-valid-targets") throw;
        if (original_has_halos)
        {
            postprocessing_metrics.halo_median_unavailable_halos=original.atom_ids.size()-target_count;
            snapshot.evidence.push_back({std::string(JointHaloMedianReconstructionEvidenceName),
                JointCheckStatus::Unavailable,JointEvidenceScope::AssembledGlobal,{}, {},error.what()});
            halo_reconstruction_unavailable=true;
            LogJointWarning("[Joint] Halo median reconstruction unavailable: " + std::string(error.what()));
        }
    }
    data_internal::ApplyJointStageEstimates(model, snapshot, joint_run_id);
    model.EditAnalysis().SetJointResult(std::move(snapshot));
    if (halo_reconstruction_unavailable)
    {
        finish_without_downstream();
        return;
    }
    const auto & joint_snapshot = *model.GetAnalysisView().GetJointResult();
    auto phase_start = std::chrono::steady_clock::now();
    auto contributor_universe=BuildJointPeelingContributorUniverse(model,problem,&joint_snapshot);
    const JointPostprocessingContext postprocessing(map, model, problem, &joint_snapshot,
        &postprocessing_metrics,halo_median_build ? &*halo_median_build : nullptr,
        std::move(contributor_universe));
    postprocessing_metrics.context_construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start).count();
    postprocessing_metrics.rss_after_context_bytes = joint_component::CurrentRssBytes();
    phase_start = std::chrono::steady_clock::now();
    auto peeling_updates = postprocessing.BuildPostFitPeelingBatch(original.selection_domain->target_indices);
    model.EditAnalysis().ApplyPostFitPeelingsBatch(std::move(peeling_updates));
    if (postprocessing_metrics.mapping_mismatch_samples || postprocessing_metrics.missing_contributor_samples)
    {
        std::ostringstream warning;
        warning << "[Joint] Peeling data integrity issue | mapping_mismatch="
            << postprocessing_metrics.mapping_mismatch_samples
            << " | missing_contributors=" << postprocessing_metrics.missing_contributor_samples;
        LogJointWarning(warning.str());
    }
    std::vector<const AtomObject *> joint_targets;
    for (const auto * atom : model.GetSelectedAtoms())
    {
        const auto local = AtomLocalPotentialView::For(*atom);
        const auto & stage = local.GetStageEstimate(FittingStage::Second);
        if (stage.source.method != EstimateMethod::JointComponents || stage.source.role != FittingRole::Target)
            continue;
        joint_targets.push_back(atom);
        const auto assessment = local.AssessLocalFittingPeelingRatio(1.0, 2.0);
        ++postprocessing_metrics.tail_target_count;
        postprocessing_metrics.tail_raw_sample_count += assessment.total_samples;
        postprocessing_metrics.tail_valid_sample_count += assessment.valid_samples;
        postprocessing_metrics.tail_unavailable_sample_count +=
            assessment.total_samples - assessment.valid_samples;
        postprocessing_metrics.tail_raw_response_sum += assessment.raw_response_sum;
        postprocessing_metrics.tail_peeling_response_sum += assessment.peeled_response_sum;
        if (assessment.ratio) ++postprocessing_metrics.tail_ratio_available_targets;
        if (!assessment.reason.empty()) ++postprocessing_metrics.tail_unavailable_reasons[assessment.reason];
        if (!assessment.total_samples)
        {
            ++postprocessing_metrics.tail_no_sample_targets;
            ++postprocessing_metrics.tail_coverage_distribution["no-tail-samples"];
        }
        else if (assessment.valid_samples == assessment.total_samples)
        {
            ++postprocessing_metrics.tail_full_coverage_targets;
            ++postprocessing_metrics.tail_coverage_distribution["complete"];
        }
        else if (assessment.valid_samples == 0)
        {
            ++postprocessing_metrics.tail_zero_coverage_targets;
            ++postprocessing_metrics.tail_coverage_distribution["zero"];
        }
        else
        {
            ++postprocessing_metrics.tail_partial_coverage_targets;
            const double coverage = static_cast<double>(assessment.valid_samples) /
                static_cast<double>(assessment.total_samples);
            const auto bucket = coverage <= 0.25 ? "(0,0.25]" : coverage <= 0.5 ? "(0.25,0.5]" :
                coverage <= 0.75 ? "(0.5,0.75]" : "(0.75,1)";
            ++postprocessing_metrics.tail_coverage_distribution[bucket];
        }
    }
    postprocessing_metrics.peeling_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start).count();
    postprocessing_metrics.rss_after_peeling_bytes = joint_component::CurrentRssBytes();
    std::vector<std::pair<int, PostJointLocalResult>> local_updates;
    for (const auto * atom : joint_targets)
    {
        const auto local = AtomLocalPotentialView::For(*atom);
        const auto & stage = local.GetStageEstimate(FittingStage::Second);
        if (stage.source.method != EstimateMethod::JointComponents || stage.source.role != FittingRole::Target ||
            stage.source.component_id.empty()) continue;
        const auto & peeling = local.GetPostFitPeeling();
        PostJointLocalResult result;
        if (peeling) result = FitPostJointLocalTarget(
            postprocessing, model, atom->GetSerialID(), *peeling, options, postprocessing_metrics);
        else
        {
            result.source = stage.source;
            result.status = EvidenceStatus::Unavailable;
            result.reason = "peeling-unavailable";
            result.sample_count = local.GetRawSamplingEntriesView().size();
            if (stage.point)
                result.fixed_charge = {stage.point->GetOffset(), stage.point->GetWidth()};
        }
        local_updates.emplace_back(atom->GetSerialID(), std::move(result));
    }
    model.EditAnalysis().ApplyPostJointLocalResultsBatch(std::move(local_updates));
    postprocessing_metrics.estimated_peak_postprocessing_payload_bytes +=
        joint_targets.capacity() * sizeof(const AtomObject *);
    postprocessing_metrics.rss_after_local_bytes = joint_component::CurrentRssBytes();
    if (!options.quiet_mode) LogJointDebug(BuildSecondStageSpotSummary(model));
    RunJointGroupPotentialFitting(model, options, &postprocessing_metrics);
    postprocessing_metrics.rss_after_group_bytes = joint_component::CurrentRssBytes();
    if (!options.quiet_mode)
    {
        std::ostringstream message;
        message << "[Joint] Halo reconstruction | donors="
            << postprocessing_metrics.halo_median_valid_target_count
            << " | reconstructed=" << postprocessing_metrics.halo_median_geometry_local_halos
            << " | unavailable=" << postprocessing_metrics.halo_median_unavailable_halos;
        LogJointInfo(message.str());
        std::ostringstream tail_message;
        tail_message << "[Joint] Peeling | targets="
            << postprocessing_metrics.tail_target_count
            << " | valid_samples=" << postprocessing_metrics.tail_valid_sample_count << '/'
            << postprocessing_metrics.tail_raw_sample_count
            << " | unavailable=" << postprocessing_metrics.tail_unavailable_sample_count
            << " | coverage=" << postprocessing_metrics.tail_full_coverage_targets << " complete/"
            << postprocessing_metrics.tail_partial_coverage_targets << " partial/"
            << postprocessing_metrics.tail_zero_coverage_targets << " zero/"
            << postprocessing_metrics.tail_no_sample_targets << " no-samples";
        LogJointInfo(tail_message.str());
        std::ostringstream tail_detail;
        tail_detail << "[Joint] Peeling diagnostics | ratio_available="
            << postprocessing_metrics.tail_ratio_available_targets
            << " | raw_sum=" << postprocessing_metrics.tail_raw_response_sum
            << " | peeled_sum=" << postprocessing_metrics.tail_peeling_response_sum
            << " | map_boundary_samples=" << postprocessing_metrics.map_geometry_boundary_samples
            << " | original_support_boundary_samples=" << postprocessing_metrics.original_support_boundary_samples
            << " | fitting_mask_only_samples=" << postprocessing_metrics.fitting_mask_only_intersection_samples
            << " | mapping_mismatch_samples=" << postprocessing_metrics.mapping_mismatch_samples
            << " | missing_contributor_samples=" << postprocessing_metrics.missing_contributor_samples;
        for (const auto & [bucket, count] : postprocessing_metrics.tail_coverage_distribution)
            tail_detail << " | coverage_" << bucket << '=' << count;
        for (const auto & [reason, count] : postprocessing_metrics.tail_unavailable_reasons)
            tail_detail << " | " << reason << '=' << count;
        LogJointDebug(tail_detail.str());
        std::ostringstream reconstruction_message;
        reconstruction_message << "[Joint] Reconstruction-only lookup: atoms="
            << postprocessing_metrics.reconstruction_lookup_atom_count
            << " | node_evaluations=" << postprocessing_metrics.reconstruction_only_node_evaluation_count
            << " | cache_hits=" << postprocessing_metrics.reconstruction_node_cache_hit_count
            << " | cache_peak_entries=" << postprocessing_metrics.reconstruction_node_cache_peak_entries
            << " | estimated_workspace_bytes=" << postprocessing_metrics.estimated_reconstruction_workspace_bytes;
        LogJointDebug(reconstruction_message.str());
    }
    postprocessing_metrics.total_postprocessing_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - postprocessing_start).count();
    if (!options.quiet_mode) LogJointInfo("[Joint] Completed.");
}

} // namespace rhbm_gem::core::detail

namespace rhbm_gem::core {

void RunPotentialFittingWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    model.ApplyElementSelection(Element::HYDROGEN, true);
    if (model.GetSelectedAtoms().empty())
        throw std::invalid_argument("Potential fitting requires at least one selected non-hydrogen target.");

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
