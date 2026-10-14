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
#include "detail/second_stage/IterationProcess.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <chrono>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem::core::detail {
namespace {
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
    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model));
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
