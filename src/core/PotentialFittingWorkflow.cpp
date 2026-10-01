#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include "detail/FirstStageInitialization.hpp"
#include "detail/FittingWorkset.hpp"
#include "detail/PotentialFittingWorkflow.hpp"
#include "detail/StageSummary.hpp"
#include "detail/PostFitPeeling.hpp"
#include "detail/JointUncertainty.hpp"
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
#include <stdexcept>
#include <vector>

namespace rhbm_gem::core::detail {

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
    RunPotentialSamplingWorkflow(map, model, first_stage_atoms, options.sampling_method, options.thread_size);
    RunTwoStageFromPreparedSamples(model, workset.contributors, options);
}

void RunJointComponentWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    if (options.sampling_method != SphereSamplingMethod::FibonacciDeterministic)
        throw std::invalid_argument("Joint initialization requires Fibonacci sampling.");
    const auto construction_start = std::chrono::steady_clock::now();
    const auto problem = BuildJointProblem(map, model);
    const auto construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - construction_start).count();
    const auto initialization_start = std::chrono::steady_clock::now();
    const auto workset = MakeJointFittingWorkset(model, problem);
    model.EditAnalysis().InitializeFromSelection();
    const auto first_stage_atoms = CollectFirstStageAtoms(workset);
    RunPotentialSamplingWorkflow(map, model, first_stage_atoms, options.sampling_method, 1);
    const auto initialization = RunJointFirstStageInitializationFromPreparedSamples(model, workset, options);
    const auto initialization_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - initialization_start).count();
    auto snapshot = [&] {
        auto fit = FitJointComponents(problem, initialization.b);
        fit.costs.construction_seconds = construction_seconds;
        fit.costs.initialization_seconds = initialization_seconds;
        fit.initialization = initialization;
        return CaptureJointAnalysisResult(fit);
    }();
    data_internal::ApplyJointStageEstimates(model, snapshot, boost::uuids::to_string(boost::uuids::random_generator()()));
    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model));
    for (auto & [id, peeling] : BuildPostFitPeelingSamples(map, model, problem, problem.Input().selection_domain->target_indices, &snapshot))
        model.EditAnalysis().SetAtomPostFitPeeling(*model.FindAtomPtr(id), std::move(peeling));
    for (auto & [id, uncertainty] : ComputeJointUncertainty(problem, snapshot, problem.Input().selection_domain->target_indices))
    {
        auto & atom = *model.FindAtomPtr(id);
        auto stage = AtomLocalPotentialView::For(atom).GetStageEstimate(FittingStage::Second);
        stage.uncertainty = std::move(uncertainty);
        model.EditAnalysis().SetAtomStageEstimate(FittingStage::Second, atom, stage);
        model.EditAnalysis().SetAtomGroupEvidence(atom, BuildJointParameterEvidence(stage));
    }
    model.EditAnalysis().SetJointResult(std::move(snapshot));
    RunGroupPotentialFitting(model, options);
}

} // namespace rhbm_gem::core::detail

namespace rhbm_gem::core {

void RunPotentialFittingWorkflow(ModelObject & model, const FitOptions & options)
{
    if (options.estimator != PotentialEstimator::TWO_STAGE)
        throw std::invalid_argument("Joint fitting requires the map-aware workflow.");
    const auto workset = detail::MakeTwoStageFittingWorkset(model);
    detail::RunTwoStageFromPreparedSamples(model, workset.contributors, options);
}

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
