#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include "detail/FirstStageInitialization.hpp"
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

namespace rhbm_gem::core {
void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    if (options.estimator != PotentialEstimator::TWO_STAGE)
        throw std::invalid_argument("Joint fitting requires the map-aware workflow.");
    model_object.EditAnalysis().InitializeLocalFittingSeedModels();

    const detail::FittingWorkset workset{model_object.GetSelectedAtoms(),
        std::vector<bool>(model_object.GetSelectedAtomCount(), true)};
    detail::RunFirstStage(model_object, workset, options, detail::FirstStageMode::ExistingSamplesBatch);

    detail::RunSecondStageIterations(model_object, options);

    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model_object));
    detail::RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
}

void RunPotentialFittingWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    if (options.estimator == PotentialEstimator::TWO_STAGE)
    {
        model.EditAnalysis().InitializeFromSelection();
        RunPotentialSamplingWorkflow(map, model, options.sampling_method, options.thread_size);
        RunPotentialFittingWorkflow(model, options);
        return;
    }
    if (options.sampling_method != SphereSamplingMethod::FibonacciDeterministic)
        throw std::invalid_argument("Joint initialization requires Fibonacci sampling.");
    const auto construction_start = std::chrono::steady_clock::now();
    const auto problem = BuildJointProblem(map, model);
    const auto construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - construction_start).count();
    const auto initialization_start = std::chrono::steady_clock::now();
    const auto workset = detail::MakeJointFittingWorkset(model, problem);
    model.EditAnalysis().InitializeFromSelection();
    const auto initialization = detail::RunFirstStage(model, workset, options, detail::FirstStageMode::SampleContributorsIsolated, &map);
    const auto initialization_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - initialization_start).count();
    auto snapshot = [&] {
        auto fit = FitJointComponents(problem, initialization.b);
        fit.costs.construction_seconds = construction_seconds;
        fit.costs.initialization_seconds = initialization_seconds;
        fit.initialization.atoms = initialization.atoms;
        fit.initialization.data_scope = initialization.data_scope;
        return CaptureJointAnalysisResult(fit);
    }();
    data_internal::ApplyJointStageEstimates(model, snapshot, boost::uuids::to_string(boost::uuids::random_generator()()));
    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model));
    for (auto & [id, peeling] : detail::BuildPostFitPeelingSamples(map, model, problem, problem.Input().selection_domain->target_indices))
        model.EditAnalysis().SetAtomPostFitPeeling(*model.FindAtomPtr(id), std::move(peeling));
    for (auto & [id, uncertainty] : detail::ComputeJointUncertainty(problem, snapshot, problem.Input().selection_domain->target_indices))
    {
        auto & atom = *model.FindAtomPtr(id);
        auto stage = AtomLocalPotentialView::For(atom).GetStageEstimate(FittingStage::Second);
        stage.uncertainty = std::move(uncertainty);
        model.EditAnalysis().SetAtomStageEstimate(FittingStage::Second, atom, stage);
        model.EditAnalysis().SetAtomGroupEvidence(atom, detail::BuildJointParameterEvidence(stage));
    }
    model.EditAnalysis().SetJointResult(std::move(snapshot));
    RunGroupPotentialFitting(model, options);
}

} // namespace rhbm_gem::core
