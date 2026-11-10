#pragma once

#include "core/estimator/GaussianEstimator.hpp"
#include "joint_component/PostFitPeeling.hpp"

#include <string>

namespace rhbm_gem { struct JointAnalysisResult; }

namespace rhbm_gem::core {
void RunPotentialFittingWorkflow(MapObject &, ModelObject &, const FitOptions &);

std::string BuildJointTargetSpotSummary(const ModelObject &);
}

namespace rhbm_gem::core::detail {

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
// Exercises command outcomes with a qualified or rejected postprocessing snapshot.
using JointPostprocessingSnapshotHook = void (*)(JointAnalysisResult &);
JointPostprocessingSnapshotHook & JointPostprocessingSnapshotHookForTesting();
#endif

PostJointLocalResult FitPostJointLocalTarget(
    const JointPostprocessingContext &, const ModelObject &, int, const PostFitPeelingResult &,
    const FitOptions &, JointPostprocessingMetrics &);

void RunPotentialFittingWorkflowForSelection(
    MapObject & map, ModelObject & model, const FitOptions & options,
    JointPostprocessingMetrics * metrics = nullptr);

// Completes JointStage publication and runs GroupingStage only when downstream inputs are ready.
void CompletePotentialFittingWorkflow(
    MapObject &, ModelObject &, const JointProblem &, JointAnalysisResult,
    const FitOptions &, JointPostprocessingMetrics &);

} // namespace rhbm_gem::core::detail
