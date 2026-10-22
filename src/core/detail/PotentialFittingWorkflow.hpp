#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include "joint_component/PostFitPeeling.hpp"

#include <vector>

namespace rhbm_gem::core { struct JointFitResult; }

namespace rhbm_gem::core::detail {

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using JointWorkflowFitHook = void (*)(JointFitResult &);
JointWorkflowFitHook & JointWorkflowFitHookForTesting();
#endif

void RunTwoStageFromPreparedSamples(
    ModelObject & model,
    const FitOptions & options);

void RunTwoStageFromPreparedSamples(
    ModelObject & model,
    const std::vector<AtomObject *> & atoms,
    const FitOptions & options);

void RunJointComponentWorkflow(
    MapObject & map, ModelObject & model, const FitOptions & options,
    JointPostprocessingMetrics * metrics = nullptr);

// Research-only entry point for replaying local post-processing from a fixed
// Joint snapshot. It does not build or alter a Joint objective.
PostJointLocalResult FitPostJointLocalTargetForResearch(
    const JointPostprocessingContext &, const ModelObject &, int,
    const PostFitPeelingResult &, const FitOptions &, JointPostprocessingMetrics &);

} // namespace rhbm_gem::core::detail
