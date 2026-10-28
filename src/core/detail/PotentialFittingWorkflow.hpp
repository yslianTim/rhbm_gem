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

} // namespace rhbm_gem::core::detail
