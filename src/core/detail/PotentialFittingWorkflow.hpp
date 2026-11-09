#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include "joint_component/PostFitPeeling.hpp"

#include <vector>

namespace rhbm_gem { struct JointAnalysisResult; }

namespace rhbm_gem::core::detail {

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
// Exercises command outcomes with a qualified or rejected postprocessing snapshot.
using JointPostprocessingSnapshotHook = void (*)(JointAnalysisResult &);
JointPostprocessingSnapshotHook & JointPostprocessingSnapshotHookForTesting();
#endif

PostJointLocalResult FitPostJointLocalTarget(
    const JointPostprocessingContext &, const ModelObject &, int, const PostFitPeelingResult &,
    const FitOptions &, JointPostprocessingMetrics &);

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

void RunJointPostprocessing(
    MapObject &, ModelObject &, const JointProblem &, JointAnalysisResult,
    const FitOptions &, JointPostprocessingMetrics &);

} // namespace rhbm_gem::core::detail
