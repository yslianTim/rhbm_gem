#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include "joint_component/PostFitPeeling.hpp"
#include "JointSeedSampling.hpp"

#include <vector>

namespace rhbm_gem::core { struct JointFitResult; }
namespace rhbm_gem { struct JointAnalysisResult; }

namespace rhbm_gem::core::detail {

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using JointWorkflowFitHook = void (*)(JointFitResult &);
JointWorkflowFitHook & JointWorkflowFitHookForTesting();
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
    JointPostprocessingMetrics * metrics = nullptr,
    const JointSeedSamplingPolicy & policy = {}, JointSeedSamplingReport * sampling_report = nullptr);

void RunJointPostprocessing(
    MapObject &, ModelObject &, const JointProblem &, JointAnalysisResult,
    const FitOptions &, JointPostprocessingMetrics &);

} // namespace rhbm_gem::core::detail
