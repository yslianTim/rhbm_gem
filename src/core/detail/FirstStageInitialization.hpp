#pragma once
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include "joint_component/JointComponentEstimator.hpp"
#include "joint_component/JointProgress.hpp"

namespace rhbm_gem::core::detail {
struct FittingWorkset;
double EstimateFirstStageAlpha(const LocalPotentialSampleList &, const FitOptions &);
void ApplyJointSeedFallback(JointInitialization &);
LocalGaussianResult FitFirstStageAtom(const AtomObject &, const FitOptions &);
JointInitialization RunJointFirstStageInitializationFromPreparedSamples(
    ModelObject &, const FittingWorkset &, const FitOptions &);
JointInitialization PrepareJointWidthInitialization(
    MapObject &, ModelObject &, const FittingWorkset &, const FitOptions &);
JointFitResult FitJointComponentsInitializedWithProgress(
    const JointProblem &, const JointInitialization &,
    const joint_component::JointProgressObserver &);
}
