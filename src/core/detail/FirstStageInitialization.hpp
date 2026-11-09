#pragma once
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include "joint_component/JointComponentEstimator.hpp"
#include "joint_component/JointProgress.hpp"
#include <vector>

namespace rhbm_gem { class MapObject; class ModelAnalysisEditor; }

namespace rhbm_gem::core::detail {
struct FittingWorkset;
double EstimateLocalAlphaForAtom(const AtomObject &, const FitOptions &, FittingStage);
void TrainLocalAlphaForAtom(ModelAnalysisEditor &, const FitOptions &, FittingStage, AtomObject &);
void ApplyJointSeedFallback(JointInitialization &);
LocalGaussianResult FitFirstStageAtom(const AtomObject &, const FitOptions &);
void RunBatchFirstStageFromPreparedSamples(
    ModelObject &, const std::vector<AtomObject *> &, const FitOptions &);
JointInitialization RunJointFirstStageInitializationFromPreparedSamples(
    ModelObject &, const FittingWorkset &, const FitOptions &);
JointInitialization PrepareJointWidthInitialization(
    MapObject &, ModelObject &, const FittingWorkset &, const FitOptions &);
JointFitResult FitJointComponentsInitializedWithProgress(
    const JointProblem &, const JointInitialization &,
    const joint_component::JointProgressObserver &);
}
