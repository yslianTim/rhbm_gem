#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include "joint_component/JointComponentEstimator.hpp"

#include <vector>

namespace rhbm_gem::core::detail {
// Contributors and target flags share the Joint Problem's fitting-domain order.
struct InitialStageWorkset
{
    std::vector<AtomObject *> contributors;
    std::vector<bool> target_mask;
};

InitialStageWorkset BuildInitialStageWorkset(ModelObject &, const JointProblem &);

double EstimateFirstStageAlpha(const LocalPotentialSampleList &, const FitOptions &);
void ApplyJointSeedFallback(JointInitialization &);
LocalGaussianResult FitFirstStageAtom(const AtomObject &, const FitOptions &);
JointInitialization FitInitialStageFromPreparedSamples(
    ModelObject &, const InitialStageWorkset &, const FitOptions &);
JointInitialization RunInitialStage(
    MapObject &, ModelObject &, const InitialStageWorkset &, const FitOptions &);
}
