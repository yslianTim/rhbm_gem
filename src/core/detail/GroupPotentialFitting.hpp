#pragma once
#include "joint_component/PostFitPeeling.hpp"
namespace rhbm_gem {
class ModelObject;
namespace core {
struct FitOptions;
namespace detail {
void RunGroupAlphaTraining(ModelObject &, const FitOptions &);
void RunJointGroupPotentialFitting(ModelObject &, const FitOptions &, JointPostprocessingMetrics * = nullptr);
}
}
}
