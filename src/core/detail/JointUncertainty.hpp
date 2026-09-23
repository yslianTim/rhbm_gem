#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <map>
namespace rhbm_gem::core { struct FitOptions; }
namespace rhbm_gem::core::detail {
std::map<int, StageUncertainty> ComputeJointUncertainty(const JointProblem &, const JointAnalysisResult &);
GroupParameterEvidence BuildJointParameterEvidence(const LocalStageEstimate &);
void RunJointGroupPotentialFitting(ModelObject &, const FitOptions &);
}
