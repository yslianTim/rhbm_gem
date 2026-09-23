#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <map>
#include <span>
namespace rhbm_gem::core::detail {
std::map<int, StageUncertainty> ComputeJointUncertainty(const JointProblem &, const JointAnalysisResult &,
    std::optional<std::span<const std::size_t>> outputs = std::nullopt);
GroupParameterEvidence BuildJointParameterEvidence(const LocalStageEstimate &);
}
