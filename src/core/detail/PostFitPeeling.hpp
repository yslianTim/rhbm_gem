#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <map>
#include <span>
namespace rhbm_gem::core::detail {
// Reads geometry, immutable support, raw samples and fixed Second states only.
std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const MapObject &, const ModelObject &, const JointProblem &,
    std::optional<std::span<const std::size_t>> outputs = std::nullopt);
}
