#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <cstddef>
#include <string>
#include <vector>

namespace rhbm_gem::core::joint_component {
inline constexpr std::size_t FullABCHaloMinimumObservations=3;
inline constexpr std::size_t FullABCHaloMinimumDistinctRadialInformation=3;

struct FullABCHaloEligibility
{
    std::size_t support_rows{},distinct_squared_distances{};
    bool eligible{};
    std::string reason;
};

// The support vector is scratch storage and may be reordered while counting.
// Each entry must represent one unique observation row for this atom.
FullABCHaloEligibility EvaluateFullABCHaloEligibility(
    std::vector<JointSupport> & support,
    std::size_t minimum_observations=FullABCHaloMinimumObservations);
}
