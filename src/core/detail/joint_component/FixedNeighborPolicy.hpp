#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include "ConvergenceCriteria.hpp"

namespace rhbm_gem::core::joint_component {

inline constexpr std::size_t ProductionFixedNeighborCoreAtoms=12;
inline constexpr std::size_t ProductionFixedNeighborMaximumSweeps=200;

// The production search accepts its calibrated bound and stability-work budget.
// Component identity is supplied internally for failure-only logs.
struct FixedNeighborSearchPolicy
{
    std::size_t component_index{};
    std::string component_id;
    std::function<bool(int)> claim_parameter_stability_block;
    std::optional<double> eta_upper_bound;
};

} // namespace rhbm_gem::core::joint_component
