#pragma once

#include "FixedNeighborPolicy.hpp"
#include <cstddef>

namespace rhbm_gem::core::joint_component {

struct JointSolverConfiguration
{
    std::size_t fixed_neighbor_core_atoms{};
    bool bounded_width_profile{};
};

JointSolverConfiguration ResolveJointSolverConfiguration(const FixedNeighborSearchPolicy & policy);

} // namespace rhbm_gem::core::joint_component
