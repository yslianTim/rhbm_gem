#pragma once

#include "FixedNeighborPolicy.hpp"
#include <cstddef>
#include <string_view>

namespace rhbm_gem::core::joint_component {

struct JointSolverConfiguration
{
    std::size_t fixed_neighbor_core_atoms{};
};

JointSolverConfiguration ResolveJointSolverConfiguration(const FixedNeighborSearchPolicy & policy);
std::string_view FixedNeighborBlockOrderName(FixedNeighborBlockOrder order);

} // namespace rhbm_gem::core::joint_component
