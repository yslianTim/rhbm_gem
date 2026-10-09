#pragma once

#include "FixedNeighborPolicy.hpp"
#include <cstddef>
#include <string_view>

namespace rhbm_gem::core::joint_component {

enum class SparseBackend {Eigen,Spqr};
struct JointSolverConfiguration
{
    SparseBackend sparse_backend{};
    std::size_t fixed_neighbor_core_atoms{};
};

SparseBackend ActiveSparseBackend();
JointSolverConfiguration ResolveJointSolverConfiguration(const FixedNeighborSearchPolicy & policy);
std::string_view SparseBackendName(SparseBackend backend);
std::string_view FixedNeighborBlockOrderName(FixedNeighborBlockOrder order);

} // namespace rhbm_gem::core::joint_component
