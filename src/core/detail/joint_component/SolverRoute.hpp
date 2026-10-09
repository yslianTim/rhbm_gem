#pragma once

#include "FixedNeighborPolicy.hpp"
#include <optional>
#include <cstddef>
#include <string_view>

namespace rhbm_gem::core::joint_component {

enum class SparseBackend {Eigen,Spqr};
enum class SearchMethod {LegacyCompact,FixedNeighbor};
enum class FreeDesignRankBackend {Dense,SpqrBounds};
struct RankBudget
{
    double seconds{120};
    std::size_t entries{100000000},workspace_bytes{256*1024*1024};
};
struct SearchPolicy
{
    SearchMethod method{SearchMethod::LegacyCompact};
    FixedNeighborSearchPolicy fixed_neighbor{};
};

struct JointSolverRoute
{
    SparseBackend sparse_backend{};
    SearchMethod search_method{};
    std::optional<std::size_t> fixed_neighbor_core_atoms;
};

SparseBackend ActiveSparseBackend();
JointSolverRoute ResolveJointSolverRoute(const SearchPolicy & policy);
std::string_view SparseBackendName(SparseBackend backend);
std::string_view SearchMethodName(SearchMethod method);
std::string_view SearchMethodToken(SearchMethod method);
std::string_view FixedNeighborBlockOrderName(FixedNeighborBlockOrder order);
std::string_view FreeDesignRankBackendName(FreeDesignRankBackend backend);

} // namespace rhbm_gem::core::joint_component
