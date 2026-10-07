#pragma once

#include "FixedNeighborPolicy.hpp"
#include <optional>
#include <cstddef>
#include <string_view>

namespace rhbm_gem::core::joint_component {

enum class SparseBackend {Eigen,Spqr};
enum class SearchMethod {LegacyCompact,OperatorPcg,FixedNeighbor};
enum class PreconditionerKind {Identity,Diagonal,Schwarz};
enum class FreeDesignRankBackend {Dense,SpqrBounds};
enum class OperatorRankMode {Auto,Dense,SpqrBounds};
struct RankBudget
{
    double seconds{120};
    std::size_t entries{100000000},workspace_bytes{256*1024*1024};
};
struct OperatorRankPolicy
{
    OperatorRankMode mode{OperatorRankMode::Auto};
    RankBudget budget{};
};
struct SchwarzPolicy
{
    std::size_t core_atoms{128},overlap_hops{1},max_block_atoms{512};
    std::size_t storage_bytes{512ULL*1024*1024},scratch_bytes{256ULL*1024*1024};
};

struct SearchPolicy
{
    SearchMethod method{SearchMethod::LegacyCompact};
    PreconditionerKind preconditioner{PreconditionerKind::Schwarz};
    OperatorRankPolicy operator_rank{};
    SchwarzPolicy schwarz{};
    int pcg_iterations{-1},damping_trials{20};
    FixedNeighborSearchPolicy fixed_neighbor{};
};

struct JointSolverRoute
{
    SparseBackend sparse_backend{};
    SearchMethod search_method{};
    std::optional<PreconditionerKind> preconditioner;
    std::optional<std::size_t> fixed_neighbor_core_atoms;
    std::optional<FixedNeighborLocalWork> fixed_neighbor_local_work;
};

SparseBackend ActiveSparseBackend();
std::optional<FreeDesignRankBackend> ResolveOperatorRankBackend(OperatorRankMode,SparseBackend);
JointSolverRoute ResolveJointSolverRoute(const SearchPolicy & policy);
std::string_view SparseBackendName(SparseBackend backend);
std::string_view SearchMethodName(SearchMethod method);
std::string_view SearchMethodToken(SearchMethod method);
std::string_view PreconditionerName(PreconditionerKind preconditioner);
std::string_view FixedNeighborBlockOrderName(FixedNeighborBlockOrder order);
std::string_view FixedNeighborLocalWorkName(FixedNeighborLocalWork work);
std::string_view OperatorRankModeName(OperatorRankMode mode);
std::string_view FreeDesignRankBackendName(FreeDesignRankBackend backend);

} // namespace rhbm_gem::core::joint_component
