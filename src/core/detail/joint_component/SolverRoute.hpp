#pragma once

#include <optional>
#include <string_view>

namespace rhbm_gem::core::joint_component {

enum class SparseBackend {Eigen,Spqr};
enum class SearchMethod {LegacyCompact,OperatorPcg};
enum class PreconditionerKind {Identity,Diagonal,Schwarz};
enum class FreeDesignRankBackend {Dense,SpqrBounds};
enum class OperatorRankMode {Auto,Dense,SpqrBounds};

struct SearchPolicy
{
    SearchMethod method{SearchMethod::LegacyCompact};
    PreconditionerKind preconditioner{PreconditionerKind::Schwarz};
    OperatorRankMode operator_rank{OperatorRankMode::Auto};
    int pcg_iterations{-1},damping_trials{20};
};

struct JointSolverRoute
{
    SparseBackend sparse_backend{};
    SearchMethod search_method{};
    std::optional<PreconditionerKind> preconditioner;
};

SparseBackend ActiveSparseBackend();
std::optional<FreeDesignRankBackend> ResolveOperatorRankBackend(OperatorRankMode,SparseBackend);
JointSolverRoute ResolveJointSolverRoute(const SearchPolicy & policy);
std::string_view SparseBackendName(SparseBackend backend);
std::string_view SearchMethodName(SearchMethod method);
std::string_view PreconditionerName(PreconditionerKind preconditioner);
std::string_view OperatorRankModeName(OperatorRankMode mode);

} // namespace rhbm_gem::core::joint_component
