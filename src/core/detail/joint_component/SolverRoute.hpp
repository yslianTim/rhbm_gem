#pragma once

#include <optional>
#include <string_view>

namespace rhbm_gem::core::joint_component {

enum class SparseBackend {Eigen,Spqr};
enum class SearchMethod {LegacyCompact,OperatorPcg};
enum class PreconditionerKind {Identity,Diagonal,Schwarz};

struct SearchPolicy
{
    SearchMethod method{SearchMethod::LegacyCompact};
    PreconditionerKind preconditioner{PreconditionerKind::Schwarz};
    int pcg_iterations{-1},damping_trials{20};
};

struct JointSolverRoute
{
    SparseBackend sparse_backend{};
    SearchMethod search_method{};
    std::optional<PreconditionerKind> preconditioner;
};

SparseBackend ActiveSparseBackend();
JointSolverRoute ResolveJointSolverRoute(const SearchPolicy & policy);
std::string_view SparseBackendName(SparseBackend backend);
std::string_view SearchMethodName(SearchMethod method);
std::string_view PreconditionerName(PreconditionerKind preconditioner);

} // namespace rhbm_gem::core::joint_component
