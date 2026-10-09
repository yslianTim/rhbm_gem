#pragma once
#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <string_view>

namespace rhbm_gem {

struct JointMapNormalization
{
    bool requested{}, applied{};
    double divisor{1};
};
struct JointSoftwareProvenance
{
    std::string version, source_sha256, configuration_sha256, build_sha256;
};
struct JointSolverProvenance
{
    std::string search_method;
    std::optional<std::size_t> fixed_neighbor_core_atoms;
    // Legacy provenance only. Current FixedNeighbor v2 has intrinsic OneAccepted behavior.
    std::optional<std::string> fixed_neighbor_local_work;
    std::optional<std::string> contract_version;
    std::optional<std::string> sparse_backend;
    // Historical multi-route provenance only. These are not current solver controls.
    std::optional<std::string> preconditioner;
    std::optional<std::string> operator_rank_mode,operator_rank_backend;
    std::optional<int> operator_pcg_iterations,operator_damping_trials;
    std::optional<double> operator_rank_budget_seconds;
    std::optional<std::size_t> operator_rank_budget_entries,operator_rank_budget_workspace_bytes;
    std::optional<std::size_t> schwarz_core_atoms,schwarz_overlap_hops,schwarz_max_block_atoms,
        schwarz_storage_bytes,schwarz_scratch_bytes;
    std::optional<std::string> fixed_neighbor_policy_version,fixed_neighbor_order,fixed_neighbor_local_search;
    std::optional<std::size_t> fixed_neighbor_maximum_sweeps;
};
inline constexpr std::string_view JointSolverProvenanceContractVersion="joint-solver-provenance-v3";
inline constexpr std::string_view JointSolverProvenanceHistoricalContractVersion="joint-solver-provenance-v2";
inline constexpr std::string_view FixedNeighborLegacyPolicyContractVersion="fixed-neighbor-production-v1";
inline constexpr std::string_view FixedNeighborPolicyContractVersion="fixed-neighbor-production-v2";
struct JointAnalysisMetadata
{
    std::string model_path, map_path;
    std::array<int,3> grid_size{};
    std::array<double,3> grid_spacing{}, origin{};
    bool simulation{};
    // Unknown for callers that supply only an in-memory problem.
    std::optional<JointMapNormalization> map_normalization;
    std::optional<std::string> model_sha256, map_sha256;
    std::optional<JointSoftwareProvenance> software;
    std::optional<JointSolverProvenance> solver;
};
struct JointAnalysisComponent : JointComponentData
{
    JointCheckStatus runtime_convergence{JointCheckStatus::Unavailable};
    JointCheckStatus target_runtime_convergence{JointCheckStatus::NotRun};
};
// A saved outcome, not an audit/restart bundle. Convergence values are captured
// from the runtime result once and are never inferred by storage or readers.
struct JointAnalysisResult
{
    std::string parameterization_contract{"full-abc-v1"};
    std::optional<JointParameterLayout> layout;
    JointAnalysisMetadata metadata;
    std::vector<std::string> atom_ids, row_ids;
    std::optional<JointSelectionDomain> selection_domain;
    JointInitialization initialization;
    JointFitCosts costs;
    std::vector<JointAnalysisComponent> components;
    std::optional<JointState> assembled_state;
    std::optional<double> objective;
    std::vector<bool> available_row_mask;
    std::vector<JointCheck> evidence;
    std::vector<JointRankEvidence> ranks;
    bool search_completed{};
    double observation_scale{1};
    JointCheckStatus runtime_convergence{JointCheckStatus::Unavailable};
    JointCheckStatus regular_certificate{JointCheckStatus::NotRun};
    std::optional<JointTargetEvidence> target_evidence;
    JointCheckStatus target_runtime_convergence{JointCheckStatus::NotRun};
};
}
