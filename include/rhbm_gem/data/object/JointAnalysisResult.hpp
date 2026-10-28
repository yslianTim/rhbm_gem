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
    std::optional<std::string> contract_version;
    std::optional<std::string> sparse_backend;
    std::optional<std::string> fixed_neighbor_policy_version,fixed_neighbor_order,fixed_neighbor_local_search;
    std::optional<std::size_t> fixed_neighbor_maximum_sweeps;
};
inline constexpr std::string_view JointSolverProvenanceContractVersion="joint-solver-provenance-v5";
inline constexpr std::string_view FixedNeighborPolicyContractVersion="fixed-neighbor-production-v4";
inline constexpr std::string_view FixedNeighborBoundedLocalSearchName="active-set-lm";
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
// A saved outcome with optional parameter-only candidate diagnostics, not an
// audit/restart bundle. Convergence is captured and never inferred by storage.
struct JointAnalysisResult
{
    std::string parameterization_contract{JointPureHaloSupportExclusionParameterizationContract};
    JointAnalysisMetadata metadata;
    std::vector<std::string> atom_ids, row_ids;
    std::vector<std::string> fitting_row_ids;
    std::vector<std::size_t> fitting_to_original_contributors;
    std::vector<std::size_t> fitting_to_original_rows;
    std::optional<JointSelectionDomain> selection_domain;
    JointInitialization initialization;
    JointFitCosts costs;
    std::vector<JointAnalysisComponent> components;
    std::optional<JointState> assembled_state;
    std::optional<double> objective;
    std::vector<bool> available_row_mask;
    std::vector<JointCheck> evidence;
    bool search_completed{};
    double observation_scale{1};
    JointCheckStatus runtime_convergence{JointCheckStatus::Unavailable};
    std::optional<JointTargetEvidence> target_evidence;
    JointCheckStatus target_runtime_convergence{JointCheckStatus::NotRun};
    std::optional<JointHaloMedianSnapshot> halo_median_reconstruction;
};

inline constexpr std::string_view JointHaloMedianReconstructionEvidenceName =
    "halo-median-reconstruction";

inline bool IsQualifiedJointAnalysisResult(const JointAnalysisResult & result)
{
    return result.initialization.valid && result.search_completed &&
        result.runtime_convergence == JointCheckStatus::Passed &&
        result.target_runtime_convergence == JointCheckStatus::Passed;
}

inline bool HasUnavailableJointHaloMedianReconstruction(const JointAnalysisResult & result)
{
    for (const auto & check : result.evidence)
        if (check.name == JointHaloMedianReconstructionEvidenceName &&
            check.status == JointCheckStatus::Unavailable && !check.reason.empty())
            return true;
    return false;
}
}
