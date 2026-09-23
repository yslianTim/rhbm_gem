#pragma once
#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>

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
};
struct JointAnalysisComponent : JointComponentData
{
    JointCheckStatus runtime_convergence{JointCheckStatus::Unavailable};
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
};
}
