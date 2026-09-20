#pragma once
#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>

namespace rhbm_gem {

struct JointAnalysisMetadata
{
    std::string model_path, map_path;
    std::array<int,3> grid_size{};
    std::array<double,3> grid_spacing{}, origin{};
    bool map_normalization_applied{}, simulation{};
};
struct JointAnalysisComponent : JointComponentData
{
    JointCheckStatus runtime_convergence{JointCheckStatus::Unavailable};
};
// A saved outcome, not an audit/restart bundle. Convergence values are captured
// from the runtime result once and are never inferred by storage or readers.
struct JointAnalysisResult
{
    JointAnalysisMetadata metadata;
    std::vector<std::string> atom_ids, row_ids;
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
