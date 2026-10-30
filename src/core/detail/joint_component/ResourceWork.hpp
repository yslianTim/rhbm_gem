#pragma once
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <cstddef>
#include <optional>
#include <string>

namespace rhbm_gem::core::joint_component {
inline constexpr std::size_t JointMemoryHardLimitBytes=4'000'000'000;
inline constexpr std::size_t JointMemoryPlanningBudgetBytes=3'200'000'000;
enum class MemoryStage {ProblemConstruction,StructuralConstruction,FixedNeighborPreparation,GlobalStateEvaluation,
    EndpointVerification,LocalConditioning};
enum class MemoryStatus {Pass,Reject,Overflow};
struct MemoryShape
{
    std::size_t component_atoms{},component_observations{},component_support_entries{},component_blocks{},
        affected_row_entries{},component_neighbor_entries{},maximum_block_rows{},maximum_block_support_entries{},
        total_block_support_entries{},maximum_neighbor_atoms{},maximum_block_replay_entries{},maximum_block_atoms{},
        component_free_columns{},component_width_columns{},
        identity_validation_bytes{},resident_bytes{},
        maximum_tile_rows{8192};
    bool problem_scope{},affected_replay_scratch{};
};
struct MemoryEstimate
{
    MemoryShape shape;
    MemoryStage stage{MemoryStage::ProblemConstruction};
    std::size_t resident_bytes{},temporary_bytes{},estimated_peak_bytes{},budget_bytes{JointMemoryPlanningBudgetBytes};
    MemoryStatus status{MemoryStatus::Pass};
    std::string reason;
    std::size_t sparse_design_bytes{},dense_factor_bytes{},factor_workspace_bytes{},factor_rows{},factor_columns{};
};
bool CheckedAddBytes(std::size_t,std::size_t,std::size_t &);
bool CheckedMultiplyBytes(std::size_t,std::size_t,std::size_t &);
std::optional<std::size_t> EstimateSnapshotBytes(const JointProblemInput &);
std::optional<std::size_t> EstimateConstructionScratchBytes(const JointProblemInput &);
std::optional<std::size_t> CountSupportEntries(const JointProblemInput &);
MemoryEstimate PlanMemory(const MemoryShape &,MemoryStage,std::size_t=JointMemoryPlanningBudgetBytes);
const char * MemoryStageName(MemoryStage);
const char * MemoryStatusName(MemoryStatus);
std::optional<std::size_t> CurrentRssBytes();
std::string FormatMemoryEstimate(const MemoryEstimate &);
void LogMemoryEstimate(const MemoryEstimate &);
}
