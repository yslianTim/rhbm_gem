#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace rhbm_gem::core::detail {

enum class PreObjectiveFailureReason
{
    None,
    InvalidModel,
    NoCandidateWithinTrustRegion
};

struct CandidateDecisionEvidence
{
    std::optional<double> accepted_factor{};
    PreObjectiveFailureReason pre_objective_failure_reason{ PreObjectiveFailureReason::None };
    std::optional<ObjectiveBreakdown> candidate_objective{};
    std::optional<ObjectiveBreakdown> previous_objective{};
    bool rejected_by_previous{ false };
    std::size_t invalid_trial_count{ 0 };
    std::size_t guard_rejected_trial_count{ 0 };
    std::size_t objective_rejected_trial_count{ 0 };
    std::vector<StabilizationTerminalEvidence> terminal_evidence_list{};
};

struct TrustRegionRadiusUpdate
{
    std::vector<ClusterKey> changed_key_list{};
    std::vector<ClusterKey> saturated_key_list{};
};

struct ClusterCandidateDecision
{
    ClusterKey key{};
    CandidateDecisionEvidence evidence{};
};

struct PolishProgress
{
    std::size_t eligible_count{ 0 };
    std::size_t accepted_count{ 0 };
    std::size_t rejected_count{ 0 };
    std::size_t skipped_count{ 0 };
};

enum class BoundaryComponentAcceptedSource
{
    None,
    Endpoint,
    JointCorrection,
    Backtracking
};

struct BoundaryComponentDecision
{
    std::vector<ClusterKey> key_list{};
    std::optional<double> accepted_factor{};
    BoundaryComponentAcceptedSource accepted_source{ BoundaryComponentAcceptedSource::None };
    bool exhausted{ false };
};

struct CandidateSelection
{
    SuspiciousBlockActivity block_activity{};
    FitState assembled_state{};
    PolishProvenance assembled_polish_provenance{};
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    std::vector<ClusterKey> shrink_trust_region_key_list{};
    std::vector<ClusterKey> exhausted_key_list{};
    std::vector<ClusterCandidateDecision> accepted_cluster_evidence_list{};
    std::vector<ClusterCandidateDecision> rejected_cluster_evidence_list{};
    std::vector<BoundaryComponentDecision> boundary_decision_list{};
    std::optional<ObjectiveBreakdown> final_audit_objective{};
    PolishProgress polish_progress{};
};

} // namespace rhbm_gem::core::detail
