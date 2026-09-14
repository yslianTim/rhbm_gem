#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"

#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

namespace rhbm_gem::core::detail {

struct ObjectiveProgressGateEvidence
{
    bool previous_checked{ false }, best_checked{ false };
    std::string_view reason{ "objective-unavailable" };
};

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

enum class BoundaryComponentAcceptedSource
{
    None,
    Endpoint,
    JointCorrection,
    Backtracking
};

} // namespace rhbm_gem::core::detail
