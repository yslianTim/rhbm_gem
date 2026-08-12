#pragma once

#include <cstddef>
#include <optional>

#include "core/detail/LocalFittingAudit.hpp"

namespace rhbm_gem::core::detail {

enum class LocalFittingPreObjectiveFailureReason
{
    None,
    InvalidModel,
    PreviousSharedOffsetProjectionOutsideTrustRegion,
    NoCandidateWithinTrustRegion
};

struct LocalFittingObjectiveAttemptDiagnostic
{
    double effective_damping{ 1.0 };
    bool is_invalid_model{ false };
    LocalFittingPreObjectiveFailureReason pre_objective_failure_reason{
        LocalFittingPreObjectiveFailureReason::None
    };
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<double> fit_scale{};
    std::optional<double> tail_scale{};
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
    std::optional<LocalFittingObjectiveBreakdown> candidate_objective{};
    std::optional<LocalFittingObjectiveBreakdown> previous_objective{};
    std::optional<LocalFittingObjectiveBreakdown> best_objective{};
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    std::size_t backtracking_trial_count{ 0 };
    std::optional<double> accepted_backtracking_factor{};
    bool backtracking_exhausted{ false };
};

} // namespace rhbm_gem::core::detail
