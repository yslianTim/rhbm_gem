#pragma once

#include "core/detail/second_stage/CandidateState.hpp"

#include <cstddef>
#include <vector>

namespace rhbm_gem::core::detail {

enum class SecondStageStopReason
{
    None,
    Quarantine,
    Converged,
    AuditPatience,
    AllRejectedBacktrackingExhausted,
    AllRejectedAtMaximumIterations,
    MaximumIterations
};

struct IterationResult
{
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    TrustRegionRadiusUpdate trust_region_update{};
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t quarantine_atom_count{ 0 };
    PolishProgress polish_progress{};
    std::size_t suspicious_atom_count{ 0 };
    SecondStageStopReason stop_reason{ SecondStageStopReason::None };
    bool objective_domain_changed{ false };
    TransformedChange transformed_change_percentile{};
};

} // namespace rhbm_gem::core::detail
