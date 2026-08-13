#pragma once

#include <cstddef>

namespace rhbm_gem::core::detail {

enum class LocalFittingFinalStateSource
{
    BestAudit,
    LatestValidated,
    Unavailable
};

inline LocalFittingFinalStateSource SelectLocalFittingFinalStateSource(
    bool apply_best_iteration,
    bool has_latest_validated_state,
    bool has_best_audit_state)
{
    if (!has_latest_validated_state)
    {
        return LocalFittingFinalStateSource::Unavailable;
    }
    return apply_best_iteration && has_best_audit_state ?
        LocalFittingFinalStateSource::BestAudit :
        LocalFittingFinalStateSource::LatestValidated;
}

inline std::size_t AdvanceLocalFittingAuditPatience(
    std::size_t current_count,
    bool improved_best_audit,
    bool changed_rejected_trust_radius)
{
    if (improved_best_audit || changed_rejected_trust_radius) return 0;
    return current_count + 1;
}

} // namespace rhbm_gem::core::detail
