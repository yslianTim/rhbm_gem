#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace rhbm_gem::core::detail {

struct LocalFittingObjectiveBreakdown
{
    double residual_objective{ 0.0 };
    double width_prior_penalty{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    double total_objective{ 0.0 };
};

inline std::optional<LocalFittingObjectiveBreakdown>
BuildLocalFittingMeanObjectiveBreakdown(
    double residual_objective,
    double width_prior_penalty_sum,
    double offset_plausibility_penalty_sum,
    std::size_t atom_count,
    double width_prior_penalty_weight,
    double offset_plausibility_penalty_weight)
{
    if (atom_count == 0 ||
        !std::isfinite(residual_objective) ||
        !std::isfinite(width_prior_penalty_sum) ||
        !std::isfinite(offset_plausibility_penalty_sum) ||
        !std::isfinite(width_prior_penalty_weight) ||
        !std::isfinite(offset_plausibility_penalty_weight))
    {
        return std::nullopt;
    }

    const auto atom_count_double{ static_cast<double>(atom_count) };
    LocalFittingObjectiveBreakdown breakdown;
    breakdown.residual_objective = residual_objective;
    breakdown.width_prior_penalty =
        width_prior_penalty_weight * width_prior_penalty_sum / atom_count_double;
    breakdown.offset_plausibility_penalty =
        offset_plausibility_penalty_weight *
        offset_plausibility_penalty_sum / atom_count_double;
    breakdown.total_objective =
        breakdown.residual_objective +
        breakdown.width_prior_penalty +
        breakdown.offset_plausibility_penalty;
    if (!std::isfinite(breakdown.total_objective)) return std::nullopt;
    return breakdown;
}

inline bool IsBetterLocalFittingAuditObjective(
    double candidate,
    double best,
    double relative_tolerance)
{
    if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerance must be finite and non-negative.");
    }
    if (!std::isfinite(candidate)) return false;
    if (!std::isfinite(best)) return true;
    const auto scale{ std::max({ std::abs(candidate), std::abs(best), 1.0 }) };
    return candidate < best - relative_tolerance * scale;
}

inline bool IsLocalFittingAuditObjectiveAcceptableForProgress(
    const std::optional<double> & candidate,
    const std::optional<double> & previous,
    const std::optional<double> & best,
    double relative_tolerance)
{
    if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerance must be finite and non-negative.");
    }
    if (!candidate.has_value() || !previous.has_value() ||
        !std::isfinite(*candidate) || !std::isfinite(*previous))
    {
        return false;
    }
    const auto is_deteriorated = [&](double reference)
    {
        if (!std::isfinite(reference)) return true;
        const auto scale{ std::max(std::abs(reference), 1.0) };
        return *candidate > reference + relative_tolerance * scale;
    };
    return !is_deteriorated(*previous) &&
        (!best.has_value() || !is_deteriorated(*best));
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
