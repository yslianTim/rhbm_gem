#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace rhbm_gem::core::detail {

struct LocalFittingObjectiveBreakdown
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double tail_validation_penalty{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    double total_objective{ 0.0 };
};

struct LocalFittingObjectiveTolerance
{
    double absolute_tolerance{ 0.0 };
    double relative_tolerance{ 0.0 };
};

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

inline double CalculateLocalFittingClusterAtomWeight(
    std::size_t cluster_atom_count,
    std::size_t active_atom_count)
{
    if (cluster_atom_count == 0 || active_atom_count == 0 ||
        cluster_atom_count > active_atom_count)
    {
        throw std::invalid_argument(
            "Local fitting cluster atom counts are invalid.");
    }
    return static_cast<double>(cluster_atom_count) /
        static_cast<double>(active_atom_count);
}

inline void ValidateLocalFittingObjectiveTolerance(
    const LocalFittingObjectiveTolerance & tolerance)
{
    if (!std::isfinite(tolerance.absolute_tolerance) ||
        tolerance.absolute_tolerance < 0.0 ||
        !std::isfinite(tolerance.relative_tolerance) ||
        tolerance.relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerances must be finite and "
            "non-negative.");
    }
}

inline double CalculateLocalFittingObjectiveTolerance(
    double reference,
    const LocalFittingObjectiveTolerance & tolerance)
{
    ValidateLocalFittingObjectiveTolerance(tolerance);
    if (!std::isfinite(reference))
    {
        throw std::invalid_argument(
            "Local fitting audit objective reference must be finite.");
    }
    return tolerance.absolute_tolerance +
        tolerance.relative_tolerance * std::abs(reference);
}

inline std::optional<LocalFittingObjectiveBreakdown>
BuildLocalFittingObjectiveBreakdown(
    double fit_range_residual_objective,
    double tail_validation_loss,
    double offset_plausibility_penalty,
    double tail_validation_weight)
{
    if (!std::isfinite(fit_range_residual_objective) ||
        !std::isfinite(tail_validation_loss) ||
        !std::isfinite(offset_plausibility_penalty) ||
        !std::isfinite(tail_validation_weight))
    {
        return std::nullopt;
    }

    LocalFittingObjectiveBreakdown breakdown;
    breakdown.fit_range_residual_objective = fit_range_residual_objective;
    breakdown.tail_validation_loss = tail_validation_loss;
    breakdown.tail_validation_penalty =
        tail_validation_weight * tail_validation_loss;
    breakdown.offset_plausibility_penalty = offset_plausibility_penalty;
    breakdown.total_objective =
        breakdown.fit_range_residual_objective +
        breakdown.tail_validation_penalty +
        breakdown.offset_plausibility_penalty;
    if (!std::isfinite(breakdown.total_objective)) return std::nullopt;
    return breakdown;
}

inline bool IsBetterLocalFittingAuditObjective(
    double candidate,
    double best,
    const LocalFittingObjectiveTolerance & tolerance)
{
    ValidateLocalFittingObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate)) return false;
    if (!std::isfinite(best)) return true;
    return candidate < best -
        CalculateLocalFittingObjectiveTolerance(best, tolerance);
}

inline bool IsLocalFittingAuditObjectiveAcceptableForProgress(
    const std::optional<double> & candidate,
    const std::optional<double> & previous,
    const std::optional<double> & best,
    const LocalFittingObjectiveTolerance & tolerance)
{
    ValidateLocalFittingObjectiveTolerance(tolerance);
    if (!candidate.has_value() || !previous.has_value() ||
        !std::isfinite(*candidate) || !std::isfinite(*previous))
    {
        return false;
    }
    const auto is_deteriorated = [&](double reference)
    {
        if (!std::isfinite(reference)) return true;
        return *candidate > reference +
            CalculateLocalFittingObjectiveTolerance(reference, tolerance);
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
