#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rhbm_gem::core::detail {

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

} // namespace rhbm_gem::core::detail
