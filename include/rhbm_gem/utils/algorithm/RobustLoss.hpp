#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>

namespace rhbm_gem::algorithm {

namespace detail {

inline void ValidateRobustLossCutoff(double cutoff)
{
    if (!std::isfinite(cutoff) || cutoff <= 0.0)
    {
        throw std::invalid_argument("Robust loss cutoff must be finite and positive.");
    }
}

inline double CalculateLogOnePlusSquare(double value)
{
    const auto absolute_value{ std::abs(value) };
    if (!std::isfinite(absolute_value))
    {
        return std::numeric_limits<double>::infinity();
    }
    if (absolute_value > std::sqrt(std::numeric_limits<double>::max()))
    {
        return 2.0 * std::log(absolute_value);
    }
    return std::log1p(absolute_value * absolute_value);
}

} // namespace detail

inline double CalculateCauchyLoss(
    double residual,
    double cutoff)
{
    detail::ValidateRobustLossCutoff(cutoff);
    if (!std::isfinite(residual))
    {
        return std::numeric_limits<double>::infinity();
    }

    const auto normalized_residual{ residual / cutoff };
    return 0.5 * cutoff * cutoff *
        detail::CalculateLogOnePlusSquare(normalized_residual);
}

inline double CalculateCauchyWeight(
    double residual,
    double residual_scale,
    double cutoff_multiplier)
{
    detail::ValidateRobustLossCutoff(residual_scale);
    detail::ValidateRobustLossCutoff(cutoff_multiplier);
    if (!std::isfinite(residual))
    {
        return 0.0;
    }

    const auto cutoff{ cutoff_multiplier * residual_scale };
    detail::ValidateRobustLossCutoff(cutoff);
    const auto absolute_residual{ std::abs(residual) };
    const auto normalized_residual{ absolute_residual / cutoff };
    if (!std::isfinite(normalized_residual))
    {
        return 0.0;
    }
    const auto denominator{ 1.0 + normalized_residual * normalized_residual };
    return std::isfinite(denominator) ? 1.0 / denominator : 0.0;
}

} // namespace rhbm_gem::algorithm
