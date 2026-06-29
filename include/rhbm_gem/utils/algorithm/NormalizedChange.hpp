#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

namespace rhbm_gem::algorithm {

inline double CalculateNormalizedChange(
    double current,
    double previous,
    double scale_floor)
{
    if (!std::isfinite(scale_floor) || scale_floor <= 0.0)
    {
        throw std::invalid_argument("Normalized change scale floor must be positive and finite.");
    }
    const auto scale{
        std::max({ std::abs(current), std::abs(previous), scale_floor })
    };
    return std::abs(current - previous) / scale;
}

inline double CalculateMaximumNormalizedVectorChange(
    const Eigen::VectorXd & current,
    const Eigen::VectorXd & previous,
    double scale_floor)
{
    if (current.size() != previous.size())
    {
        throw std::invalid_argument("Normalized vector change input sizes are inconsistent.");
    }
    double maximum_change{ 0.0 };
    for (Eigen::Index i = 0; i < current.size(); i++)
    {
        maximum_change = std::max(
            maximum_change,
            CalculateNormalizedChange(current(i), previous(i), scale_floor));
    }
    return maximum_change;
}

} // namespace rhbm_gem::algorithm
