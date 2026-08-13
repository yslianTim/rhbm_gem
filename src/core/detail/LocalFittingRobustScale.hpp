#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

constexpr double kLocalFittingRobustScaleMultiplier{ 1.4826 };
constexpr double kLocalFittingRobustScaleMin{ 1.0e-12 };

inline double CalculateLocalFittingMedianAbsoluteDeviationScale(
    const std::vector<double> & value_list)
{
    if (value_list.empty())
    {
        return std::numeric_limits<double>::infinity();
    }
    for (const auto value : value_list)
    {
        if (!std::isfinite(value))
        {
            return std::numeric_limits<double>::infinity();
        }
    }
    const auto median{ array_helper::ComputeMedian(value_list) };
    std::vector<double> deviation_list;
    deviation_list.reserve(value_list.size());
    for (const auto value : value_list)
    {
        deviation_list.emplace_back(std::abs(value - median));
    }
    return kLocalFittingRobustScaleMultiplier *
        array_helper::ComputeMedian(deviation_list);
}

// Kept as a compatibility name while the implementation lives in the
// responsibility-specific robust-scale header.
inline double CalculateMedianAbsoluteDeviationScale(
    const std::vector<double> & value_list)
{
    return CalculateLocalFittingMedianAbsoluteDeviationScale(value_list);
}

} // namespace rhbm_gem::core::detail
