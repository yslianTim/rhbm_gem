#pragma once

#include <cmath>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

constexpr double kRobustScaleMultiplier{ 1.4826 };
constexpr double kRobustScaleMin{ 1.0e-12 };

inline double CalculateMedianAbsoluteDeviationScale(const std::vector<double> & value_list)
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
    return kRobustScaleMultiplier * array_helper::ComputeMedian(deviation_list);
}

} // namespace rhbm_gem::core::detail
