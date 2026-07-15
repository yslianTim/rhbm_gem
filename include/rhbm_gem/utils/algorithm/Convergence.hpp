#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::algorithm {

struct ParameterChange
{
    std::vector<double> value_list{};
};

struct ParameterChangeStats
{
    std::vector<double> percentile_list{};
};

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

inline ParameterChangeStats SummarizeParameterChangeStats(
    const std::vector<ParameterChange> & change_list,
    const std::vector<std::size_t> & index_list,
    double percentile)
{
    std::size_t parameter_size{ 0 };
    for (const auto index : index_list)
    {
        if (index >= change_list.size())
        {
            throw std::invalid_argument("Parameter change index is out of range.");
        }
        parameter_size = std::max(parameter_size, change_list.at(index).value_list.size());
    }

    ParameterChangeStats stats;
    stats.percentile_list.resize(parameter_size, 0.0);
    for (std::size_t parameter_index = 0; parameter_index < parameter_size; parameter_index++)
    {
        std::vector<double> parameter_change_list;
        parameter_change_list.reserve(index_list.size());
        for (const auto index : index_list)
        {
            const auto & values{ change_list.at(index).value_list };
            if (parameter_index >= values.size())
            {
                throw std::invalid_argument("Parameter change sizes are inconsistent.");
            }
            parameter_change_list.emplace_back(values.at(parameter_index));
        }
        stats.percentile_list.at(parameter_index) = array_helper::ComputePercentile(
            parameter_change_list,
            percentile);
    }
    return stats;
}

inline double GetMaximumParameterChange(const ParameterChangeStats & stats)
{
    if (stats.percentile_list.empty()) return 0.0;
    return *std::max_element(stats.percentile_list.begin(), stats.percentile_list.end());
}

} // namespace rhbm_gem::algorithm
