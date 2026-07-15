#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/algorithm/ParameterChange.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::algorithm {

struct ParameterChangeStats
{
    std::vector<double> percentile_list{};
};

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
