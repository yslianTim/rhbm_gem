#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include "core/detail/FitStateView.hpp"
#include "core/detail/TransformedGaussianModel.hpp"

namespace rhbm_gem::core::detail {

constexpr double kTransformedChangePercentile{ 0.99 };
constexpr double kTransformedChangeTolerance{ 1.0e-4 };
constexpr double kTransformedMaximumChangeTolerance{ 1.0e-3 };

struct TransformedChangeSummary
{
    algorithm::ParameterChangeStats percentile_stats{};
    std::vector<double> maximum_list{};
};

inline algorithm::ParameterChange MakeInfiniteTransformedChange()
{
    return algorithm::ParameterChange{
        std::vector<double>(
            kTransformedChangeSize,
            std::numeric_limits<double>::infinity())
    };
}

inline algorithm::ParameterChange CalculateTransformedChange(
    const GaussianModel3D & current,
    const GaussianModel3D & previous)
{
    const auto current_coordinates{
        EncodeTransformedCoordinates(current)
    };
    const auto previous_coordinates{
        EncodeTransformedCoordinates(previous)
    };
    if (!current_coordinates.has_value() || !previous_coordinates.has_value())
    {
        return MakeInfiniteTransformedChange();
    }

    algorithm::ParameterChange change{
        std::vector<double>(kTransformedChangeSize, 0.0)
    };
    for (std::size_t i = 0; i < kTransformedChangeSize; i++)
    {
        const auto eigen_index{ static_cast<Eigen::Index>(i) };
        const auto value{
            std::abs(
                (*current_coordinates)(eigen_index) -
                (*previous_coordinates)(eigen_index))
        };
        if (!std::isfinite(value))
        {
            return MakeInfiniteTransformedChange();
        }
        change.value_list.at(i) = value;
    }
    return change;
}

inline double GetMaximumTransformedChange(const std::vector<double> & value_list)
{
    if (value_list.empty()) return 0.0;
    if (value_list.size() != kTransformedChangeSize)
    {
        throw std::invalid_argument(
            "Local fitting transformed change input is inconsistent.");
    }
    return std::ranges::max(value_list);
}

inline bool IsTransformedChangeMaterial(
    const algorithm::ParameterChange & change,
    double minimum_change)
{
    if (!std::isfinite(minimum_change) || minimum_change < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting transformed change threshold is invalid.");
    }
    if (change.value_list.size() != kTransformedChangeSize)
    {
        throw std::invalid_argument(
            "Local fitting transformed change input is inconsistent.");
    }
    return std::any_of(
        change.value_list.begin(),
        change.value_list.end(),
        [minimum_change](double value)
        {
            return std::isfinite(value) && value >= minimum_change;
        });
}

inline std::vector<double> SummarizeMaximumTransformedChanges(
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & index_list)
{
    std::vector<double> maximum_list(kTransformedChangeSize, 0.0);
    for (const auto index : index_list)
    {
        if (index >= change_list.size() ||
            change_list.at(index).value_list.size() != kTransformedChangeSize)
        {
            throw std::invalid_argument(
                "Local fitting maximum transformed change input is inconsistent.");
        }
        for (std::size_t parameter_index = 0;
            parameter_index < kTransformedChangeSize;
            parameter_index++)
        {
            maximum_list.at(parameter_index) = std::max(
                maximum_list.at(parameter_index),
                change_list.at(index).value_list.at(parameter_index));
        }
    }
    return maximum_list;
}

inline bool IsTransformedChangeConverged(
    const algorithm::ParameterChangeStats & percentile_stats,
    const std::vector<double> & maximum_list)
{
    if (percentile_stats.percentile_list.size() != kTransformedChangeSize ||
        maximum_list.size() != kTransformedChangeSize)
    {
        throw std::invalid_argument(
            "Local fitting transformed convergence statistics are inconsistent.");
    }
    for (std::size_t i = 0; i < kTransformedChangeSize; i++)
    {
        if (!std::isfinite(percentile_stats.percentile_list.at(i)) ||
            percentile_stats.percentile_list.at(i) >= kTransformedChangeTolerance ||
            !std::isfinite(maximum_list.at(i)) ||
            maximum_list.at(i) >= kTransformedMaximumChangeTolerance)
        {
            return false;
        }
    }
    return true;
}

template <typename CurrentState, typename PreviousState>
inline TransformedChangeSummary SummarizeTransformedChanges(
    const CurrentState & current_state,
    const PreviousState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    std::vector<algorithm::ParameterChange> change_list;
    change_list.reserve(index_list.size());
    for (const auto i : index_list)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, i),
            GetFitModel(previous_state, i)));
    }

    std::vector<std::size_t> local_index_list(change_list.size());
    for (std::size_t i = 0; i < local_index_list.size(); i++)
    {
        local_index_list.at(i) = i;
    }
    return TransformedChangeSummary{
        algorithm::SummarizeParameterChangeStats(
            change_list,
            local_index_list,
            kTransformedChangePercentile),
        SummarizeMaximumTransformedChanges(change_list, local_index_list)
    };
}

inline double GetMaximumTransformedChange(const TransformedChangeSummary & summary)
{
    return GetMaximumTransformedChange(summary.maximum_list);
}

inline bool IsTransformedChangeConverged(const TransformedChangeSummary & summary)
{
    return IsTransformedChangeConverged(summary.percentile_stats, summary.maximum_list);
}

inline bool IsTransformedPercentileConverged(const TransformedChangeSummary & summary)
{
    return std::all_of(
        summary.percentile_stats.percentile_list.begin(),
        summary.percentile_stats.percentile_list.end(),
        [](double value)
        {
            return std::isfinite(value) && value < kTransformedChangeTolerance;
        });
}

} // namespace rhbm_gem::core::detail
