#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/algorithm/ParameterChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr std::size_t kLogPeakHeightChangeIndex{ 0 };
constexpr std::size_t kLogWidthChangeIndex{ 1 };
constexpr std::size_t kOffsetToPeakRatioChangeIndex{ 2 };
constexpr std::size_t kTransformedChangeSize{ 3 };

inline algorithm::ParameterChange MakeInfiniteLocalFittingTransformedChange()
{
    return algorithm::ParameterChange{
        std::vector<double>(
            kTransformedChangeSize,
            std::numeric_limits<double>::infinity())
    };
}

inline std::optional<std::array<double, kTransformedChangeSize>>
BuildLocalFittingTransformedCoordinates(const GaussianModel3D & model)
{
    const auto amplitude{ model.GetAmplitude() };
    const auto width{ model.GetWidth() };
    const auto offset{ model.GetOffset() };
    if (!std::isfinite(amplitude) || amplitude <= 0.0 ||
        !std::isfinite(width) || width <= 0.0 ||
        !std::isfinite(offset))
    {
        return std::nullopt;
    }

    const auto log_width{ std::log(width) };
    const auto log_peak_height{
        std::log(amplitude) -
        1.5 * std::log(Constants::two_pi) -
        3.0 * log_width
    };
    double offset_to_peak_ratio{ 0.0 };
    if (offset != 0.0)
    {
        const auto log_abs_offset_to_peak_ratio{
            std::log(std::abs(offset)) +
            0.5 * std::log(4.0 / Constants::two_pi) -
            log_width -
            log_peak_height
        };
        if (log_abs_offset_to_peak_ratio >
            std::log(std::numeric_limits<double>::max()))
        {
            return std::nullopt;
        }
        offset_to_peak_ratio = std::copysign(
            std::exp(log_abs_offset_to_peak_ratio),
            offset);
    }

    if (!std::isfinite(log_peak_height) ||
        !std::isfinite(log_width) ||
        !std::isfinite(offset_to_peak_ratio))
    {
        return std::nullopt;
    }
    return std::array<double, kTransformedChangeSize>{
        log_peak_height,
        log_width,
        offset_to_peak_ratio
    };
}

inline std::optional<Eigen::VectorXd> EncodeLocalFittingTransformedCoordinates(
    const GaussianModel3D & model)
{
    const auto coordinates{ BuildLocalFittingTransformedCoordinates(model) };
    if (!coordinates.has_value()) return std::nullopt;

    Eigen::VectorXd encoded{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(kTransformedChangeSize))
    };
    for (std::size_t i = 0; i < kTransformedChangeSize; i++)
    {
        encoded(static_cast<Eigen::Index>(i)) = coordinates->at(i);
    }
    return encoded;
}

inline std::optional<GaussianModel3D> DecodeLocalFittingTransformedCoordinates(
    const Eigen::VectorXd & coordinates)
{
    if (coordinates.size() != static_cast<Eigen::Index>(kTransformedChangeSize) ||
        !coordinates.allFinite())
    {
        return std::nullopt;
    }

    const auto log_peak_height{
        coordinates(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex))
    };
    const auto log_width{
        coordinates(static_cast<Eigen::Index>(kLogWidthChangeIndex))
    };
    const auto offset_to_peak_ratio{
        coordinates(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex))
    };
    const auto log_amplitude{
        log_peak_height + 1.5 * std::log(Constants::two_pi) + 3.0 * log_width
    };
    const auto amplitude{ std::exp(log_amplitude) };
    const auto width{ std::exp(log_width) };
    double offset{ 0.0 };
    if (offset_to_peak_ratio != 0.0)
    {
        const auto log_abs_offset{
            std::log(std::abs(offset_to_peak_ratio)) +
            log_peak_height + log_width -
            0.5 * std::log(4.0 / Constants::two_pi)
        };
        offset = std::copysign(std::exp(log_abs_offset), offset_to_peak_ratio);
    }

    const GaussianModel3D model{ amplitude, width, offset };
    if (!std::isfinite(amplitude) || amplitude <= 0.0 ||
        !std::isfinite(width) || width <= 0.0 ||
        !std::isfinite(offset))
    {
        return std::nullopt;
    }
    return model;
}

inline algorithm::ParameterChange CalculateLocalFittingTransformedChange(
    const GaussianModel3D & current,
    const GaussianModel3D & previous)
{
    const auto current_coordinates{ BuildLocalFittingTransformedCoordinates(current) };
    const auto previous_coordinates{ BuildLocalFittingTransformedCoordinates(previous) };
    if (!current_coordinates.has_value() || !previous_coordinates.has_value())
    {
        return MakeInfiniteLocalFittingTransformedChange();
    }

    algorithm::ParameterChange change{
        std::vector<double>(kTransformedChangeSize, 0.0)
    };
    for (std::size_t i = 0; i < kTransformedChangeSize; i++)
    {
        const auto value{
            std::abs(current_coordinates->at(i) - previous_coordinates->at(i))
        };
        if (!std::isfinite(value))
        {
            return MakeInfiniteLocalFittingTransformedChange();
        }
        change.value_list.at(i) = value;
    }
    return change;
}

inline std::vector<double> SummarizeLocalFittingMaximumTransformedChanges(
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

inline bool IsLocalFittingTransformedChangeConverged(
    const algorithm::ParameterChangeStats & percentile_stats,
    const std::vector<double> & maximum_list,
    double percentile_tolerance,
    double maximum_tolerance)
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
            percentile_stats.percentile_list.at(i) >= percentile_tolerance ||
            !std::isfinite(maximum_list.at(i)) ||
            maximum_list.at(i) >= maximum_tolerance)
        {
            return false;
        }
    }
    return true;
}

} // namespace rhbm_gem::core::detail
