#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

#include <rhbm_gem/utils/algorithm/ParameterChange.hpp>
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

} // namespace rhbm_gem::core::detail
