#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr std::size_t kLogPeakHeightChangeIndex{ 0 };
constexpr std::size_t kLogWidthChangeIndex{ 1 };
constexpr std::size_t kOffsetToPeakRatioChangeIndex{ 2 };
constexpr std::size_t kTransformedChangeSize{ 3 };

inline std::optional<Eigen::Vector3d> EncodeTransformedCoordinates(
    const GaussianModel3D & model)
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

    return Eigen::Vector3d{
        log_peak_height,
        log_width,
        offset_to_peak_ratio
    };
}

inline std::optional<GaussianModel3D> DecodeTransformedCoordinates(
    const Eigen::Vector3d & coordinates)
{
    if (!coordinates.allFinite())
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
        offset = std::copysign(
            std::exp(log_abs_offset),
            offset_to_peak_ratio);
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

inline bool IsValidSecondStageGaussianModel(const GaussianModel3D & model)
{
    return std::isfinite(model.GetAmplitude()) && model.GetAmplitude() > 0.0 &&
        std::isfinite(model.GetWidth()) && model.GetWidth() > 0.0 &&
        std::isfinite(model.GetOffset()) &&
        EncodeTransformedCoordinates(model).has_value();
}

} // namespace rhbm_gem::core::detail
