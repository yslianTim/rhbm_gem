#pragma once

#include "core/detail/LocalFittingTransformedChange.hpp"

#include <cmath>
#include <optional>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

struct LocalFittingTransformedResponse
{
    double response{ 0.0 };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
};

inline std::optional<LocalFittingTransformedResponse>
EvaluateLocalFittingTransformedResponse(
    const GaussianModel3D & model,
    double distance)
{
    const auto transformed{ EncodeLocalFittingTransformedCoordinates(model) };
    if (!transformed.has_value() || !std::isfinite(distance) || distance < 0.0)
    {
        return std::nullopt;
    }

    const auto width{ model.GetWidth() };
    const auto peak_height{ std::exp((*transformed)(kLogPeakHeightChangeIndex)) };
    const auto signal{ model.SignalAtDistance(distance) };
    const auto offset_basis{ model.OffsetBasisAtDistance(distance) };
    const auto offset_response{ model.GetOffset() * offset_basis };
    const auto response{ signal + offset_response };
    if (!std::isfinite(peak_height) || !std::isfinite(signal) ||
        !std::isfinite(offset_basis) || !std::isfinite(response))
    {
        return std::nullopt;
    }

    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
    jacobian(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) = response;
    jacobian(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) =
        peak_height * width * offset_basis / center_offset_basis_scale;

    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        signal * normalized_distance * normalized_distance
    };
    if (distance >= 1.0e-5)
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative += model.GetOffset() * (
            offset_basis - center_offset_basis_scale / width * std::exp(exponent));
    }
    jacobian(static_cast<Eigen::Index>(kLogWidthChangeIndex)) =
        log_width_derivative;
    if (!jacobian.allFinite()) return std::nullopt;

    return LocalFittingTransformedResponse{ response, jacobian };
}

} // namespace rhbm_gem::core::detail
