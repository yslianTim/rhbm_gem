#pragma once

#include <cmath>
#include <optional>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/TransformedGaussianModel.hpp"

namespace rhbm_gem::core::detail {

struct SharedOffsetResponse
{
    double response{ 0.0 };
    Eigen::Vector2d shape_jacobian{ Eigen::Vector2d::Zero() };
    double offset_jacobian{ 0.0 };
};

inline std::optional<SharedOffsetResponse>
EvaluateSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance)
{
    if (!EncodeTransformedCoordinates(model).has_value() ||
        !std::isfinite(distance) || distance < 0.0)
    {
        return std::nullopt;
    }

    const auto width{ model.GetWidth() };
    const auto signal{ model.SignalAtDistance(distance) };
    const auto offset_basis{ model.OffsetBasisAtDistance(distance) };
    const auto response{ signal + model.GetOffset() * offset_basis };
    if (!std::isfinite(signal) || !std::isfinite(offset_basis) ||
        !std::isfinite(response))
    {
        return std::nullopt;
    }

    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -=
            model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative -= model.GetOffset() *
            center_offset_basis_scale / width * std::exp(exponent);
    }

    Eigen::Vector2d shape_jacobian{
        signal,
        log_width_derivative
    };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        response,
        shape_jacobian,
        offset_basis
    };
}

struct TransformedResponse
{
    double response{ 0.0 };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
};

struct TransformedModelInvariants
{
    GaussianModel3D model{};
    Eigen::Vector3d transformed{ Eigen::Vector3d::Zero() };
    double peak_height{ 0.0 };
};

inline std::optional<TransformedModelInvariants>
BuildTransformedModelInvariants(
    const GaussianModel3D & model)
{
    const auto transformed{ EncodeTransformedCoordinates(model) };
    if (!transformed.has_value()) return std::nullopt;

    const auto peak_height{
        std::exp((*transformed)(static_cast<Eigen::Index>(
            kLogPeakHeightChangeIndex)))
    };
    if (!std::isfinite(peak_height)) return std::nullopt;

    return TransformedModelInvariants{
        model,
        *transformed,
        peak_height
    };
}

inline std::optional<SharedOffsetResponse>
EvaluateSharedOffsetResponse(
    const TransformedModelInvariants & invariants,
    double distance)
{
    if (!std::isfinite(distance) || distance < 0.0) return std::nullopt;

    const auto evaluation{ invariants.model.EvaluateAtDistance(distance) };
    if (!std::isfinite(evaluation.signal) ||
        !std::isfinite(evaluation.offset_basis) ||
        !std::isfinite(evaluation.response))
    {
        return std::nullopt;
    }

    const auto width{ invariants.model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        evaluation.signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -=
            invariants.model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative -= invariants.model.GetOffset() *
            center_offset_basis_scale / width * std::exp(exponent);
    }

    const Eigen::Vector2d shape_jacobian{
        evaluation.signal,
        log_width_derivative
    };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        evaluation.response,
        shape_jacobian,
        evaluation.offset_basis
    };
}

inline std::optional<TransformedResponse>
EvaluateTransformedResponse(
    const TransformedModelInvariants & invariants,
    double distance)
{
    const auto shared_offset_evaluation{
        EvaluateSharedOffsetResponse(invariants, distance)
    };
    if (!shared_offset_evaluation.has_value())
    {
        return std::nullopt;
    }

    const auto & model{ invariants.model };
    const auto width{ model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
    jacobian(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(0) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kLogWidthChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(1) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) =
        invariants.peak_height * width *
        shared_offset_evaluation->offset_jacobian /
        center_offset_basis_scale;
    if (!jacobian.allFinite())
    {
        return std::nullopt;
    }

    return TransformedResponse{
        shared_offset_evaluation->response,
        jacobian
    };
}

inline std::optional<TransformedResponse>
EvaluateTransformedResponse(
    const GaussianModel3D & model,
    double distance)
{
    const auto invariants{ BuildTransformedModelInvariants(model) };
    if (!invariants.has_value()) return std::nullopt;
    return EvaluateTransformedResponse(*invariants, distance);
}

} // namespace rhbm_gem::core::detail
