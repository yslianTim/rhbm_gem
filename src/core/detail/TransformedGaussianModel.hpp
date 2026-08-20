#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr std::size_t kLogPeakHeightChangeIndex{ 0 };
constexpr std::size_t kLogWidthChangeIndex{ 1 };
constexpr std::size_t kOffsetToPeakRatioChangeIndex{ 2 };
constexpr std::size_t kTransformedChangeSize{ 3 };

inline std::optional<Eigen::Vector3d> EncodeTransformedCoordinates(const GaussianModel3D & model)
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
        std::log(amplitude) - 1.5 * std::log(Constants::two_pi) - 3.0 * log_width
    };
    double offset_to_peak_ratio{ 0.0 };
    if (offset != 0.0)
    {
        const auto log_abs_offset_to_peak_ratio{
            std::log(std::abs(offset)) +
            0.5 * std::log(4.0 / Constants::two_pi) - log_width - log_peak_height
        };
        if (log_abs_offset_to_peak_ratio > std::log(std::numeric_limits<double>::max()))
        {
            return std::nullopt;
        }
        offset_to_peak_ratio = std::copysign(std::exp(log_abs_offset_to_peak_ratio), offset);
    }

    if (!std::isfinite(log_peak_height) ||
        !std::isfinite(log_width) ||
        !std::isfinite(offset_to_peak_ratio))
    {
        return std::nullopt;
    }

    return Eigen::Vector3d{ log_peak_height, log_width, offset_to_peak_ratio };
}

inline std::optional<GaussianModel3D> DecodeTransformedCoordinates(const Eigen::Vector3d & coordinates)
{
    if (!coordinates.allFinite()) return std::nullopt;

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
            log_peak_height + log_width - 0.5 * std::log(4.0 / Constants::two_pi)
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

inline bool IsValidSecondStageGaussianModel(const GaussianModel3D & model)
{
    return EncodeTransformedCoordinates(model).has_value();
}

inline std::optional<GaussianModel3D> BuildGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list)
{
    std::vector<double> amplitude_list;
    std::vector<double> width_list;
    std::vector<double> offset_list;
    amplitude_list.reserve(model_list.size());
    width_list.reserve(model_list.size());
    offset_list.reserve(model_list.size());
    for (const auto & model : model_list)
    {
        if (!IsValidSecondStageGaussianModel(model)) continue;
        amplitude_list.emplace_back(model.GetAmplitude());
        width_list.emplace_back(model.GetWidth());
        offset_list.emplace_back(model.GetOffset());
    }
    if (amplitude_list.empty()) return std::nullopt;

    const GaussianModel3D median_model{
        array_helper::ComputeMedian(amplitude_list),
        array_helper::ComputeMedian(width_list),
        array_helper::ComputeMedian(offset_list)
    };
    return IsValidSecondStageGaussianModel(median_model) ?
        std::optional<GaussianModel3D>{ median_model } : std::nullopt;
}

inline std::vector<GaussianModel3D> BuildGroupMedianModelList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    if (group_id_by_atom_position.size() != model_list.size())
    {
        throw std::invalid_argument("Local fitting group median inputs are inconsistent.");
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> atom_position_list_by_group;
    atom_position_list_by_group.reserve(model_list.size());
    for (std::size_t atom_position = 0; atom_position < model_list.size(); atom_position++)
    {
        atom_position_list_by_group[group_id_by_atom_position.at(atom_position)]
            .emplace_back(atom_position);
    }

    auto group_median_model_list{ model_list };
    std::vector<GaussianModel3D> group_model_list;
    group_model_list.reserve(model_list.size());
    for (const auto & entry : atom_position_list_by_group)
    {
        const auto & atom_position_list{ entry.second };
        group_model_list.clear();
        for (const auto atom_position : atom_position_list)
        {
            group_model_list.emplace_back(model_list.at(atom_position));
        }
        const auto median_model{ BuildGaussianParameterMedian(group_model_list) };
        if (!median_model.has_value()) continue;
        for (const auto atom_position : atom_position_list)
        {
            group_median_model_list.at(atom_position) = *median_model;
        }
    }
    return group_median_model_list;
}

inline std::vector<double> BuildGroupMedianOffsetList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    if (group_id_by_atom_position.size() != model_list.size())
    {
        throw std::invalid_argument("Local fitting group median inputs are inconsistent.");
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> atom_position_list_by_group;
    atom_position_list_by_group.reserve(model_list.size());
    std::vector<double> offset_list;
    offset_list.reserve(model_list.size());
    for (std::size_t atom_position = 0;
        atom_position < model_list.size();
        atom_position++)
    {
        atom_position_list_by_group[group_id_by_atom_position.at(atom_position)]
            .emplace_back(atom_position);
        offset_list.emplace_back(model_list.at(atom_position).GetOffset());
    }

    std::vector<double> valid_offset_list;
    valid_offset_list.reserve(model_list.size());
    for (const auto & entry : atom_position_list_by_group)
    {
        valid_offset_list.clear();
        for (const auto atom_position : entry.second)
        {
            const auto & model{ model_list.at(atom_position) };
            if (IsValidSecondStageGaussianModel(model))
            {
                valid_offset_list.emplace_back(model.GetOffset());
            }
        }
        if (valid_offset_list.empty()) continue;
        const auto median{ array_helper::ComputeMedian(valid_offset_list) };
        if (!std::isfinite(median)) continue;
        for (const auto atom_position : entry.second)
        {
            offset_list.at(atom_position) = median;
        }
    }
    return offset_list;
}

inline bool TryBuildSharedOffsetDampedModelList(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & raw_model_list,
    const std::vector<double> & previous_shared_offset_list,
    const std::vector<double> & raw_shared_offset_list,
    double damping,
    std::vector<GaussianModel3D> & candidate_model_list)
{
    candidate_model_list.clear();
    if (raw_model_list.size() != previous_model_list.size() ||
        previous_shared_offset_list.size() != previous_model_list.size() ||
        raw_shared_offset_list.size() != previous_model_list.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        return false;
    }

    candidate_model_list.reserve(previous_model_list.size());
    for (std::size_t atom_position = 0; atom_position < previous_model_list.size(); atom_position++)
    {
        const auto previous_coordinates{
            EncodeTransformedCoordinates(previous_model_list.at(atom_position))
        };
        const auto raw_coordinates{
            EncodeTransformedCoordinates(raw_model_list.at(atom_position))
        };
        if (!previous_coordinates.has_value() || !raw_coordinates.has_value())
        {
            candidate_model_list.clear();
            return false;
        }

        Eigen::Vector3d shape_coordinates{
            (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) +
                damping * (
                    (*raw_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) -
                    (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex))),
            (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) +
                damping * (
                    (*raw_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) -
                    (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex))),
            0.0
        };
        const auto shape_model{ DecodeTransformedCoordinates(shape_coordinates) };
        if (!shape_model.has_value())
        {
            candidate_model_list.clear();
            return false;
        }

        const auto previous_shared_offset{ previous_shared_offset_list.at(atom_position) };
        const auto raw_shared_offset{ raw_shared_offset_list.at(atom_position) };
        const auto candidate_model{
            shape_model->WithOffset(previous_shared_offset + damping * (raw_shared_offset - previous_shared_offset))
        };
        if (!IsValidSecondStageGaussianModel(candidate_model))
        {
            candidate_model_list.clear();
            return false;
        }
        candidate_model_list.emplace_back(candidate_model);
    }
    return true;
}

struct SharedOffsetResponse
{
    double response{ 0.0 };
    Eigen::Vector2d shape_jacobian{ Eigen::Vector2d::Zero() };
    double offset_jacobian{ 0.0 };
};

inline std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance)
{
    if (!EncodeTransformedCoordinates(model).has_value() || !std::isfinite(distance) || distance < 0.0)
    {
        return std::nullopt;
    }

    const auto width{ model.GetWidth() };
    const auto signal{ model.SignalAtDistance(distance) };
    const auto offset_basis{ model.OffsetBasisAtDistance(distance) };
    const auto response{ signal + model.GetOffset() * offset_basis };
    if (!std::isfinite(signal) || !std::isfinite(offset_basis) || !std::isfinite(response))
    {
        return std::nullopt;
    }

    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{ signal * normalized_distance * normalized_distance };
    if (distance < 1.0e-5)
    {
        log_width_derivative -= model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{ -0.5 * normalized_distance * normalized_distance };
        log_width_derivative -= model.GetOffset() * center_offset_basis_scale / width * std::exp(exponent);
    }

    Eigen::Vector2d shape_jacobian{ signal, log_width_derivative };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        response,
        shape_jacobian,
        offset_basis
    };
}

struct TransformedModelInvariants
{
    GaussianModel3D model{};
    double peak_height{ 0.0 };
};

inline std::optional<TransformedModelInvariants> BuildTransformedModelInvariants(const GaussianModel3D & model)
{
    const auto transformed{ EncodeTransformedCoordinates(model) };
    if (!transformed.has_value()) return std::nullopt;

    const auto peak_height{
        std::exp((*transformed)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)))
    };
    if (!std::isfinite(peak_height)) return std::nullopt;

    return TransformedModelInvariants{ model, peak_height };
}

inline std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(
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
    auto log_width_derivative{ evaluation.signal * normalized_distance * normalized_distance };
    if (distance < 1.0e-5)
    {
        log_width_derivative -= invariants.model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{ -0.5 * normalized_distance * normalized_distance };
        log_width_derivative -= invariants.model.GetOffset() * center_offset_basis_scale / width * std::exp(exponent);
    }

    const Eigen::Vector2d shape_jacobian{ evaluation.signal, log_width_derivative };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        evaluation.response,
        shape_jacobian,
        evaluation.offset_basis
    };
}

inline std::optional<Eigen::Vector3d> EvaluateTransformedJacobian(
    const TransformedModelInvariants & invariants,
    double distance)
{
    const auto shared_offset_evaluation{ EvaluateSharedOffsetResponse(invariants, distance) };
    if (!shared_offset_evaluation.has_value()) return std::nullopt;

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
        shared_offset_evaluation->offset_jacobian / center_offset_basis_scale;
    if (!jacobian.allFinite()) return std::nullopt;

    return jacobian;
}

} // namespace rhbm_gem::core::detail
