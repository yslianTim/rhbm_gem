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
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
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
    return std::isfinite(model.GetAmplitude()) && model.GetAmplitude() > 0.0 &&
        std::isfinite(model.GetWidth()) && model.GetWidth() > 0.0 &&
        std::isfinite(model.GetOffset()) &&
        EncodeTransformedCoordinates(model).has_value();
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
    const std::vector<GroupKey> & group_key_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    if (group_key_by_atom_position.size() != model_list.size())
    {
        throw std::invalid_argument("Local fitting group median inputs are inconsistent.");
    }

    std::unordered_map<GroupKey, std::vector<GaussianModel3D>> model_list_by_group;
    model_list_by_group.reserve(model_list.size());
    for (std::size_t atom_position = 0; atom_position < model_list.size(); atom_position++)
    {
        model_list_by_group[group_key_by_atom_position.at(atom_position)]
            .emplace_back(model_list.at(atom_position));
    }

    std::unordered_map<GroupKey, GaussianModel3D> median_model_by_group;
    median_model_by_group.reserve(model_list_by_group.size());
    for (const auto & [group_key, group_model_list] : model_list_by_group)
    {
        const auto median_model{ BuildGaussianParameterMedian(group_model_list) };
        if (median_model.has_value())
        {
            median_model_by_group.emplace(group_key, *median_model);
        }
    }

    std::vector<GaussianModel3D> group_median_model_list;
    group_median_model_list.reserve(model_list.size());
    for (std::size_t atom_position = 0; atom_position < model_list.size(); atom_position++)
    {
        const auto median_iter{ median_model_by_group.find(
            group_key_by_atom_position.at(atom_position)) };
        group_median_model_list.emplace_back(
            median_iter != median_model_by_group.end() ?
                median_iter->second : model_list.at(atom_position));
    }
    return group_median_model_list;
}

inline std::vector<double> BuildGroupMedianOffsetList(
    const std::vector<GroupKey> & group_key_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    const auto median_model_list{ BuildGroupMedianModelList(group_key_by_atom_position, model_list) };
    std::vector<double> offset_list;
    offset_list.reserve(median_model_list.size());
    for (const auto & model : median_model_list)
    {
        offset_list.emplace_back(model.GetOffset());
    }
    return offset_list;
}

inline std::optional<std::vector<GaussianModel3D>> BuildSharedOffsetDampedModelList(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & raw_model_list,
    const std::vector<double> & previous_shared_offset_list,
    const std::vector<double> & raw_shared_offset_list,
    double damping)
{
    if (raw_model_list.size() != previous_model_list.size() ||
        previous_shared_offset_list.size() != previous_model_list.size() ||
        raw_shared_offset_list.size() != previous_model_list.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        return std::nullopt;
    }

    std::vector<GaussianModel3D> candidate_model_list;
    candidate_model_list.reserve(previous_model_list.size());
    for (std::size_t atom_position = 0; atom_position < previous_model_list.size(); atom_position++)
    {
        const auto previous_coordinates{
            EncodeTransformedCoordinates(previous_model_list.at(atom_position))
        };
        const auto raw_coordinates{
            EncodeTransformedCoordinates(raw_model_list.at(atom_position))
        };
        if (!previous_coordinates.has_value() || !raw_coordinates.has_value()) return std::nullopt;

        Eigen::Vector3d shape_coordinates{
            (*previous_coordinates + damping * (*raw_coordinates - *previous_coordinates)).eval()
        };
        shape_coordinates(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) = 0.0;
        const auto shape_model{ DecodeTransformedCoordinates(shape_coordinates) };
        if (!shape_model.has_value()) return std::nullopt;

        const auto previous_shared_offset{ previous_shared_offset_list.at(atom_position) };
        const auto raw_shared_offset{ raw_shared_offset_list.at(atom_position) };
        const auto candidate_model{
            shape_model->WithOffset(previous_shared_offset + damping * (raw_shared_offset - previous_shared_offset))
        };
        if (!IsValidSecondStageGaussianModel(candidate_model)) return std::nullopt;
        candidate_model_list.emplace_back(candidate_model);
    }
    return candidate_model_list;
}

} // namespace rhbm_gem::core::detail
