#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/TransformedChange.hpp"

namespace rhbm_gem::core::detail {

inline std::optional<GaussianModel3D> BuildLocalFittingGaussianParameterMedian(
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

inline std::vector<GaussianModel3D> BuildLocalFittingGroupMedianModelList(
    const std::vector<GroupKey> & group_key_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    if (group_key_by_atom_position.size() != model_list.size())
    {
        throw std::invalid_argument(
            "Local fitting group median inputs are inconsistent.");
    }

    std::unordered_map<GroupKey, std::vector<GaussianModel3D>> model_list_by_group;
    model_list_by_group.reserve(model_list.size());
    for (std::size_t atom_position = 0;
        atom_position < model_list.size();
        atom_position++)
    {
        model_list_by_group[group_key_by_atom_position.at(atom_position)]
            .emplace_back(model_list.at(atom_position));
    }

    std::unordered_map<GroupKey, GaussianModel3D> median_model_by_group;
    median_model_by_group.reserve(model_list_by_group.size());
    for (const auto & [group_key, group_model_list] : model_list_by_group)
    {
        const auto median_model{
            BuildLocalFittingGaussianParameterMedian(group_model_list)
        };
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

inline std::optional<std::vector<GaussianModel3D>> BuildLocalFittingSharedOffsetDampedModelList(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & raw_model_list,
    const std::vector<GaussianModel3D> & previous_shared_offset_model_list,
    const std::vector<GaussianModel3D> & raw_shared_offset_model_list,
    double damping)
{
    if (raw_model_list.size() != previous_model_list.size() ||
        previous_shared_offset_model_list.size() != previous_model_list.size() ||
        raw_shared_offset_model_list.size() != previous_model_list.size() ||
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
        if (!previous_coordinates.has_value() || !raw_coordinates.has_value())
        {
            return std::nullopt;
        }

        Eigen::Vector3d shape_coordinates{
            (*previous_coordinates + damping * (*raw_coordinates - *previous_coordinates)).eval()
        };
        shape_coordinates(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) = 0.0;
        const auto shape_model{
            DecodeTransformedCoordinates(shape_coordinates)
        };
        if (!shape_model.has_value()) return std::nullopt;

        const auto previous_shared_offset{
            previous_shared_offset_model_list.at(atom_position).GetOffset()
        };
        const auto raw_shared_offset{
            raw_shared_offset_model_list.at(atom_position).GetOffset()
        };
        const auto candidate_model{ shape_model->WithOffset(
            previous_shared_offset + damping * (raw_shared_offset - previous_shared_offset)) };
        if (!IsValidSecondStageGaussianModel(candidate_model))
        {
            return std::nullopt;
        }
        candidate_model_list.emplace_back(candidate_model);
    }
    return candidate_model_list;
}

} // namespace rhbm_gem::core::detail
