#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/SamplingPointFilter.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace detail {

inline float MakeInterpolationInMapObject(
    const MapObject & data_object, const std::array<float, 3> & position)
{
    auto index{ data_object.GetIndexFromPosition(position) };
    auto origin{ data_object.GetOrigin() };
    auto grid_spacing{ data_object.GetGridSpacing() };

    std::array<float, 3> local;
    for (size_t i = 0; i < 3; i++)
    {
        local.at(i) =
            (position.at(i) - origin.at(i) - static_cast<float>(index.at(i)) * grid_spacing.at(i))
            / grid_spacing.at(i);
    }

    auto cubic_interpolate = [](float p0, float p1, float p2, float p3, float t)
    {
        float a0{ p1 };
        float a1{ 0.5f * (p2 - p0) };
        float a2{ 0.5f * (-p3 + 4.0f * p2 - 5.0f * p1 + 2.0f * p0) };
        float a3{ 0.5f * (p3 - 3.0f * p2 + 3.0f * p1 - p0) };
        return a3 * t * t * t + a2 * t * t + a1 * t + a0;
    };

    std::array<std::array<std::array<float, 4>, 4>, 4> values;
    const auto grid_size{ data_object.GetGridSize() };
    for (int i = -1; i < 3; i++)
    {
        for (int j = -1; j < 3; j++)
        {
            for (int k = -1; k < 3; k++)
            {
                size_t i_next{ static_cast<size_t>(i + 1) };
                size_t j_next{ static_cast<size_t>(j + 1) };
                size_t k_next{ static_cast<size_t>(k + 1) };
                int xi{ std::clamp(index.at(0) + i, 0, grid_size.at(0) - 1) };
                int yj{ std::clamp(index.at(1) + j, 0, grid_size.at(1) - 1) };
                int zk{ std::clamp(index.at(2) + k, 0, grid_size.at(2) - 1) };
                values[i_next][j_next][k_next] = data_object.GetMapValue(xi, yj, zk);
            }
        }
    }

    std::array<std::array<float, 4>, 4> temp_y;
    for (size_t j = 0; j < 4; j++)
    {
        for (size_t k = 0; k < 4; k++)
        {
            temp_y[j][k] = cubic_interpolate(
                values[0][j][k], values[1][j][k], values[2][j][k], values[3][j][k], local.at(0));
        }
    }

    std::array<float, 4> temp_z;
    for (size_t k = 0; k < 4; k++)
    {
        temp_z[k] = cubic_interpolate(
            temp_y[0][k], temp_y[1][k], temp_y[2][k], temp_y[3][k], local.at(1));
    }

    return cubic_interpolate(temp_z[0], temp_z[1], temp_z[2], temp_z[3], local.at(2));
}

inline void ValidateNeighborRadius(double neighbor_radius)
{
    if (!std::isfinite(neighbor_radius) || neighbor_radius < 0.0)
    {
        throw std::invalid_argument(
            "SampleMapValues neighbor radius must be finite and non-negative.");
    }
}

inline std::vector<Eigen::VectorXd> BuildRejectDirectionList(
    const AtomObject & atom,
    const std::vector<AtomObject *> & neighbor_atom_list)
{
    std::vector<Eigen::VectorXd> reject_direction_list;
    reject_direction_list.reserve(neighbor_atom_list.size());

    const auto center_position{ atom.GetPosition() };
    for (const auto * neighbor_atom : neighbor_atom_list)
    {
        const auto neighbor_position{ neighbor_atom->GetPosition() };
        Eigen::VectorXd direction{ Eigen::VectorXd::Zero(3) };
        direction(0) = static_cast<double>(neighbor_position[0] - center_position[0]);
        direction(1) = static_cast<double>(neighbor_position[1] - center_position[1]);
        direction(2) = static_cast<double>(neighbor_position[2] - center_position[2]);
        reject_direction_list.emplace_back(direction);
    }

    return reject_direction_list;
}

inline SamplingPointList TranslateSamplingPoints(
    const SamplingPointList & point_list,
    const std::array<float, 3> & reference_position,
    float scale)
{
    SamplingPointList translated_point_list;
    translated_point_list.reserve(point_list.size());
    for (const auto & point : point_list)
    {
        translated_point_list.emplace_back(SamplingPoint{
            point.distance,
            {
                point.position[0] + scale * reference_position[0],
                point.position[1] + scale * reference_position[1],
                point.position[2] + scale * reference_position[2]
            }
        });
    }
    return translated_point_list;
}

LocalPotentialSampleList BuildLocalPotentialSampleList(
    const MapObject & map_object,
    const SamplingPointList & sampling_points,
    const std::vector<float> * sampling_scores = nullptr)
{
    if (sampling_scores != nullptr && sampling_scores->size() != sampling_points.size())
    {
        throw std::invalid_argument(
            "BuildLocalPotentialSampleList sampling_scores must align with sampling_points.");
    }

    LocalPotentialSampleList sampling_data_list;
    sampling_data_list.reserve(sampling_points.size());
    for (size_t i = 0; i < sampling_points.size(); i++)
    {
        const auto & sampling_point{ sampling_points.at(i) };
        const auto sampling_score{
            (sampling_scores == nullptr) ? 1.0f : sampling_scores->at(i)
        };

        auto map_value{
            MakeInterpolationInMapObject(map_object, sampling_point.position)
        };
        sampling_data_list.emplace_back(LocalPotentialSample{
            sampling_point.distance,
            map_value,
            sampling_score,
            sampling_point.position
        });
    }
    return sampling_data_list;
}

} // namespace detail

// Prefer overloads that match the sampler's real inputs instead of forcing every sampler to
// accept the same parameter list. Add new overloads here only when a sampler introduces a new
// input shape that SampleMapValues needs to support.
template <typename Sampler>
LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const Sampler & sampler,
    const std::array<float, 3> & position)
{
    const auto sampling_points{ sampler.GenerateSamplingPoints(position) };
    return detail::BuildLocalPotentialSampleList(map_object, sampling_points);
}

template <typename Sampler>
LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const Sampler & sampler,
    const std::array<float, 3> & position,
    const std::array<float, 3> & direction_like_input)
{
    const auto sampling_points{ sampler.GenerateSamplingPoints(position, direction_like_input) };
    return detail::BuildLocalPotentialSampleList(map_object, sampling_points);
}

template <typename Sampler>
LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const Sampler & sampler,
    const AtomObject & atom,
    double neighbor_radius,
    double angle = 0.0)
{
    detail::ValidateNeighborRadius(neighbor_radius);
    detail::ValidateSamplingPointFilterAngle(angle);

    const auto position{ atom.GetPosition() };
    const auto sampling_points{ sampler.GenerateSamplingPoints(position) };
    if (angle <= 0.0)
    {
        return detail::BuildLocalPotentialSampleList(map_object, sampling_points);
    }

    const auto neighbor_atom_list{ atom.FindNeighborAtoms(neighbor_radius, false) };
    const auto reject_direction_list{
        detail::BuildRejectDirectionList(atom, neighbor_atom_list)
    };
    const auto local_sampling_points{
        detail::TranslateSamplingPoints(sampling_points, position, -1.0f)
    };
    const auto sampling_scores{
        BuildSamplingPointScoreList(local_sampling_points, reject_direction_list, angle)
    };
    return detail::BuildLocalPotentialSampleList(map_object, sampling_points, &sampling_scores);
}

} // namespace rhbm_gem
