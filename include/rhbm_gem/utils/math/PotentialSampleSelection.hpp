#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace {

inline bool IsOwnedByNeighbor(
    const std::array<float, 3> & sample_position,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & valid_neighbor_list)
{
    bool is_bad_point{ false };
    const auto to_local_distance{
        array_helper::ComputeNorm(sample_position, local_position)
    };
    for (const auto & neighbor_position : valid_neighbor_list)
    {
        const auto to_neighbor_distance{
            array_helper::ComputeNorm(sample_position, neighbor_position)
        };
        if (to_local_distance > to_neighbor_distance) is_bad_point = true;
    }
    return is_bad_point;
}

inline bool IsInsideNeighborCone(
    const std::array<float, 3> & sample_position,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & valid_neighbor_list,
    double angle)
{
    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };
    const auto sampling_unit_vector{
        array_helper::ComputeVector(local_position, sample_position, true)
    };

    for (const auto & neighbor_position : valid_neighbor_list)
    {
        const auto neighbor_unit_vector{
            array_helper::ComputeVector(local_position, neighbor_position, true)
        };
        const auto cos_theta{
            array_helper::ComputeDotProduct(sampling_unit_vector, neighbor_unit_vector)
        };
        if (cos_theta > cos_threshold) return true;
    }
    return false;
}

inline LocalPotentialSampleList KeepLowestResponseDecileByDistance(LocalPotentialSampleList sample_list)
{
    constexpr double kRetainedRatio{ 0.1 };

    std::map<float, LocalPotentialSampleList> samples_by_distance;
    for (auto & sample : sample_list)
    {
        samples_by_distance[sample.point.distance].emplace_back(std::move(sample));
    }

    LocalPotentialSampleList retained_samples;
    for (auto & distance_entry : samples_by_distance)
    {
        auto & distance_samples{ distance_entry.second };
        std::stable_sort(
            distance_samples.begin(),
            distance_samples.end(),
            [](const LocalPotentialSample & lhs, const LocalPotentialSample & rhs)
            {
                return lhs.response < rhs.response;
            });

        const auto keep_count{
            std::max<std::size_t>(
                1,
                static_cast<std::size_t>(
                    std::ceil(static_cast<double>(distance_samples.size()) * kRetainedRatio)))
        };
        retained_samples.insert(
            retained_samples.end(),
            std::make_move_iterator(distance_samples.begin()),
            std::make_move_iterator(
                distance_samples.begin() + static_cast<std::ptrdiff_t>(keep_count)));
    }

    return retained_samples;
}

} // namespace

inline std::vector<std::size_t> BuildSelectedLocalPotentialSampleIndexList(
    const SamplingPointList & sample_point_list,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
    double angle = 0.0)
{
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    if (reject_position_list.empty())
    {
        std::vector<std::size_t> selected_indices;
        selected_indices.resize(sample_point_list.size());
        std::iota(selected_indices.begin(), selected_indices.end(), std::size_t{ 0 });
        return selected_indices;
    }

    std::vector<std::array<float, 3>> valid_neighbor_list;
    valid_neighbor_list.reserve(reject_position_list.size());
    for (const auto & reject_position : reject_position_list)
    {
        numeric_validation::RequireAllFinite(reject_position, "reject positions");
        if (!numeric_validation::IsNonEqual(local_position, reject_position)) continue;
        valid_neighbor_list.emplace_back(reject_position);
    }

    std::vector<std::size_t> selected_indices;
    selected_indices.reserve(sample_point_list.size());
    for (std::size_t i = 0; i < sample_point_list.size(); i++)
    {
        const auto & sample_position{ sample_point_list.at(i).position };
        if (IsOwnedByNeighbor(sample_position, local_position, valid_neighbor_list) ||
            IsInsideNeighborCone(sample_position, local_position, valid_neighbor_list, angle))
        {
            continue;
        }
        selected_indices.emplace_back(i);
    }

    return selected_indices;
}

inline LocalPotentialSampleList FilterLocalPotentialSampleList(LocalPotentialSampleList sample_list)
{
    return KeepLowestResponseDecileByDistance(std::move(sample_list));
}

} // namespace rhbm_gem
