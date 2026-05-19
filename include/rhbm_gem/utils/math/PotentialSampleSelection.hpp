#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace detail {

inline Eigen::Vector3d ToVector3d(const std::array<float, 3> & position)
{
    return Eigen::Vector3d{
        static_cast<double>(position[0]),
        static_cast<double>(position[1]),
        static_cast<double>(position[2])
    };
}

inline Eigen::Vector3d BuildLocalPotentialDirection(
    const std::array<float, 3> & position,
    const std::array<float, 3> & reference_position)
{
    return ToVector3d(position) - ToVector3d(reference_position);
}

inline std::vector<Eigen::Vector3d> BuildLocalPotentialNormalizedRejectDirectionList(
    const std::vector<std::array<float, 3>> & reject_position_list,
    const std::array<float, 3> & reference_position)
{
    std::vector<Eigen::Vector3d> normalized_reject_direction_list;
    normalized_reject_direction_list.reserve(reject_position_list.size());

    for (const auto & reject_position : reject_position_list)
    {
        numeric_validation::RequireAllFinite(
            reject_position,
            "PotentialSampleSelection reject positions");
        const auto direction{
            BuildLocalPotentialDirection(reject_position, reference_position)
        };
        const auto norm{ direction.norm() };
        if (norm <= 0.0)
        {
            continue;
        }

        normalized_reject_direction_list.emplace_back(direction / norm);
    }

    return normalized_reject_direction_list;
}

inline bool ShouldRejectLocalPotentialSamplingPoint(
    const SamplingPoint & point,
    const std::array<float, 3> & reference_position,
    const std::vector<Eigen::Vector3d> & normalized_reject_direction_list,
    double cos_threshold)
{
    const Eigen::Vector3d sampling_vector{
        BuildLocalPotentialDirection(point.position, reference_position)
    };
    const auto sampling_norm{ sampling_vector.norm() };
    if (sampling_norm <= 0.0)
    {
        return false;
    }

    const Eigen::Vector3d normalized_sampling_direction{ sampling_vector / sampling_norm };
    for (const auto & reject_direction : normalized_reject_direction_list)
    {
        if (normalized_sampling_direction.dot(reject_direction) > cos_threshold)
        {
            return true;
        }
    }

    return false;
}

inline std::vector<std::size_t> BuildAllLocalPotentialSampleIndexList(std::size_t sample_count)
{
    std::vector<std::size_t> selected_indices;
    selected_indices.reserve(sample_count);
    for (std::size_t i = 0; i < sample_count; i++)
    {
        selected_indices.emplace_back(i);
    }
    return selected_indices;
}

} // namespace detail

inline std::vector<std::size_t> BuildSelectedLocalPotentialSampleIndexList(
    const SamplingPointList & point_list,
    const std::array<float, 3> & reference_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
    double angle = 0.0)
{
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    if (angle == 0.0)
    {
        return detail::BuildAllLocalPotentialSampleIndexList(point_list.size());
    }

    const auto normalized_reject_direction_list{
        detail::BuildLocalPotentialNormalizedRejectDirectionList(
            reject_position_list,
            reference_position)
    };
    if (normalized_reject_direction_list.empty())
    {
        return detail::BuildAllLocalPotentialSampleIndexList(point_list.size());
    }

    std::vector<std::size_t> selected_indices;
    selected_indices.reserve(point_list.size());
    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };
    for (std::size_t i = 0; i < point_list.size(); i++)
    {
        if (!detail::ShouldRejectLocalPotentialSamplingPoint(
            point_list.at(i),
            reference_position,
            normalized_reject_direction_list,
            cos_threshold))
        {
            selected_indices.emplace_back(i);
        }
    }

    return selected_indices;
}

inline LocalPotentialSampleList KeepLowestResponseDecileByDistance(
    LocalPotentialSampleList sample_list,
    double retained_ratio = 0.1)
{
    numeric_validation::RequireFiniteExclusiveInclusiveRange(
        retained_ratio, 0.0, 1.0, "retained_ratio");

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
                    std::ceil(static_cast<double>(distance_samples.size()) * retained_ratio)))
        };
        retained_samples.insert(
            retained_samples.end(),
            std::make_move_iterator(distance_samples.begin()),
            std::make_move_iterator(
                distance_samples.begin() + static_cast<std::ptrdiff_t>(keep_count)));
    }

    return retained_samples;
}

inline LocalPotentialSampleList GetMedianResponseByDistance(
    LocalPotentialSampleList sample_list)
{
    std::map<float, LocalPotentialSampleList> samples_by_distance;
    for (auto & sample : sample_list)
    {
        samples_by_distance[sample.point.distance].emplace_back(std::move(sample));
    }

    LocalPotentialSampleList median_samples;
    median_samples.reserve(samples_by_distance.size());
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

        const auto lower_index{ (distance_samples.size() - 1) / 2 };
        const auto upper_index{ distance_samples.size() / 2 };
        const auto median_response{
            lower_index == upper_index
                ? distance_samples.at(lower_index).response
                : (distance_samples.at(lower_index).response
                    + distance_samples.at(upper_index).response) / 2.0f
        };
        auto median_sample{ std::move(distance_samples.at(lower_index)) };
        median_sample.response = median_response;
        median_samples.emplace_back(std::move(median_sample));
    }

    return median_samples;
}

} // namespace rhbm_gem
