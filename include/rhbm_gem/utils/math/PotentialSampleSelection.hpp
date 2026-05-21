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
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace detail {

inline Eigen::Vector3d BuildLocalPotentialDirection(
    const std::array<float, 3> & position,
    const std::array<float, 3> & reference_position)
{
    return Eigen::Vector3d{ eigen_helper::ToEigenVector(position) } -
        Eigen::Vector3d{ eigen_helper::ToEigenVector(reference_position) };
}

inline std::vector<Eigen::Vector3d> BuildLocalPotentialNormalizedRejectDirectionList(
    const std::vector<Eigen::Vector3d> & reject_position_list,
    const Eigen::Vector3d & reference_position)
{
    std::vector<Eigen::Vector3d> normalized_reject_direction_list;
    normalized_reject_direction_list.reserve(reject_position_list.size());

    for (const auto & reject_position : reject_position_list)
    {
        const auto direction{ reject_position - reference_position };
        const auto norm{ direction.norm() };
        if (norm <= 0.0)
        {
            continue;
        }

        normalized_reject_direction_list.emplace_back(direction / norm);
    }

    return normalized_reject_direction_list;
}

inline std::vector<Eigen::Vector3d> BuildLocalPotentialNeighborPositionList(
    const std::vector<std::array<float, 3>> & reject_position_list,
    const std::array<float, 3> & reference_position)
{
    const Eigen::Vector3d reference_vector{
        eigen_helper::ToEigenVector(reference_position)
    };

    std::vector<Eigen::Vector3d> neighbor_position_list;
    neighbor_position_list.reserve(reject_position_list.size());
    for (const auto & reject_position : reject_position_list)
    {
        numeric_validation::RequireAllFinite(
            reject_position,
            "PotentialSampleSelection reject positions");
        const Eigen::Vector3d neighbor_position{
            eigen_helper::ToEigenVector(reject_position)
        };
        if ((neighbor_position - reference_vector).squaredNorm() <= 0.0)
        {
            continue;
        }

        neighbor_position_list.emplace_back(neighbor_position);
    }

    return neighbor_position_list;
}

inline bool IsLocalPotentialSamplingPointOwnedByReference(
    const SamplingPoint & point,
    const Eigen::Vector3d & reference_position,
    const std::vector<Eigen::Vector3d> & neighbor_position_list)
{
    const Eigen::Vector3d sampling_position{
        eigen_helper::ToEigenVector(point.position)
    };
    const auto reference_distance_squared{
        (sampling_position - reference_position).squaredNorm()
    };

    for (const auto & neighbor_position : neighbor_position_list)
    {
        if (!(reference_distance_squared
            < (sampling_position - neighbor_position).squaredNorm()))
        {
            return false;
        }
    }

    return true;
}

inline bool ShouldRejectLocalPotentialSamplingPointByNeighborDistance(
    const SamplingPoint & point,
    const Eigen::Vector3d & reference_position,
    const std::vector<Eigen::Vector3d> & neighbor_position_list)
{
    const Eigen::Vector3d sampling_position{
        eigen_helper::ToEigenVector(point.position)
    };
    const auto reference_distance_squared{
        (sampling_position - reference_position).squaredNorm()
    };

    for (const auto & neighbor_position : neighbor_position_list)
    {
        if ((sampling_position - neighbor_position).squaredNorm()
            < reference_distance_squared)
        {
            return true;
        }
    }

    return false;
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

} // namespace detail

inline std::vector<std::size_t> BuildSelectedLocalPotentialSampleIndexList(
    const SamplingPointList & point_list,
    const std::array<float, 3> & reference_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
    double angle = 0.0)
{
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    const Eigen::Vector3d reference_vector{
        eigen_helper::ToEigenVector(reference_position)
    };
    const auto neighbor_position_list{
        detail::BuildLocalPotentialNeighborPositionList(
            reject_position_list,
            reference_position)
    };
    if (neighbor_position_list.empty())
    {
        return detail::BuildAllLocalPotentialSampleIndexList(point_list.size());
    }

    const auto normalized_reject_direction_list{
        angle == 0.0
            ? std::vector<Eigen::Vector3d>{}
            : detail::BuildLocalPotentialNormalizedRejectDirectionList(
                neighbor_position_list,
                reference_vector)
    };

    std::vector<std::size_t> selected_indices;
    selected_indices.reserve(point_list.size());
    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };
    for (std::size_t i = 0; i < point_list.size(); i++)
    {
        const auto & point{ point_list.at(i) };
        if (detail::ShouldRejectLocalPotentialSamplingPointByNeighborDistance(
            point,
            reference_vector,
            neighbor_position_list))
        {
            continue;
        }

        if (!detail::IsLocalPotentialSamplingPointOwnedByReference(
            point,
            reference_vector,
            neighbor_position_list))
        {
            continue;
        }

        if (angle == 0.0
            || !detail::ShouldRejectLocalPotentialSamplingPoint(
                point,
                reference_position,
                normalized_reject_direction_list,
                cos_threshold))
        {
            selected_indices.emplace_back(i);
        }
    }

    return selected_indices;
}

inline LocalPotentialSampleList FilterLocalPotentialSampleList(
    LocalPotentialSampleList sample_list)
{
    return detail::KeepLowestResponseDecileByDistance(std::move(sample_list));
}

} // namespace rhbm_gem
