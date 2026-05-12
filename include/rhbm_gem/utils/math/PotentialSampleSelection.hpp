#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
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

inline Eigen::Vector3d ValidateAndConvertLocalPotentialRejectDirection(
    const Eigen::VectorXd & reject_direction)
{
    try
    {
        return eigen_validation::RequireVector3d(
            reject_direction,
            "PotentialSampleSelection reject directions");
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument(
            reject_direction.size() != 3 ?
                "PotentialSampleSelection reject directions must have dimension 3." :
                "PotentialSampleSelection reject directions must contain finite values.");
    }
}

inline std::vector<Eigen::Vector3d> BuildLocalPotentialRejectPositionList(
    const std::vector<Eigen::VectorXd> & reject_direction_list)
{
    std::vector<Eigen::Vector3d> reject_position_list;
    reject_position_list.reserve(reject_direction_list.size());

    for (const auto & reject_direction : reject_direction_list)
    {
        reject_position_list.emplace_back(
            ValidateAndConvertLocalPotentialRejectDirection(reject_direction));
    }

    return reject_position_list;
}

inline std::vector<Eigen::Vector3d> BuildLocalPotentialNormalizedRejectDirectionList(
    const std::vector<Eigen::VectorXd> & reject_direction_list)
{
    std::vector<Eigen::Vector3d> normalized_reject_direction_list;
    normalized_reject_direction_list.reserve(reject_direction_list.size());

    for (const auto & reject_direction : reject_direction_list)
    {
        const auto direction{
            ValidateAndConvertLocalPotentialRejectDirection(reject_direction)
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
    const std::vector<Eigen::Vector3d> & normalized_reject_direction_list,
    double cos_threshold)
{
    const Eigen::Vector3d sampling_vector{ ToVector3d(point.position) };
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
    const std::vector<Eigen::VectorXd> & reject_direction_list,
    double angle = 0.0)
{
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    if (angle == 0.0)
    {
        return detail::BuildAllLocalPotentialSampleIndexList(point_list.size());
    }

    const auto normalized_reject_direction_list{
        detail::BuildLocalPotentialNormalizedRejectDirectionList(reject_direction_list)
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
            normalized_reject_direction_list,
            cos_threshold))
        {
            selected_indices.emplace_back(i);
        }
    }

    return selected_indices;
}

inline LocalPotentialSampleList KeepLowestResponseDecileByDistance(
    LocalPotentialSampleList sample_list)
{
    std::map<float, LocalPotentialSampleList> samples_by_distance;
    for (auto & sample : sample_list)
    {
        samples_by_distance[sample.distance].emplace_back(std::move(sample));
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
                    std::ceil(static_cast<double>(distance_samples.size()) * 0.1)))
        };
        retained_samples.insert(
            retained_samples.end(),
            std::make_move_iterator(distance_samples.begin()),
            std::make_move_iterator(
                distance_samples.begin() + static_cast<std::ptrdiff_t>(keep_count)));
    }

    return retained_samples;
}

} // namespace rhbm_gem
