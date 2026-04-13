#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace detail {

inline void ValidateSamplingPointFilterAngle(double angle)
{
    if (!std::isfinite(angle) || angle < 0.0 || angle > 180.0)
    {
        throw std::invalid_argument(
            "SelectSamplingPoint angle must be finite and within [0, 180] degrees.");
    }
}

inline Eigen::Vector3d ToVector3d(const std::array<float, 3> & position)
{
    return Eigen::Vector3d{
        static_cast<double>(position[0]),
        static_cast<double>(position[1]),
        static_cast<double>(position[2])
    };
}

inline std::vector<Eigen::Vector3d> BuildNormalizedRejectDirections(
    const std::vector<Eigen::VectorXd> & reject_direction_list)
{
    std::vector<Eigen::Vector3d> normalized_reject_directions;
    normalized_reject_directions.reserve(reject_direction_list.size());

    for (const auto & reject_direction : reject_direction_list)
    {
        if (reject_direction.size() != 3)
        {
            throw std::invalid_argument(
                "SelectSamplingPoint reject directions must have dimension 3.");
        }

        if (!reject_direction.allFinite())
        {
            throw std::invalid_argument(
                "SelectSamplingPoint reject directions must contain finite values.");
        }

        const Eigen::Vector3d direction{
            reject_direction(0), reject_direction(1), reject_direction(2)
        };
        const auto norm{ direction.norm() };
        if (norm <= 0.0)
        {
            continue;
        }

        normalized_reject_directions.emplace_back(direction / norm);
    }

    return normalized_reject_directions;
}

} // namespace detail

inline SamplingPointList SelectSamplingPoint(
    const SamplingPointList & point_list,
    const std::vector<Eigen::VectorXd> & reject_direction_list,
    double angle = 0.0)
{
    detail::ValidateSamplingPointFilterAngle(angle);
    if (angle == 0.0)
    {
        return point_list;
    }

    const auto normalized_reject_directions{
        detail::BuildNormalizedRejectDirections(reject_direction_list)
    };
    if (normalized_reject_directions.empty())
    {
        return point_list;
    }

    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };

    SamplingPointList filtered_point_list;
    filtered_point_list.reserve(point_list.size());
    for (const auto & point : point_list)
    {
        const Eigen::Vector3d sampling_vector{ detail::ToVector3d(point.position) };
        const auto sampling_norm{ sampling_vector.norm() };
        if (sampling_norm <= 0.0)
        {
            filtered_point_list.emplace_back(point);
            continue;
        }

        const Eigen::Vector3d normalized_sampling_direction{ sampling_vector / sampling_norm };
        bool should_reject{ false };
        for (const auto & reject_direction : normalized_reject_directions)
        {
            if (normalized_sampling_direction.dot(reject_direction) > cos_threshold)
            {
                should_reject = true;
                break;
            }
        }

        if (!should_reject)
        {
            filtered_point_list.emplace_back(point);
        }
    }

    return filtered_point_list;
}

} // namespace rhbm_gem
