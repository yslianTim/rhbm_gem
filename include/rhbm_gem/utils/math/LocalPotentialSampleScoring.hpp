#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
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
            "LocalPotentialSampleScoring reject directions");
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument(
            reject_direction.size() != 3 ?
                "LocalPotentialSampleScoring reject directions must have dimension 3." :
                "LocalPotentialSampleScoring reject directions must contain finite values.");
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

} // namespace detail

inline std::vector<float> BuildDefaultLocalPotentialSampleScoreList(std::size_t sample_count)
{
    return std::vector<float>(sample_count, 1.0f);
}

inline std::vector<float> BuildLocalPotentialSampleScoreList(
    const SamplingPointList & point_list,
    const std::vector<Eigen::VectorXd> & reject_direction_list,
    double angle = 0.0)
{
    // Callers provide score-policy inputs as plain data. This layer stays agnostic
    // to domain objects such as atoms and only owns the score calculation itself.
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    std::vector<float> score_list(point_list.size(), 1.0f);
    if (angle == 0.0)
    {
        return score_list;
    }

    const auto normalized_reject_direction_list{
        detail::BuildLocalPotentialNormalizedRejectDirectionList(reject_direction_list)
    };
    if (normalized_reject_direction_list.empty())
    {
        return score_list;
    }

    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };
    for (std::size_t i = 0; i < point_list.size(); i++)
    {
        if (detail::ShouldRejectLocalPotentialSamplingPoint(
                point_list.at(i),
                normalized_reject_direction_list,
                cos_threshold))
        {
            score_list.at(i) = 0.0f;
        }
    }

    return score_list;
}

} // namespace rhbm_gem
