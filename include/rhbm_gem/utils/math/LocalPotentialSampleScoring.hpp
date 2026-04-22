#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/SamplingPointAcceptanceMask.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace detail {

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

inline float ComputeLocalPotentialCleanScore(
    const Eigen::Vector3d & sampling_position,
    const std::vector<Eigen::Vector3d> & reject_position_list)
{
    constexpr double clean_score_sigma{ 0.5 };
    constexpr double two_sigma_square{ 2.0 * clean_score_sigma * clean_score_sigma };

    double contamination_score{ 0.0 };
    for (const auto & reject_position : reject_position_list)
    {
        contamination_score += std::exp(
            -(sampling_position - reject_position).squaredNorm() / two_sigma_square);
    }

    return static_cast<float>(1.0 / (1.0 + contamination_score));
}

} // namespace detail

inline std::vector<float> BuildDefaultLocalPotentialSampleScoreList(
    std::size_t sample_count)
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
    return BuildSamplingPointAcceptanceMask(point_list, reject_direction_list, angle);
}

inline std::vector<float> BuildLocalPotentialCleanSampleScoreList(
    const SamplingPointList & point_list,
    const std::vector<Eigen::VectorXd> & reject_direction_list)
{
    std::vector<float> score_list(point_list.size(), 1.0f);
    if (reject_direction_list.empty())
    {
        return score_list;
    }

    const auto reject_position_list{
        detail::BuildLocalPotentialRejectPositionList(reject_direction_list)
    };
    for (size_t i = 0; i < point_list.size(); ++i)
    {
        score_list.at(i) = detail::ComputeLocalPotentialCleanScore(
            detail::ToVector3d(point_list.at(i).position),
            reject_position_list);
    }

    return score_list;
}

} // namespace rhbm_gem
