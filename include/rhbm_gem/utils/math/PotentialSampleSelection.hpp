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

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {
namespace detail {

struct LocalPotentialRejectContext
{
    std::vector<Eigen::VectorXd> neighbor_position_list;
    std::vector<Eigen::VectorXd> normalized_reject_direction_list;
};

inline LocalPotentialRejectContext BuildLocalPotentialRejectContext(
    const std::vector<std::array<float, 3>> & reject_position_list,
    const Eigen::VectorXd & reference_position,
    bool build_normalized_directions)
{
    LocalPotentialRejectContext context;
    context.neighbor_position_list.reserve(reject_position_list.size());
    if (build_normalized_directions)
    {
        context.normalized_reject_direction_list.reserve(reject_position_list.size());
    }

    for (const auto & reject_position : reject_position_list)
    {
        numeric_validation::RequireAllFinite(reject_position, "PotentialSampleSelection reject positions");
        const Eigen::VectorXd neighbor_position{ eigen_helper::ToEigenVector(reject_position) };
        const auto direction{ neighbor_position - reference_position };
        const auto direction_norm{ direction.norm() };
        if (direction_norm <= 0.0)
        {
            continue;
        }

        context.neighbor_position_list.emplace_back(neighbor_position);
        if (build_normalized_directions)
        {
            context.normalized_reject_direction_list.emplace_back(direction / direction_norm);
        }
    }

    return context;
}

inline bool IsLocalPotentialSamplingPointOwnedByReference(
    const Eigen::VectorXd & sampling_position,
    const Eigen::VectorXd & reference_position,
    const std::vector<Eigen::VectorXd> & neighbor_position_list)
{
    const auto distance_squared{ (sampling_position - reference_position).squaredNorm() };

    for (const auto & neighbor_position : neighbor_position_list)
    {
        if (!(distance_squared < (sampling_position - neighbor_position).squaredNorm()))
        {
            return false;
        }
    }

    return true;
}

inline bool ShouldRejectLocalPotentialSamplingPoint(
    const Eigen::VectorXd & sampling_position,
    const Eigen::VectorXd & reference_position,
    const std::vector<Eigen::VectorXd> & normalized_reject_direction_list,
    double cos_threshold)
{
    const auto sampling_vector{ sampling_position - reference_position };
    const auto sampling_norm{ sampling_vector.norm() };
    if (sampling_norm <= 0.0)
    {
        return false;
    }

    const Eigen::VectorXd normalized_sampling_direction{ sampling_vector / sampling_norm };
    for (const auto & reject_direction : normalized_reject_direction_list)
    {
        if (normalized_sampling_direction.dot(reject_direction) > cos_threshold)
        {
            return true;
        }
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

} // namespace detail

inline std::vector<std::size_t> BuildSelectedLocalPotentialSampleIndexList(
    const SamplingPointList & point_list,
    const std::array<float, 3> & reference_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
    double angle = 0.0)
{
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    const Eigen::VectorXd reference_vector{ eigen_helper::ToEigenVector(reference_position) };
    const auto reject_context{
        detail::BuildLocalPotentialRejectContext(
            reject_position_list, reference_vector, angle > 0.0)
    };
    if (reject_context.neighbor_position_list.empty())
    {
        std::vector<std::size_t> selected_indices;
        selected_indices.resize(point_list.size());
        std::iota(selected_indices.begin(), selected_indices.end(), std::size_t{ 0 });
        return selected_indices;
    }

    std::vector<std::size_t> selected_indices;
    selected_indices.reserve(point_list.size());
    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };
    for (std::size_t i = 0; i < point_list.size(); i++)
    {
        const Eigen::VectorXd sampling_position{
            eigen_helper::ToEigenVector(point_list.at(i).position)
        };
        if (!detail::IsLocalPotentialSamplingPointOwnedByReference(
            sampling_position,
            reference_vector,
            reject_context.neighbor_position_list))
        {
            continue;
        }

        if (angle > 0.0
            && detail::ShouldRejectLocalPotentialSamplingPoint(
                sampling_position,
                reference_vector,
                reject_context.normalized_reject_direction_list,
                cos_threshold))
        {
            continue;
        }

        selected_indices.emplace_back(i);
    }

    return selected_indices;
}

inline LocalPotentialSampleList FilterLocalPotentialSampleList(LocalPotentialSampleList sample_list)
{
    return detail::KeepLowestResponseDecileByDistance(std::move(sample_list));
}

} // namespace rhbm_gem
