#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "core/detail/LocalFittingStateView.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"

namespace rhbm_gem::core::detail {

constexpr std::array<double, kTransformedChangeSize>
kLocalFittingTrustRegionParameterScale{ 0.50, 0.35, 1.0 };

using LocalFittingTrustRegionClusterKey = std::vector<std::size_t>;

struct LocalFittingTrustRegionOptions
{
    double initial_radius{ 1.0 };
    double minimum_radius{ 0.0625 };
    double maximum_radius{ 4.0 };
    double shrink_factor{ 0.5 };
    double growth_factor{ 2.0 };
};

struct LocalFittingTrustRegionRadiusUpdate
{
    std::vector<LocalFittingTrustRegionClusterKey> changed_key_list{};
    std::vector<LocalFittingTrustRegionClusterKey> saturated_key_list{};
};

struct LocalFittingRejectedClusterPartition
{
    std::vector<LocalFittingTrustRegionClusterKey> exhausted_key_list{};
    std::vector<LocalFittingTrustRegionClusterKey> retryable_key_list{};
};

enum class LocalFittingAllRejectedResolution
{
    Retry,
    MaximumIterations,
    BacktrackingExhausted,
    MinimumRadius,
    NoRetryProgress
};

inline LocalFittingRejectedClusterPartition
PartitionLocalFittingRejectedClusters(
    const std::vector<LocalFittingTrustRegionClusterKey> & rejected_key_list,
    const std::vector<LocalFittingTrustRegionClusterKey> & exhausted_key_list)
{
    LocalFittingRejectedClusterPartition partition;
    for (const auto & key : rejected_key_list)
    {
        if (std::find(
                exhausted_key_list.begin(),
                exhausted_key_list.end(),
                key) != exhausted_key_list.end())
        {
            partition.exhausted_key_list.emplace_back(key);
        }
        else
        {
            partition.retryable_key_list.emplace_back(key);
        }
    }
    return partition;
}

inline LocalFittingAllRejectedResolution ResolveLocalFittingAllRejected(
    bool maximum_iterations_reached,
    const LocalFittingRejectedClusterPartition & partition,
    const LocalFittingTrustRegionRadiusUpdate & radius_update)
{
    if (maximum_iterations_reached)
    {
        return LocalFittingAllRejectedResolution::MaximumIterations;
    }
    if (partition.exhausted_key_list.empty() &&
        partition.retryable_key_list.empty())
    {
        throw std::invalid_argument(
            "All-rejected local fitting resolution requires rejected clusters.");
    }
    if (partition.retryable_key_list.empty())
    {
        return LocalFittingAllRejectedResolution::BacktrackingExhausted;
    }
    if (!radius_update.changed_key_list.empty())
    {
        return LocalFittingAllRejectedResolution::Retry;
    }

    const auto all_retryable_saturated{
        std::all_of(
            partition.retryable_key_list.begin(),
            partition.retryable_key_list.end(),
            [&](const LocalFittingTrustRegionClusterKey & key)
            {
                return std::find(
                    radius_update.saturated_key_list.begin(),
                    radius_update.saturated_key_list.end(),
                    key) != radius_update.saturated_key_list.end();
            })
    };
    if (!all_retryable_saturated)
    {
        throw std::logic_error(
            "Trust-region shrink did not classify every retryable rejection.");
    }
    return partition.exhausted_key_list.empty() ?
        LocalFittingAllRejectedResolution::MinimumRadius :
        LocalFittingAllRejectedResolution::NoRetryProgress;
}

class LocalFittingTrustRegionStateSet
{
private:
    LocalFittingTrustRegionOptions m_options{};
    std::map<LocalFittingTrustRegionClusterKey, double> m_radius_by_key{};

    static void ValidateOptions(const LocalFittingTrustRegionOptions & options)
    {
        if (!std::isfinite(options.initial_radius) ||
            !std::isfinite(options.minimum_radius) ||
            !std::isfinite(options.maximum_radius) ||
            options.minimum_radius <= 0.0 ||
            options.initial_radius < options.minimum_radius ||
            options.maximum_radius < options.initial_radius ||
            !std::isfinite(options.shrink_factor) ||
            options.shrink_factor <= 0.0 ||
            options.shrink_factor >= 1.0 ||
            !std::isfinite(options.growth_factor) ||
            options.growth_factor <= 1.0)
        {
            throw std::invalid_argument("Local fitting trust-region options are invalid.");
        }
    }

public:
    explicit LocalFittingTrustRegionStateSet(
        LocalFittingTrustRegionOptions options = {})
        : m_options{ options }
    {
        ValidateOptions(m_options);
    }

    void Reconcile(const std::vector<LocalFittingTrustRegionClusterKey> & key_list)
    {
        std::map<LocalFittingTrustRegionClusterKey, double> next_radius_by_key;
        for (const auto & key : key_list)
        {
            const auto iter{ m_radius_by_key.find(key) };
            next_radius_by_key.emplace(
                key,
                iter == m_radius_by_key.end() ?
                    m_options.initial_radius : iter->second);
        }
        m_radius_by_key = std::move(next_radius_by_key);
    }

    double GetRadius(const LocalFittingTrustRegionClusterKey & key) const
    {
        const auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument("Local fitting trust-region state is missing.");
        }
        return iter->second;
    }

    LocalFittingTrustRegionRadiusUpdate Shrink(
        const std::vector<LocalFittingTrustRegionClusterKey> & key_list)
    {
        LocalFittingTrustRegionRadiusUpdate update;
        for (const auto & key : key_list)
        {
            auto iter{ m_radius_by_key.find(key) };
            if (iter == m_radius_by_key.end())
            {
                throw std::invalid_argument("Local fitting trust-region state is missing.");
            }
            if (iter->second <= m_options.minimum_radius)
            {
                update.saturated_key_list.emplace_back(key);
                continue;
            }
            iter->second = std::max(
                m_options.minimum_radius,
                iter->second * m_options.shrink_factor);
            update.changed_key_list.emplace_back(key);
        }
        return update;
    }

    void Grow(const std::vector<LocalFittingTrustRegionClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_radius_by_key.find(key) };
            if (iter == m_radius_by_key.end())
            {
                throw std::invalid_argument("Local fitting trust-region state is missing.");
            }
            iter->second = std::min(
                m_options.maximum_radius,
                iter->second * m_options.growth_factor);
        }
    }

};

struct LocalFittingTrustRegionDamping
{
    double effective_damping{ 1.0 };
    double step_norm{ 0.0 };
};

inline void ValidateLocalFittingTrustRegionInputs(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    const std::array<double, 3> & parameter_scale,
    double requested_damping,
    double radius)
{
    if (candidate_estimation_list.size() != previous_estimation_list.size() ||
        !std::isfinite(requested_damping) ||
        requested_damping <= 0.0 ||
        requested_damping > 1.0 ||
        !std::isfinite(radius) ||
        radius <= 0.0)
    {
        throw std::invalid_argument("Local fitting trust-region inputs are invalid.");
    }
    for (const auto scale : parameter_scale)
    {
        if (!std::isfinite(scale) || scale <= 0.0)
        {
            throw std::invalid_argument("Local fitting trust-region scale is invalid.");
        }
    }
    for (std::size_t i = 0; i < previous_estimation_list.size(); i++)
    {
        if (!previous_estimation_list.at(i).allFinite() ||
            !candidate_estimation_list.at(i).allFinite())
        {
            throw std::invalid_argument("Local fitting trust-region estimation is invalid.");
        }
    }
}

inline LocalFittingTrustRegionDamping LimitLocalFittingTrustRegionDamping(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    const std::array<double, 3> & parameter_scale,
    double requested_damping,
    double radius)
{
    ValidateLocalFittingTrustRegionInputs(
        previous_estimation_list,
        candidate_estimation_list,
        parameter_scale,
        requested_damping,
        radius);

    double undamped_step_norm{ 0.0 };
    for (std::size_t i = 0; i < previous_estimation_list.size(); i++)
    {
        for (std::size_t parameter_index = 0; parameter_index < 3; parameter_index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
            undamped_step_norm = std::max(
                undamped_step_norm,
                std::abs(
                    candidate_estimation_list.at(i)(eigen_index) -
                    previous_estimation_list.at(i)(eigen_index)) /
                    parameter_scale.at(parameter_index));
        }
    }

    const auto effective_damping{
        undamped_step_norm > 0.0 ?
            std::min(requested_damping, radius / undamped_step_norm) :
            requested_damping
    };
    return LocalFittingTrustRegionDamping{
        effective_damping,
        effective_damping * undamped_step_norm
    };
}

inline LocalFittingTrustRegionDamping LimitLocalFittingTrustRegionSubstepDamping(
    const std::vector<Eigen::Vector3d> & outer_previous_estimation_list,
    const std::vector<Eigen::Vector3d> & substep_previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    const std::array<double, 3> & parameter_scale,
    double requested_damping,
    double radius)
{
    ValidateLocalFittingTrustRegionInputs(
        outer_previous_estimation_list,
        substep_previous_estimation_list,
        parameter_scale,
        requested_damping,
        radius);
    ValidateLocalFittingTrustRegionInputs(
        substep_previous_estimation_list,
        candidate_estimation_list,
        parameter_scale,
        requested_damping,
        radius);

    double maximum_damping{ requested_damping };
    constexpr double tolerance{ 1.0e-12 };
    for (std::size_t i = 0; i < outer_previous_estimation_list.size(); i++)
    {
        for (std::size_t parameter_index = 0; parameter_index < 3; parameter_index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
            const auto limit{ radius * parameter_scale.at(parameter_index) };
            const auto base_step{
                substep_previous_estimation_list.at(i)(eigen_index) -
                outer_previous_estimation_list.at(i)(eigen_index)
            };
            if (std::abs(base_step) > limit + tolerance)
            {
                throw std::invalid_argument(
                    "Local fitting trust-region substep starts outside the radius.");
            }
            const auto direction{
                candidate_estimation_list.at(i)(eigen_index) -
                substep_previous_estimation_list.at(i)(eigen_index)
            };
            if (direction > 0.0)
            {
                maximum_damping = std::min(
                    maximum_damping,
                    std::max(0.0, (limit - base_step) / direction));
            }
            else if (direction < 0.0)
            {
                maximum_damping = std::min(
                    maximum_damping,
                    std::max(0.0, (-limit - base_step) / direction));
            }
        }
    }

    double step_norm{ 0.0 };
    for (std::size_t i = 0; i < outer_previous_estimation_list.size(); i++)
    {
        for (std::size_t parameter_index = 0; parameter_index < 3; parameter_index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
            const auto base_step{
                substep_previous_estimation_list.at(i)(eigen_index) -
                outer_previous_estimation_list.at(i)(eigen_index)
            };
            const auto direction{
                candidate_estimation_list.at(i)(eigen_index) -
                substep_previous_estimation_list.at(i)(eigen_index)
            };
            step_norm = std::max(
                step_norm,
                std::abs(base_step + maximum_damping * direction) /
                    parameter_scale.at(parameter_index));
        }
    }
    return LocalFittingTrustRegionDamping{
        maximum_damping,
        step_norm
    };
}

inline std::optional<double> CalculateLocalFittingClusterModelTrustRegionStepNorm(
    const LocalFittingState & outer_previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<GaussianModel3D> & candidate_model_list)
{
    if (candidate_model_list.size() != key.size()) return std::nullopt;
    double step_norm{ 0.0 };
    for (std::size_t atom_position = 0;
        atom_position < key.size();
        atom_position++)
    {
        const auto atom_index{ key.at(atom_position) };
        const auto previous{
            EncodeLocalFittingTransformedCoordinates(
                outer_previous_state.at(atom_index).mdpde.GetModel())
        };
        const auto candidate{
            EncodeLocalFittingTransformedCoordinates(
                candidate_model_list.at(atom_position))
        };
        if (!previous.has_value() || !candidate.has_value())
        {
            return std::nullopt;
        }
        for (std::size_t parameter_index = 0;
            parameter_index < kTransformedChangeSize;
            parameter_index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
            step_norm = std::max(
                step_norm,
                std::abs((*candidate)(eigen_index) - (*previous)(eigen_index)) /
                    kLocalFittingTrustRegionParameterScale.at(parameter_index));
        }
    }
    return std::isfinite(step_norm) ?
        std::optional<double>{ step_norm } : std::nullopt;
}

} // namespace rhbm_gem::core::detail
