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

#include "core/detail/FitStateView.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"

namespace rhbm_gem::core::detail {

constexpr std::array<double, kTransformedChangeSize> kTrustRegionParameterScale{ 0.50, 0.35, 1.0 };
constexpr double kTrustRegionBoundaryTolerance{ 1.0e-12 };
constexpr double kTrustRegionGrowthBoundaryRatio{ 0.8 };

using TrustRegionClusterKey = std::vector<std::size_t>;

struct TrustRegionOptions
{
    double initial_radius{ 1.0 };
    double minimum_radius{ 0.0625 };
    double maximum_radius{ 4.0 };
    double shrink_factor{ 0.5 };
    double growth_factor{ 2.0 };
};

struct TrustRegionRadiusUpdate
{
    std::vector<TrustRegionClusterKey> changed_key_list{};
    std::vector<TrustRegionClusterKey> saturated_key_list{};
};

struct RejectedClusterPartition
{
    std::vector<TrustRegionClusterKey> exhausted_key_list{};
    std::vector<TrustRegionClusterKey> retryable_key_list{};
};

struct TrustRegionIterationUpdate
{
    RejectedClusterPartition rejected_cluster_partition{};
    TrustRegionRadiusUpdate radius_update{};
};

enum class AllRejectedResolution
{
    Retry,
    MaximumIterations,
    BacktrackingExhausted,
    MinimumRadius,
    NoRetryProgress
};

inline bool IsTrustRegionStepWithinRadius(double step_norm, double radius)
{
    return std::isfinite(step_norm) &&
        std::isfinite(radius) &&
        radius > 0.0 &&
        step_norm <= radius + kTrustRegionBoundaryTolerance;
}

inline bool IsTrustRegionStepAtGrowthBoundary(double step_norm, double radius)
{
    return std::isfinite(step_norm) &&
        std::isfinite(radius) &&
        radius > 0.0 &&
        step_norm >= kTrustRegionGrowthBoundaryRatio * radius;
}

inline bool IsTrustRegionGrowthEligible(
    double step_norm,
    double radius,
    bool objective_improved)
{
    return objective_improved && IsTrustRegionStepAtGrowthBoundary(step_norm, radius);
}

inline RejectedClusterPartition PartitionRejectedClusters(
    const std::vector<TrustRegionClusterKey> & rejected_key_list,
    const std::vector<TrustRegionClusterKey> & exhausted_key_list)
{
    RejectedClusterPartition partition;
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

inline AllRejectedResolution ResolveAllRejected(
    bool maximum_iterations_reached,
    const RejectedClusterPartition & partition,
    const TrustRegionRadiusUpdate & radius_update)
{
    if (maximum_iterations_reached)
    {
        return AllRejectedResolution::MaximumIterations;
    }
    if (partition.exhausted_key_list.empty() && partition.retryable_key_list.empty())
    {
        throw std::invalid_argument(
            "All-rejected local fitting resolution requires rejected clusters.");
    }
    if (partition.retryable_key_list.empty())
    {
        return AllRejectedResolution::BacktrackingExhausted;
    }
    if (!radius_update.changed_key_list.empty())
    {
        return AllRejectedResolution::Retry;
    }

    const auto all_retryable_saturated{
        std::all_of(
            partition.retryable_key_list.begin(),
            partition.retryable_key_list.end(),
            [&](const TrustRegionClusterKey & key)
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
        AllRejectedResolution::MinimumRadius : AllRejectedResolution::NoRetryProgress;
}

class TrustRegionStateSet
{
private:
    TrustRegionOptions m_options{};
    std::map<TrustRegionClusterKey, double> m_radius_by_key{};

    static void ValidateOptions(const TrustRegionOptions & options)
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
    explicit TrustRegionStateSet(TrustRegionOptions options = {})
        : m_options{ options }
    {
        ValidateOptions(m_options);
    }

    void Reconcile(const std::vector<TrustRegionClusterKey> & key_list)
    {
        std::map<TrustRegionClusterKey, double> next_radius_by_key;
        for (const auto & key : key_list)
        {
            const auto iter{ m_radius_by_key.find(key) };
            next_radius_by_key.emplace(
                key, iter == m_radius_by_key.end() ? m_options.initial_radius : iter->second);
        }
        m_radius_by_key = std::move(next_radius_by_key);
    }

    double GetRadius(const TrustRegionClusterKey & key) const
    {
        const auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument("Local fitting trust-region state is missing.");
        }
        return iter->second;
    }

    TrustRegionRadiusUpdate Shrink(const std::vector<TrustRegionClusterKey> & key_list)
    {
        TrustRegionRadiusUpdate update;
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
            iter->second = std::max(m_options.minimum_radius, iter->second * m_options.shrink_factor);
            update.changed_key_list.emplace_back(key);
        }
        return update;
    }

    void Grow(const std::vector<TrustRegionClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_radius_by_key.find(key) };
            if (iter == m_radius_by_key.end())
            {
                throw std::invalid_argument("Local fitting trust-region state is missing.");
            }
            iter->second = std::min(m_options.maximum_radius, iter->second * m_options.growth_factor);
        }
    }

    TrustRegionIterationUpdate UpdateAfterIteration(
        const std::vector<TrustRegionClusterKey> & grow_key_list,
        const std::vector<TrustRegionClusterKey> & rejected_key_list,
        const std::vector<TrustRegionClusterKey> & backtracking_exhausted_key_list)
    {
        Grow(grow_key_list);
        auto rejected_cluster_partition{
            PartitionRejectedClusters(
                rejected_key_list,
                backtracking_exhausted_key_list)
        };
        auto radius_update{
            Shrink(rejected_cluster_partition.retryable_key_list)
        };
        return TrustRegionIterationUpdate{
            std::move(rejected_cluster_partition),
            std::move(radius_update)
        };
    }

};

struct TrustRegionDamping
{
    double effective_damping{ 1.0 };
    double step_norm{ 0.0 };
};

inline void ValidateTrustRegionInputs(
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

inline TrustRegionDamping LimitTrustRegionDamping(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    const std::array<double, 3> & parameter_scale,
    double requested_damping,
    double radius)
{
    ValidateTrustRegionInputs(
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
    return TrustRegionDamping{
        effective_damping,
        effective_damping * undamped_step_norm
    };
}

inline TrustRegionDamping LimitTrustRegionSubstepDamping(
    const std::vector<Eigen::Vector3d> & outer_previous_estimation_list,
    const std::vector<Eigen::Vector3d> & substep_previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    const std::array<double, 3> & parameter_scale,
    double requested_damping,
    double radius)
{
    ValidateTrustRegionInputs(
        outer_previous_estimation_list,
        substep_previous_estimation_list,
        parameter_scale,
        requested_damping,
        radius);
    ValidateTrustRegionInputs(
        substep_previous_estimation_list,
        candidate_estimation_list,
        parameter_scale,
        requested_damping,
        radius);

    double maximum_damping{ requested_damping };
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
            if (std::abs(base_step) > limit + kTrustRegionBoundaryTolerance)
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
                std::abs(base_step + maximum_damping * direction) / parameter_scale.at(parameter_index));
        }
    }
    return TrustRegionDamping{
        maximum_damping,
        step_norm
    };
}

inline std::optional<double> CalculateClusterModelTrustRegionStepNorm(
    const FitState & outer_previous_state,
    const ClusterKey & key,
    const std::vector<GaussianModel3D> & candidate_model_list)
{
    if (candidate_model_list.size() != key.size()) return std::nullopt;
    double step_norm{ 0.0 };
    for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
    {
        const auto atom_index{ key.at(atom_position) };
        const auto previous{
            EncodeLocalFittingTransformedCoordinates(outer_previous_state.at(atom_index).mdpde.GetModel())
        };
        const auto candidate{
            EncodeLocalFittingTransformedCoordinates(candidate_model_list.at(atom_position))
        };
        if (!previous.has_value() || !candidate.has_value())
        {
            return std::nullopt;
        }
        for (std::size_t index = 0; index < kTransformedChangeSize; index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(index) };
            step_norm = std::max(
                step_norm,
                std::abs((*candidate)(eigen_index) - (*previous)(eigen_index)) /
                    kTrustRegionParameterScale.at(index));
        }
    }
    return std::isfinite(step_norm) ? std::optional<double>{ step_norm } : std::nullopt;
}

} // namespace rhbm_gem::core::detail
