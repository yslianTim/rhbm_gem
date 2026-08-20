#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "core/detail/SecondStageContext.hpp"
#include "core/detail/TransformedGaussianModel.hpp"

namespace rhbm_gem::core::detail {

constexpr std::array<double, kTransformedChangeSize> kTrustRegionParameterScale{ 0.50, 0.35, 1.0 };
constexpr double kTrustRegionBoundaryTolerance{ 1.0e-12 };
constexpr double kTrustRegionGrowthBoundaryRatio{ 0.8 };

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
    std::vector<ClusterKey> changed_key_list{};
    std::vector<ClusterKey> saturated_key_list{};
};

struct RejectedClusterPartition
{
    std::vector<ClusterKey> exhausted_key_list{};
    std::vector<ClusterKey> retryable_key_list{};
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

inline RejectedClusterPartition PartitionRejectedClusters(
    const std::vector<ClusterKey> & rejected_key_list,
    const std::vector<ClusterKey> & exhausted_key_list)
{
    RejectedClusterPartition partition;
    for (const auto & key : rejected_key_list)
    {
        if (std::ranges::find(exhausted_key_list, key) != exhausted_key_list.end())
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
            [&](const ClusterKey & key)
            {
                return std::ranges::find(radius_update.saturated_key_list, key) !=
                    radius_update.saturated_key_list.end();
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
    TrustRegionOptions m_options{};
    std::map<ClusterKey, double> m_radius_by_key{};

public:
    explicit TrustRegionStateSet(TrustRegionOptions options = {})
        : m_options{ options }
    {
        if (!std::isfinite(m_options.initial_radius) ||
            !std::isfinite(m_options.minimum_radius) ||
            !std::isfinite(m_options.maximum_radius) ||
            m_options.minimum_radius <= 0.0 ||
            m_options.initial_radius < m_options.minimum_radius ||
            m_options.maximum_radius < m_options.initial_radius ||
            !std::isfinite(m_options.shrink_factor) ||
            m_options.shrink_factor <= 0.0 ||
            m_options.shrink_factor >= 1.0 ||
            !std::isfinite(m_options.growth_factor) ||
            m_options.growth_factor <= 1.0)
        {
            throw std::invalid_argument("Local fitting trust-region options are invalid.");
        }
    }

    void Reconcile(const std::vector<ClusterKey> & key_list)
    {
        std::map<ClusterKey, double> next_radius_by_key;
        for (const auto & key : key_list)
        {
            const auto iter{ m_radius_by_key.find(key) };
            next_radius_by_key.emplace(
                key, iter == m_radius_by_key.end() ? m_options.initial_radius : iter->second);
        }
        m_radius_by_key = std::move(next_radius_by_key);
    }

    double GetRadius(const ClusterKey & key) const
    {
        const auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument("Local fitting trust-region state is missing.");
        }
        return iter->second;
    }

    TrustRegionRadiusUpdate Shrink(const std::vector<ClusterKey> & key_list)
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

    void Grow(const std::vector<ClusterKey> & key_list)
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
        const std::vector<ClusterKey> & grow_key_list,
        const std::vector<ClusterKey> & rejected_key_list,
        const std::vector<ClusterKey> & backtracking_exhausted_key_list)
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

inline TrustRegionDamping LimitTrustRegionDamping(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
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
    for (std::size_t i = 0; i < previous_estimation_list.size(); i++)
    {
        if (!previous_estimation_list.at(i).allFinite() ||
            !candidate_estimation_list.at(i).allFinite())
        {
            throw std::invalid_argument("Local fitting trust-region estimation is invalid.");
        }
    }

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
                    kTrustRegionParameterScale.at(parameter_index));
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

inline std::optional<double> CalculateModelTrustRegionStepNorm(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & candidate_model_list)
{
    if (candidate_model_list.size() != previous_model_list.size())
    {
        return std::nullopt;
    }
    double step_norm{ 0.0 };
    for (std::size_t atom_position = 0; atom_position < previous_model_list.size(); atom_position++)
    {
        const auto previous{
            EncodeTransformedCoordinates(previous_model_list.at(atom_position))
        };
        const auto candidate{
            EncodeTransformedCoordinates(candidate_model_list.at(atom_position))
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
