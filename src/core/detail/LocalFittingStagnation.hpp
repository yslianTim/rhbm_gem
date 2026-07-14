#pragma once

#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>

namespace rhbm_gem::core::detail {

inline bool IsLocalFittingForcedFixedPointObjectiveImproved(
    const std::optional<double> & candidate,
    const std::optional<double> & previous,
    double relative_tolerance)
{
    return candidate.has_value() && previous.has_value() &&
        IsBetterLocalFittingAuditObjective(
            *candidate,
            *previous,
            relative_tolerance);
}

inline bool AreAllLocalFittingActiveClustersStalled(
    const std::vector<algorithm::ClusterKey> & active_key_list,
    const std::vector<algorithm::ClusterKey> & stalled_key_list)
{
    return !active_key_list.empty() && active_key_list == stalled_key_list;
}

struct LocalFittingStagnationOptions
{
    std::size_t mismatch_iteration_limit{ 3 };
    double percentile_tolerance{ 1.0e-4 };
    double maximum_tolerance{ 1.0e-3 };
};

class LocalFittingStagnationTracker
{
    struct ClusterState
    {
        std::size_t mismatch_iteration_count{ 0 };
        bool force_fixed_point{ false };
    };

    LocalFittingStagnationOptions m_options{};
    std::map<algorithm::ClusterKey, ClusterState> m_state_by_key{};

    static void ValidateOptions(const LocalFittingStagnationOptions & options)
    {
        if (options.mismatch_iteration_limit == 0 ||
            !std::isfinite(options.percentile_tolerance) ||
            options.percentile_tolerance <= 0.0 ||
            !std::isfinite(options.maximum_tolerance) ||
            options.maximum_tolerance <= 0.0)
        {
            throw std::invalid_argument("Local fitting stagnation options are invalid.");
        }
    }

public:
    explicit LocalFittingStagnationTracker(
        LocalFittingStagnationOptions options = {})
        : m_options{ options }
    {
        ValidateOptions(m_options);
    }

    void Reconcile(const std::vector<algorithm::ClusterKey> & key_list)
    {
        std::map<algorithm::ClusterKey, ClusterState> next_state_by_key;
        for (const auto & key : key_list)
        {
            const auto iter{ m_state_by_key.find(key) };
            next_state_by_key.emplace(
                key,
                iter == m_state_by_key.end() ? ClusterState{} : iter->second);
        }
        m_state_by_key = std::move(next_state_by_key);
    }

    void Update(
        const algorithm::ClusterKey & key,
        bool eligible,
        const algorithm::ParameterChangeStats & accepted_stats,
        const std::vector<double> & accepted_maximum_list,
        const algorithm::ParameterChangeStats & raw_stats,
        const std::vector<double> & raw_maximum_list)
    {
        const auto iter{ m_state_by_key.find(key) };
        if (iter == m_state_by_key.end())
        {
            throw std::invalid_argument("Local fitting stagnation cluster state is missing.");
        }
        auto & state{ iter->second };
        if (state.force_fixed_point) return;
        if (!eligible)
        {
            state.mismatch_iteration_count = 0;
            return;
        }

        const auto accepted_converged{
            IsLocalFittingTransformedChangeConverged(
                accepted_stats,
                accepted_maximum_list,
                m_options.percentile_tolerance,
                m_options.maximum_tolerance)
        };
        const auto raw_converged{
            IsLocalFittingTransformedChangeConverged(
                raw_stats,
                raw_maximum_list,
                m_options.percentile_tolerance,
                m_options.maximum_tolerance)
        };
        if (!accepted_converged || raw_converged)
        {
            state.mismatch_iteration_count = 0;
            return;
        }

        state.mismatch_iteration_count++;
        if (state.mismatch_iteration_count >= m_options.mismatch_iteration_limit)
        {
            state.force_fixed_point = true;
        }
    }

    void MarkRecovered(const std::vector<algorithm::ClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            const auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end()) continue;
            iter->second = ClusterState{};
        }
    }

    void MarkIneligible(const std::vector<algorithm::ClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            const auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end() || iter->second.force_fixed_point) continue;
            iter->second.mismatch_iteration_count = 0;
        }
    }

    std::vector<algorithm::ClusterKey> BuildForcedFixedPointKeyList() const
    {
        std::vector<algorithm::ClusterKey> key_list;
        for (const auto & [key, state] : m_state_by_key)
        {
            if (state.force_fixed_point) key_list.emplace_back(key);
        }
        return key_list;
    }
};

} // namespace rhbm_gem::core::detail
