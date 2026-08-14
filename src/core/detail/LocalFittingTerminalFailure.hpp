#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <ostream>
#include <variant>
#include <vector>

#include "core/detail/JointOffset.hpp"
#include "core/detail/ClusterHealth.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/TransformedChange.hpp"
#include "core/detail/SuspiciousUpdate.hpp"

namespace rhbm_gem::core::detail {

constexpr std::size_t kPersistentTerminalFailureIterationLimit{ 5 };

using PersistentSuspiciousRollbackReason = std::vector<std::size_t>;
using PersistentTerminalFailureReason = std::variant<PersistentSuspiciousRollbackReason, JointOffsetSolveStatus>;

struct PersistentTerminalFailureState
{
    PersistentTerminalFailureReason reason{};
    std::size_t stable_iteration_count{ 0 };
};

using PersistentTerminalFailureStateMap = std::map<ClusterKey, PersistentTerminalFailureState>;
using TerminalPersistentFailureMap = std::map<ClusterKey, PersistentTerminalFailureReason>;

struct LocalFittingTerminalSummary
{
    std::size_t suspicious_cluster_count{ 0 };
    std::size_t suspicious_atom_count{ 0 };
    std::size_t joint_offset_failure_cluster_count{ 0 };
    std::size_t joint_offset_failure_atom_count{ 0 };
    std::map<JointOffsetSolveStatus, std::size_t> joint_offset_failure_status_count{};

    std::size_t AtomCount() const
    {
        return suspicious_atom_count + joint_offset_failure_atom_count;
    }

    bool HasFailures() const
    {
        return AtomCount() > 0;
    }
};

inline void AppendLocalFittingTerminalSummary(
    std::ostream & stream,
    const LocalFittingTerminalSummary & summary)
{
    if (summary.suspicious_atom_count > 0)
    {
        stream << "; terminal suspicious rollback fallback clusters/atoms = "
            << summary.suspicious_cluster_count
            << "/" << summary.suspicious_atom_count;
    }
    if (summary.joint_offset_failure_atom_count > 0)
    {
        stream << "; terminal joint-offset failure fallback clusters/atoms = "
            << summary.joint_offset_failure_cluster_count
            << "/" << summary.joint_offset_failure_atom_count;
        if (!summary.joint_offset_failure_status_count.empty())
        {
            stream << ", statuses = ";
            bool is_first_status{ true };
            for (const auto & [status, count] :
                summary.joint_offset_failure_status_count)
            {
                if (!is_first_status) stream << ",";
                stream << GetJointOffsetSolveStatusText(status) << ":" << count;
                is_first_status = false;
            }
        }
    }
}

inline std::vector<ClusterKey> AccumulateTerminalFailureSummary(
    const TerminalPersistentFailureMap & terminal_failure_by_key,
    LocalFittingTerminalSummary & terminal_summary)
{
    std::vector<ClusterKey> terminal_key_list;
    terminal_key_list.reserve(terminal_failure_by_key.size());
    for (const auto & [key, reason] : terminal_failure_by_key)
    {
        terminal_key_list.emplace_back(key);
        if (std::holds_alternative<PersistentSuspiciousRollbackReason>(reason))
        {
            terminal_summary.suspicious_cluster_count++;
            terminal_summary.suspicious_atom_count += key.size();
            continue;
        }
        const auto status{ std::get<JointOffsetSolveStatus>(reason) };
        terminal_summary.joint_offset_failure_cluster_count++;
        terminal_summary.joint_offset_failure_atom_count += key.size();
        terminal_summary.joint_offset_failure_status_count[status]++;
    }
    return terminal_key_list;
}

inline TerminalPersistentFailureMap UpdatePersistentTerminalFailureState(
    const std::vector<ClusterKey> & accepted_key_list,
    const SuspiciousUpdateMask & suspicious_atom_mask,
    const ClusterHealthMap & health_by_key,
    const FitState & assembled_state,
    const FitState & previous_state,
    PersistentTerminalFailureStateMap & state_by_key)
{
    PersistentTerminalFailureStateMap next_state_by_key;
    TerminalPersistentFailureMap terminal_failure_by_key;
    for (const auto & [key, health] : health_by_key)
    {
        if (std::find(accepted_key_list.begin(), accepted_key_list.end(), key) ==
            accepted_key_list.end())
        {
            continue;
        }

        auto cluster_suspicious_atom_index_list{
            CollectSuspiciousAtomIndices(key, suspicious_atom_mask)
        };
        PersistentTerminalFailureReason reason;
        if (!cluster_suspicious_atom_index_list.empty())
        {
            reason = std::move(cluster_suspicious_atom_index_list);
        }
        else
        {
            const auto status{ health.joint_offset_status };
            if (!IsJointOffsetSolveHardFailure(status)) continue;
            reason = status;
        }

        const auto transformed_change_summary{
            SummarizeTransformedChanges(assembled_state, previous_state, key)
        };
        if (!IsTransformedPercentileConverged(transformed_change_summary))
        {
            continue;
        }

        PersistentTerminalFailureState next_state{ std::move(reason), 1 };
        const auto previous_iter{ state_by_key.find(key) };
        if (previous_iter != state_by_key.end() && previous_iter->second.reason == next_state.reason)
        {
            next_state.stable_iteration_count = previous_iter->second.stable_iteration_count + 1;
        }

        if (next_state.stable_iteration_count >= kPersistentTerminalFailureIterationLimit)
        {
            terminal_failure_by_key.emplace(key, std::move(next_state.reason));
            continue;
        }
        next_state_by_key.emplace(key, std::move(next_state));
    }
    state_by_key = std::move(next_state_by_key);
    return terminal_failure_by_key;
}

inline void ApplyTerminalFallbackClusters(
    const std::vector<ClusterKey> & terminal_key_list,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    std::vector<char> & terminal_atom_mask,
    FitState & assembled_state,
    PolishProvenance & assembled_polish_provenance)
{
    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            terminal_atom_mask.at(atom_index) = 1;
            assembled_state.at(atom_index) = previous_state.at(atom_index);
            assembled_polish_provenance.at(atom_index) = previous_polish_provenance.at(atom_index);
        }
    }
}

inline std::vector<std::size_t> BuildEligibleLocalFittingActiveIndexList(const std::vector<char> & terminal_atom_mask)
{
    const auto atom_size{ terminal_atom_mask.size() };
    std::vector<std::size_t> active_index_list;
    active_index_list.reserve(atom_size);
    for (std::size_t atom_index = 0; atom_index < atom_size; atom_index++)
    {
        if (terminal_atom_mask.at(atom_index) == 0)
        {
            active_index_list.emplace_back(atom_index);
        }
    }
    return active_index_list;
}

struct LocalFittingTerminalFailureState
{
    PersistentTerminalFailureStateMap persistent_state_by_key{};
    std::vector<char> terminal_atom_mask{};
    LocalFittingTerminalSummary terminal_summary{};

    LocalFittingTerminalFailureState() = default;

    explicit LocalFittingTerminalFailureState(std::size_t atom_count)
        : terminal_atom_mask(atom_count, 0)
    {
    }

    const LocalFittingTerminalSummary & Summary() const
    {
        return terminal_summary;
    }

    bool HasFailures() const
    {
        return terminal_summary.HasFailures();
    }

    std::size_t AtomCount() const
    {
        return terminal_summary.AtomCount();
    }

    std::vector<std::size_t> BuildEligibleActiveIndexList() const
    {
        return BuildEligibleLocalFittingActiveIndexList(terminal_atom_mask);
    }

    std::vector<ClusterKey> IsolatePersistentFailures(
        const std::vector<ClusterKey> & accepted_key_list,
        SuspiciousUpdateMask & suspicious_atom_mask,
        const ClusterHealthMap & health_by_key,
        FitState & assembled_state,
        const FitState & previous_state,
        const PolishProvenance & previous_polish_provenance,
        PolishProvenance & assembled_polish_provenance)
    {
        const auto terminal_failure_by_key{
            UpdatePersistentTerminalFailureState(
                accepted_key_list,
                suspicious_atom_mask,
                health_by_key,
                assembled_state,
                previous_state,
                persistent_state_by_key)
        };
        const auto terminal_key_list{
            AccumulateTerminalFailureSummary(
                terminal_failure_by_key,
                terminal_summary)
        };
        ApplyTerminalFallbackClusters(
            terminal_key_list,
            previous_state,
            previous_polish_provenance,
            terminal_atom_mask,
            assembled_state,
            assembled_polish_provenance);
        if (!terminal_key_list.empty())
        {
            ClearSuspiciousUpdateMaskForClusters(
                terminal_key_list,
                suspicious_atom_mask);
        }
        return terminal_key_list;
    }
};

} // namespace rhbm_gem::core::detail
