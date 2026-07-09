#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/algorithm/AndersonAcceleration.hpp>

namespace rhbm_gem::algorithm {

using ClusterKey = std::vector<std::size_t>;

struct ClusteredAndersonCandidate
{
    std::vector<Eigen::VectorXd> state_list{};
    std::vector<ClusterKey> used_cluster_key_list{};
};

class ClusteredAndersonAccelerationHistorySet
{
    struct ClusterState
    {
        AndersonAccelerationHistory history{};
        bool suppress_anderson{ false };
    };

    AndersonAccelerationOptions m_options{};
    std::map<ClusterKey, ClusterState> m_state_by_key{};

    static bool ContainsAtom(const ClusterKey & key, std::size_t atom_index)
    {
        return std::binary_search(key.begin(), key.end(), atom_index);
    }

public:
    explicit ClusteredAndersonAccelerationHistorySet(AndersonAccelerationOptions options)
        : m_options{ options }
    {
    }

    void Reconcile(const std::vector<ClusterKey> & key_list)
    {
        std::map<ClusterKey, ClusterState> next_state_by_key;
        for (const auto & key : key_list)
        {
            const auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end())
            {
                next_state_by_key.emplace(
                    key,
                    ClusterState{
                        AndersonAccelerationHistory{ m_options },
                        false
                    });
                continue;
            }
            next_state_by_key.emplace(key, std::move(iter->second));
        }
        m_state_by_key = std::move(next_state_by_key);
    }

    std::optional<ClusteredAndersonCandidate> BuildCandidate(
        const std::vector<ClusterKey> & key_list,
        const std::vector<Eigen::VectorXd> & previous_state_list,
        const std::vector<Eigen::VectorXd> & raw_state_list)
    {
        ClusteredAndersonCandidate candidate;
        candidate.state_list = raw_state_list;
        for (const auto & key : key_list)
        {
            auto state_iter{ m_state_by_key.find(key) };
            if (state_iter == m_state_by_key.end() ||
                state_iter->second.suppress_anderson)
            {
                continue;
            }

            auto & history{ state_iter->second.history };
            if (!history.HasCompatibleActiveIndexList(key))
            {
                history.Clear();
                continue;
            }

            const auto cluster_candidate{
                history.BuildCandidate(
                    key,
                    previous_state_list,
                    raw_state_list)
            };
            if (!cluster_candidate.has_value()) continue;

            bool valid_cluster_candidate_structure{ true };
            for (const auto active_index : key)
            {
                if (active_index >= cluster_candidate->size())
                {
                    valid_cluster_candidate_structure = false;
                    break;
                }
            }
            if (!valid_cluster_candidate_structure) continue;

            for (const auto active_index : key)
            {
                candidate.state_list.at(active_index) = cluster_candidate->at(active_index);
            }
            candidate.used_cluster_key_list.emplace_back(key);
        }

        if (candidate.used_cluster_key_list.empty())
        {
            return std::nullopt;
        }
        return candidate;
    }

    void ClearAndSuppress(const std::vector<ClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end()) continue;

            iter->second.history.Clear();
            iter->second.suppress_anderson = true;
        }
    }

    void ClearAndSuppressContaining(const std::vector<std::size_t> & atom_index_list)
    {
        for (auto & entry : m_state_by_key)
        {
            const auto & key{ entry.first };
            auto & state{ entry.second };
            const auto has_affected_atom{
                std::any_of(
                    atom_index_list.begin(),
                    atom_index_list.end(),
                    [&](std::size_t atom_index)
                    {
                        return ContainsAtom(key, atom_index);
                    })
            };
            if (!has_affected_atom) continue;

            state.history.Clear();
            state.suppress_anderson = true;
        }
    }

    void ReleaseSuppression(const std::vector<ClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end()) continue;

            iter->second.suppress_anderson = false;
        }
    }

    void Commit(
        const std::vector<ClusterKey> & key_list,
        const std::vector<Eigen::VectorXd> & input_list,
        const std::vector<Eigen::VectorXd> & output_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end() || iter->second.suppress_anderson) continue;

            iter->second.history.Commit(key, input_list, output_list);
        }
    }
};

} // namespace rhbm_gem::algorithm
