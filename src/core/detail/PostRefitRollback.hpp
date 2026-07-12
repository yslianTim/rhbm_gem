#pragma once

#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace rhbm_gem::core::detail {

using PostRefitRollbackClusterKey = std::vector<std::size_t>;

inline std::vector<std::size_t> ExpandPostRefitRollbackClusters(
    const std::vector<std::size_t> & active_index_list,
    const std::vector<PostRefitRollbackClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & seed_active_position_list,
    std::vector<char> & suspicious_mask)
{
    if (suspicious_mask.size() != active_index_list.size())
    {
        throw std::invalid_argument(
            "Post-refit rollback mask size is inconsistent.");
    }
    if (seed_active_position_list.empty()) return {};

    std::unordered_map<std::size_t, std::size_t> active_position_by_atom_index;
    active_position_by_atom_index.reserve(active_index_list.size());
    for (std::size_t active_position = 0;
        active_position < active_index_list.size();
        active_position++)
    {
        if (!active_position_by_atom_index.emplace(
                active_index_list.at(active_position),
                active_position).second)
        {
            throw std::invalid_argument(
                "Post-refit rollback active atom index is duplicated.");
        }
    }

    std::unordered_map<std::size_t, std::size_t> cluster_position_by_atom_index;
    cluster_position_by_atom_index.reserve(active_index_list.size());
    for (std::size_t cluster_position = 0;
        cluster_position < cluster_key_list.size();
        cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        if (key.empty())
        {
            throw std::invalid_argument(
                "Post-refit rollback cluster key must not be empty.");
        }
        for (const auto atom_index : key)
        {
            if (active_position_by_atom_index.find(atom_index) ==
                active_position_by_atom_index.end())
            {
                throw std::invalid_argument(
                    "Post-refit rollback cluster atom is not active.");
            }
            if (!cluster_position_by_atom_index.emplace(
                    atom_index,
                    cluster_position).second)
            {
                throw std::invalid_argument(
                    "Post-refit rollback cluster atom is duplicated.");
            }
        }
    }
    if (cluster_position_by_atom_index.size() != active_index_list.size())
    {
        throw std::invalid_argument(
            "Post-refit rollback clusters do not cover all active atoms.");
    }

    std::vector<char> affected_cluster_mask(cluster_key_list.size(), 0);
    for (const auto seed_active_position : seed_active_position_list)
    {
        if (seed_active_position >= active_index_list.size())
        {
            throw std::invalid_argument(
                "Post-refit rollback seed position is out of range.");
        }
        const auto atom_index{ active_index_list.at(seed_active_position) };
        affected_cluster_mask.at(
            cluster_position_by_atom_index.at(atom_index)) = 1;
    }

    std::vector<std::size_t> affected_active_position_list;
    for (std::size_t active_position = 0;
        active_position < active_index_list.size();
        active_position++)
    {
        const auto atom_index{ active_index_list.at(active_position) };
        const auto cluster_position{
            cluster_position_by_atom_index.at(atom_index)
        };
        if (affected_cluster_mask.at(cluster_position) == 0) continue;

        suspicious_mask.at(active_position) = 1;
        affected_active_position_list.emplace_back(active_position);
    }
    return affected_active_position_list;
}

} // namespace rhbm_gem::core::detail
