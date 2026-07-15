#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace rhbm_gem::core::detail {

using PostRefitRollbackClusterKey = std::vector<std::size_t>;

inline void ExpandPostRefitRollbackClusters(
    const std::vector<PostRefitRollbackClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & seed_atom_index_list,
    std::vector<char> & suspicious_mask)
{
    if (seed_atom_index_list.empty()) return;

    for (const auto & key : cluster_key_list)
    {
        const auto is_affected{
            std::any_of(
                key.begin(),
                key.end(),
                [&](std::size_t atom_index)
                {
                    return std::find(
                        seed_atom_index_list.begin(),
                        seed_atom_index_list.end(),
                        atom_index) != seed_atom_index_list.end();
                })
        };
        if (!is_affected) continue;
        for (const auto atom_index : key)
        {
            suspicious_mask.at(atom_index) = 1;
        }
    }
}

} // namespace rhbm_gem::core::detail
