#pragma once

#include "core/detail/joint_component/SnapshotViews.hpp"

namespace second_stage_test::joint_domain {

inline std::vector<std::vector<rhbm_gem::core::joint_component::Support>> CopySupport(
    const rhbm_gem::core::joint_component::Domain & domain)
{
    std::vector<std::vector<rhbm_gem::core::joint_component::Support>> support(domain.atoms.size());
    for (std::size_t atom = 0; atom < domain.atoms.size(); ++atom)
        for (const auto entry : domain.atoms[atom])
            support[atom].push_back(entry);
    return support;
}

}
