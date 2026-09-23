#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <optional>
#include <span>

namespace rhbm_gem::core::detail {
// Null selects recorded targets; problems without selection metadata expose all atoms.
inline std::vector<bool> JointOutputMask(const JointProblemInput & input,
    std::optional<std::span<const std::size_t>> requested)
{
    if (!requested && input.selection_domain) requested = input.selection_domain->target_indices;
    std::vector<bool> mask(input.atom_ids.size(), !requested);
    if (requested) for (auto atom : *requested) mask.at(atom) = true;
    return mask;
}
}
