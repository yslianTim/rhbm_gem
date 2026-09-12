#pragma once

#include "core/detail/SecondStageFitting.hpp"

#include <optional>
#include <span>

namespace rhbm_gem::core::detail {

// Null entries are unselected components. Order and exclusion positions belong
// to the caller; assembly never changes component membership or acceptance.
void ApplyComponentPatches(
    FitState & state,
    std::span<const FitStatePatch * const> patches,
    std::optional<std::size_t> excluded_position = std::nullopt);

FitState AssembleComponentState(
    const FitState & base_state,
    std::span<const FitStatePatch * const> patches,
    std::optional<std::size_t> excluded_position = std::nullopt);

// Policies own removal ordering and audit semantics. Applying a removal returns
// its audit, allowing a policy to reuse an already evaluated trial objective.
template<typename Evaluate, typename Passed, typename SelectRemoval, typename Remove>
auto AuditAndSalvageComponents(
    Evaluate evaluate, Passed passed, SelectRemoval select_removal, Remove remove)
{
    auto objective{ evaluate() };
    while (!passed(objective))
    {
        auto removal{ select_removal(objective) };
        if (!removal.has_value()) break;
        objective = remove(*removal);
    }
    return objective;
}

} // namespace rhbm_gem::core::detail
