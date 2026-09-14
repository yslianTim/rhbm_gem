#pragma once

#include "core/detail/second_stage/SecondStageState.hpp"

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

} // namespace rhbm_gem::core::detail
