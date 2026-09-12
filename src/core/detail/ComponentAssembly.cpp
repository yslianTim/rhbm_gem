#include "core/detail/ComponentAssembly.hpp"

namespace rhbm_gem::core::detail {

void ApplyComponentPatches(
    FitState & state,
    std::span<const FitStatePatch * const> patches,
    std::optional<std::size_t> excluded_position)
{
    for (std::size_t position = 0; position < patches.size(); position++)
    {
        if (position != excluded_position && patches[position] != nullptr)
        {
            patches[position]->ApplyTo(state);
        }
    }
}

FitState AssembleComponentState(
    const FitState & base_state,
    std::span<const FitStatePatch * const> patches,
    std::optional<std::size_t> excluded_position)
{
    auto state{ base_state };
    ApplyComponentPatches(state, patches, excluded_position);
    return state;
}

} // namespace rhbm_gem::core::detail
