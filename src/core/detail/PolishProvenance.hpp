#pragma once

#include <functional>
#include <ranges>
#include <vector>

namespace rhbm_gem::core::detail {

using PolishProvenance = std::vector<char>;

inline bool UsesPolish(const PolishProvenance & provenance)
{
    return std::ranges::any_of(provenance, std::identity{});
}

} // namespace rhbm_gem::core::detail
