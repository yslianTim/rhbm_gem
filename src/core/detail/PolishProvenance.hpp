#pragma once

#include <algorithm>
#include <vector>

namespace rhbm_gem::core::detail {

using PolishProvenance = std::vector<char>;

inline bool UsesPolish(const PolishProvenance & provenance)
{
    return std::any_of(
        provenance.begin(),
        provenance.end(),
        [](char is_polished)
        {
            return is_polished != 0;
        });
}

} // namespace rhbm_gem::core::detail
