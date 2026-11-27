#pragma once

#include <array>
#include <cstddef>
#include <string_view>

namespace rhbm_gem::core::classifier {

inline constexpr std::size_t kLocalFittingFeatureCount{ 10 };

inline constexpr std::array<std::string_view, kLocalFittingFeatureCount>
    kLocalFittingFeatureNames{
        "neighbor distance sum in 2A",
        "distance to closest neighbor",
        "signal peeling ratio",
        "tail peeling ratio",
        "amplitude joint",
        "width joint",
        "offset joint",
        "joint model center height",
        "joint signal gaussian integral",
        "joint signal charge integral",
    };

} // namespace rhbm_gem::core::classifier
