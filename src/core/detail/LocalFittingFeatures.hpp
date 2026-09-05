#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {

class ModelObject;

} // namespace rhbm_gem

namespace rhbm_gem::core::detail {

inline constexpr std::size_t kLocalFittingIdentifierColumnCount{ 3 };
inline constexpr std::size_t kLocalFittingFeatureCount{ 22 };
inline constexpr std::size_t kLocalFittingColumnCount{
    kLocalFittingIdentifierColumnCount + kLocalFittingFeatureCount
};

inline constexpr std::array<std::string_view, kLocalFittingIdentifierColumnCount>
    kLocalFittingIdentifierNames{
        "serial id",
        "residue",
        "spot",
    };

inline constexpr std::array<std::string_view, kLocalFittingFeatureCount>
    kLocalFittingFeatureNames{
        "neighbor count for peeling",
        "neighbor count in 2A",
        "signal peeling ratio",
        "tail peeling ratio",
        "amplitude 1st",
        "amplitude 2nd",
        "amplitude 3rd",
        "width 1st",
        "width 2nd",
        "width 3rd",
        "offset 1st",
        "offset 2nd",
        "offset 3rd",
        "amplitude rank 1st",
        "amplitude rank 2nd",
        "amplitude rank 3rd",
        "width rank 1st",
        "width rank 2nd",
        "width rank 3rd",
        "offset rank 1st",
        "offset rank 2nd",
        "offset rank 3rd",
    };

inline constexpr std::array<bool, kLocalFittingFeatureCount>
    kLocalFittingFeatureIsIntegral{
        true, true,
        false, false,
        false, false, false,
        false, false, false,
        false, false, false,
        true, true, true,
        true, true, true,
        true, true, true,
    };

struct LocalFittingFeatureRow
{
    int serial_id{ 0 };
    std::string residue{};
    std::string spot{};
    std::array<double, kLocalFittingFeatureCount> features{};
};

std::string BuildLocalFittingCsvHeader();

std::vector<LocalFittingFeatureRow> BuildLocalFittingFeatureRows(
    const ModelObject & model_object,
    bool peeling_applied);

} // namespace rhbm_gem::core::detail
