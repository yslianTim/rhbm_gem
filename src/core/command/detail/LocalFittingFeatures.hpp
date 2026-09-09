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
inline constexpr std::size_t kLocalFittingFeatureCount{ 14 };
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
        "neighbor count in 1.5A",
        "neighbor distance sum in 2A",
        "neighbor distance sum in 1.5A",
        "signal peeling ratio",
        "tail peeling ratio",
        "amplitude 2nd",
        "width 2nd",
        "offset 2nd",
        "amplitude rank 2nd",
        "width rank 2nd",
        "offset rank 2nd",
        "distance to closest neighbor",
    };

inline constexpr std::array<bool, kLocalFittingFeatureCount>
    kLocalFittingFeatureIsIntegral{
        true, true, true,
        false, false,
        false, false,
        false, false, false,
        true, true, true,
        false,
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
