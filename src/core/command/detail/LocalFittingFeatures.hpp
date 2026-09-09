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
inline constexpr std::size_t kLocalFittingFeatureCount{ 10 };
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
        "neighbor distance sum in 2A",
        "distance to closest neighbor",
        "signal peeling ratio",
        "tail peeling ratio",
        "amplitude 2nd",
        "width 2nd",
        "offset 2nd",
        "amplitude rank 2nd",
        "width rank 2nd",
        "offset rank 2nd",
    };

inline constexpr std::array<bool, kLocalFittingFeatureCount>
    kLocalFittingFeatureIsIntegral{
        false,
        false,
        false, false,
        false, false, false,
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
    const ModelObject & model_object);

} // namespace rhbm_gem::core::detail
