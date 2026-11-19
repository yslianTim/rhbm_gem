#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

namespace rhbm_gem::core::classifier {

inline constexpr std::size_t kOutputDimensionCount{ 2 };

struct UmapFeatureMatrix
{
    std::vector<double> values;
    std::size_t column_count{};
};

struct UmapBin
{
    std::vector<std::size_t> points;
    std::size_t number{};
};

struct UmapBinning
{
    std::map<std::array<std::int64_t, kOutputDimensionCount>, UmapBin> bins;
    std::vector<std::size_t> point_bin_numbers;
};

} // namespace rhbm_gem::core::classifier
