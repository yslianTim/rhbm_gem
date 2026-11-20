#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <vector>

namespace rhbm_gem::core::classifier {

inline constexpr std::size_t kOutputDimensionCount{ 2 };
inline constexpr double kUmapBinWidthRatio{ 0.5 };

struct UmapFeatureMatrix
{
    std::vector<double> values;
    std::size_t column_count{};
};

struct UmapBin
{
    std::vector<std::size_t> points;
    std::size_t number{};
    std::array<std::optional<double>, 3> joint_medians{}; // Amplitude, width, offset.
};

struct UmapSpatialScale
{
    std::array<double, kOutputDimensionCount> origin{};
    double diagonal{};
    double radius{};
    double scaled_diagonal{};
};

struct UmapMeanShiftResult
{
    std::vector<double> coordinates;
    UmapSpatialScale scale;
};

struct UmapBinning
{
    std::map<std::array<std::int64_t, kOutputDimensionCount>, UmapBin> bins;
    std::vector<std::size_t> point_bin_numbers;
    UmapSpatialScale scale;
};

struct UmapEmbeddingResult
{
    std::vector<double> coordinates;
    std::vector<double> mean_shift_coordinates;
    UmapBinning binning;
};

} // namespace rhbm_gem::core::classifier
