#include <gtest/gtest.h>

#include "core/classifier/UmapEmbeddingAlgorithm.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rhbm_gem::core::classifier {

TEST(UmapEmbeddingAlgorithmTest, HandlesEmptyIdenticalAndIsolatedPoints)
{
    EXPECT_TRUE(RunUmapMeanShift({}, {}, 1).coordinates.empty());

    const double max{ std::numeric_limits<double>::max() };
    const std::vector<double> identical{ max, -max, max, -max, max, -max };
    const std::vector<double> isolated{ -3.0, 5.0, 4.0, -5.0, 0.0, 0.0 };
    const std::vector<LocalFittingFeatureRow> rows(3);
    EXPECT_EQ(RunUmapMeanShift(identical, rows, 2).coordinates, identical);
    EXPECT_EQ(RunUmapMeanShift(isolated, rows, 4).coordinates, isolated);
}

TEST(UmapEmbeddingAlgorithmTest, KeepsQueryPositionForAnEmptyNeighborhoodOrZeroRadius)
{
    const std::array sources{ 0.0, 0.0, 1.0, 0.0 };
    const std::array outside{ 100.0, 0.0 };
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, outside, 98.0), outside);
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, outside, 0.0), outside);
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, { 1.0, 0.0 }, 0.0), (std::array{ 1.0, 0.0 }));
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, { 1.5, 0.0 }, 0.5), (std::array{ 1.0, 0.0 }));
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(std::array{ 0.0, 0.0 }, outside, 100.0),
        (std::array{ 0.0, 0.0 }));
}

TEST(UmapEmbeddingAlgorithmTest, UsesGaussianWeightsAtRelativeDistances)
{
    constexpr double radius{ 4.0 };
    const std::array query{ 0.0, 0.0 };
    // The source at the query has weight exp(0) = 1. Each second source tests a weight relative to it.
    const std::array cases{
        std::pair{ 0.25, std::exp(-0.25) },
        std::pair{ 0.5, std::exp(-1.0) },
        std::pair{ 1.0, std::exp(-4.0) },
    };
    for (const auto & [relative_distance, weight] : cases)
    {
        SCOPED_TRACE(relative_distance);
        const double distance{ relative_distance * radius };
        const std::array sources{ 0.0, 0.0, distance, 0.0, std::nextafter(radius, 5.0), 0.0 };
        const auto shifted{ RunUmapMeanShiftStepForTesting(sources, query, radius) };
        EXPECT_NEAR(shifted[0], distance * weight / (1.0 + weight), 1e-15);
        EXPECT_EQ(shifted[1], 0.0);
    }
}

TEST(UmapEmbeddingAlgorithmTest, HandlesSubnormalAndVeryLargePositiveRadii)
{
    for (const double scale : { 0x1p-1060, 1e300 })
    {
        SCOPED_TRACE(scale);
        const std::vector<double> embedding{ 0.0, 0.0, scale, 0.0, 20.0 * scale, 0.0 };
        const std::vector<LocalFittingFeatureRow> rows(3);
        const auto step{ RunUmapMeanShiftStepForTesting(embedding, { 0.0, 0.0 }, 2.0 * scale) };
        const auto shifted{ RunUmapMeanShift(embedding, rows, 2).coordinates };
        // Subnormal coordinates here are quantized in units of scale / 16384.
        const double tolerance{ scale < std::numeric_limits<double>::min() ? 2e-4 : 1e-12 };
        EXPECT_NEAR(step[0] / scale, 1.0 / (1.0 + std::exp(1.0)), tolerance);
        EXPECT_NEAR(shifted[0] / scale, 0.5, tolerance);
        EXPECT_NEAR(shifted[2] / scale, 0.5, tolerance);
        EXPECT_EQ(shifted[4], embedding[4]);
    }
}

TEST(UmapEmbeddingAlgorithmTest, ValidatesCoordinateShapeAndReportsNonFiniteMath)
{
    const double infinity{ std::numeric_limits<double>::infinity() };
    const double nan{ std::numeric_limits<double>::quiet_NaN() };
    const double max{ std::numeric_limits<double>::max() };
    const std::vector<LocalFittingFeatureRow> row(1);
    const std::vector<LocalFittingFeatureRow> two_rows(2);
    EXPECT_THROW(RunUmapMeanShift({ 0.0 }, row, 1), std::invalid_argument);
    EXPECT_THROW(RunUmapMeanShift({ 0.0, 1.0 }, two_rows, 1), std::invalid_argument);
    EXPECT_THROW(RunUmapMeanShift({ infinity, 0.0 }, row, 1), std::runtime_error);
    EXPECT_THROW(RunUmapMeanShift({ 0.0, nan }, row, 1), std::runtime_error);
    EXPECT_THROW(RunUmapMeanShift({ -max, 0.0, max, 0.0 }, two_rows, 1), std::runtime_error);

    // The unscaled diagonal overflows, while the radius is finite.
    const std::vector<double> wide{ 0.0, 0.0, max, max };
    EXPECT_EQ(RunUmapMeanShift(wide, two_rows, 2).coordinates, wide);
}

TEST(UmapEmbeddingAlgorithmTest, ComputesFiniteMeanWhenWeightedSumWouldOverflow)
{
    const double max{ std::numeric_limits<double>::max() };
    const double nearby{ 0.95 * max };
    const std::vector<double> embedding{ nearby, -nearby, max, -max, max, -max, 0.0, 0.0 };
    const std::vector<LocalFittingFeatureRow> rows(4);
    const auto step{ RunUmapMeanShiftStepForTesting(embedding, { nearby, -nearby }, 0.1 * max) };
    const double weight{ std::exp(-2.0) };
    const double expected{ (0.95 + 2.0 * weight) / (1.0 + 2.0 * weight) };
    EXPECT_NEAR(step[0] / max, expected, 1e-15);
    EXPECT_NEAR(step[1] / max, -expected, 1e-15);

    const auto shifted{ RunUmapMeanShift(embedding, rows, 4).coordinates };
    EXPECT_TRUE(HasFiniteCoordinates(shifted));
    EXPECT_GT(shifted[0], nearby);
    EXPECT_EQ(shifted[0], -shifted[1]);
    EXPECT_EQ(shifted[2], shifted[4]);
    EXPECT_EQ(shifted[6], 0.0);
}

TEST(UmapEmbeddingAlgorithmTest, PreservesTrajectoriesUnderUniformScalingTranslationAndRowReversal)
{
    const std::vector<double> embedding{
        -1.25, 0.5, -0.75, 1.25, 0.0, 0.75, 0.0, 0.75,
        4.0, 2.0, 4.6, 2.2, 17.0, 6.0,
    };
    const std::vector<LocalFittingFeatureRow> rows(embedding.size() / 2);
    const auto original{ RunUmapMeanShift(embedding, rows, 1).coordinates };
    EXPECT_GT(std::abs(original[0] - embedding[0]), 0.1);
    for (const double scale : { 1.0, 10.0, 0.1 })
    {
        SCOPED_TRACE(scale);
        const std::array offset{ 13.0, -7.0 };
        std::vector<double> transformed;
        for (std::size_t i = embedding.size(); i > 0; i -= 2)
        {
            transformed.push_back(scale * embedding[i - 2] + offset[0]);
            transformed.push_back(scale * embedding[i - 1] + offset[1]);
        }
        const auto shifted{ RunUmapMeanShift(transformed, rows, 1).coordinates };
        for (std::size_t i = 0; i < embedding.size(); i += 2)
        {
            EXPECT_NEAR(shifted[i], scale * original[embedding.size() - i - 2] + offset[0], 1e-11);
            EXPECT_NEAR(shifted[i + 1], scale * original[embedding.size() - i - 1] + offset[1], 1e-11);
        }
    }
}

TEST(UmapEmbeddingAlgorithmTest, BinsAtAdaptiveBoundariesAndNumbersInCoordinateOrder)
{
    // The original 3-by-4 box gives R = 0.5 and cells of width 0.25.
    const std::array original{ 0.0, -2.0, 3.0, 2.0 };
    const auto scale{ GetUmapSpatialScale(original) };
    EXPECT_EQ(scale.origin, (std::array{ 0.0, -2.0 }));
    EXPECT_DOUBLE_EQ(scale.diagonal, 5.0);
    EXPECT_DOUBLE_EQ(scale.radius, 0.5);
    const std::vector<LocalFittingFeatureRow> rows(6);
    const std::vector<double> shifted{
        0.25, -1.75,
        std::nextafter(0.25, 0.0), -1.75,
        0.0, -2.0,
        0.25, -1.75,
        std::nextafter(0.25, 1.0), -1.8,
        0.26, -1.7,
    };

    const auto binning{ BinMeanShiftEmbedding(rows, shifted, scale) };
    EXPECT_EQ(binning.point_bin_numbers, (std::vector<std::size_t>{ 4, 2, 1, 4, 3, 4 }));
    ASSERT_EQ(binning.bins.size(), 4u);
    const auto & shared_bin{ binning.bins.at({ 1, 1 }) };
    EXPECT_EQ(shared_bin.number, 4u);
    EXPECT_EQ(shared_bin.points, (std::vector<std::size_t>{ 0, 3, 5 }));
}

TEST(UmapEmbeddingAlgorithmTest, AdaptiveIndicesAreInvariantUnderUniformScalingAndTranslation)
{
    const std::vector<double> original{ -2.0, -1.0, 2.0, 2.0, 0.0, 0.0 };
    const std::vector<double> shifted{ -1.875, -0.875, 1.875, 1.875, 0.125, 0.125 };
    const std::vector<LocalFittingFeatureRow> rows(3);
    for (const double factor : { 0.25, 1.0, 16.0 })
    {
        for (const auto & translation : { std::array{ 0.0, 0.0 }, std::array{ -10.0, 12.0 } })
        {
            auto transformed_original{ original };
            auto transformed_shifted{ shifted };
            for (std::size_t i = 0; i < original.size(); ++i)
            {
                transformed_original[i] = factor * original[i] + translation[i % 2];
                transformed_shifted[i] = factor * shifted[i] + translation[i % 2];
            }
            const auto scale{ GetUmapSpatialScale(transformed_original) };
            EXPECT_DOUBLE_EQ(scale.radius, 0.5 * factor);
            const auto binning{ BinMeanShiftEmbedding(rows, transformed_shifted, scale) };
            EXPECT_EQ(binning.point_bin_numbers, (std::vector<std::size_t>{ 1, 3, 2 }));
            ASSERT_EQ(binning.bins.size(), 3u);
            EXPECT_EQ(binning.bins.at({ 0, 0 }).number, 1u);
            EXPECT_EQ(binning.bins.at({ 8, 4 }).number, 2u);
            EXPECT_EQ(binning.bins.at({ 15, 11 }).number, 3u);
        }
    }
}

TEST(UmapEmbeddingAlgorithmTest, BinsZeroSubnormalAndExtremeRanges)
{
    EXPECT_TRUE(BinMeanShiftEmbedding({}, {}, GetUmapSpatialScale({})).bins.empty());
    const double max{ std::numeric_limits<double>::max() };
    const double tiny{ std::numeric_limits<double>::denorm_min() };
    const std::vector<LocalFittingFeatureRow> rows(2);
    const std::vector<double> identical{ max, -max, max, -max };
    const auto zero{ GetUmapSpatialScale(identical) };
    EXPECT_EQ(zero.radius, 0.0);
    EXPECT_EQ(zero.diagonal, 0.0);
    const auto binning{ BinMeanShiftEmbedding(rows, identical, zero) };
    ASSERT_EQ(binning.bins.size(), 1u);
    EXPECT_EQ(binning.bins.at({ 0, 0 }).points, (std::vector<std::size_t>{ 0, 1 }));

    // The radius is either zero or the smallest positive double. In both cases h underflows.
    for (const double diagonal : { tiny, 7.0 * tiny })
    {
        const std::vector<double> original{ -diagonal, tiny, 0.0, tiny };
        const auto shifted{ RunUmapMeanShift(original, rows, 1) };
        EXPECT_EQ(shifted.scale.diagonal, diagonal);
        EXPECT_EQ(shifted.scale.radius, diagonal == tiny ? 0.0 : tiny);
        EXPECT_EQ(shifted.coordinates, original);
        const auto bins{ BinMeanShiftEmbedding(rows, shifted.coordinates, shifted.scale) };
        ASSERT_EQ(bins.bins.size(), 2u);
        EXPECT_EQ(bins.bins.at({ 20, 0 }).number, 2u);
    }

    const std::vector<double> tiny_diagonal{ 0.0, 0.0, tiny, tiny };
    const auto tiny_bins{ BinMeanShiftEmbedding(rows, tiny_diagonal, GetUmapSpatialScale(tiny_diagonal)) };
    EXPECT_EQ(tiny_bins.bins.at({ 14, 14 }).number, 2u);

    const std::vector<double> wide{ -max, 0.0, 0.0, max };
    const auto scale{ GetUmapSpatialScale(wide) };
    EXPECT_TRUE(std::isinf(scale.diagonal));
    EXPECT_NEAR(scale.radius / max, std::sqrt(2.0) / 10.0, 1e-15);
    const auto bins{ BinMeanShiftEmbedding(rows, wide, scale) };
    ASSERT_EQ(bins.bins.size(), 2u);
    EXPECT_EQ(bins.bins.at({ 14, 14 }).points, (std::vector<std::size_t>{ 1 }));
}

TEST(UmapEmbeddingAlgorithmTest, DetectsNonFiniteCoordinatesAndOutOfRangeBinIndices)
{
    const auto infinity{ std::numeric_limits<double>::infinity() };
    const auto nan{ std::numeric_limits<double>::quiet_NaN() };
    EXPECT_TRUE(HasFiniteCoordinates(std::array{ -0.5, 0.5 }));
    EXPECT_FALSE(HasFiniteCoordinates(std::array{ 0.0, infinity }));
    EXPECT_FALSE(HasFiniteCoordinates(std::array{ nan, 0.0 }));

    const std::vector<LocalFittingFeatureRow> rows(1);
    const auto max{ std::numeric_limits<double>::max() };
    const auto scale{ GetUmapSpatialScale(std::array{ 0.0, 0.0, 1.0, 1.0 }) };
    EXPECT_THROW(BinMeanShiftEmbedding(rows, { max, 0.0 }, scale), std::runtime_error);
    EXPECT_THROW(BinMeanShiftEmbedding(rows, { nan, 0.0 }, scale), std::runtime_error);
    EXPECT_THROW(BinMeanShiftEmbedding(rows, { 0.0 }, scale), std::invalid_argument);
}

} // namespace rhbm_gem::core::classifier
