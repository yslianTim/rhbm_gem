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
    EXPECT_TRUE(RunUmapMeanShift({}, {}, 1).empty());

    const double max{ std::numeric_limits<double>::max() };
    const std::vector<double> identical{ max, -max, max, -max, max, -max };
    const std::vector<double> isolated{ -3.0, 5.0, 4.0, -5.0, 0.0, 0.0 };
    const std::vector<LocalFittingFeatureRow> rows(3);
    EXPECT_EQ(RunUmapMeanShift(identical, rows, 2), identical);
    EXPECT_EQ(RunUmapMeanShift(isolated, rows, 4), isolated);
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
        const auto shifted{ RunUmapMeanShift(embedding, rows, 2) };
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
    EXPECT_EQ(RunUmapMeanShift(wide, two_rows, 2), wide);
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

    const auto shifted{ RunUmapMeanShift(embedding, rows, 4) };
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
    const auto original{ RunUmapMeanShift(embedding, rows, 1) };
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
        const auto shifted{ RunUmapMeanShift(transformed, rows, 1) };
        for (std::size_t i = 0; i < embedding.size(); i += 2)
        {
            EXPECT_NEAR(shifted[i], scale * original[embedding.size() - i - 2] + offset[0], 1e-11);
            EXPECT_NEAR(shifted[i + 1], scale * original[embedding.size() - i - 1] + offset[1], 1e-11);
        }
    }
}

TEST(UmapEmbeddingAlgorithmTest, BinsAtHalfUnitBoundariesAndNumbersInCoordinateOrder)
{
    const std::vector<LocalFittingFeatureRow> rows(6);
    const std::vector<double> coordinates{
        -0.5, 0.5,
        0.0, std::nextafter(0.5, 0.0),
        std::nextafter(-0.5, -1.0), -0.01,
        -0.5, 0.99,
        0.5, -0.5,
        -0.001, -0.5,
    };

    const auto binning{ BinMeanShiftEmbedding(rows, coordinates) };
    EXPECT_EQ(binning.point_bin_numbers, (std::vector<std::size_t>{ 3, 4, 1, 3, 5, 2 }));
    ASSERT_EQ(binning.bins.size(), 5u);
    const auto & shared_bin{ binning.bins.at({ -1, 1 }) };
    EXPECT_EQ(shared_bin.number, 3u);
    EXPECT_EQ(shared_bin.points, (std::vector<std::size_t>{ 0, 3 }));
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
    EXPECT_THROW(BinMeanShiftEmbedding(rows, { max, 0.0 }), std::runtime_error);
}

} // namespace rhbm_gem::core::classifier
