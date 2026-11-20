#include <gtest/gtest.h>

#include "core/classifier/UmapEmbeddingAlgorithm.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
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

TEST(UmapEmbeddingAlgorithmTest, KeepsQueryPositionForAnEmptyNeighborhoodOrUnderflowedWeights)
{
    const std::array sources{ 0.0, 0.0, 1.0, 0.0 };
    const std::array outside{ 100.0, 0.0 };
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, outside, 98.0), outside);
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, outside, 100.0), outside);
    EXPECT_EQ(RunUmapMeanShiftStepForTesting(sources, { 1.5, 0.0 }, 0.5), (std::array{ 1.0, 0.0 }));
}

TEST(UmapEmbeddingAlgorithmTest, ValidatesCoordinateShapeAndReportsNonFiniteMath)
{
    const double infinity{ std::numeric_limits<double>::infinity() };
    const double nan{ std::numeric_limits<double>::quiet_NaN() };
    const double max{ std::numeric_limits<double>::max() };
    const std::vector<LocalFittingFeatureRow> row(1);
    const std::vector<LocalFittingFeatureRow> two_rows(2);
    const std::vector<LocalFittingFeatureRow> four_rows(4);
    EXPECT_THROW(RunUmapMeanShift({ 0.0 }, row, 1), std::invalid_argument);
    EXPECT_THROW(RunUmapMeanShift({ 0.0, 1.0 }, two_rows, 1), std::invalid_argument);
    EXPECT_THROW(RunUmapMeanShift({ infinity, 0.0 }, row, 1), std::runtime_error);
    EXPECT_THROW(RunUmapMeanShift({ 0.0, nan }, row, 1), std::runtime_error);
    EXPECT_THROW(RunUmapMeanShift({ -max, 0.0, max, 0.0 }, two_rows, 1), std::runtime_error);

    // The unscaled diagonal overflows, while the radius is finite.
    const std::vector<double> wide{ 0.0, 0.0, max, max };
    EXPECT_EQ(RunUmapMeanShift(wide, two_rows, 2), wide);

    // Each neighbor is finite, but adding three equally weighted coordinates overflows.
    const double half{ max / 2.0 };
    EXPECT_THROW(RunUmapMeanShift({ half, 0.0, half, 0.0, half, 0.0, 0.0, 0.0 }, four_rows, 4),
        std::runtime_error);
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
