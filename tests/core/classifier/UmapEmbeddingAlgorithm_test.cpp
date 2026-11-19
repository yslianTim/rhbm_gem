#include <gtest/gtest.h>

#include "core/classifier/UmapEmbeddingAlgorithm.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::core::classifier {

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
