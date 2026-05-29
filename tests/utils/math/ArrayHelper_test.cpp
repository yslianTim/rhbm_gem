#include <gtest/gtest.h>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

TEST(ArrayHelperTest, ComputeMin)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    EXPECT_DOUBLE_EQ(1.0, rhbm_gem::array_helper::ComputeMin(data.data(), data.size()));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMin<double>(nullptr, 0)));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMin(data.data(), 0)));
}

TEST(ArrayHelperTest, ComputeMax)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    EXPECT_DOUBLE_EQ(5.0, rhbm_gem::array_helper::ComputeMax(data.data(), data.size()));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMax<double>(nullptr, 0)));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMax(data.data(), 0)));
}

TEST(ArrayHelperTest, ComputeMinMaxNegativeValues)
{
    const std::array<double, 6> data{ -3.0, -1.0, -7.0, 2.0, 5.0, -2.0 };
    EXPECT_DOUBLE_EQ(-7.0, rhbm_gem::array_helper::ComputeMin(data.data(), data.size()));
    EXPECT_DOUBLE_EQ(5.0, rhbm_gem::array_helper::ComputeMax(data.data(), data.size()));
}

TEST(ArrayHelperTest, ComputeMean)
{
    const std::array<double, 4> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(2.5, rhbm_gem::array_helper::ComputeMean(data.data(), data.size()));
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputeMean<double>(nullptr, 0));
}

TEST(ArrayHelperTest, ComputeStandardDeviation)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    const double mean{ rhbm_gem::array_helper::ComputeMean(data.data(), data.size()) };
    const double std_dev{ rhbm_gem::array_helper::ComputeStandardDeviation(data.data(), data.size(), mean) };

    EXPECT_NEAR(1.5811388300841898, std_dev, 1e-9);

    const double single_value{ 2.0 };
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputeStandardDeviation(&single_value, 1, single_value));
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputeStandardDeviation<double>(nullptr, 0, 0.0));
}

TEST(ArrayHelperTest, ComputePercentileBounds)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputePercentile(data, 0.0));
    EXPECT_DOUBLE_EQ(4.0, rhbm_gem::array_helper::ComputePercentile(data, 1.0));
}

TEST(ArrayHelperTest, ComputePercentileInterpolation)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(2.5, rhbm_gem::array_helper::ComputePercentile(data, 0.5));
}

TEST(ArrayHelperTest, ComputePercentileEdgeCases)
{
    // Empty data returns 0 regardless of percentile within (0, 1)
    std::vector<double> empty_data{};
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputePercentile(empty_data, 0.5));

    // Percentile less than 0 returns 0, greater than 1 returns max value
    std::vector<double> data{ 1.0, 2.0, 3.0 };
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputePercentile(data, -0.1));
    EXPECT_DOUBLE_EQ(3.0, rhbm_gem::array_helper::ComputePercentile(data, 1.5));

    // Unsorted data uses interpolation for percentile
    std::vector<double> unsorted{ 3.0, 1.0, 4.0, 2.0 };
    EXPECT_DOUBLE_EQ(1.75, rhbm_gem::array_helper::ComputePercentile(unsorted, 0.25));
}

TEST(ArrayHelperTest, ComputePercentileSingleElement)
{
    std::vector<double> data{ 42.0 };
    for (double p : {0.1, 0.5, 0.9})
    {
        EXPECT_DOUBLE_EQ(42.0, rhbm_gem::array_helper::ComputePercentile(data, p));
    }
}

TEST(ArrayHelperTest, ComputePercentileWithDuplicates)
{
    std::vector<double> data{ 1.0, 2.0, 2.0, 3.0 };
    EXPECT_DOUBLE_EQ(2.0, rhbm_gem::array_helper::ComputePercentile(data, 0.5));
    EXPECT_DOUBLE_EQ(1.75, rhbm_gem::array_helper::ComputePercentile(data, 0.25));
    EXPECT_DOUBLE_EQ(2.25, rhbm_gem::array_helper::ComputePercentile(data, 0.75));
}

TEST(ArrayHelperTest, ComputeMedianOddSizedDataset)
{
    std::vector<double> data{ 3.0, 1.0, 2.0 };
    EXPECT_DOUBLE_EQ(2.0, rhbm_gem::array_helper::ComputeMedian(data));
}

TEST(ArrayHelperTest, ComputeMedianEvenSizedDataset)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(2.5, rhbm_gem::array_helper::ComputeMedian(data));
}

TEST(ArrayHelperTest, ComputeMedianEmptyInput)
{
    std::vector<double> data{};
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::array_helper::ComputeMedian(data));
}

TEST(ArrayHelperTest, ComputeSmallestProportionValuesReturnsSortedSubset)
{
    const std::vector<double> data{ 4.0, 1.0, 5.0, 2.0, 3.0 };
    const std::vector<double> expected{ 1.0, 2.0 };

    const auto actual{ rhbm_gem::array_helper::ComputeSmallestProportionValues(data, 0.4) };

    EXPECT_EQ(expected, actual);
    EXPECT_EQ((std::vector<double>{ 4.0, 1.0, 5.0, 2.0, 3.0 }), data);
}

TEST(ArrayHelperTest, ComputeLargestProportionValuesReturnsSortedSubset)
{
    const std::vector<double> data{ 4.0, 1.0, 5.0, 2.0, 3.0 };
    const std::vector<double> expected{ 5.0, 4.0 };

    const auto actual{ rhbm_gem::array_helper::ComputeLargestProportionValues(data, 0.4) };

    EXPECT_EQ(expected, actual);
    EXPECT_EQ((std::vector<double>{ 4.0, 1.0, 5.0, 2.0, 3.0 }), data);
}

TEST(ArrayHelperTest, ComputeProportionValuesReturnEmptyForInvalidRatioOrInput)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0 };
    const std::vector<double> empty_data{};

    EXPECT_TRUE(rhbm_gem::array_helper::ComputeSmallestProportionValues(data, 0.0).empty());
    EXPECT_TRUE(rhbm_gem::array_helper::ComputeLargestProportionValues(data, -0.1).empty());
    EXPECT_TRUE(rhbm_gem::array_helper::ComputeSmallestProportionValues(
        data, std::numeric_limits<double>::quiet_NaN()).empty());
    EXPECT_TRUE(rhbm_gem::array_helper::ComputeLargestProportionValues(empty_data, 0.5).empty());
}

TEST(ArrayHelperTest, ComputeProportionValuesClampRatioAboveOne)
{
    const std::vector<double> data{ 4.0, 1.0, 3.0 };
    const std::vector<double> expected_smallest{ 1.0, 3.0, 4.0 };
    const std::vector<double> expected_largest{ 4.0, 3.0, 1.0 };

    EXPECT_EQ(
        expected_smallest,
        rhbm_gem::array_helper::ComputeSmallestProportionValues(data, 1.5));
    EXPECT_EQ(
        expected_largest,
        rhbm_gem::array_helper::ComputeLargestProportionValues(data, 1.5));
}

TEST(ArrayHelperTest, ComputeProportionValuesPreserveDuplicateCount)
{
    const std::vector<double> data{ 2.0, 4.0, 4.0, 1.0, 3.0 };
    const std::vector<double> expected_smallest{ 1.0, 2.0, 3.0 };
    const std::vector<double> expected_largest{ 4.0, 4.0, 3.0 };

    EXPECT_EQ(
        expected_smallest,
        rhbm_gem::array_helper::ComputeSmallestProportionValues(data, 0.6));
    EXPECT_EQ(
        expected_largest,
        rhbm_gem::array_helper::ComputeLargestProportionValues(data, 0.6));
}

TEST(ArrayHelperTest, ComputePercentileRangeTupleReturnsExpectedPercentiles)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    auto [min_value, max_value]{ rhbm_gem::array_helper::ComputePercentileRangeTuple(data, 0.2, 0.8) };
    EXPECT_NEAR(1.8, min_value, 1e-9);
    EXPECT_NEAR(4.2, max_value, 1e-9);
}

TEST(ArrayHelperTest, ComputeScalingPercentileRangeTupleExtendsSymmetrically)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingPercentileRangeTuple(data, 0.5, 0.2, 0.8) };
    EXPECT_NEAR(0.6, std::get<0>(scaled), 1e-9);
    EXPECT_NEAR(5.4, std::get<1>(scaled), 1e-9);

    auto base{ rhbm_gem::array_helper::ComputePercentileRangeTuple(data, 0.2, 0.8) };
    EXPECT_NEAR(std::get<0>(base) - std::get<0>(scaled),
                std::get<1>(scaled) - std::get<1>(base), 1e-9);
}

TEST(ArrayHelperTest, ComputeScalingPercentileRangeTupleContractsSymmetrically)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingPercentileRangeTuple(data, -0.25, 0.2, 0.8) };
    EXPECT_NEAR(2.4, std::get<0>(scaled), 1e-9);
    EXPECT_NEAR(3.6, std::get<1>(scaled), 1e-9);

    auto base{ rhbm_gem::array_helper::ComputePercentileRangeTuple(data, 0.2, 0.8) };
    EXPECT_NEAR(std::get<0>(scaled) - std::get<0>(base),
                std::get<1>(base) - std::get<1>(scaled), 1e-9);
    EXPECT_LE(std::get<0>(scaled), std::get<1>(scaled));
}

TEST(ArrayHelperTest, ComputeRangeTupleReturnsExpectedMinMax)
{
    std::vector<double> data{ 3.0, 1.0, 4.0, -1.0, 5.0 };
    auto range{ rhbm_gem::array_helper::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(-1.0, std::get<0>(range));
    EXPECT_DOUBLE_EQ(5.0, std::get<1>(range));
}

TEST(ArrayHelperTest, ComputeScalingRangeTupleExtendsRangeSymmetrically)
{
    std::vector<double> data{ 3.0, 1.0, 4.0, -1.0, 5.0 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.5) };
    EXPECT_DOUBLE_EQ(-4.0, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(8.0, std::get<1>(scaled));

    auto base{ rhbm_gem::array_helper::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(base) - std::get<0>(scaled),
                     std::get<1>(scaled) - std::get<1>(base));
}

TEST(ArrayHelperTest, ComputeScalingRangeTupleContractsRangeSymmetrically)
{
    std::vector<double> data{ 3.0, 1.0, 4.0, -1.0, 5.0 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, -0.25) };
    EXPECT_DOUBLE_EQ(0.5, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(3.5, std::get<1>(scaled));

    auto base{ rhbm_gem::array_helper::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(scaled) - std::get<0>(base),
                     std::get<1>(base) - std::get<1>(scaled));
    EXPECT_LE(std::get<0>(scaled), std::get<1>(scaled));
}

TEST(ArrayHelperTest, ComputeScalingRangeTupleKeepsRangeAboveMinRange)
{
    std::vector<double> data{ 1.0, 3.0 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.0, 0.1) };
    EXPECT_DOUBLE_EQ(1.0, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(3.0, std::get<1>(scaled));
}

TEST(ArrayHelperTest, ComputeScalingRangeTupleExpandsSmallRangeAroundCenter)
{
    std::vector<double> data{ 1.0, 1.05 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.0, 0.1) };
    EXPECT_DOUBLE_EQ(0.925, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(1.125, std::get<1>(scaled));
}

TEST(ArrayHelperTest, ComputeScalingRangeTupleDefaultMinRangeExpandsFlatData)
{
    std::vector<double> data{ 2.0, 2.0 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.0) };
    EXPECT_DOUBLE_EQ(1.9, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(2.1, std::get<1>(scaled));
}

TEST(ArrayHelperTest, ComputeScalingRangeTupleCustomMinRange)
{
    std::vector<double> data{ 2.0, 2.02 };
    auto scaled{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.0, 0.25) };
    EXPECT_DOUBLE_EQ(1.76, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(2.26, std::get<1>(scaled));
}

TEST(ArrayHelperTest, ScalingRangeDefaultAndEmpty)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };

    auto scaled_percentile{ rhbm_gem::array_helper::ComputeScalingPercentileRangeTuple(data, 0.0) };
    auto base_percentile{ rhbm_gem::array_helper::ComputePercentileRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(base_percentile), std::get<0>(scaled_percentile));
    EXPECT_DOUBLE_EQ(std::get<1>(base_percentile), std::get<1>(scaled_percentile));

    auto scaled_range{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.0) };
    auto base_range{ rhbm_gem::array_helper::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(base_range), std::get<0>(scaled_range));
    EXPECT_DOUBLE_EQ(std::get<1>(base_range), std::get<1>(scaled_range));

    std::vector<double> empty{};
    auto empty_percentile{ rhbm_gem::array_helper::ComputeScalingPercentileRangeTuple(empty, 0.0) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(empty_percentile));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(empty_percentile));

    auto empty_range{ rhbm_gem::array_helper::ComputeScalingRangeTuple(empty, 0.0) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(empty_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(empty_range));
}

TEST(ArrayHelperTest, RangeFunctionsEmptyInput)
{
    std::vector<double> data{};

    auto percentile_range{ rhbm_gem::array_helper::ComputePercentileRangeTuple(data) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(percentile_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(percentile_range));

    auto scaling_percentile_range{
        rhbm_gem::array_helper::ComputeScalingPercentileRangeTuple(data, 0.5)
    };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(scaling_percentile_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(scaling_percentile_range));

    auto range{ rhbm_gem::array_helper::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(range));

    auto scaling_range{ rhbm_gem::array_helper::ComputeScalingRangeTuple(data, 0.5) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(scaling_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(scaling_range));
}

TEST(ArrayHelperTest, ComputeNormEuclideanDistance)
{
    std::array<double, 3> v1{ 0.0, 0.0, 0.0 };
    std::array<double, 3> v2{ 3.0, 4.0, 0.0 };
    auto distance{ rhbm_gem::array_helper::ComputeNorm(v1, v2) };
    EXPECT_DOUBLE_EQ(5.0, distance);
}

TEST(ArrayHelperTest, ComputeNormWithNegativeComponents)
{
    std::array<double, 3> v1{ -1.0, -2.0, -3.0 };
    std::array<double, 3> v2{ -4.0, -6.0, -3.0 };
    auto distance{ rhbm_gem::array_helper::ComputeNorm(v1, v2) };
    EXPECT_DOUBLE_EQ(5.0, distance);
}

TEST(ArrayHelperTest, ComputeNormZeroDistanceForIdenticalVectors)
{
    std::array<double, 3> v1{ 1.0, 2.0, 3.0 };
    auto distance{ rhbm_gem::array_helper::ComputeNorm(v1, v1) };
    EXPECT_DOUBLE_EQ(0.0, distance);
}

TEST(ArrayHelperTest, ComputeNormSupportsTwoDimensionalVector)
{
    std::array<double, 2> vector{ 3.0, 4.0 };
    auto norm{ rhbm_gem::array_helper::ComputeNorm(vector) };
    EXPECT_DOUBLE_EQ(5.0, norm);
}

TEST(ArrayHelperTest, ComputeNormSupportsFourDimensionalDistance)
{
    std::array<double, 4> v1{ 0.0, 0.0, 0.0, 0.0 };
    std::array<double, 4> v2{ 1.0, 2.0, 2.0, 0.0 };
    auto distance{ rhbm_gem::array_helper::ComputeNorm(v1, v2) };
    EXPECT_DOUBLE_EQ(3.0, distance);
}

TEST(ArrayHelperTest, ComputeVectorReturnsP1ToP2Vector)
{
    std::array<double, 3> p1{ 1.0, 2.0, 3.0 };
    std::array<double, 3> p2{ 4.0, 0.0, 8.0 };
    auto vector{ rhbm_gem::array_helper::ComputeVector(p1, p2, false) };

    EXPECT_DOUBLE_EQ(3.0, vector[0]);
    EXPECT_DOUBLE_EQ(-2.0, vector[1]);
    EXPECT_DOUBLE_EQ(5.0, vector[2]);
}

TEST(ArrayHelperTest, ComputeVectorNormalizesWhenRequested)
{
    std::array<double, 3> p1{ 0.0, 0.0, 0.0 };
    std::array<double, 3> p2{ 3.0, 4.0, 0.0 };
    auto vector{ rhbm_gem::array_helper::ComputeVector(p1, p2, true) };

    EXPECT_DOUBLE_EQ(0.6, vector[0]);
    EXPECT_DOUBLE_EQ(0.8, vector[1]);
    EXPECT_DOUBLE_EQ(0.0, vector[2]);
}

TEST(ArrayHelperTest, ComputeVectorReturnsZeroVectorWhenNormalizingZeroVector)
{
    std::array<double, 3> p1{ 1.0, 2.0, 3.0 };
    auto vector{ rhbm_gem::array_helper::ComputeVector(p1, p1, true) };

    EXPECT_DOUBLE_EQ(0.0, vector[0]);
    EXPECT_DOUBLE_EQ(0.0, vector[1]);
    EXPECT_DOUBLE_EQ(0.0, vector[2]);
}

TEST(ArrayHelperTest, ComputeVectorSupportsFloatCoordinates)
{
    std::array<float, 3> p1{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> p2{ 1.0f, 4.0f, 5.0f };
    auto vector{ rhbm_gem::array_helper::ComputeVector(p1, p2, true) };

    EXPECT_NEAR(0.0f, vector[0], 1e-6f);
    EXPECT_NEAR(0.6f, vector[1], 1e-6f);
    EXPECT_NEAR(0.8f, vector[2], 1e-6f);
}

TEST(ArrayHelperTest, ComputeVectorSupportsTwoDimensionalCoordinates)
{
    std::array<double, 2> p1{ 1.0, 2.0 };
    std::array<double, 2> p2{ 4.0, -1.0 };
    auto vector{ rhbm_gem::array_helper::ComputeVector(p1, p2, false) };

    EXPECT_DOUBLE_EQ(3.0, vector[0]);
    EXPECT_DOUBLE_EQ(-3.0, vector[1]);
}

TEST(ArrayHelperTest, ComputeVectorNormalizesFourDimensionalCoordinates)
{
    std::array<double, 4> p1{ 0.0, 0.0, 0.0, 0.0 };
    std::array<double, 4> p2{ 1.0, 2.0, 2.0, 0.0 };
    auto vector{ rhbm_gem::array_helper::ComputeVector(p1, p2, true) };

    EXPECT_DOUBLE_EQ(1.0 / 3.0, vector[0]);
    EXPECT_DOUBLE_EQ(2.0 / 3.0, vector[1]);
    EXPECT_DOUBLE_EQ(2.0 / 3.0, vector[2]);
    EXPECT_DOUBLE_EQ(0.0, vector[3]);
}

TEST(ArrayHelperTest, ComputeRankReturnsExpectedPosition)
{
    std::array<double, 3> values{ 3.0, 1.0, 2.0 };
    auto rank{ rhbm_gem::array_helper::ComputeRank(values, 1) };
    EXPECT_EQ(1, rank);
}

TEST(ArrayHelperTest, ComputeRankWithNegativeValues)
{
    std::array<double, 3> values{ -3.0, 2.0, -1.0 };
    EXPECT_EQ(1, rhbm_gem::array_helper::ComputeRank(values, 0));
    EXPECT_EQ(3, rhbm_gem::array_helper::ComputeRank(values, 1));
    EXPECT_EQ(2, rhbm_gem::array_helper::ComputeRank(values, 2));
}

TEST(ArrayHelperTest, ComputeRankWithDuplicates)
{
    std::array<double, 4> values{ 2.0, 3.0, 3.0, 4.0 };
    auto rank_index1{ rhbm_gem::array_helper::ComputeRank(values, 1) };
    auto rank_index2{ rhbm_gem::array_helper::ComputeRank(values, 2) };

    EXPECT_EQ(rank_index1, rank_index2);
    EXPECT_EQ(2, rank_index1);
    EXPECT_EQ(1, rhbm_gem::array_helper::ComputeRank(values, 0));
    EXPECT_EQ(4, rhbm_gem::array_helper::ComputeRank(values, 3));
}

TEST(ArrayHelperTest, ComputeRankIndexOutOfRangeThrows)
{
    std::array<double, 3> values{ 3.0, 1.0, 2.0 };
    EXPECT_THROW(rhbm_gem::array_helper::ComputeRank(values, 3), std::out_of_range);
}
