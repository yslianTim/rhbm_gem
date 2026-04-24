#include <gtest/gtest.h>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <stdexcept>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

TEST(ArrayHelperTest, ComputeMin)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    EXPECT_DOUBLE_EQ(1.0, rhbm_gem::array_helper::ComputeMin(data.data(), data.size()));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMin<double>(nullptr, 0)));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMin(data.data(), 0)));
}

TEST(ArrayHelperTest, ComputeMinWithThreadSize)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    auto min_single{ rhbm_gem::array_helper::ComputeMin(data.data(), data.size(), 1) };
    auto min_multi{ rhbm_gem::array_helper::ComputeMin(data.data(), data.size(), 4) };
    EXPECT_DOUBLE_EQ(min_single, min_multi);
}

TEST(ArrayHelperTest, ComputeMax)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    EXPECT_DOUBLE_EQ(5.0, rhbm_gem::array_helper::ComputeMax(data.data(), data.size()));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMax<double>(nullptr, 0)));
    EXPECT_TRUE(std::isnan(rhbm_gem::array_helper::ComputeMax(data.data(), 0)));
}

TEST(ArrayHelperTest, ComputeMaxWithThreadSize)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    auto max_single{ rhbm_gem::array_helper::ComputeMax(data.data(), data.size(), 1) };
    auto max_multi{ rhbm_gem::array_helper::ComputeMax(data.data(), data.size(), 4) };
    EXPECT_DOUBLE_EQ(max_single, max_multi);
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

TEST(ArrayHelperTest, ComputeMeanWithThreadSize)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    auto mean_single{ rhbm_gem::array_helper::ComputeMean(data.data(), data.size(), 1) };
    auto mean_multi{ rhbm_gem::array_helper::ComputeMean(data.data(), data.size(), 4) };
    EXPECT_DOUBLE_EQ(mean_single, mean_multi);
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

TEST(ArrayHelperTest, ComputeStandardDeviationWithThreadSize)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    auto mean_single{ rhbm_gem::array_helper::ComputeMean(data.data(), data.size(), 1) };
    auto mean_multi{ rhbm_gem::array_helper::ComputeMean(data.data(), data.size(), 4) };
    auto std_single{ rhbm_gem::array_helper::ComputeStandardDeviation(
        data.data(), data.size(), mean_single, 1)
    };
    auto std_multi{ rhbm_gem::array_helper::ComputeStandardDeviation(
        data.data(), data.size(), mean_multi, 4)
    };
    EXPECT_DOUBLE_EQ(std_single, std_multi);
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

TEST(ArrayHelperTest, ComputeRangeTupleWithThreadSize)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    auto range_single{ rhbm_gem::array_helper::ComputeRangeTuple(data, 1) };
    auto range_multi{ rhbm_gem::array_helper::ComputeRangeTuple(data, 4) };
    EXPECT_DOUBLE_EQ(std::get<0>(range_single), std::get<0>(range_multi));
    EXPECT_DOUBLE_EQ(std::get<1>(range_single), std::get<1>(range_multi));
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
