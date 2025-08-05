#include <gtest/gtest.h>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <stdexcept>

#include "ArrayStats.hpp"

TEST(ArrayStatsTest, ComputeMinMax)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };

    EXPECT_DOUBLE_EQ(1.0, ArrayStats<double>::ComputeMin(data.data(), data.size()));
    EXPECT_DOUBLE_EQ(5.0, ArrayStats<double>::ComputeMax(data.data(), data.size()));

    EXPECT_TRUE(std::isnan(ArrayStats<double>::ComputeMin(nullptr, 0)));
    EXPECT_TRUE(std::isnan(ArrayStats<double>::ComputeMax(nullptr, 0)));
}

TEST(ArrayStatsTest, ComputeMean)
{
    const std::array<double, 4> data{ 1.0, 2.0, 3.0, 4.0 };

    EXPECT_DOUBLE_EQ(2.5, ArrayStats<double>::ComputeMean(data.data(), data.size()));
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputeMean(nullptr, 0));
}

TEST(ArrayStatsTest, ComputeStandardDeviation)
{
    const std::array<double, 5> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    const double mean{ ArrayStats<double>::ComputeMean(data.data(), data.size()) };
    const double std_dev{ ArrayStats<double>::ComputeStandardDeviation(data.data(), data.size(), mean) };

    EXPECT_NEAR(1.5811388300841898, std_dev, 1e-9);

    const double single_value{ 2.0 };
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputeStandardDeviation(&single_value, 1, single_value));
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputeStandardDeviation(nullptr, 0, 0.0));
}

TEST(ArrayStatsTest, ComputePercentileBounds)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputePercentile(data, 0.0));
    EXPECT_DOUBLE_EQ(4.0, ArrayStats<double>::ComputePercentile(data, 1.0));
}

TEST(ArrayStatsTest, ComputePercentileInterpolation)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(2.5, ArrayStats<double>::ComputePercentile(data, 0.5));
}

TEST(ArrayStatsTest, ComputePercentileEdgeCases)
{
    // Empty data returns 0 regardless of percentile within (0, 1)
    std::vector<double> empty_data{};
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputePercentile(empty_data, 0.5));

    // Percentile less than 0 returns 0, greater than 1 returns max value
    std::vector<double> data{ 1.0, 2.0, 3.0 };
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputePercentile(data, -0.1));
    EXPECT_DOUBLE_EQ(3.0, ArrayStats<double>::ComputePercentile(data, 1.5));

    // Unsorted data uses interpolation for percentile
    std::vector<double> unsorted{ 3.0, 1.0, 4.0, 2.0 };
    EXPECT_DOUBLE_EQ(1.75, ArrayStats<double>::ComputePercentile(unsorted, 0.25));
}

TEST(ArrayStatsTest, ComputeMedianOddSizedDataset)
{
    std::vector<double> data{ 3.0, 1.0, 2.0 };
    EXPECT_DOUBLE_EQ(2.0, ArrayStats<double>::ComputeMedian(data));
}

TEST(ArrayStatsTest, ComputeMedianEvenSizedDataset)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };
    EXPECT_DOUBLE_EQ(2.5, ArrayStats<double>::ComputeMedian(data));
}

TEST(ArrayStatsTest, ComputeMedianEmptyInput)
{
    std::vector<double> data{};
    EXPECT_DOUBLE_EQ(0.0, ArrayStats<double>::ComputeMedian(data));
}

TEST(ArrayStatsTest, ComputePercentileRangeTupleReturnsExpectedPercentiles)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    auto [min_value, max_value] = ArrayStats<double>::ComputePercentileRangeTuple(data, 0.2, 0.8);
    EXPECT_NEAR(1.8, min_value, 1e-9);
    EXPECT_NEAR(4.2, max_value, 1e-9);
}

TEST(ArrayStatsTest, ComputeScalingPercentileRangeTupleExtendsSymmetrically)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    auto scaled = ArrayStats<double>::ComputeScalingPercentileRangeTuple(data, 0.5, 0.2, 0.8);
    EXPECT_NEAR(0.6, std::get<0>(scaled), 1e-9);
    EXPECT_NEAR(5.4, std::get<1>(scaled), 1e-9);

    auto base{ ArrayStats<double>::ComputePercentileRangeTuple(data, 0.2, 0.8) };
    EXPECT_NEAR(std::get<0>(base) - std::get<0>(scaled),
                std::get<1>(scaled) - std::get<1>(base), 1e-9);
}

TEST(ArrayStatsTest, ComputeRangeTupleReturnsExpectedMinMax)
{
    std::vector<double> data{ 3.0, 1.0, 4.0, -1.0, 5.0 };
    auto range{ ArrayStats<double>::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(-1.0, std::get<0>(range));
    EXPECT_DOUBLE_EQ(5.0, std::get<1>(range));
}

TEST(ArrayStatsTest, ComputeScalingRangeTupleExtendsRangeSymmetrically)
{
    std::vector<double> data{ 3.0, 1.0, 4.0, -1.0, 5.0 };
    auto scaled{ ArrayStats<double>::ComputeScalingRangeTuple(data, 0.5) };
    EXPECT_DOUBLE_EQ(-4.0, std::get<0>(scaled));
    EXPECT_DOUBLE_EQ(8.0, std::get<1>(scaled));

    auto base{ ArrayStats<double>::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(base) - std::get<0>(scaled),
                     std::get<1>(scaled) - std::get<1>(base));
}

TEST(ArrayStatsTest, ScalingRangeDefaultAndEmpty)
{
    std::vector<double> data{ 1.0, 2.0, 3.0, 4.0, 5.0 };

    auto scaled_percentile{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(data, 0.0) };
    auto base_percentile{ ArrayStats<double>::ComputePercentileRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(base_percentile), std::get<0>(scaled_percentile));
    EXPECT_DOUBLE_EQ(std::get<1>(base_percentile), std::get<1>(scaled_percentile));

    auto scaled_range{ ArrayStats<double>::ComputeScalingRangeTuple(data, 0.0) };
    auto base_range{ ArrayStats<double>::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(std::get<0>(base_range), std::get<0>(scaled_range));
    EXPECT_DOUBLE_EQ(std::get<1>(base_range), std::get<1>(scaled_range));

    std::vector<double> empty{};
    auto empty_percentile{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(empty, 0.0) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(empty_percentile));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(empty_percentile));

    auto empty_range{ ArrayStats<double>::ComputeScalingRangeTuple(empty, 0.0) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(empty_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(empty_range));
}

TEST(ArrayStatsTest, RangeFunctionsEmptyInput)
{
    std::vector<double> data{};

    auto percentile_range{ ArrayStats<double>::ComputePercentileRangeTuple(data) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(percentile_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(percentile_range));

    auto scaling_percentile_range{
        ArrayStats<double>::ComputeScalingPercentileRangeTuple(data, 0.5)
    };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(scaling_percentile_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(scaling_percentile_range));

    auto range{ ArrayStats<double>::ComputeRangeTuple(data) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(range));

    auto scaling_range{ ArrayStats<double>::ComputeScalingRangeTuple(data, 0.5) };
    EXPECT_DOUBLE_EQ(0.0, std::get<0>(scaling_range));
    EXPECT_DOUBLE_EQ(0.0, std::get<1>(scaling_range));
}

TEST(ArrayStatsTest, ComputeNormEuclideanDistance)
{
    std::array<double, 3> v1{ 0.0, 0.0, 0.0 };
    std::array<double, 3> v2{ 3.0, 4.0, 0.0 };
    auto distance{ ArrayStats<double>::ComputeNorm(v1, v2) };
    EXPECT_DOUBLE_EQ(5.0, distance);
}

TEST(ArrayStatsTest, ComputeRankReturnsExpectedPosition)
{
    std::array<double, 3> values{ 3.0, 1.0, 2.0 };
    auto rank{ ArrayStats<double>::ComputeRank(values, 1) };
    EXPECT_EQ(1, rank);
}

TEST(ArrayStatsTest, ComputeRankWithDuplicates)
{
    std::array<double, 4> values{ 2.0, 3.0, 3.0, 4.0 };
    auto rank_index1{ ArrayStats<double>::ComputeRank(values, 1) };
    auto rank_index2{ ArrayStats<double>::ComputeRank(values, 2) };

    EXPECT_EQ(rank_index1, rank_index2);
    EXPECT_EQ(2, rank_index1);
    EXPECT_EQ(1, ArrayStats<double>::ComputeRank(values, 0));
    EXPECT_EQ(4, ArrayStats<double>::ComputeRank(values, 3));
}

TEST(ArrayStatsTest, ComputeRankIndexOutOfRangeThrows)
{
    std::array<double, 3> values{ 3.0, 1.0, 2.0 };
    EXPECT_THROW(ArrayStats<double>::ComputeRank(values, 3), std::out_of_range);
}

TEST(ArrayStatsTest, ComputeWithThreadSize)
{
    const std::vector<double> data{ 1.0, 2.0, 3.0, 4.0 };

    auto min_single{ ArrayStats<double>::ComputeMin(data.data(), data.size(), 1) };
    auto min_multi{ ArrayStats<double>::ComputeMin(data.data(), data.size(), 4) };
    EXPECT_DOUBLE_EQ(min_single, min_multi);

    auto max_single{ ArrayStats<double>::ComputeMax(data.data(), data.size(), 1) };
    auto max_multi{ ArrayStats<double>::ComputeMax(data.data(), data.size(), 4) };
    EXPECT_DOUBLE_EQ(max_single, max_multi);

    auto mean_single{ ArrayStats<double>::ComputeMean(data.data(), data.size(), 1) };
    auto mean_multi{ ArrayStats<double>::ComputeMean(data.data(), data.size(), 4) };
    EXPECT_DOUBLE_EQ(mean_single, mean_multi);

    auto std_single{ ArrayStats<double>::ComputeStandardDeviation(
        data.data(), data.size(), mean_single, 1)
    };
    auto std_multi{ ArrayStats<double>::ComputeStandardDeviation(
        data.data(), data.size(), mean_multi, 4)
    };
    EXPECT_DOUBLE_EQ(std_single, std_multi);

    auto range_single{ ArrayStats<double>::ComputeRangeTuple(data, 1) };
    auto range_multi{ ArrayStats<double>::ComputeRangeTuple(data, 4) };
    EXPECT_DOUBLE_EQ(std::get<0>(range_single), std::get<0>(range_multi));
    EXPECT_DOUBLE_EQ(std::get<1>(range_single), std::get<1>(range_multi));
}
