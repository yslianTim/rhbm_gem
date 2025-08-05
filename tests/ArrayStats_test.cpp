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

TEST(ArrayStatsTest, ComputeRankIndexOutOfRangeThrows)
{
    std::array<double, 3> values{ 3.0, 1.0, 2.0 };
    EXPECT_THROW(ArrayStats<double>::ComputeRank(values, 3), std::out_of_range);
}
