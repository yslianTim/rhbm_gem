#include <gtest/gtest.h>

#include <cmath>
#include <tuple>
#include <vector>

#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

#include "detail/PotentialSeriesOps.hpp"

namespace rg = rhbm_gem;

TEST(PotentialSeriesOpsTest, RangeAndBinningStayStableForLocalEntry)
{
    rg::LocalPotentialEntry entry;
    entry.SetDistanceAndMapValueList({
        {0.0f, 2.0f},
        {0.2f, 4.0f},
        {0.7f, 8.0f},
        {1.0f, -1.0f},
    });
    const auto distance_range{ rg::series_ops::ComputeDistanceRange(entry, 0.0) };
    const auto map_value_range{ rg::series_ops::ComputeMapValueRange(entry, 0.0) };
    const auto binned{ rg::series_ops::BuildBinnedDistanceAndMapValueList(entry, 2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(distance_range), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(distance_range), 1.0f);
    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), -1.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 8.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(0)), 3.0f);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(1)), 8.0f);
}

TEST(PotentialSeriesOpsTest, LinearModelTransformKeepsPositiveSamplesOnly)
{
    rg::LocalPotentialEntry entry;
    entry.SetDistanceAndMapValueList({
        {0.1f, 4.0f},
        {0.2f, -2.0f},
        {0.3f, 8.0f},
    });
    const auto transformed{ rg::series_ops::BuildLinearModelDistanceAndMapValueList(entry) };
    ASSERT_EQ(transformed.size(), 2U);

    Eigen::VectorXd init{ Eigen::VectorXd::Zero(3) };
    const auto expected0{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.1f, 4.0f, init) };
    const auto expected1{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.3f, 8.0f, init) };

    EXPECT_NEAR(std::get<0>(transformed.at(0)), expected0(1), 1e-6);
    EXPECT_NEAR(std::get<1>(transformed.at(0)), expected0(2), 1e-6);
    EXPECT_NEAR(std::get<0>(transformed.at(1)), expected1(1), 1e-6);
    EXPECT_NEAR(std::get<1>(transformed.at(1)), expected1(2), 1e-6);
}
