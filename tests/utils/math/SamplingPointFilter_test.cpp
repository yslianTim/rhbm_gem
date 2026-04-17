#include <gtest/gtest.h>

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/SamplingPointFilter.hpp>

namespace rg = rhbm_gem;

namespace {

Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (double value : values)
    {
        result(index++) = value;
    }
    return result;
}

SamplingPointList MakePointList()
{
    return {
        SamplingPoint{ 0.0f, { 0.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 1.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 0.0f, 1.0f, 0.0f } },
        SamplingPoint{ 1.0f, { -1.0f, 0.0f, 0.0f } }
    };
}

} // namespace

TEST(SamplingPointFilterTest, ZeroAngleKeepsAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto point_weights{
        rg::BuildSamplingPointWeightList(point_list, { MakeVector({ 1.0, 0.0, 0.0 }) }, 0.0)
    };

    const auto filtered_point_list{
        rg::SelectSamplingPoint(point_list, { MakeVector({ 1.0, 0.0, 0.0 }) }, 0.0)
    };

    ASSERT_EQ(point_weights.size(), point_list.size());
    EXPECT_TRUE(std::all_of(point_weights.begin(), point_weights.end(), [](float weight)
    {
        return weight == 1.0f;
    }));
    ASSERT_EQ(filtered_point_list.size(), point_list.size());
    EXPECT_EQ(filtered_point_list.front().position, point_list.front().position);
    EXPECT_EQ(filtered_point_list.back().position, point_list.back().position);
}

TEST(SamplingPointFilterTest, RejectsPointsWithinAngleThresholdOfRejectDirections)
{
    const auto point_weights{
        rg::BuildSamplingPointWeightList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            30.0)
    };
    const auto filtered_point_list{
        rg::SelectSamplingPoint(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            30.0)
    };

    ASSERT_EQ(point_weights.size(), 4u);
    EXPECT_FLOAT_EQ(point_weights.at(0), 1.0f);
    EXPECT_FLOAT_EQ(point_weights.at(1), 0.0f);
    EXPECT_FLOAT_EQ(point_weights.at(2), 1.0f);
    EXPECT_FLOAT_EQ(point_weights.at(3), 1.0f);
    ASSERT_EQ(filtered_point_list.size(), 3u);
    EXPECT_EQ(filtered_point_list.at(0).position, (std::array<float, 3>{ 0.0f, 0.0f, 0.0f }));
    EXPECT_EQ(filtered_point_list.at(1).position, (std::array<float, 3>{ 0.0f, 1.0f, 0.0f }));
    EXPECT_EQ(filtered_point_list.at(2).position, (std::array<float, 3>{ -1.0f, 0.0f, 0.0f }));
}

TEST(SamplingPointFilterTest, KeepsPerpendicularAndOppositeDirections)
{
    const auto filtered_point_list{
        rg::SelectSamplingPoint(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            89.0)
    };

    ASSERT_EQ(filtered_point_list.size(), 3u);
    EXPECT_EQ(filtered_point_list.at(1).position, (std::array<float, 3>{ 0.0f, 1.0f, 0.0f }));
    EXPECT_EQ(filtered_point_list.at(2).position, (std::array<float, 3>{ -1.0f, 0.0f, 0.0f }));
}

TEST(SamplingPointFilterTest, IgnoresZeroLengthRejectDirections)
{
    const auto point_list{ MakePointList() };
    const auto point_weights{
        rg::BuildSamplingPointWeightList(point_list, { MakeVector({ 0.0, 0.0, 0.0 }) }, 45.0)
    };

    const auto filtered_point_list{
        rg::SelectSamplingPoint(point_list, { MakeVector({ 0.0, 0.0, 0.0 }) }, 45.0)
    };

    ASSERT_EQ(point_weights.size(), point_list.size());
    EXPECT_TRUE(std::all_of(point_weights.begin(), point_weights.end(), [](float weight)
    {
        return weight == 1.0f;
    }));
    ASSERT_EQ(filtered_point_list.size(), point_list.size());
}

TEST(SamplingPointFilterTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    const auto point_weights{
        rg::BuildSamplingPointWeightList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            45.0)
    };
    const auto filtered_point_list{
        rg::SelectSamplingPoint(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_EQ(point_weights.size(), 4u);
    EXPECT_FLOAT_EQ(point_weights.front(), 1.0f);
    ASSERT_FALSE(filtered_point_list.empty());
    EXPECT_EQ(filtered_point_list.front().position, (std::array<float, 3>{ 0.0f, 0.0f, 0.0f }));
}

TEST(SamplingPointFilterTest, RejectsInvalidAngles)
{
    const auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)rg::SelectSamplingPoint(point_list, {}, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::SelectSamplingPoint(
            point_list,
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::SelectSamplingPoint(point_list, {}, 181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointWeightList(point_list, {}, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointWeightList(
            point_list,
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointWeightList(point_list, {}, 181.0),
        std::invalid_argument);
}

TEST(SamplingPointFilterTest, RejectsRejectDirectionsWithUnexpectedDimension)
{
    EXPECT_THROW(
        (void)rg::SelectSamplingPoint(
            MakePointList(),
            { MakeVector({ 1.0, 0.0 }) },
            45.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointWeightList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0 }) },
            45.0),
        std::invalid_argument);
}

TEST(SamplingPointFilterTest, RejectsRejectDirectionsWithNonFiniteValues)
{
    EXPECT_THROW(
        (void)rg::BuildSamplingPointWeightList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, std::numeric_limits<double>::quiet_NaN() }) },
            45.0),
        std::invalid_argument);
}
