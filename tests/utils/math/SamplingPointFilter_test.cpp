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
    const auto point_scores{
        rg::BuildSamplingPointScoreList(point_list, { MakeVector({ 1.0, 0.0, 0.0 }) }, 0.0)
    };

    ASSERT_EQ(point_scores.size(), point_list.size());
    EXPECT_TRUE(std::all_of(point_scores.begin(), point_scores.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(SamplingPointFilterTest, RejectsPointsWithinAngleThresholdOfRejectDirections)
{
    const auto point_scores{
        rg::BuildSamplingPointScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            30.0)
    };

    ASSERT_EQ(point_scores.size(), 4u);
    EXPECT_FLOAT_EQ(point_scores.at(0), 1.0f);
    EXPECT_FLOAT_EQ(point_scores.at(1), 0.0f);
    EXPECT_FLOAT_EQ(point_scores.at(2), 1.0f);
    EXPECT_FLOAT_EQ(point_scores.at(3), 1.0f);
}

TEST(SamplingPointFilterTest, KeepsPerpendicularAndOppositeDirections)
{
    const auto point_scores{
        rg::BuildSamplingPointScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            89.0)
    };

    ASSERT_EQ(point_scores.size(), 4u);
    EXPECT_FLOAT_EQ(point_scores.at(0), 1.0f);
    EXPECT_FLOAT_EQ(point_scores.at(1), 0.0f);
    EXPECT_FLOAT_EQ(point_scores.at(2), 1.0f);
    EXPECT_FLOAT_EQ(point_scores.at(3), 1.0f);
}

TEST(SamplingPointFilterTest, IgnoresZeroLengthRejectDirections)
{
    const auto point_list{ MakePointList() };
    const auto point_scores{
        rg::BuildSamplingPointScoreList(point_list, { MakeVector({ 0.0, 0.0, 0.0 }) }, 45.0)
    };

    ASSERT_EQ(point_scores.size(), point_list.size());
    EXPECT_TRUE(std::all_of(point_scores.begin(), point_scores.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(SamplingPointFilterTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    const auto point_scores{
        rg::BuildSamplingPointScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_EQ(point_scores.size(), 4u);
    EXPECT_FLOAT_EQ(point_scores.front(), 1.0f);
}

TEST(SamplingPointFilterTest, RejectsInvalidAngles)
{
    const auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)rg::BuildSamplingPointScoreList(point_list, {}, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointScoreList(
            point_list,
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointScoreList(point_list, {}, 181.0),
        std::invalid_argument);
}

TEST(SamplingPointFilterTest, RejectsRejectDirectionsWithUnexpectedDimension)
{
    EXPECT_THROW(
        (void)rg::BuildSamplingPointScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0 }) },
            45.0),
        std::invalid_argument);
}

TEST(SamplingPointFilterTest, RejectsRejectDirectionsWithNonFiniteValues)
{
    EXPECT_THROW(
        (void)rg::BuildSamplingPointScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, std::numeric_limits<double>::quiet_NaN() }) },
            45.0),
        std::invalid_argument);
}
