#include <gtest/gtest.h>

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/LocalPotentialSampleScoring.hpp>

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

TEST(LocalPotentialSampleScoringTest, DefaultScoreListKeepsAllSamplingPoints)
{
    const auto score_list{ rg::BuildDefaultLocalPotentialSampleScoreList(4) };

    ASSERT_EQ(score_list.size(), 4u);
    EXPECT_TRUE(std::all_of(score_list.begin(), score_list.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(LocalPotentialSampleScoringTest, ZeroAngleKeepsAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto score_list{
        rg::BuildLocalPotentialSampleScoreList(
            point_list,
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            0.0)
    };

    ASSERT_EQ(score_list.size(), point_list.size());
    EXPECT_TRUE(std::all_of(score_list.begin(), score_list.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(LocalPotentialSampleScoringTest, RejectsPointsWithinAngleThresholdOfRejectDirections)
{
    const auto score_list{
        rg::BuildLocalPotentialSampleScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            30.0)
    };

    ASSERT_EQ(score_list.size(), 4u);
    EXPECT_FLOAT_EQ(score_list.at(0), 1.0f);
    EXPECT_FLOAT_EQ(score_list.at(1), 0.0f);
    EXPECT_FLOAT_EQ(score_list.at(2), 1.0f);
    EXPECT_FLOAT_EQ(score_list.at(3), 1.0f);
}

TEST(LocalPotentialSampleScoringTest, KeepsPerpendicularAndOppositeDirections)
{
    const auto score_list{
        rg::BuildLocalPotentialSampleScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            89.0)
    };

    ASSERT_EQ(score_list.size(), 4u);
    EXPECT_FLOAT_EQ(score_list.at(0), 1.0f);
    EXPECT_FLOAT_EQ(score_list.at(1), 0.0f);
    EXPECT_FLOAT_EQ(score_list.at(2), 1.0f);
    EXPECT_FLOAT_EQ(score_list.at(3), 1.0f);
}

TEST(LocalPotentialSampleScoringTest, IgnoresZeroLengthRejectDirections)
{
    const auto point_list{ MakePointList() };
    const auto score_list{
        rg::BuildLocalPotentialSampleScoreList(
            point_list,
            { MakeVector({ 0.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_EQ(score_list.size(), point_list.size());
    EXPECT_TRUE(std::all_of(score_list.begin(), score_list.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(LocalPotentialSampleScoringTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    const auto score_list{
        rg::BuildLocalPotentialSampleScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_EQ(score_list.size(), 4u);
    EXPECT_FLOAT_EQ(score_list.front(), 1.0f);
}

TEST(LocalPotentialSampleScoringTest, RejectsInvalidAngles)
{
    const auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)rg::BuildLocalPotentialSampleScoreList(point_list, {}, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildLocalPotentialSampleScoreList(
            point_list,
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildLocalPotentialSampleScoreList(point_list, {}, 181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildLocalPotentialSampleScoreList(
            point_list,
            {},
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(LocalPotentialSampleScoringTest, RejectsRejectDirectionsWithUnexpectedDimension)
{
    EXPECT_THROW(
        (void)rg::BuildLocalPotentialSampleScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0 }) },
            45.0),
        std::invalid_argument);
}

TEST(LocalPotentialSampleScoringTest, RejectsRejectDirectionsWithNonFiniteValues)
{
    EXPECT_THROW(
        (void)rg::BuildLocalPotentialSampleScoreList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, std::numeric_limits<double>::quiet_NaN() }) },
            45.0),
        std::invalid_argument);
}

TEST(LocalPotentialSampleScoringTest, CleanScoreListWithEmptyRejectDirectionsKeepsAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto score_list{ rg::BuildLocalPotentialSampleCleanScoreList(point_list, {}) };

    ASSERT_EQ(score_list.size(), point_list.size());
    EXPECT_TRUE(std::all_of(score_list.begin(), score_list.end(), [](float score)
    {
        return score == 1.0f;
    }));
}
