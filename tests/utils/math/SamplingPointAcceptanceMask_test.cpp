#include <gtest/gtest.h>

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/SamplingPointAcceptanceMask.hpp>

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

TEST(SamplingPointAcceptanceMaskTest, ZeroAngleKeepsAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto acceptance_mask{
        rg::BuildSamplingPointAcceptanceMask(
            point_list,
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            0.0)
    };

    ASSERT_EQ(acceptance_mask.size(), point_list.size());
    EXPECT_TRUE(std::all_of(acceptance_mask.begin(), acceptance_mask.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(SamplingPointAcceptanceMaskTest, RejectsPointsWithinAngleThresholdOfRejectDirections)
{
    const auto acceptance_mask{
        rg::BuildSamplingPointAcceptanceMask(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            30.0)
    };

    ASSERT_EQ(acceptance_mask.size(), 4u);
    EXPECT_FLOAT_EQ(acceptance_mask.at(0), 1.0f);
    EXPECT_FLOAT_EQ(acceptance_mask.at(1), 0.0f);
    EXPECT_FLOAT_EQ(acceptance_mask.at(2), 1.0f);
    EXPECT_FLOAT_EQ(acceptance_mask.at(3), 1.0f);
}

TEST(SamplingPointAcceptanceMaskTest, KeepsPerpendicularAndOppositeDirections)
{
    const auto acceptance_mask{
        rg::BuildSamplingPointAcceptanceMask(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            89.0)
    };

    ASSERT_EQ(acceptance_mask.size(), 4u);
    EXPECT_FLOAT_EQ(acceptance_mask.at(0), 1.0f);
    EXPECT_FLOAT_EQ(acceptance_mask.at(1), 0.0f);
    EXPECT_FLOAT_EQ(acceptance_mask.at(2), 1.0f);
    EXPECT_FLOAT_EQ(acceptance_mask.at(3), 1.0f);
}

TEST(SamplingPointAcceptanceMaskTest, IgnoresZeroLengthRejectDirections)
{
    const auto point_list{ MakePointList() };
    const auto acceptance_mask{
        rg::BuildSamplingPointAcceptanceMask(
            point_list,
            { MakeVector({ 0.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_EQ(acceptance_mask.size(), point_list.size());
    EXPECT_TRUE(std::all_of(acceptance_mask.begin(), acceptance_mask.end(), [](float score)
    {
        return score == 1.0f;
    }));
}

TEST(SamplingPointAcceptanceMaskTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    const auto acceptance_mask{
        rg::BuildSamplingPointAcceptanceMask(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_EQ(acceptance_mask.size(), 4u);
    EXPECT_FLOAT_EQ(acceptance_mask.front(), 1.0f);
}

TEST(SamplingPointAcceptanceMaskTest, RejectsInvalidAngles)
{
    const auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)rg::BuildSamplingPointAcceptanceMask(point_list, {}, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointAcceptanceMask(
            point_list,
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointAcceptanceMask(point_list, {}, 181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSamplingPointAcceptanceMask(
            point_list,
            {},
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(SamplingPointAcceptanceMaskTest, RejectsRejectDirectionsWithUnexpectedDimension)
{
    EXPECT_THROW(
        (void)rg::BuildSamplingPointAcceptanceMask(
            MakePointList(),
            { MakeVector({ 1.0, 0.0 }) },
            45.0),
        std::invalid_argument);
}

TEST(SamplingPointAcceptanceMaskTest, RejectsRejectDirectionsWithNonFiniteValues)
{
    EXPECT_THROW(
        (void)rg::BuildSamplingPointAcceptanceMask(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, std::numeric_limits<double>::quiet_NaN() }) },
            45.0),
        std::invalid_argument);
}
