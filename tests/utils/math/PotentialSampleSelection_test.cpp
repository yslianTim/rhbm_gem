#include <gtest/gtest.h>

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/PotentialSampleSelection.hpp>

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

TEST(PotentialSampleSelectionTest, EmptyRejectDirectionsSelectAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(point_list, {})
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 1u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, ZeroAngleSelectsAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            0.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 1u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, RejectsPointsWithinAngleThresholdOfRejectDirections)
{
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            30.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, KeepsPerpendicularAndOppositeDirections)
{
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            89.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, IgnoresZeroLengthRejectDirections)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { MakeVector({ 0.0, 0.0, 0.0 }) },
            45.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 1u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, 0.0 }) },
            45.0)
    };

    ASSERT_FALSE(selected_indices.empty());
    EXPECT_EQ(selected_indices.front(), 0u);
}

TEST(PotentialSampleSelectionTest, RejectsInvalidAngles)
{
    const auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(point_list, {}, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(point_list, {}, 181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            {},
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(PotentialSampleSelectionTest, RejectsRejectDirectionsWithUnexpectedDimension)
{
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0 }) },
            45.0),
        std::invalid_argument);
}

TEST(PotentialSampleSelectionTest, RejectsRejectDirectionsWithNonFiniteValues)
{
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { MakeVector({ 1.0, 0.0, std::numeric_limits<double>::quiet_NaN() }) },
            45.0),
        std::invalid_argument);
}
