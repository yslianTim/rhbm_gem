#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/math/PotentialSampleSelection.hpp>

namespace rg = rhbm_gem;

namespace {

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

TEST(PotentialSampleSelectionTest, EmptyRejectPositionsSelectAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {})
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 1u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, ZeroAngleSelectsAllSamplingPoints)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            0.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 1u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, RejectsPointsWithinAngleThresholdOfRejectPositions)
{
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            30.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, KeepsPerpendicularAndOppositeDirections)
{
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            89.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, IgnoresRejectPositionsAtReferencePosition)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 0.0f, 0.0f, 0.0f } },
            45.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 1u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            45.0)
    };

    ASSERT_FALSE(selected_indices.empty());
    EXPECT_EQ(selected_indices.front(), 0u);
}

TEST(PotentialSampleSelectionTest, RejectsInvalidAngles)
{
    const auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(PotentialSampleSelectionTest, RejectsRejectPositionsWithNonFiniteValues)
{
    EXPECT_THROW(
        (void)rg::BuildSelectedLocalPotentialSampleIndexList(
            MakePointList(),
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{
                1.0f,
                0.0f,
                std::numeric_limits<float>::quiet_NaN()
            } },
            45.0),
        std::invalid_argument);
}

TEST(PotentialSampleSelectionTest, UsesReferencePositionForWorldSpaceSamplingPoints)
{
    const SamplingPointList point_list{
        SamplingPoint{ 1.0f, { 10.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 11.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 10.0f, 1.0f, 0.0f } }
    };

    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 10.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 11.0f, 0.0f, 0.0f } },
            30.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u }));
}
