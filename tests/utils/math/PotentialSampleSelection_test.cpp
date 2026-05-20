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

TEST(PotentialSampleSelectionTest, ZeroAngleAppliesVoronoiOwnership)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            0.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u, 3u }));
}

TEST(PotentialSampleSelectionTest, VoronoiOwnershipRejectsTiePoints)
{
    const SamplingPointList point_list{
        SamplingPoint{ 0.0f, { 0.0f, 0.0f, 0.0f } },
        SamplingPoint{ 0.5f, { 0.5f, 0.0f, 0.0f } },
        SamplingPoint{ 0.49f, { 0.49f, 0.0f, 0.0f } }
    };

    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            0.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u }));
}

TEST(PotentialSampleSelectionTest, VoronoiOwnershipRejectsNeighborOwnedPoints)
{
    const SamplingPointList point_list{
        SamplingPoint{ 0.25f, { 0.25f, 0.0f, 0.0f } },
        SamplingPoint{ 0.75f, { 0.75f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 0.0f, 1.0f, 0.0f } }
    };

    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
            0.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u, 2u }));
}

TEST(PotentialSampleSelectionTest, RejectsPointsCloserToNearestNeighborThanReference)
{
    const SamplingPointList point_list{
        SamplingPoint{ 0.05f, { 0.05f, 0.0f, 0.0f } },
        SamplingPoint{ 0.15f, { 0.15f, 0.0f, 0.0f } },
        SamplingPoint{ 0.25f, { 0.25f, 0.0f, 0.0f } }
    };

    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {
                std::array<float, 3>{ 10.0f, 0.0f, 0.0f },
                std::array<float, 3>{ 0.2f, 0.0f, 0.0f }
            },
            0.0)
    };

    EXPECT_EQ(selected_indices, std::vector<std::size_t>({ 0u }));
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

TEST(PotentialSampleSelectionTest, IgnoresRejectPositionsAtReferencePositionWithZeroAngle)
{
    const auto point_list{ MakePointList() };
    const auto selected_indices{
        rg::BuildSelectedLocalPotentialSampleIndexList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{ 0.0f, 0.0f, 0.0f } },
            0.0)
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

TEST(PotentialSampleSelectionTest, RejectsRejectPositionsWithNonFiniteValuesAtZeroAngle)
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
            0.0),
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

TEST(PotentialSampleSelectionTest, LowestResponseDecileUsesDefaultRatio)
{
    LocalPotentialSampleList sample_list;
    for (std::size_t i = 0; i < 11u; i++)
    {
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(10u - i),
            SamplingPoint{ 1.0f }
        });
    }

    const auto retained_samples{
        rg::KeepLowestResponseDecileByDistance(std::move(sample_list))
    };

    ASSERT_EQ(retained_samples.size(), 2u);
    EXPECT_FLOAT_EQ(retained_samples.at(0).response, 0.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).response, 1.0f);
}

TEST(PotentialSampleSelectionTest, LowestResponseDecileUsesCustomRatioByDistance)
{
    const LocalPotentialSampleList sample_list{
        LocalPotentialSample{ 9.0f, SamplingPoint{ 2.0f } },
        LocalPotentialSample{ 4.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 1.0f, SamplingPoint{ 2.0f } },
        LocalPotentialSample{ 2.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 3.0f, SamplingPoint{ 2.0f } },
        LocalPotentialSample{ 8.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 5.0f, SamplingPoint{ 2.0f } }
    };

    const auto retained_samples{
        rg::KeepLowestResponseDecileByDistance(sample_list, 0.5)
    };

    ASSERT_EQ(retained_samples.size(), 4u);
    EXPECT_FLOAT_EQ(retained_samples.at(0).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(0).response, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).response, 4.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(2).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(2).response, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(3).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(3).response, 3.0f);
}

TEST(PotentialSampleSelectionTest, LowestResponseDecileAllowsKeepingAllSamples)
{
    const LocalPotentialSampleList sample_list{
        LocalPotentialSample{ 9.0f, SamplingPoint{ 2.0f } },
        LocalPotentialSample{ 4.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 1.0f, SamplingPoint{ 2.0f } },
        LocalPotentialSample{ 2.0f, SamplingPoint{ 1.0f } }
    };

    const auto retained_samples{
        rg::KeepLowestResponseDecileByDistance(sample_list, 1.0)
    };

    ASSERT_EQ(retained_samples.size(), 4u);
    EXPECT_FLOAT_EQ(retained_samples.at(0).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(0).response, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).response, 4.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(2).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(2).response, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(3).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(3).response, 9.0f);
}

TEST(PotentialSampleSelectionTest, LowestResponseDecileRejectsInvalidRatios)
{
    const LocalPotentialSampleList sample_list{
        LocalPotentialSample{ 1.0f, SamplingPoint{ 1.0f } }
    };

    EXPECT_THROW(
        (void)rg::KeepLowestResponseDecileByDistance(sample_list, 0.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::KeepLowestResponseDecileByDistance(sample_list, -0.1),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::KeepLowestResponseDecileByDistance(sample_list, 1.1),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::KeepLowestResponseDecileByDistance(
            sample_list,
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::KeepLowestResponseDecileByDistance(
            sample_list,
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
}

TEST(PotentialSampleSelectionTest, MedianResponseByDistanceHandlesEmptySampleList)
{
    EXPECT_TRUE(rg::GetMedianResponseByDistance({}).empty());
}

TEST(PotentialSampleSelectionTest, MedianResponseByDistanceReturnsOneSamplePerDistance)
{
    const LocalPotentialSampleList sample_list{
        LocalPotentialSample{ 10.0f, SamplingPoint{ 2.0f, { 20.0f, 0.0f, 0.0f } } },
        LocalPotentialSample{ 3.0f, SamplingPoint{ 1.0f, { 10.0f, 0.0f, 0.0f } } },
        LocalPotentialSample{ 7.0f, SamplingPoint{ 1.0f, { 11.0f, 0.0f, 0.0f } } },
        LocalPotentialSample{ 5.0f, SamplingPoint{ 1.0f, { 12.0f, 0.0f, 0.0f } } },
        LocalPotentialSample{ 2.0f, SamplingPoint{ 2.0f, { 21.0f, 0.0f, 0.0f } } }
    };

    const auto median_samples{ rg::GetMedianResponseByDistance(sample_list) };

    ASSERT_EQ(median_samples.size(), 2u);
    EXPECT_FLOAT_EQ(median_samples.at(0).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(median_samples.at(0).response, 5.0f);
    EXPECT_EQ(median_samples.at(0).point.position, (std::array<float, 3>{ 12.0f, 0.0f, 0.0f }));
    EXPECT_FLOAT_EQ(median_samples.at(1).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(median_samples.at(1).response, 6.0f);
    EXPECT_EQ(median_samples.at(1).point.position, (std::array<float, 3>{ 21.0f, 0.0f, 0.0f }));
}
