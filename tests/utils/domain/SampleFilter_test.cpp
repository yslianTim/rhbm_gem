#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/domain/SampleFilter.hpp>

namespace sf = rhbm_gem::sample_filter;

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

std::vector<std::array<float, 3>> GetPositions(const SamplingPointList & point_list)
{
    std::vector<std::array<float, 3>> positions;
    positions.reserve(point_list.size());
    for (const auto & point : point_list)
    {
        positions.emplace_back(point.position);
    }
    return positions;
}

std::vector<std::array<float, 3>> GetSelectedPositions(const SamplingPointList & point_list)
{
    std::vector<std::array<float, 3>> positions;
    for (const auto & point : point_list)
    {
        if (point.is_selected)
        {
            positions.emplace_back(point.position);
        }
    }
    return positions;
}

bool AllSelected(const SamplingPointList & point_list)
{
    for (const auto & point : point_list)
    {
        if (!point.is_selected) return false;
    }
    return true;
}

} // namespace

TEST(SampleFilterTest, BuildsMedianResponseSampleEntriesByRadius)
{
    const LocalPotentialSampleList sample_entries{
        LocalPotentialSample{ 100.0f, SamplingPoint{ 1.2f } },
        LocalPotentialSample{ 2.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 1.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 3.0f, SamplingPoint{ 1.0f } },
        LocalPotentialSample{ 10.0f, SamplingPoint{ 1.2f } },
        LocalPotentialSample{ 11.0f, SamplingPoint{ 1.2f } }
    };

    const auto actual{ sf::BuildMedianResponseSampleEntriesByRadius(sample_entries) };

    ASSERT_EQ(actual.size(), 2u);
    EXPECT_FLOAT_EQ(actual.at(0).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(actual.at(0).response, 2.0f);
    EXPECT_FLOAT_EQ(actual.at(1).point.distance, 1.2f);
    EXPECT_FLOAT_EQ(actual.at(1).response, 11.0f);
}

TEST(SampleFilterTest, BuildsEmptyMedianResponseSampleEntriesByRadius)
{
    const LocalPotentialSampleList sample_entries;

    const auto actual{ sf::BuildMedianResponseSampleEntriesByRadius(sample_entries) };

    EXPECT_TRUE(actual.empty());
}

TEST(SampleFilterTest, BuildsResponseShiftedSampleEntries)
{
    const LocalPotentialSampleList sample_entries{
        LocalPotentialSample{
            3.5f,
            SamplingPoint{ 1.0f, { 1.0f, 2.0f, 3.0f }, false }
        },
        LocalPotentialSample{
            -1.0f,
            SamplingPoint{ 2.0f, { 4.0f, 5.0f, 6.0f }, true }
        }
    };

    const auto actual{ sf::BuildResponseShiftedSampleEntries(sample_entries, 0.5) };

    ASSERT_EQ(actual.size(), 2u);
    EXPECT_FLOAT_EQ(actual.at(0).response, 3.0f);
    EXPECT_FLOAT_EQ(actual.at(1).response, -1.5f);
    EXPECT_EQ(actual.at(0).point.position, sample_entries.at(0).point.position);
    EXPECT_EQ(actual.at(1).point.position, sample_entries.at(1).point.position);
    EXPECT_EQ(actual.at(0).point.is_selected, sample_entries.at(0).point.is_selected);
    EXPECT_EQ(actual.at(1).point.is_selected, sample_entries.at(1).point.is_selected);
}

TEST(SampleFilterTest, EmptyRejectPositionsSelectAllSamplingPoints)
{
    auto point_list{ MakePointList() };
    point_list.at(1).is_selected = false;

    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        {});

    EXPECT_EQ(GetSelectedPositions(point_list), GetPositions(point_list));
    EXPECT_TRUE(AllSelected(point_list));
}

TEST(SampleFilterTest, RejectsPointsCloserToNearestNeighborThanReference)
{
    SamplingPointList point_list{
        SamplingPoint{ 0.05f, { 0.05f, 0.0f, 0.0f } },
        SamplingPoint{ 0.15f, { 0.15f, 0.0f, 0.0f } },
        SamplingPoint{ 0.25f, { 0.25f, 0.0f, 0.0f } }
    };

    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        {
            std::array<float, 3>{ 10.0f, 0.0f, 0.0f },
            std::array<float, 3>{ 0.2f, 0.0f, 0.0f }
        },
        0.0);

    const std::vector<std::array<float, 3>> expected_positions{
        { 0.05f, 0.0f, 0.0f }
    };
    EXPECT_EQ(point_list.size(), 3u);
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
    EXPECT_FALSE(point_list.at(1).is_selected);
    EXPECT_FALSE(point_list.at(2).is_selected);
}

TEST(SampleFilterTest, UsesThirtyDegreeDefaultRejectAngle)
{
    SamplingPointList point_list{
        SamplingPoint{ 0.0f, { 0.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 1.0f, 0.2f, 0.0f } },
        SamplingPoint{ 1.0f, { 1.0f, 0.3f, 0.0f } }
    };

    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 2.0f, 0.0f, 0.0f } });

    const std::vector<std::array<float, 3>> expected_positions{
        { 0.0f, 0.0f, 0.0f }
    };
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
    EXPECT_FALSE(point_list.at(1).is_selected);
    EXPECT_FALSE(point_list.at(2).is_selected);
}

TEST(SampleFilterTest, RejectsPointsWithinAngleThresholdOfRejectPositions)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
        30.0);

    const std::vector<std::array<float, 3>> expected_positions{
        { 0.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { -1.0f, 0.0f, 0.0f }
    };
    EXPECT_EQ(point_list.size(), 4u);
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
    EXPECT_FALSE(point_list.at(1).is_selected);
}

TEST(SampleFilterTest, KeepsPerpendicularAndOppositeDirections)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
        89.0);

    const std::vector<std::array<float, 3>> expected_positions{
        { 0.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { -1.0f, 0.0f, 0.0f }
    };
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
}

TEST(SampleFilterTest, IgnoresRejectPositionsAtReferencePosition)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 0.0f, 0.0f, 0.0f } },
        45.0);

    EXPECT_EQ(GetSelectedPositions(point_list), GetPositions(point_list));
}

TEST(SampleFilterTest, IgnoresRejectPositionsAtReferencePositionWithZeroAngle)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 0.0f, 0.0f, 0.0f } },
        0.0);

    EXPECT_EQ(GetSelectedPositions(point_list), GetPositions(point_list));
}

TEST(SampleFilterTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 1.0f, 0.0f, 0.0f } },
        45.0);

    ASSERT_FALSE(point_list.empty());
    EXPECT_EQ(point_list.front().position, (std::array<float, 3>{ 0.0f, 0.0f, 0.0f }));
    EXPECT_TRUE(point_list.front().is_selected);
}

TEST(SampleFilterTest, RejectsInvalidAngles)
{
    auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            {},
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(SampleFilterTest, RejectsRejectPositionsWithNonFiniteValues)
{
    auto point_list{ MakePointList() };
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{
                1.0f,
                0.0f,
                std::numeric_limits<float>::quiet_NaN()
            } },
            45.0),
        std::invalid_argument);
}

TEST(SampleFilterTest, RejectsRejectPositionsWithNonFiniteValuesAtZeroAngle)
{
    auto point_list{ MakePointList() };
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0f, 0.0f, 0.0f },
            { std::array<float, 3>{
                1.0f,
                0.0f,
                std::numeric_limits<float>::quiet_NaN()
            } },
            0.0),
        std::invalid_argument);
}

TEST(SampleFilterTest, UsesReferencePositionForWorldSpaceSamplingPoints)
{
    SamplingPointList point_list{
        SamplingPoint{ 1.0f, { 10.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 11.0f, 0.0f, 0.0f } },
        SamplingPoint{ 1.0f, { 10.0f, 1.0f, 0.0f } }
    };

    sf::FilterSamplingPointList(
        point_list,
        { 10.0f, 0.0f, 0.0f },
        { std::array<float, 3>{ 11.0f, 0.0f, 0.0f } },
        30.0);

    const std::vector<std::array<float, 3>> expected_positions{
        { 10.0f, 0.0f, 0.0f },
        { 10.0f, 1.0f, 0.0f }
    };
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
}
