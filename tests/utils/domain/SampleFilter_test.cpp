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
        SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 } },
        SamplingPoint{ 1.0, { 1.0, 0.0, 0.0 } },
        SamplingPoint{ 1.0, { 0.0, 1.0, 0.0 } },
        SamplingPoint{ 1.0, { -1.0, 0.0, 0.0 } }
    };
}

std::vector<std::array<double, 3>> GetPositions(const SamplingPointList & point_list)
{
    std::vector<std::array<double, 3>> positions;
    positions.reserve(point_list.size());
    for (const auto & point : point_list)
    {
        positions.emplace_back(point.position);
    }
    return positions;
}

std::vector<std::array<double, 3>> GetSelectedPositions(const SamplingPointList & point_list)
{
    std::vector<std::array<double, 3>> positions;
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
        LocalPotentialSample{ 100.0, SamplingPoint{ 1.2 } },
        LocalPotentialSample{ 2.0, SamplingPoint{ 1.0 } },
        LocalPotentialSample{ 1.0, SamplingPoint{ 1.0 } },
        LocalPotentialSample{ 3.0, SamplingPoint{ 1.0 } },
        LocalPotentialSample{ 10.0, SamplingPoint{ 1.2 } },
        LocalPotentialSample{ 11.0, SamplingPoint{ 1.2 } }
    };

    const auto actual{ sf::BuildMedianResponseSampleEntriesByRadius(sample_entries) };

    ASSERT_EQ(actual.size(), 2u);
    EXPECT_DOUBLE_EQ(actual.at(0).point.distance, 1.0);
    EXPECT_DOUBLE_EQ(actual.at(0).response, 2.0);
    EXPECT_DOUBLE_EQ(actual.at(1).point.distance, 1.2);
    EXPECT_DOUBLE_EQ(actual.at(1).response, 11.0);
}

TEST(SampleFilterTest, BuildsEmptyMedianResponseSampleEntriesByRadius)
{
    const LocalPotentialSampleList sample_entries;

    const auto actual{ sf::BuildMedianResponseSampleEntriesByRadius(sample_entries) };

    EXPECT_TRUE(actual.empty());
}

TEST(SampleFilterTest, EmptyRejectPositionsSelectAllSamplingPoints)
{
    auto point_list{ MakePointList() };
    point_list.at(1).is_selected = false;

    sf::FilterSamplingPointList(
        point_list,
        { 0.0, 0.0, 0.0 },
        {});

    EXPECT_EQ(GetSelectedPositions(point_list), GetPositions(point_list));
    EXPECT_TRUE(AllSelected(point_list));
}

TEST(SampleFilterTest, RejectsPointsCloserToNearestNeighborThanReference)
{
    SamplingPointList point_list{
        SamplingPoint{ 0.05, { 0.05, 0.0, 0.0 } },
        SamplingPoint{ 0.15, { 0.15, 0.0, 0.0 } },
        SamplingPoint{ 0.25, { 0.25, 0.0, 0.0 } }
    };

    sf::FilterSamplingPointList(
        point_list,
        { 0.0, 0.0, 0.0 },
        {
            std::array<double, 3>{ 10.0, 0.0, 0.0 },
            std::array<double, 3>{ 0.2, 0.0, 0.0 }
        },
        0.0);

    const std::vector<std::array<double, 3>> expected_positions{
        { 0.05, 0.0, 0.0 }
    };
    EXPECT_EQ(point_list.size(), 3u);
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
    EXPECT_FALSE(point_list.at(1).is_selected);
    EXPECT_FALSE(point_list.at(2).is_selected);
}

TEST(SampleFilterTest, UsesThirtyDegreeDefaultRejectAngle)
{
    SamplingPointList point_list{
        SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 } },
        SamplingPoint{ 1.0, { 1.0, 0.2, 0.0 } },
        SamplingPoint{ 1.0, { 1.0, 0.3, 0.0 } }
    };

    sf::FilterSamplingPointList(
        point_list,
        { 0.0, 0.0, 0.0 },
        { std::array<double, 3>{ 2.0, 0.0, 0.0 } });

    const std::vector<std::array<double, 3>> expected_positions{
        { 0.0, 0.0, 0.0 }
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
        { 0.0, 0.0, 0.0 },
        { std::array<double, 3>{ 1.0, 0.0, 0.0 } },
        30.0);

    const std::vector<std::array<double, 3>> expected_positions{
        { 0.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { -1.0, 0.0, 0.0 }
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
        { 0.0, 0.0, 0.0 },
        { std::array<double, 3>{ 1.0, 0.0, 0.0 } },
        89.0);

    const std::vector<std::array<double, 3>> expected_positions{
        { 0.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { -1.0, 0.0, 0.0 }
    };
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
}

TEST(SampleFilterTest, IgnoresRejectPositionsAtReferencePosition)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0, 0.0, 0.0 },
        { std::array<double, 3>{ 0.0, 0.0, 0.0 } },
        45.0);

    EXPECT_EQ(GetSelectedPositions(point_list), GetPositions(point_list));
}

TEST(SampleFilterTest, IgnoresRejectPositionsAtReferencePositionWithZeroAngle)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0, 0.0, 0.0 },
        { std::array<double, 3>{ 0.0, 0.0, 0.0 } },
        0.0);

    EXPECT_EQ(GetSelectedPositions(point_list), GetPositions(point_list));
}

TEST(SampleFilterTest, KeepsOriginSamplingPointEvenWhenFilteringEnabled)
{
    auto point_list{ MakePointList() };
    sf::FilterSamplingPointList(
        point_list,
        { 0.0, 0.0, 0.0 },
        { std::array<double, 3>{ 1.0, 0.0, 0.0 } },
        45.0);

    ASSERT_FALSE(point_list.empty());
    EXPECT_EQ(point_list.front().position, (std::array<double, 3>{ 0.0, 0.0, 0.0 }));
    EXPECT_TRUE(point_list.front().is_selected);
}

TEST(SampleFilterTest, RejectsInvalidAngles)
{
    auto point_list{ MakePointList() };

    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0, 0.0, 0.0 },
            {},
            -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0, 0.0, 0.0 },
            {},
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0, 0.0, 0.0 },
            {},
            181.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)sf::FilterSamplingPointList(
            point_list,
            { 0.0, 0.0, 0.0 },
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
            { 0.0, 0.0, 0.0 },
            { std::array<double, 3>{
                1.0,
                0.0,
                std::numeric_limits<double>::quiet_NaN()
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
            { 0.0, 0.0, 0.0 },
            { std::array<double, 3>{
                1.0,
                0.0,
                std::numeric_limits<double>::quiet_NaN()
            } },
            0.0),
        std::invalid_argument);
}

TEST(SampleFilterTest, UsesReferencePositionForWorldSpaceSamplingPoints)
{
    SamplingPointList point_list{
        SamplingPoint{ 1.0, { 10.0, 0.0, 0.0 } },
        SamplingPoint{ 1.0, { 11.0, 0.0, 0.0 } },
        SamplingPoint{ 1.0, { 10.0, 1.0, 0.0 } }
    };

    sf::FilterSamplingPointList(
        point_list,
        { 10.0, 0.0, 0.0 },
        { std::array<double, 3>{ 11.0, 0.0, 0.0 } },
        30.0);

    const std::vector<std::array<double, 3>> expected_positions{
        { 10.0, 0.0, 0.0 },
        { 10.0, 1.0, 0.0 }
    };
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
}
