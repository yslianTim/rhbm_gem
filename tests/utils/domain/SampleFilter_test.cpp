#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <limits>
#include <utility>
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

TEST(SampleFilterTest, UsesFifteenDegreeDefaultRejectAngle)
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
        { 0.0f, 0.0f, 0.0f },
        { 1.0f, 0.3f, 0.0f }
    };
    EXPECT_EQ(GetSelectedPositions(point_list), expected_positions);
    EXPECT_FALSE(point_list.at(1).is_selected);
    EXPECT_TRUE(point_list.at(2).is_selected);
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

TEST(SampleFilterTest, FilterLocalPotentialSampleListHandlesEmptySampleList)
{
    EXPECT_TRUE(sf::FilterLocalPotentialSampleList({}).empty());
}

TEST(SampleFilterTest, FilterLocalPotentialSampleListUsesDefaultRetainedRatio)
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
        sf::FilterLocalPotentialSampleList(std::move(sample_list))
    };

    ASSERT_EQ(retained_samples.size(), 2u);
    EXPECT_FLOAT_EQ(retained_samples.at(0).response, 0.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).response, 1.0f);
}

TEST(SampleFilterTest, FilterLocalPotentialSampleListFiltersByDistance)
{
    LocalPotentialSampleList sample_list;
    for (std::size_t i = 0; i < 11u; i++)
    {
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(10u - i),
            SamplingPoint{ 1.0f }
        });
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(20u - i),
            SamplingPoint{ 2.0f }
        });
    }

    const auto retained_samples{
        sf::FilterLocalPotentialSampleList(std::move(sample_list))
    };

    ASSERT_EQ(retained_samples.size(), 4u);
    EXPECT_FLOAT_EQ(retained_samples.at(0).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(0).response, 0.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).point.distance, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(1).response, 1.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(2).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(2).response, 10.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(3).point.distance, 2.0f);
    EXPECT_FLOAT_EQ(retained_samples.at(3).response, 11.0f);
}
