#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "data/detail/MapSpatialIndex.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectMapBehaviorTest, NormalizeMapObjectNormalizesMapValues)
{
    auto map{ data_test::MakeMapObject() };
    const auto original_value{ map.GetMapValue(0) };
    const auto original_sd{ map.GetMapValueSD() };
    ASSERT_GT(original_sd, 0.0f);

    map.MapValueArrayNormalization();

    EXPECT_NEAR(map.GetMapValue(0), original_value / original_sd, 1.0e-5f);
    EXPECT_NEAR(map.GetMapValueSD(), 1.0f, 1.0e-5f);
}

TEST(DataObjectMapBehaviorTest, SetMapValueArrayRefreshesStatisticsAndPreservesSpatialQueries)
{
    auto map{ data_test::MakeMapObject() };
    std::vector<size_t> initial_hits;
    rg::MapSpatialIndex(map).CollectGridIndicesInRange(
        map.GetGridPosition(0),
        0.1f,
        initial_hits);
    ASSERT_FALSE(initial_hits.empty());

    auto replacement_values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        replacement_values[i] = static_cast<float>(10 + i);
    }

    map.SetMapValueArray(std::move(replacement_values));

    std::vector<size_t> refreshed_hits;
    rg::MapSpatialIndex(map).CollectGridIndicesInRange(
        map.GetGridPosition(0),
        0.1f,
        refreshed_hits);
    EXPECT_EQ(refreshed_hits, initial_hits);
    EXPECT_FLOAT_EQ(map.GetMapValueMin(), 10.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMax(), 17.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMean(), 13.5f);
}
