#include <gtest/gtest.h>

#include <memory>

#include <rhbm_gem/data/object/MapObject.hpp>
#include "support/DataObjectTestSupport.hpp"

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

TEST(DataObjectMapBehaviorTest, SetMapValueArrayRefreshesStatistics)
{
    auto map{ data_test::MakeMapObject() };
    auto replacement_values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        replacement_values[i] = static_cast<float>(10 + i);
    }

    map.SetMapValueArray(std::move(replacement_values));

    EXPECT_FLOAT_EQ(map.GetMapValueMin(), 10.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMax(), 17.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMean(), 13.5f);
}

TEST(DataObjectMapBehaviorTest, ClearMapValueArrayResetsValuesAndStatistics)
{
    auto map{ data_test::MakeMapObject() };
    ASSERT_GT(map.GetMapValueMax(), 0.0f);

    map.ClearMapValueArray();

    for (size_t i = 0; i < map.GetMapValueArraySize(); ++i)
    {
        EXPECT_FLOAT_EQ(map.GetMapValue(i), 0.0f);
    }
    EXPECT_FLOAT_EQ(map.GetMapValueMin(), 0.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMax(), 0.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMean(), 0.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueSD(), 0.0f);
}
