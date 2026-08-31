#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/math/GridSampler.hpp>

TEST(GridSamplerTest, PointCountTracksGridResolution)
{
    GridSampler sampler;
    sampler.SetGridResolution(4);
    sampler.SetWindowSize(6.0);

    const auto samples{
        sampler.GenerateSamplingPoints({ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 })
    };

    EXPECT_EQ(4u, sampler.GetGridResolution());
    EXPECT_EQ(16u, sampler.GetPointCount());
    ASSERT_EQ(16u, samples.size());
}

TEST(GridSamplerTest, ThrowsWhenPlaneNormalIsZero)
{
    GridSampler sampler;

    try
    {
        (void)sampler.GenerateSamplingPoints({ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 });
        FAIL() << "Expected std::invalid_argument";
    }
    catch (const std::invalid_argument & ex)
    {
        EXPECT_STREQ("GridSampler: plane normal cannot be zero.", ex.what());
    }
}

TEST(GridSamplerTest, ThrowsWhenWindowSizeIsNotPositive)
{
    GridSampler sampler;
    sampler.SetWindowSize(0.0);

    EXPECT_THROW(
        sampler.GenerateSamplingPoints({ 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }),
        std::invalid_argument);
}
