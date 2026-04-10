#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/math/GridSampler.hpp>

TEST(GridSamplerTest, PointCountTracksGridResolution)
{
    GridSampler sampler;
    sampler.SetGridResolution(4);
    sampler.SetWindowSize(6.0f);

    const auto samples{
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 1.0f })
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
        (void)sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f });
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
    sampler.SetWindowSize(0.0f);

    EXPECT_THROW(
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }),
        std::invalid_argument);
}
