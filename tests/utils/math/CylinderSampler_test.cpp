#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/math/CylinderSampler.hpp>

TEST(CylinderSamplerTest, ProducesConfiguredNumberOfSamples)
{
    CylinderSampler sampler;
    sampler.SetSampleCount(12);
    sampler.SetDistanceRangeMinimum(1.0);
    sampler.SetDistanceRangeMaximum(2.0);
    sampler.SetHeight(3.0);

    const auto samples{
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 1.0f })
    };

    ASSERT_EQ(12u, samples.size());
    EXPECT_EQ(12u, sampler.GetSampleCount());
    EXPECT_EQ(12u, sampler.GetPointCount());
}

TEST(CylinderSamplerTest, ThrowsWhenAxisVectorIsZero)
{
    CylinderSampler sampler;
    EXPECT_THROW(
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }),
        std::invalid_argument);
}

TEST(CylinderSamplerTest, ThrowsWhenRadiusRangeIsInvalid)
{
    CylinderSampler sampler;
    sampler.SetDistanceRangeMinimum(2.0);
    sampler.SetDistanceRangeMaximum(1.0);

    EXPECT_THROW(
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }),
        std::invalid_argument);
}
