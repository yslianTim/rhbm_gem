#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/math/CylinderSampler.hpp>

TEST(CylinderSamplerTest, ProducesConfiguredNumberOfSamples)
{
    CylinderSampler sampler;
    sampler.SetSampleCount(12);
    sampler.SetDistanceRange(1.0, 2.0);
    sampler.SetHeight(3.0);

    const auto samples{
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 1.0f })
    };

    ASSERT_EQ(12u, samples.size());
    EXPECT_EQ(12u, sampler.GetSampleCount());
}

TEST(CylinderSamplerTest, ThrowsWhenAxisVectorIsZero)
{
    CylinderSampler sampler;
    EXPECT_THROW(
        sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }),
        std::invalid_argument);
}

TEST(CylinderSamplerTest, SetDistanceRangeThrowsWhenRadiusRangeIsInvalid)
{
    CylinderSampler sampler;
    EXPECT_THROW(sampler.SetDistanceRange(2.0, 1.0), std::invalid_argument);
}

TEST(CylinderSamplerTest, SetDistanceRangeThrowsWhenRadiusRangeIsNegative)
{
    CylinderSampler sampler;
    EXPECT_THROW(sampler.SetDistanceRange(-1.0, 1.0), std::invalid_argument);
}
