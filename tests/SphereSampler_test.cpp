#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <stdexcept>

#include "SphereSampler.hpp"

TEST(SphereSamplerTest, DefaultConfiguration)
{
    SphereSampler sampler;
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_GE(radius, 0.0f);
        EXPECT_LE(radius, 1.0f);
    }
}

TEST(SphereSamplerTest, PrintOutputsDefaultConfiguration)
{
    SphereSampler sampler;

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("Thread size: 1"), std::string::npos);
    EXPECT_NE(output.find("Sampling size: 10"), std::string::npos);
    EXPECT_NE(output.find("Distance range: [0.0, 1.0] Angstrom."), std::string::npos);
}

TEST(SphereSamplerTest, PrintOutputsConfiguration)
{
    SphereSampler sampler;
    sampler.SetThreadSize(4);
    sampler.SetSamplingSize(20);
    sampler.SetDistanceRangeMinimum(0.5);
    sampler.SetDistanceRangeMaximum(2.0);

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("Thread size: 4"), std::string::npos);
    EXPECT_NE(output.find("Sampling size: 20"), std::string::npos);
    EXPECT_NE(output.find("Distance range: [0.5, 2.0] Angstrom."), std::string::npos);
}

TEST(SphereSamplerTest, Setters)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(5);
    sampler.SetDistanceRangeMinimum(1.0);
    sampler.SetDistanceRangeMaximum(1.0);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(5u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_NEAR(1.0f, radius, 1e-5f);
    }
}

TEST(SphereSamplerTest, RadiusWithinConfiguredRange)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(1.0);
    sampler.SetDistanceRangeMaximum(2.0);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_GE(radius, 1.0f);
        EXPECT_LE(radius, 2.0f);
    }
}

TEST(SphereSamplerTest, PositionMath)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(5);
    sampler.SetDistanceRangeMinimum(1.0);
    sampler.SetDistanceRangeMaximum(1.0);
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    auto samples{ sampler.GenerateSamplingPoints(center) };
    for (const auto & sample : samples)
    {
        float radius{std::get<0>(sample)};
        const auto & pos{ std::get<1>(sample) };
        float dx{pos[0] - center[0]};
        float dy{pos[1] - center[1]};
        float dz{pos[2] - center[2]};
        float distance{ std::sqrt(dx * dx + dy * dy + dz * dz) };
        EXPECT_NEAR(radius, distance, 1e-5f);
    }
}

TEST(SphereSamplerTest, ValidThreadSize)
{
    SphereSampler sampler;
    sampler.SetThreadSize(4);
    sampler.SetSamplingSize(5);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(5u, samples.size());
}

TEST(SphereSamplerTest, ZeroThreadSizeClampedToOne)
{
    SphereSampler sampler;
    sampler.SetThreadSize(0);

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };
    EXPECT_NE(output.find("Thread size: 1"), std::string::npos);

    sampler.SetSamplingSize(5);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(5u, samples.size());
}

TEST(SphereSamplerTest, ZeroSamplingSizeReturnsEmptyVector)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(0);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    EXPECT_TRUE(samples.empty());
}

TEST(SphereSamplerTest, ZeroDistanceRange)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(0.0);
    sampler.SetDistanceRangeMaximum(0.0);
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    auto samples{ sampler.GenerateSamplingPoints(center) };
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        const auto & pos{ std::get<1>(sample) };
        EXPECT_FLOAT_EQ(0.0f, radius);
        EXPECT_FLOAT_EQ(center[0], pos[0]);
        EXPECT_FLOAT_EQ(center[1], pos[1]);
        EXPECT_FLOAT_EQ(center[2], pos[2]);
    }
}

TEST(SphereSamplerTest, ThreadSizeDoesNotAffectOutput)
{
    SphereSampler sampler;
    sampler.SetThreadSize(2);
    sampler.SetSamplingSize(5);
    sampler.SetDistanceRangeMinimum(0.5);
    sampler.SetDistanceRangeMaximum(1.0);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(5u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_GE(radius, 0.5f);
        EXPECT_LE(radius, 1.0f);
    }
}

TEST(SphereSamplerTest, ThrowsWhenMinGreaterThanMax)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(2.0);
    sampler.SetDistanceRangeMaximum(1.0);
    EXPECT_THROW(sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}), std::invalid_argument);
}

TEST(SphereSamplerTest, ThrowsWhenMinNegative)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(-0.1);
    sampler.SetDistanceRangeMaximum(1.0);
    EXPECT_THROW(sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}), std::invalid_argument);
}

TEST(SphereSamplerTest, ValidRangeProducesSamples)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(0.2);
    sampler.SetDistanceRangeMaximum(0.4);
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_GE(radius, 0.2f);
        EXPECT_LE(radius, 0.4f);
    }
}
