#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <stdexcept>

#include "SphereSampler.hpp"

TEST(SphereSamplerTest, DefaultConfiguration)
{
    SphereSampler sampler;
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples);
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

    EXPECT_NE(output.find("Sampling size: 10"), std::string::npos);
    EXPECT_NE(output.find("Distance range: [0.0, 1.0] Angstrom."), std::string::npos);
}

TEST(SphereSamplerTest, PrintOutputsHeader)
{
    SphereSampler sampler;

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("SphereSampler Configuration:"), std::string::npos);
}

TEST(SphereSamplerTest, PrintOutputsConfiguration)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(20);
    sampler.SetDistanceRangeMinimum(0.5);
    sampler.SetDistanceRangeMaximum(2.0);

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("Sampling size: 20"), std::string::npos);
    EXPECT_NE(output.find("Distance range: [0.5, 2.0] Angstrom."), std::string::npos);
}

TEST(SphereSamplerTest, Setters)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(5);
    sampler.SetDistanceRangeMinimum(1.0);
    sampler.SetDistanceRangeMaximum(1.0);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples);
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
    const std::array<float, 3> center{ 0.f, 0.f, 0.f };
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints(center, samples);
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        const auto & pos{ std::get<1>(sample) };
        EXPECT_GE(radius, 1.0f);
        EXPECT_LE(radius, 2.0f);
        if (radius > 0.0f)
        {
            float dx{ pos[0] - center[0] };
            float dy{ pos[1] - center[1] };
            float dz{ pos[2] - center[2] };
            float norm{ std::sqrt((dx / radius) * (dx / radius) +
                                  (dy / radius) * (dy / radius) +
                                  (dz / radius) * (dz / radius)) };
            EXPECT_NEAR(1.0f, norm, 1e-5f);
        }
    }
}

TEST(SphereSamplerTest, PositionMath)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(5);
    sampler.SetDistanceRangeMinimum(1.0);
    sampler.SetDistanceRangeMaximum(1.0);
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints(center, samples);
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

TEST(SphereSamplerTest, ZeroSamplingSizeReturnsEmptyVector)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(0);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples);
    EXPECT_TRUE(samples.empty());
}

TEST(SphereSamplerTest, ZeroDistanceRange)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(0.0);
    sampler.SetDistanceRangeMaximum(0.0);
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints(center, samples);
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

TEST(SphereSamplerTest, ThrowsWhenMinGreaterThanMax)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(2.0);
    sampler.SetDistanceRangeMaximum(1.0);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    EXPECT_THROW(sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples), std::invalid_argument);
}

TEST(SphereSamplerTest, ZeroSamplingSizeWithInvalidRangeThrows)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(0);
    sampler.SetDistanceRangeMinimum(2.0);
    sampler.SetDistanceRangeMaximum(1.0);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    EXPECT_THROW(sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples), std::invalid_argument);
}

TEST(SphereSamplerTest, ThrowsWhenMinNegative)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(-0.1);
    sampler.SetDistanceRangeMaximum(1.0);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    EXPECT_THROW(sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples), std::invalid_argument);
}

TEST(SphereSamplerTest, ThrowsWhenMaxNegative)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(0.0);
    sampler.SetDistanceRangeMaximum(-1.0);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    EXPECT_THROW(sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples), std::invalid_argument);
}

TEST(SphereSamplerTest, ValidRangeProducesSamples)
{
    SphereSampler sampler;
    sampler.SetDistanceRangeMinimum(0.2);
    sampler.SetDistanceRangeMaximum(0.4);
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}, samples);
    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_GE(radius, 0.2f);
        EXPECT_LE(radius, 0.4f);
    }
}

TEST(SphereSamplerTest, PositionMatchesRadiusForVariableRange)
{
    SphereSampler sampler;
    sampler.SetSamplingSize(20);
    sampler.SetDistanceRangeMinimum(0.5);
    sampler.SetDistanceRangeMaximum(1.5);
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    sampler.GenerateSamplingPoints(center, samples);
    ASSERT_EQ(20u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        const auto & pos{ std::get<1>(sample) };
        float dx{ pos[0] - center[0] };
        float dy{ pos[1] - center[1] };
        float dz{ pos[2] - center[2] };
        float distance{ std::sqrt(dx * dx + dy * dy + dz * dz) };
        EXPECT_NEAR(radius, distance, 1e-5f);
    }
}
