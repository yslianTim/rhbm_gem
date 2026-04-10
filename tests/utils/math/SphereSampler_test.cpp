#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <stdexcept>

#include <rhbm_gem/utils/math/SphereSampler.hpp>

TEST(SphereSamplerTest, DefaultConfiguration)
{
    SphereSampler sampler;
    const auto & profile{ sampler.GetSamplingProfile() };
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    EXPECT_EQ(SphereSamplingMethod::RadiusUniformRandom, profile.method);
    EXPECT_DOUBLE_EQ(0.0, profile.distance_range.min);
    EXPECT_DOUBLE_EQ(1.0, profile.distance_range.max);
    EXPECT_EQ(10u, std::get<SphereRandomSamplingConfig>(profile.method_config).sample_count);
    EXPECT_EQ(10u, sampler.GetExpectedPointCount());
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

    EXPECT_NE(output.find("Sampling method: RadiusUniformRandom"), std::string::npos);
    EXPECT_NE(output.find("Sample count: 10"), std::string::npos);
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
    sampler.SetSamplingProfile(
        SphereSamplingProfile::VolumeUniformRandom(
            SphereDistanceRange{ 0.5, 2.0 },
            20));

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("Sampling method: VolumeUniformRandom"), std::string::npos);
    EXPECT_NE(output.find("Sample count: 20"), std::string::npos);
    EXPECT_NE(output.find("Distance range: [0.5, 2.0] Angstrom."), std::string::npos);
}

TEST(SphereSamplerTest, RadiusUniformRandomProfileProducesConfiguredSampleCount)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 1.0, 1.0 },
            5));
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    EXPECT_EQ(5u, sampler.GetExpectedPointCount());
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
    sampler.SetSamplingProfile(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 1.0, 2.0 },
            10));
    const std::array<float, 3> center{ 0.f, 0.f, 0.f };
    auto samples{ sampler.GenerateSamplingPoints(center) };
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
    sampler.SetSamplingProfile(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 1.0, 1.0 },
            5));
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

TEST(SphereSamplerTest, ZeroSampleCountReturnsEmptyVector)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 0.0, 1.0 },
            0));
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    EXPECT_EQ(0u, sampler.GetExpectedPointCount());
    EXPECT_TRUE(samples.empty());
}

TEST(SphereSamplerTest, ZeroDistanceRange)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 0.0, 0.0 },
            10));
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

TEST(SphereSamplerTest, RadiusUniformRandomProfileThrowsWhenMinGreaterThanMax)
{
    EXPECT_THROW(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 2.0, 1.0 },
            10),
        std::invalid_argument);
}

TEST(SphereSamplerTest, SetSamplingProfileThrowsWhenRangeIsInvalid)
{
    SphereSampler sampler;
    SphereSamplingProfile invalid_profile;
    invalid_profile.method = SphereSamplingMethod::RadiusUniformRandom;
    invalid_profile.distance_range = SphereDistanceRange{ 2.0, 1.0 };
    invalid_profile.method_config = SphereRandomSamplingConfig{ 0 };

    EXPECT_THROW(sampler.SetSamplingProfile(invalid_profile), std::invalid_argument);
}

TEST(SphereSamplerTest, RadiusUniformRandomProfileThrowsWhenMinNegative)
{
    EXPECT_THROW(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ -0.1, 1.0 },
            10),
        std::invalid_argument);
}

TEST(SphereSamplerTest, RadiusUniformRandomProfileThrowsWhenMaxNegative)
{
    EXPECT_THROW(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 0.0, -1.0 },
            10),
        std::invalid_argument);
}

TEST(SphereSamplerTest, RadiusUniformRandomProfileProducesSamples)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::RadiusUniformRandom(
            SphereDistanceRange{ 0.2, 0.4 },
            10));
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        float radius{ std::get<0>(sample) };
        EXPECT_GE(radius, 0.2f);
        EXPECT_LE(radius, 0.4f);
    }
}

TEST(SphereSamplerTest, VolumeUniformRandomProfileProducesConfiguredSampleCount)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::VolumeUniformRandom(
            SphereDistanceRange{ 0.5, 1.5 },
            20));
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    auto samples{ sampler.GenerateSamplingPoints(center) };

    EXPECT_EQ(20u, sampler.GetExpectedPointCount());
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
