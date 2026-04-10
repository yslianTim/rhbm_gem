#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace {

float ComputeDistance(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    const float dx{ lhs[0] - rhs[0] };
    const float dy{ lhs[1] - rhs[1] };
    const float dz{ lhs[2] - rhs[2] };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void ExpectSamplesEqual(
    const SamplingPointList & lhs,
    const SamplingPointList & rhs)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        EXPECT_FLOAT_EQ(std::get<0>(lhs[i]), std::get<0>(rhs[i]));
        const auto & lhs_position{ std::get<1>(lhs[i]) };
        const auto & rhs_position{ std::get<1>(rhs[i]) };
        EXPECT_FLOAT_EQ(lhs_position[0], rhs_position[0]);
        EXPECT_FLOAT_EQ(lhs_position[1], rhs_position[1]);
        EXPECT_FLOAT_EQ(lhs_position[2], rhs_position[2]);
    }
}

} // namespace

TEST(SphereSamplerTest, DefaultConfiguration)
{
    SphereSampler sampler;
    const auto & profile{ sampler.GetSamplingProfile() };
    auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    EXPECT_EQ(SphereSamplingMethod::RadiusUniformRandom, profile.GetMethod());
    EXPECT_DOUBLE_EQ(0.0, profile.GetDistanceRange().min);
    EXPECT_DOUBLE_EQ(1.0, profile.GetDistanceRange().max);
    EXPECT_EQ(10u, profile.GetRandomConfig().sample_count);
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

TEST(SphereSamplerTest, PrintOutputsFibonacciConfiguration)
{
    SphereSampler sampler;
    constexpr double radius_bin_size{ 0.3 };
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            radius_bin_size,
            4));

    testing::internal::CaptureStdout();
    sampler.Print();
    std::string output{ testing::internal::GetCapturedStdout() };
    const auto expected_radius_bin_size{
        StringHelper::ToStringWithPrecision<double>(radius_bin_size, 2) };

    EXPECT_NE(output.find("Sampling method: FibonacciDeterministic"), std::string::npos);
    EXPECT_NE(
        output.find("Radius bin size: " + expected_radius_bin_size + " Angstrom."),
        std::string::npos);
    EXPECT_NE(output.find("Samples per radius: 4"), std::string::npos);
    EXPECT_NE(output.find("Expected point count: 12"), std::string::npos);
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

TEST(SphereSamplerTest, SetSamplingProfileAcceptsFactoryCreatedProfile)
{
    SphereSampler sampler;
    const auto profile{
        SphereSamplingProfile::VolumeUniformRandom(
            SphereDistanceRange{ 0.5, 2.0 },
            12)
    };

    sampler.SetSamplingProfile(profile);

    EXPECT_EQ(SphereSamplingMethod::VolumeUniformRandom, sampler.GetSamplingProfile().GetMethod());
    EXPECT_DOUBLE_EQ(0.5, sampler.GetSamplingProfile().GetDistanceRange().min);
    EXPECT_DOUBLE_EQ(2.0, sampler.GetSamplingProfile().GetDistanceRange().max);
    EXPECT_EQ(12u, sampler.GetExpectedPointCount());
}

TEST(SphereSamplerTest, SetSamplingProfileAcceptsFibonacciProfile)
{
    SphereSampler sampler;
    const auto profile{
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.25,
            7)
    };

    sampler.SetSamplingProfile(profile);

    EXPECT_EQ(
        SphereSamplingMethod::FibonacciDeterministic,
        sampler.GetSamplingProfile().GetMethod());
    EXPECT_DOUBLE_EQ(0.5, sampler.GetSamplingProfile().GetDistanceRange().min);
    EXPECT_DOUBLE_EQ(1.0, sampler.GetSamplingProfile().GetDistanceRange().max);
    EXPECT_DOUBLE_EQ(0.25, sampler.GetSamplingProfile().GetFibonacciConfig().radius_bin_size);
    EXPECT_EQ(7u, sampler.GetSamplingProfile().GetFibonacciConfig().samples_per_radius);
    EXPECT_EQ(21u, sampler.GetExpectedPointCount());
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

TEST(SphereSamplerTest, FibonacciDeterministicProfileProducesDeterministicSamples)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.3,
            4));
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };

    const auto first_samples{ sampler.GenerateSamplingPoints(center) };
    const auto second_samples{ sampler.GenerateSamplingPoints(center) };

    ExpectSamplesEqual(first_samples, second_samples);
}

TEST(SphereSamplerTest, FibonacciDeterministicSingleShellProducesConfiguredPointCount)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 1.25, 1.25 },
            0.4,
            6));
    const std::array<float, 3> center{ 0.f, 0.f, 0.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    EXPECT_EQ(6u, sampler.GetExpectedPointCount());
    ASSERT_EQ(6u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_FLOAT_EQ(1.25f, std::get<0>(sample));
        EXPECT_NEAR(1.25f, ComputeDistance(center, std::get<1>(sample)), 1e-5f);
    }
}

TEST(SphereSamplerTest, FibonacciDeterministicProducesExpectedShellPointCount)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.3,
            4));
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };
    const std::vector<float> expected_shells{ 0.5f, 0.8f, 1.0f };

    EXPECT_EQ(12u, sampler.GetExpectedPointCount());
    ASSERT_EQ(12u, samples.size());
    for (std::size_t shell_index = 0; shell_index < expected_shells.size(); ++shell_index)
    {
        for (std::size_t point_index = 0; point_index < 4; ++point_index)
        {
            const auto & sample{ samples[shell_index * 4 + point_index] };
            const float radius{ std::get<0>(sample) };
            EXPECT_NEAR(expected_shells[shell_index], radius, 1e-5f);
            EXPECT_NEAR(radius, ComputeDistance(center, std::get<1>(sample)), 1e-5f);
        }
    }
}

TEST(SphereSamplerTest, FibonacciDeterministicZeroRadiusShellReturnsCenterPoints)
{
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 0.0 },
            0.5,
            5));
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    EXPECT_EQ(5u, sampler.GetExpectedPointCount());
    ASSERT_EQ(5u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_FLOAT_EQ(0.0f, std::get<0>(sample));
        const auto & position{ std::get<1>(sample) };
        EXPECT_FLOAT_EQ(center[0], position[0]);
        EXPECT_FLOAT_EQ(center[1], position[1]);
        EXPECT_FLOAT_EQ(center[2], position[2]);
    }
}

TEST(SphereSamplerTest, FibonacciDeterministicProfileThrowsWhenRadiusBinSizeNonPositive)
{
    EXPECT_THROW(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 1.0 },
            0.0,
            4),
        std::invalid_argument);
}

TEST(SphereSamplerTest, FibonacciDeterministicProfileThrowsWhenSamplesPerRadiusZero)
{
    EXPECT_THROW(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 1.0 },
            0.5,
            0),
        std::invalid_argument);
}

TEST(SphereSamplerTest, FibonacciDeterministicProfileThrowsWhenRangeMinNegative)
{
    EXPECT_THROW(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ -0.1, 1.0 },
            0.5,
            4),
        std::invalid_argument);
}

TEST(SphereSamplerTest, FibonacciDeterministicProfileThrowsWhenRangeMinGreaterThanMax)
{
    EXPECT_THROW(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ 2.0, 1.0 },
            0.5,
            4),
        std::invalid_argument);
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
