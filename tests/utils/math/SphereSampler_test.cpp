#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

namespace {

float ComputeDistance(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    return rhbm_gem::array_helper::ComputeNorm(lhs, rhs);
}

void ExpectSamplesEqual(
    const SamplingPointList & lhs,
    const SamplingPointList & rhs)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        EXPECT_FLOAT_EQ(lhs[i].distance, rhs[i].distance);
        EXPECT_EQ(lhs[i].position, rhs[i].position);
    }
}

std::size_t CountSamplesAtRadius(
    const SamplingPointList & samples,
    float expected_radius)
{
    std::size_t count{ 0 };
    for (const auto & sample : samples)
    {
        if (std::abs(sample.distance - expected_radius) <= 1e-5f)
        {
            count++;
        }
    }

    return count;
}

void ExpectAnalysisSampleDistances(
    const SamplingPointList & samples,
    const std::array<float, 3> & center)
{
    for (const auto & sample : samples)
    {
        EXPECT_GE(sample.distance, 0.0f);
        EXPECT_LE(sample.distance, 1.5f);
        EXPECT_NEAR(sample.distance, ComputeDistance(center, sample.position), 1e-5f);
    }
}

} // namespace

TEST(SphereSamplerTest, ConstructorRadiusUniformProducesFixedAnalysisSamples)
{
    const SphereSampler sampler{ SphereSamplingMethod::RadiusUniformRandom };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(50u, samples.size());
    ExpectAnalysisSampleDistances(samples, center);
}

TEST(SphereSamplerTest, ConstructorVolumeUniformProducesFixedAnalysisSamples)
{
    const SphereSampler sampler{ SphereSamplingMethod::VolumeUniformRandom };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(50u, samples.size());
    ExpectAnalysisSampleDistances(samples, center);
}

TEST(SphereSamplerTest, ConstructorFibonacciProducesFixedDeterministicShells)
{
    const SphereSampler sampler{ SphereSamplingMethod::FibonacciDeterministic };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };

    const auto first_samples{ sampler.GenerateSamplingPoints(center) };
    const auto second_samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(750u, first_samples.size());
    EXPECT_EQ(0u, CountSamplesAtRadius(first_samples, 0.0f));
    EXPECT_EQ(50u, CountSamplesAtRadius(first_samples, 0.05f));
    EXPECT_EQ(50u, CountSamplesAtRadius(first_samples, 1.45f));
    ExpectAnalysisSampleDistances(first_samples, center);
    ExpectSamplesEqual(first_samples, second_samples);
}

TEST(SphereSamplerTest, ConstructorThrowsForUnsupportedMethod)
{
    EXPECT_THROW(
        SphereSampler{ static_cast<SphereSamplingMethod>(99) },
        std::invalid_argument);
}
