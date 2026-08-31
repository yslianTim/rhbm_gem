#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

namespace {

double ComputeDistance(
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs)
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
        EXPECT_DOUBLE_EQ(lhs[i].distance, rhs[i].distance);
        EXPECT_EQ(lhs[i].position, rhs[i].position);
    }
}

std::size_t CountSamplesAtRadius(
    const SamplingPointList & samples,
    double expected_radius)
{
    std::size_t count{ 0 };
    for (const auto & sample : samples)
    {
        if (std::abs(sample.distance - expected_radius) <= 1e-5)
        {
            count++;
        }
    }

    return count;
}

void ExpectAnalysisSampleDistances(
    const SamplingPointList & samples,
    const std::array<double, 3> & center)
{
    for (const auto & sample : samples)
    {
        EXPECT_GE(sample.distance, 0.0);
        EXPECT_LE(sample.distance, 2.0);
        EXPECT_NEAR(sample.distance, ComputeDistance(center, sample.position), 1e-5);
    }
}

} // namespace

TEST(SphereSamplerTest, RadiusUniformProducesFixedAnalysisSamples)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };
    const auto samples{ rhbm_gem::sphere_sampler::GenerateRadiusUniformRandom(center) };

    ASSERT_EQ(10u, samples.size());
    ExpectAnalysisSampleDistances(samples, center);
}

TEST(SphereSamplerTest, VolumeUniformProducesFixedAnalysisSamples)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };
    const auto samples{ rhbm_gem::sphere_sampler::GenerateVolumeUniformRandom(center) };

    ASSERT_EQ(10u, samples.size());
    ExpectAnalysisSampleDistances(samples, center);
}

TEST(SphereSamplerTest, FibonacciProducesFixedDeterministicShells)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };

    const auto first_samples{ rhbm_gem::sphere_sampler::GenerateFibonacciDeterministic(center) };
    const auto second_samples{ rhbm_gem::sphere_sampler::GenerateFibonacciDeterministic(center) };

    ASSERT_EQ(200u, first_samples.size());
    EXPECT_EQ(0u, CountSamplesAtRadius(first_samples, 0.0));
    EXPECT_EQ(10u, CountSamplesAtRadius(first_samples, 0.05));
    EXPECT_EQ(10u, CountSamplesAtRadius(first_samples, 1.95));
    ExpectAnalysisSampleDistances(first_samples, center);
    ExpectSamplesEqual(first_samples, second_samples);
}

TEST(SphereSamplerTest, DispatchThrowsForUnsupportedMethod)
{
    EXPECT_THROW(
        rhbm_gem::sphere_sampler::GenerateSamplingPointList(
            { 0.0, 0.0, 0.0 },
            static_cast<SphereSamplingMethod>(99)),
        std::invalid_argument);
}

TEST(SphereSamplerTest, DispatchUsesRequestedMethod)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };
    const auto samples{
        rhbm_gem::sphere_sampler::GenerateSamplingPointList(
            center,
            SphereSamplingMethod::FibonacciDeterministic)
    };

    ASSERT_EQ(200u, samples.size());
    EXPECT_EQ(10u, CountSamplesAtRadius(samples, 0.05));
    ExpectAnalysisSampleDistances(samples, center);
}
