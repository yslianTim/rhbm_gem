#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
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
        const auto & lhs_position{ lhs[i].position };
        const auto & rhs_position{ rhs[i].position };
        EXPECT_FLOAT_EQ(lhs_position[0], rhs_position[0]);
        EXPECT_FLOAT_EQ(lhs_position[1], rhs_position[1]);
        EXPECT_FLOAT_EQ(lhs_position[2], rhs_position[2]);
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

void ExpectDistanceMatchesPosition(
    const SamplingPointList & samples,
    const std::array<float, 3> & center)
{
    for (const auto & sample : samples)
    {
        EXPECT_NEAR(sample.distance, ComputeDistance(center, sample.position), 1e-5f);
    }
}

} // namespace

TEST(SphereSamplerTest, DefaultConfigurationProducesRadiusUniformSamples)
{
    SphereSampler sampler;
    const auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_GE(sample.distance, 0.0f);
        EXPECT_LE(sample.distance, 1.0f);
    }
}

TEST(SphereSamplerTest, AnalysisDefaultUsesCommandDistanceRange)
{
    const auto sampler{
        SphereSampler::VolumeUniformRandom(
            SphereDistanceRange{ 0.0, 1.5 },
            12)
    };
    const auto analysis_sampler{
        SphereSampler::AnalysisDefault(
            SphereSamplingMethod::VolumeUniformRandom,
            12)
    };

    const auto samples{ analysis_sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    ASSERT_EQ(12u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_GE(sample.distance, 0.0f);
        EXPECT_LE(sample.distance, 1.5f);
    }

    const auto comparison_samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };
    EXPECT_EQ(comparison_samples.size(), samples.size());
}

TEST(SphereSamplerTest, AnalysisDefaultFibonacciKeepsCommandSamplingSizePerRadius)
{
    const auto sampler{
        SphereSampler::AnalysisDefault(
            SphereSamplingMethod::FibonacciDeterministic,
            7)
    };
    const auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    ASSERT_EQ(105u, samples.size());
    EXPECT_EQ(7u, CountSamplesAtRadius(samples, 0.05f));
    EXPECT_EQ(7u, CountSamplesAtRadius(samples, 1.45f));
}

TEST(SphereSamplerTest, RadiusUniformRandomProducesConfiguredSampleCount)
{
    const auto sampler{
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 1.0, 1.0 },
            5)
    };
    const auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    ASSERT_EQ(5u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_NEAR(1.0f, sample.distance, 1e-5f);
    }
}

TEST(SphereSamplerTest, RadiusUniformRandomKeepsRadiusWithinConfiguredRange)
{
    const auto sampler{
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 1.0, 2.0 },
            10)
    };
    const std::array<float, 3> center{ 0.f, 0.f, 0.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_GE(sample.distance, 1.0f);
        EXPECT_LE(sample.distance, 2.0f);
    }
    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, RadiusUniformRandomUsesReferencePosition)
{
    const auto sampler{
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 1.0, 1.0 },
            5)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, ZeroSampleCountReturnsEmptyVector)
{
    const auto sampler{
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 0.0, 1.0 },
            0)
    };
    const auto samples{ sampler.GenerateSamplingPoints({0.f, 0.f, 0.f}) };

    EXPECT_TRUE(samples.empty());
}

TEST(SphereSamplerTest, ZeroDistanceRangeReturnsCenterPoints)
{
    const auto sampler{
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 0.0, 0.0 },
            10)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(10u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_FLOAT_EQ(0.0f, sample.distance);
        EXPECT_FLOAT_EQ(center[0], sample.position[0]);
        EXPECT_FLOAT_EQ(center[1], sample.position[1]);
        EXPECT_FLOAT_EQ(center[2], sample.position[2]);
    }
}

TEST(SphereSamplerTest, RadiusUniformRandomThrowsWhenRangeIsInvalid)
{
    EXPECT_THROW(
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 2.0, 1.0 },
            10),
        std::invalid_argument);
    EXPECT_THROW(
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ -0.1, 1.0 },
            10),
        std::invalid_argument);
    EXPECT_THROW(
        SphereSampler::RadiusUniformRandom(
            SphereDistanceRange{ 0.0, -1.0 },
            10),
        std::invalid_argument);
}

TEST(SphereSamplerTest, FibonacciDeterministicProducesDeterministicSamples)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.3,
            4)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };

    const auto first_samples{ sampler.GenerateSamplingPoints(center) };
    const auto second_samples{ sampler.GenerateSamplingPoints(center) };

    ExpectSamplesEqual(first_samples, second_samples);
}

TEST(SphereSamplerTest, FibonacciDeterministicVaryingRadiusProducesDeterministicSamples)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.3,
            4,
            true)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };

    const auto first_samples{ sampler.GenerateSamplingPoints(center) };
    const auto second_samples{ sampler.GenerateSamplingPoints(center) };

    ExpectSamplesEqual(first_samples, second_samples);
}

TEST(SphereSamplerTest, FibonacciDeterministicSingleShellProducesConfiguredPointCount)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 1.25, 1.25 },
            0.4,
            6)
    };
    const std::array<float, 3> center{ 0.f, 0.f, 0.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(6u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_FLOAT_EQ(1.25f, sample.distance);
    }
    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, FibonacciDeterministicProducesExpectedShellPointCount)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.3,
            4)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(8u, samples.size());
    EXPECT_EQ(4u, CountSamplesAtRadius(samples, 0.65f));
    EXPECT_EQ(4u, CountSamplesAtRadius(samples, 0.95f));
    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, FibonacciDeterministicZeroRadiusShellReturnsCenterPoints)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 0.0 },
            0.5,
            5)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(5u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_FLOAT_EQ(0.0f, sample.distance);
        EXPECT_FLOAT_EQ(center[0], sample.position[0]);
        EXPECT_FLOAT_EQ(center[1], sample.position[1]);
        EXPECT_FLOAT_EQ(center[2], sample.position[2]);
    }
}

TEST(SphereSamplerTest, FibonacciDeterministicRangeStartingAtZeroUsesBinCenters)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 1.0 },
            0.5,
            4)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(8u, samples.size());
    EXPECT_EQ(0u, CountSamplesAtRadius(samples, 0.0f));
    EXPECT_EQ(4u, CountSamplesAtRadius(samples, 0.25f));
    EXPECT_EQ(4u, CountSamplesAtRadius(samples, 0.75f));
    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, FibonacciDeterministicNarrowRangeUsesRangeMidpointShell)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 0.2 },
            0.5,
            4)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(4u, samples.size());
    EXPECT_EQ(4u, CountSamplesAtRadius(samples, 0.1f));
    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, FibonacciDeterministicVaryingRadiusScalesShellPointCountByArea)
{
    const auto sampler{
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.5, 1.0 },
            0.3,
            4,
            true)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(6u, samples.size());
    EXPECT_EQ(2u, CountSamplesAtRadius(samples, 0.65f));
    EXPECT_EQ(4u, CountSamplesAtRadius(samples, 0.95f));
    ExpectDistanceMatchesPosition(samples, center);
}

TEST(SphereSamplerTest, FibonacciDeterministicThrowsWhenConfigIsInvalid)
{
    EXPECT_THROW(
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 1.0 },
            0.0,
            4),
        std::invalid_argument);
    EXPECT_THROW(
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 1.0 },
            0.5,
            0),
        std::invalid_argument);
    EXPECT_THROW(
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 0.0, 1.0 },
            std::numeric_limits<double>::infinity(),
            4),
        std::invalid_argument);
    EXPECT_THROW(
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ -0.1, 1.0 },
            0.5,
            4),
        std::invalid_argument);
    EXPECT_THROW(
        SphereSampler::FibonacciDeterministic(
            SphereDistanceRange{ 2.0, 1.0 },
            0.5,
            4),
        std::invalid_argument);
}

TEST(SphereSamplerTest, VolumeUniformRandomProducesConfiguredSampleCount)
{
    const auto sampler{
        SphereSampler::VolumeUniformRandom(
            SphereDistanceRange{ 0.5, 1.5 },
            20)
    };
    const std::array<float, 3> center{ 1.f, 2.f, 3.f };
    const auto samples{ sampler.GenerateSamplingPoints(center) };

    ASSERT_EQ(20u, samples.size());
    for (const auto & sample : samples)
    {
        EXPECT_GE(sample.distance, 0.5f);
        EXPECT_LE(sample.distance, 1.5f);
    }
    ExpectDistanceMatchesPosition(samples, center);
}
