#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

namespace {

using namespace rhbm_gem::sphere_sampler;

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
        EXPECT_EQ(lhs[i].distance, rhs[i].distance);
        EXPECT_EQ(lhs[i].position, rhs[i].position);
        EXPECT_EQ(lhs[i].is_selected, rhs[i].is_selected);
    }
}

void ExpectAnalysisSampleDistances(
    const SamplingPointList & samples,
    const std::array<double, 3> & center)
{
    for (const auto & sample : samples)
    {
        EXPECT_GT(sample.distance, 0.0);
        EXPECT_LE(sample.distance, 2.0);
        EXPECT_NEAR(sample.distance, ComputeDistance(center, sample.position), 1e-5);
        EXPECT_TRUE(sample.is_selected);
    }
}

} // namespace

TEST(SphereSamplerTest, FibonacciProducesFixedDeterministicShells)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };

    const auto first_samples{ GenerateFibonacciDeterministic(center) };
    const auto second_samples{ GenerateFibonacciDeterministic(center) };

    ASSERT_EQ(200u, first_samples.size());
    for (std::size_t i = 0; i < first_samples.size(); ++i)
    {
        EXPECT_NEAR(first_samples[i].distance, 0.05 + 0.1 * static_cast<double>(i / 10), 1e-14);
    }
    ExpectAnalysisSampleDistances(first_samples, center);
    ExpectSamplesEqual(first_samples, second_samples);
}

TEST(SphereSamplerTest, DispatchThrowsForUnsupportedMethod)
{
    for (const int value : { 0, 1, 99 })
    {
        EXPECT_THROW(
            GenerateSamplingPointList({ 0.0, 0.0, 0.0 }, static_cast<SphereSamplingMethod>(value)),
            std::invalid_argument);
    }
}

TEST(SphereSamplerTest, DispatchUsesRequestedMethod)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };
    const auto fibonacci{ GenerateSamplingPointList(center, SphereSamplingMethod::FibonacciDeterministic) };
    ExpectSamplesEqual(fibonacci, GenerateSamplingPointList(center));
    ExpectSamplesEqual(fibonacci, GenerateFibonacciDeterministic(center));

    const auto spiral{ GenerateSamplingPointList(center, SphereSamplingMethod::SpiralDeterministic) };
    ASSERT_EQ(200u, spiral.size());
    ExpectAnalysisSampleDistances(spiral, center);
    ExpectSamplesEqual(spiral, GenerateSamplingPointList(center, SphereSamplingMethod::SpiralDeterministic));
    for (std::size_t offset = 0; offset < spiral.size(); offset += 10)
    {
        const auto shell{ GenerateSpiralShell(center, fibonacci[offset].distance, 10) };
        ExpectSamplesEqual(shell, SamplingPointList(
            spiral.begin() + static_cast<std::ptrdiff_t>(offset),
            spiral.begin() + static_cast<std::ptrdiff_t>(offset + 10)));
    }
    EXPECT_NE(spiral.front().position, fibonacci.front().position);
}

TEST(SphereSamplerTest, SpiralShellPreservesPolarOrderAndInteriorRecurrence)
{
    const std::array<double, 3> center{ 1.0, 2.0, 3.0 };
    const auto points{ GenerateSpiralShell(center, 1.0, 4) };
    ASSERT_EQ(4u, points.size());
    EXPECT_NEAR(points.front().position[0], 1.0, 1e-15);
    EXPECT_EQ(points.front().position[1], 2.0);
    EXPECT_EQ(points.front().position[2], 2.0);
    EXPECT_EQ(points.back().position, (std::array<double, 3>{ 1.0, 2.0, 4.0 }));

    const double height{ 1.0 / 3.0 };
    const double radial_xy{ std::sqrt(1.0 - height * height) };
    const double theta_step{ 3.6 / (2.0 * radial_xy) };
    for (std::size_t i = 1; i <= 2; ++i)
    {
        const double theta{ static_cast<double>(i) * theta_step };
        EXPECT_NEAR(points[i].position[0], 1.0 + radial_xy * std::cos(theta), 1e-15);
        EXPECT_NEAR(points[i].position[1], 2.0 + radial_xy * std::sin(theta), 1e-15);
        EXPECT_NEAR(points[i].position[2], 3.0 + (i == 1 ? -height : height), 1e-15);
    }
}

TEST(SphereSamplerTest, SpiralShellSupportsVariableCountsAndTranslation)
{
    const std::array<double, 3> origin{};
    const std::array<double, 3> center{ 1.0, -2.0, 3.0 };
    for (const std::size_t count : { 2u, 7u, 106u })
    {
        const auto points{ GenerateSpiralShell(origin, 1.25, count) };
        const auto translated{ GenerateSpiralShell(center, 1.25, count) };
        ASSERT_EQ(count, points.size());
        ASSERT_EQ(count, translated.size());
        ExpectAnalysisSampleDistances(translated, center);
        for (std::size_t i = 0; i < count; ++i)
        {
            EXPECT_EQ(translated[i].distance, 1.25);
            for (std::size_t axis = 0; axis < center.size(); ++axis)
            {
                EXPECT_NEAR(translated[i].position[axis] - center[axis], points[i].position[axis], 1e-15);
            }
        }
    }
}

TEST(SphereSamplerTest, SpiralShellRejectsInvalidRadiusAndCount)
{
    const std::array<double, 3> center{};
    for (const double radius : { 0.0, -1.0, std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::quiet_NaN() })
    {
        EXPECT_THROW(GenerateSpiralShell(center, radius, 2), std::invalid_argument);
    }
    EXPECT_THROW(GenerateSpiralShell(center, 1.0, 0), std::invalid_argument);
    EXPECT_THROW(GenerateSpiralShell(center, 1.0, 1), std::invalid_argument);
}
