#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::sphere_sampler {

namespace {

constexpr double kAnalysisDistanceMin{ 0.0 };
constexpr double kAnalysisDistanceMax{ 2.0 };
constexpr double kAnalysisFibonacciRadiusBinSize{ 0.1 };
constexpr unsigned int kAnalysisSampleCount{ 10 };

static_assert(kAnalysisDistanceMin >= 0.0);
static_assert(kAnalysisDistanceMax >= kAnalysisDistanceMin);
static_assert(kAnalysisFibonacciRadiusBinSize > 0.0);
static_assert(kAnalysisSampleCount > 0);

std::vector<double> BuildFibonacciShellCenters()
{
    std::vector<double> shell_centers;
    constexpr double epsilon{ 1e-9 };
    if (kAnalysisDistanceMin == kAnalysisDistanceMax)
    {
        shell_centers.push_back(kAnalysisDistanceMin);
        return shell_centers;
    }

    const double distance_span{ kAnalysisDistanceMax - kAnalysisDistanceMin };
    if (distance_span + epsilon < kAnalysisFibonacciRadiusBinSize)
    {
        shell_centers.push_back((kAnalysisDistanceMin + kAnalysisDistanceMax) * 0.5);
        return shell_centers;
    }

    for (double shell_center{ kAnalysisDistanceMin + kAnalysisFibonacciRadiusBinSize * 0.5 };
         shell_center <= kAnalysisDistanceMax + epsilon;
         shell_center += kAnalysisFibonacciRadiusBinSize)
    {
        shell_centers.push_back(shell_center);
    }
    return shell_centers;
}

std::array<double, 3> GenerateRandomUnitDirection(
    std::mt19937 & engine,
    std::uniform_real_distribution<double> & dist_phi,
    std::uniform_real_distribution<double> & dist_cos_theta)
{
    const double phi{ dist_phi(engine) };
    const double cos_theta{ dist_cos_theta(engine) };
    const double sin_theta{ std::sqrt(1.0 - cos_theta * cos_theta) };

    return {
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta
    };
}

SamplingPoint BuildSamplingPoint(
    const std::array<double, 3> & center_position,
    double radius,
    const std::array<double, 3> & unit_vector)
{
    return SamplingPoint{
        radius,
        {
            center_position.at(0) + radius * unit_vector[0],
            center_position.at(1) + radius * unit_vector[1],
            center_position.at(2) + radius * unit_vector[2]
        }
    };
}

} // namespace

SamplingPointList GenerateVolumeUniformRandom(const std::array<double, 3> & center_position)
{
    SamplingPointList out;
    out.resize(kAnalysisSampleCount);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<double> dist_unit(0.0, 1.0);
    std::uniform_real_distribution<double> dist_phi(0.0, Constants::two_pi);
    std::uniform_real_distribution<double> dist_cos_theta(-1.0, 1.0);
    const double min_radius_cube{
        kAnalysisDistanceMin * kAnalysisDistanceMin * kAnalysisDistanceMin
    };
    const double max_radius_cube{
        kAnalysisDistanceMax * kAnalysisDistanceMax * kAnalysisDistanceMax
    };

    for (unsigned int i = 0; i < kAnalysisSampleCount; i++)
    {
        const double radius_unit{ dist_unit(engine) };
        const double radius{
            std::cbrt(min_radius_cube + radius_unit * (max_radius_cube - min_radius_cube))
        };

        out[i] = BuildSamplingPoint(
            center_position,
            radius,
            GenerateRandomUnitDirection(engine, dist_phi, dist_cos_theta));
    }
    return out;
}

SamplingPointList GenerateFibonacciDeterministic(const std::array<double, 3> & center_position)
{
    const auto shell_radii{ BuildFibonacciShellCenters() };
    SamplingPointList out;
    out.reserve(shell_radii.size() * kAnalysisSampleCount);

    const double golden_angle{ Constants::pi * (3.0 - std::sqrt(5.0)) };

    // Shell radii are bin centers within the configured range.
    // Each shell uses the same deterministic Fibonacci sphere pattern and reports its shell radius.
    for (const double radius : shell_radii)
    {
        for (std::size_t i = 0; i < kAnalysisSampleCount; i++)
        {
            const double z{
                1.0 - 2.0 * (static_cast<double>(i) + 0.5)
                    / static_cast<double>(kAnalysisSampleCount)
            };
            const double radial_xy{ std::sqrt(std::max(0.0, 1.0 - z * z)) };
            const double theta{ golden_angle * static_cast<double>(i) };
            const double x{ radial_xy * std::cos(theta) };
            const double y{ radial_xy * std::sin(theta) };

            const std::array<double, 3> unit_vector{ x, y, z };
            out.emplace_back(BuildSamplingPoint(
                center_position,
                radius,
                unit_vector));
        }
    }

    return out;
}

SamplingPointList GenerateRadiusUniformRandom(const std::array<double, 3> & center_position)
{
    SamplingPointList out;
    out.resize(kAnalysisSampleCount);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<double> dist_radius(
        kAnalysisDistanceMin, kAnalysisDistanceMax
    );
    std::uniform_real_distribution<double> dist_phi(0.0, Constants::two_pi);
    std::uniform_real_distribution<double> dist_cos_theta(-1.0, 1.0);

    for (unsigned int i = 0; i < kAnalysisSampleCount; i++)
    {
        const double radius{ dist_radius(engine) };
        out[i] = BuildSamplingPoint(
            center_position,
            radius,
            GenerateRandomUnitDirection(engine, dist_phi, dist_cos_theta));
    }
    return out;
}

SamplingPointList GenerateSamplingPointList(
    const std::array<double, 3> & center_position,
    SphereSamplingMethod method)
{
    switch (method)
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            return GenerateRadiusUniformRandom(center_position);
        case SphereSamplingMethod::VolumeUniformRandom:
            return GenerateVolumeUniformRandom(center_position);
        case SphereSamplingMethod::FibonacciDeterministic:
            return GenerateFibonacciDeterministic(center_position);
    }
    throw std::invalid_argument("Unsupported SphereSamplingMethod.");
}

} // namespace rhbm_gem::sphere_sampler
