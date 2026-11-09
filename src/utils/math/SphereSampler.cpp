#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::sphere_sampler {

namespace {

constexpr double kAnalysisDistanceMin{ 0.0 };
constexpr double kAnalysisDistanceMax{ 2.0 };
constexpr double kAnalysisRadiusBinSize{ 0.1 };
constexpr unsigned int kAnalysisSampleCount{ 10 };
constexpr double kSpiralSpacing{ 3.6 };

static_assert(kAnalysisDistanceMin >= 0.0);
static_assert(kAnalysisDistanceMax >= kAnalysisDistanceMin);
static_assert(kAnalysisRadiusBinSize > 0.0);
static_assert(kAnalysisSampleCount > 0);

std::vector<double> BuildAnalysisShellCenters()
{
    std::vector<double> shell_centers;
    constexpr double epsilon{ 1e-9 };
    if (kAnalysisDistanceMin == kAnalysisDistanceMax)
    {
        shell_centers.push_back(kAnalysisDistanceMin);
        return shell_centers;
    }

    const double distance_span{ kAnalysisDistanceMax - kAnalysisDistanceMin };
    if (distance_span + epsilon < kAnalysisRadiusBinSize)
    {
        shell_centers.push_back((kAnalysisDistanceMin + kAnalysisDistanceMax) * 0.5);
        return shell_centers;
    }

    for (double shell_center{ kAnalysisDistanceMin + kAnalysisRadiusBinSize * 0.5 };
         shell_center <= kAnalysisDistanceMax + epsilon;
         shell_center += kAnalysisRadiusBinSize)
    {
        shell_centers.push_back(shell_center);
    }
    return shell_centers;
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

void AppendSpiralShell(
    SamplingPointList & points,
    const std::array<double, 3> & center,
    double radius,
    std::size_t point_count)
{
    double previous_theta{ 0.0 };
    for (std::size_t i = 0; i < point_count; ++i)
    {
        const auto height{
            -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(point_count - 1)
        };
        const auto phi{ std::acos(height) };
        double theta{ 0.0 };
        if (i > 0 && i + 1 < point_count)
        {
            theta = std::fmod(
                previous_theta +
                    kSpiralSpacing /
                    std::sqrt(static_cast<double>(point_count) * (1.0 - height * height)),
                Constants::two_pi);
        }
        previous_theta = theta;

        const auto sin_phi{ std::sin(phi) };
        points.emplace_back(SamplingPoint{
            radius,
            {
                center.at(0) + radius * sin_phi * std::cos(theta),
                center.at(1) + radius * sin_phi * std::sin(theta),
                center.at(2) + radius * std::cos(phi)
            }
        });
    }
}

} // namespace

SamplingPointList GenerateFibonacciDeterministic(const std::array<double, 3> & center_position)
{
    const auto shell_radii{ BuildAnalysisShellCenters() };
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

SamplingPointList GenerateSpiralShell(
    const std::array<double, 3> & center_position,
    double radius,
    std::size_t point_count)
{
    numeric_validation::RequireFinitePositive(radius, "GenerateSpiralShell radius");
    numeric_validation::RequireAtLeast(point_count, std::size_t{ 2 }, "GenerateSpiralShell point_count");

    SamplingPointList points;
    points.reserve(point_count);
    AppendSpiralShell(points, center_position, radius, point_count);
    return points;
}

SamplingPointList GenerateSamplingPointList(
    const std::array<double, 3> & center_position,
    SphereSamplingMethod method)
{
    switch (method)
    {
        case SphereSamplingMethod::FibonacciDeterministic:
            return GenerateFibonacciDeterministic(center_position);
        case SphereSamplingMethod::SpiralDeterministic:
        {
            const auto shell_radii{ BuildAnalysisShellCenters() };
            SamplingPointList points;
            points.reserve(shell_radii.size() * kAnalysisSampleCount);
            for (const double radius : shell_radii)
            {
                AppendSpiralShell(points, center_position, radius, kAnalysisSampleCount);
            }
            return points;
        }
    }
    throw std::invalid_argument("Unsupported SphereSamplingMethod.");
}

} // namespace rhbm_gem::sphere_sampler
