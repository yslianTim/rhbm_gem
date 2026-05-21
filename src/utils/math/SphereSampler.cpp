#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kAnalysisDistanceMin{ 0.0 };
constexpr double kAnalysisDistanceMax{ 1.5 };
constexpr double kAnalysisFibonacciRadiusBinSize{ 0.1 };
constexpr unsigned int kAnalysisSampleCount{ 50 };

std::vector<double> BuildFibonacciShellCenters(
    double distance_min,
    double distance_max,
    double radius_bin_size)
{
    std::vector<double> shell_centers;

    constexpr double epsilon{ 1e-9 };
    if (distance_min == distance_max)
    {
        shell_centers.push_back(distance_min);
        return shell_centers;
    }

    const double distance_span{ distance_max - distance_min };
    if (distance_span + epsilon < radius_bin_size)
    {
        shell_centers.push_back((distance_min + distance_max) * 0.5);
        return shell_centers;
    }

    for (double shell_center{ distance_min + radius_bin_size * 0.5 };
         shell_center <= distance_max + epsilon;
         shell_center += radius_bin_size)
    {
        shell_centers.push_back(shell_center);
    }

    return shell_centers;
}

std::size_t GetFibonacciSampleCountForRadius(
    double radius,
    double distance_max,
    unsigned int sample_count,
    bool vary_with_radius)
{
    constexpr double epsilon{ 1e-9 };

    if (!vary_with_radius)
    {
        return sample_count;
    }

    if (distance_max <= epsilon)
    {
        return 1;
    }

    const double normalized_radius{ radius / distance_max };
    const double scaled_sample_count{
        static_cast<double>(sample_count) * normalized_radius * normalized_radius
    };
    return std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(scaled_sample_count)));
}

std::array<float, 3> GenerateRandomUnitDirection(
    std::mt19937 & engine,
    std::uniform_real_distribution<float> & dist_phi,
    std::uniform_real_distribution<float> & dist_cos_theta)
{
    const float phi{ dist_phi(engine) };
    const float cos_theta{ dist_cos_theta(engine) };
    const float sin_theta{ std::sqrt(1.0f - cos_theta * cos_theta) };

    return {
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta
    };
}

SamplingPoint MakeSamplingPoint(
    const std::array<float, 3> & center_position,
    float radius,
    const std::array<float, 3> & unit_vector)
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

SphereSampler::SphereSampler(
    SphereSamplingMethod method,
    double distance_min,
    double distance_max,
    unsigned int sample_count,
    double radius_bin_size,
    bool vary_with_radius) :
    m_method{ method },
    m_distance_min{ distance_min },
    m_distance_max{ distance_max },
    m_sample_count{ sample_count },
    m_radius_bin_size{ radius_bin_size },
    m_vary_with_radius{ vary_with_radius }
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegativeRange(
        distance_min,
        distance_max,
        "SphereSampler distance range");

    if (method == SphereSamplingMethod::FibonacciDeterministic)
    {
        rhbm_gem::numeric_validation::RequireFinitePositive(
            radius_bin_size,
            "SphereSampler Fibonacci radius bin size");
        rhbm_gem::numeric_validation::RequirePositive(
            sample_count,
            "SphereSampler analysis sample count");
    }
}

SphereSampler SphereSampler::AnalysisDefault(SphereSamplingMethod method)
{
    switch (method)
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            return SphereSampler(
                method,
                kAnalysisDistanceMin,
                kAnalysisDistanceMax,
                kAnalysisSampleCount,
                kAnalysisFibonacciRadiusBinSize,
                false);
        case SphereSamplingMethod::VolumeUniformRandom:
            return SphereSampler(
                method,
                kAnalysisDistanceMin,
                kAnalysisDistanceMax,
                kAnalysisSampleCount,
                kAnalysisFibonacciRadiusBinSize,
                false);
        case SphereSamplingMethod::FibonacciDeterministic:
            return SphereSampler(
                method,
                kAnalysisDistanceMin,
                kAnalysisDistanceMax,
                kAnalysisSampleCount,
                kAnalysisFibonacciRadiusBinSize,
                false);
    }

    throw std::invalid_argument("Unsupported SphereSamplingMethod.");
}

SamplingPointList SphereSampler::GenerateSamplingPoints(
    const std::array<float, 3> & center_position) const
{
    SamplingPointList out;
    switch (m_method)
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            GenerateRadiusUniformRandom(center_position, out);
            break;
        case SphereSamplingMethod::VolumeUniformRandom:
            GenerateVolumeUniformRandom(center_position, out);
            break;
        case SphereSamplingMethod::FibonacciDeterministic:
            GenerateFibonacciDeterministic(center_position, out);
            break;
    }

    return out;
}

void SphereSampler::GenerateVolumeUniformRandom(
    const std::array<float, 3> & center_position,
    SamplingPointList & out) const
{
    out.resize(m_sample_count);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_unit(0.0f, 1.0f);
    std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_cos_theta(-1.0f, 1.0f);
    const float min_radius_cube{
        static_cast<float>(m_distance_min * m_distance_min * m_distance_min)
    };
    const float max_radius_cube{
        static_cast<float>(m_distance_max * m_distance_max * m_distance_max)
    };

    for (unsigned int i = 0; i < m_sample_count; i++)
    {
        const float radius_unit{ dist_unit(engine) };
        const float radius{
            std::cbrt(min_radius_cube + radius_unit * (max_radius_cube - min_radius_cube))
        };

        out[i] = MakeSamplingPoint(
            center_position,
            radius,
            GenerateRandomUnitDirection(engine, dist_phi, dist_cos_theta));
    }
}

void SphereSampler::GenerateFibonacciDeterministic(
    const std::array<float, 3> & center_position,
    SamplingPointList & out) const
{
    const auto shell_radii{
        BuildFibonacciShellCenters(m_distance_min, m_distance_max, m_radius_bin_size)
    };
    out.clear();

    std::size_t total_point_count{ 0 };
    for (const double radius : shell_radii)
    {
        total_point_count += GetFibonacciSampleCountForRadius(
            radius,
            m_distance_max,
            m_sample_count,
            m_vary_with_radius);
    }
    out.reserve(total_point_count);

    const double golden_angle{ Constants::pi * (3.0 - std::sqrt(5.0)) };

    // Shell radii are bin centers within the configured range.
    // Each shell uses the same deterministic Fibonacci sphere pattern and reports its shell radius.
    for (const double radius : shell_radii)
    {
        const std::size_t shell_sample_count{
            GetFibonacciSampleCountForRadius(
                radius,
                m_distance_max,
                m_sample_count,
                m_vary_with_radius)
        };

        for (std::size_t i = 0; i < shell_sample_count; i++)
        {
            const double z{
                1.0 - 2.0 * (static_cast<double>(i) + 0.5)
                    / static_cast<double>(shell_sample_count)
            };
            const double radial_xy{ std::sqrt(std::max(0.0, 1.0 - z * z)) };
            const double theta{ golden_angle * static_cast<double>(i) };
            const double x{ radial_xy * std::cos(theta) };
            const double y{ radial_xy * std::sin(theta) };

            const std::array<float, 3> unit_vector{
                static_cast<float>(x),
                static_cast<float>(y),
                static_cast<float>(z)
            };
            out.emplace_back(MakeSamplingPoint(
                center_position,
                static_cast<float>(radius),
                unit_vector));
        }
    }
}

void SphereSampler::GenerateRadiusUniformRandom(
    const std::array<float, 3> & center_position,
    SamplingPointList & out) const
{
    out.resize(m_sample_count);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_radius(
        static_cast<float>(m_distance_min), static_cast<float>(m_distance_max)
    );
    std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_cos_theta(-1.0f, 1.0f);

    for (unsigned int i = 0; i < m_sample_count; i++)
    {
        const float radius{ dist_radius(engine) };
        out[i] = MakeSamplingPoint(
            center_position,
            radius,
            GenerateRandomUnitDirection(engine, dist_phi, dist_cos_theta));
    }
}
