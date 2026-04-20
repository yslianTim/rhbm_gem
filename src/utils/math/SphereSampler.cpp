#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <string_view>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace {

void ValidateSphereRandomSamplingInputs(
    SphereDistanceRange range,
    SphereRandomSamplingConfig random_config)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegativeRange(
        range.min,
        range.max,
        "SphereSampler distance range");
    (void)random_config;
}

void ValidateSphereFibonacciSamplingInputs(
    SphereDistanceRange range,
    SphereDeterministicSamplingConfig fibonacci_config)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegativeRange(
        range.min,
        range.max,
        "SphereSampler distance range");
    rhbm_gem::NumericValidation::RequireFinitePositive(
        fibonacci_config.radius_bin_size,
        "SphereSampler Fibonacci radius bin size");
    rhbm_gem::NumericValidation::RequirePositive(
        fibonacci_config.samples_per_radius,
        "SphereSampler Fibonacci samples per radius");
}

std::string_view GetSphereSamplingMethodName(SphereSamplingMethod method)
{
    switch (method)
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            return "RadiusUniformRandom";
        case SphereSamplingMethod::VolumeUniformRandom:
            return "VolumeUniformRandom";
        case SphereSamplingMethod::FibonacciDeterministic:
            return "FibonacciDeterministic";
    }

    return "Unknown";
}

std::vector<double> BuildFibonacciShellCenters(
    const SphereDistanceRange & distance_range,
    double radius_bin_size)
{
    std::vector<double> shell_centers;

    constexpr double epsilon{ 1e-9 };
    if (distance_range.min == distance_range.max)
    {
        shell_centers.push_back(distance_range.min);
        return shell_centers;
    }

    const double distance_span{ distance_range.max - distance_range.min };
    if (distance_span + epsilon < radius_bin_size)
    {
        shell_centers.push_back((distance_range.min + distance_range.max) * 0.5);
        return shell_centers;
    }

    for (double shell_center{ distance_range.min + radius_bin_size * 0.5 };
         shell_center <= distance_range.max + epsilon;
         shell_center += radius_bin_size)
    {
        shell_centers.push_back(shell_center);
    }

    return shell_centers;
}

std::size_t GetFibonacciSampleCountForRadius(
    double radius,
    const SphereDistanceRange & distance_range,
    const SphereDeterministicSamplingConfig & config)
{
    constexpr double epsilon{ 1e-9 };

    if (!config.vary_with_radius)
    {
        return config.samples_per_radius;
    }

    if (distance_range.max <= epsilon)
    {
        return 1;
    }

    const double normalized_radius{ radius / distance_range.max };
    const double scaled_sample_count{
        static_cast<double>(config.samples_per_radius) * normalized_radius * normalized_radius
    };
    return std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(scaled_sample_count)));
}

} // namespace

SphereSamplingProfile::SphereSamplingProfile(
    SphereSamplingMethod method,
    SphereDistanceRange range,
    SphereRandomSamplingConfig random_config) :
    m_method{ method },
    m_distance_range{ range },
    m_method_config{ random_config }
{
    ValidateSphereRandomSamplingInputs(range, random_config);
}

SphereSamplingProfile::SphereSamplingProfile(
    SphereSamplingMethod method,
    SphereDistanceRange range,
    SphereDeterministicSamplingConfig fibonacci_config) :
    m_method{ method },
    m_distance_range{ range },
    m_method_config{ fibonacci_config }
{
    ValidateSphereFibonacciSamplingInputs(range, fibonacci_config);
}

SphereSamplingProfile SphereSamplingProfile::RadiusUniformRandom(
    SphereDistanceRange range,
    unsigned int sample_count)
{
    return SphereSamplingProfile(
        SphereSamplingMethod::RadiusUniformRandom,
        range,
        SphereRandomSamplingConfig{ sample_count });
}

SphereSamplingProfile SphereSamplingProfile::VolumeUniformRandom(
    SphereDistanceRange range,
    unsigned int sample_count)
{
    return SphereSamplingProfile(
        SphereSamplingMethod::VolumeUniformRandom,
        range,
        SphereRandomSamplingConfig{ sample_count });
}

SphereSamplingProfile SphereSamplingProfile::FibonacciDeterministic(
    SphereDistanceRange range,
    double radius_bin_size,
    unsigned int samples_per_radius,
    bool vary_with_radius)
{
    return SphereSamplingProfile(
        SphereSamplingMethod::FibonacciDeterministic,
        range,
        SphereDeterministicSamplingConfig{
            radius_bin_size,
            samples_per_radius,
            vary_with_radius
        });
}

const SphereRandomSamplingConfig & SphereSamplingProfile::GetRandomConfig() const
{
    return std::get<SphereRandomSamplingConfig>(m_method_config);
}

const SphereDeterministicSamplingConfig & SphereSamplingProfile::GetFibonacciConfig() const
{
    return std::get<SphereDeterministicSamplingConfig>(m_method_config);
}

SphereSampler::SphereSampler() :
    m_profile{ SphereSamplingProfile::RadiusUniformRandom(SphereDistanceRange{ 0.0, 1.0 }, 10) }
{
}

void SphereSampler::Print() const
{
    const auto & distance_range{ m_profile.GetDistanceRange() };

    std::ostringstream oss;
    oss << "SphereSampler Configuration:\n"
        << " - Sampling method: " << GetSphereSamplingMethodName(m_profile.GetMethod()) << '\n'
        << " - Distance range: ["
        << StringHelper::ToStringWithPrecision<double>(distance_range.min, 1) << ", "
        << StringHelper::ToStringWithPrecision<double>(distance_range.max, 1)
        << "] Angstrom.\n";

    switch (m_profile.GetMethod())
    {
        case SphereSamplingMethod::RadiusUniformRandom:
        case SphereSamplingMethod::VolumeUniformRandom:
        {
            const auto & config{ m_profile.GetRandomConfig() };
            oss << " - Sample count: " << config.sample_count;
            break;
        }
        case SphereSamplingMethod::FibonacciDeterministic:
        {
            const auto & config{ m_profile.GetFibonacciConfig() };
            oss << " - Radius bin size: "
                << StringHelper::ToStringWithPrecision<double>(config.radius_bin_size, 2)
                << " Angstrom.\n"
                << " - Samples per radius: " << config.samples_per_radius << '\n'
                << " - Vary with radius: " << std::boolalpha << config.vary_with_radius << '\n'
                << " - Expected point count: " << GetExpectedPointCount();
            break;
        }
    }

    Logger::Log(LogLevel::Info, oss.str());
}

void SphereSampler::SetSamplingProfile(const SphereSamplingProfile & profile)
{
    m_profile = profile;
}

std::size_t SphereSampler::GetExpectedPointCount() const
{
    switch (m_profile.GetMethod())
    {
        case SphereSamplingMethod::RadiusUniformRandom:
        case SphereSamplingMethod::VolumeUniformRandom:
            return m_profile.GetRandomConfig().sample_count;
        case SphereSamplingMethod::FibonacciDeterministic:
        {
            const auto & config{ m_profile.GetFibonacciConfig() };
            const auto shell_radii{
                BuildFibonacciShellCenters(
                    m_profile.GetDistanceRange(),
                    config.radius_bin_size)
            };

            std::size_t expected_point_count{ 0 };
            for (const double radius : shell_radii)
            {
                expected_point_count += GetFibonacciSampleCountForRadius(
                    radius,
                    m_profile.GetDistanceRange(),
                    config);
            }

            return expected_point_count;
        }
    }

    return 0;
}

SamplingPointList SphereSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position) const
{
    SamplingPointList out;
    const auto & distance_range{ m_profile.GetDistanceRange() };
    switch (m_profile.GetMethod())
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            GenerateRadiusUniformRandom(
                reference_position,
                distance_range,
                m_profile.GetRandomConfig(),
                out);
            break;
        case SphereSamplingMethod::VolumeUniformRandom:
            GenerateVolumeUniformRandom(
                reference_position,
                distance_range,
                m_profile.GetRandomConfig(),
                out);
            break;
        case SphereSamplingMethod::FibonacciDeterministic:
            GenerateFibonacciDeterministic(
                reference_position,
                distance_range,
                m_profile.GetFibonacciConfig(),
                out);
            break;
    }

    return out;
}

void SphereSampler::GenerateVolumeUniformRandom(
    const std::array<float, 3> & reference_position,
    const SphereDistanceRange & distance_range,
    const SphereRandomSamplingConfig & config,
    SamplingPointList & out) const
{
    out.resize(config.sample_count);
    
    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_unit(0.0f, 1.0f);
    std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_cos_theta(-1.0f, 1.0f);
    const float min_radius_cube{
        static_cast<float>(distance_range.min * distance_range.min * distance_range.min)
    };
    const float max_radius_cube{
        static_cast<float>(distance_range.max * distance_range.max * distance_range.max)
    };

    for (unsigned int i = 0; i < config.sample_count; i++)
    {
        const float radius_unit{ dist_unit(engine) };
        const float radius{
            std::cbrt(min_radius_cube + radius_unit * (max_radius_cube - min_radius_cube))
        };
        float phi{ dist_phi(engine) };
        float cos_theta{ dist_cos_theta(engine) };
        float sin_theta{ std::sqrt(1.0f - cos_theta * cos_theta) };

        std::array<float, 3> direction_vector{
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta
        };

        std::array<float, 3> sampling_position{
            reference_position.at(0) + radius * direction_vector[0],
            reference_position.at(1) + radius * direction_vector[1],
            reference_position.at(2) + radius * direction_vector[2]
        };

        out[i] = SamplingPoint{ radius, sampling_position };
    }
}

void SphereSampler::GenerateFibonacciDeterministic(
    const std::array<float, 3> & reference_position,
    const SphereDistanceRange & distance_range,
    const SphereDeterministicSamplingConfig & config,
    SamplingPointList & out) const
{
    const auto shell_radii{ BuildFibonacciShellCenters(distance_range, config.radius_bin_size) };
    out.clear();

    std::size_t total_point_count{ 0 };
    for (const double radius : shell_radii)
    {
        total_point_count += GetFibonacciSampleCountForRadius(radius, distance_range, config);
    }
    out.reserve(total_point_count);

    const double golden_angle{ Constants::pi * (3.0 - std::sqrt(5.0)) };

    // Shell radii are bin centers within the configured range.
    // Each shell uses the same deterministic Fibonacci sphere pattern and reports its shell radius.
    for (const double radius : shell_radii)
    {
        const std::size_t shell_sample_count{
            GetFibonacciSampleCountForRadius(radius, distance_range, config)
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

            std::array<float, 3> sampling_position{
                reference_position.at(0) + static_cast<float>(radius * x),
                reference_position.at(1) + static_cast<float>(radius * y),
                reference_position.at(2) + static_cast<float>(radius * z)
            };

            out.emplace_back(SamplingPoint{
                static_cast<float>(radius),
                sampling_position
            });
        }
    }
}

void SphereSampler::GenerateRadiusUniformRandom(
    const std::array<float, 3> & reference_position,
    const SphereDistanceRange & distance_range,
    const SphereRandomSamplingConfig & config,
    SamplingPointList & out) const
{
    out.resize(config.sample_count);
    
    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_radius(
        static_cast<float>(distance_range.min), static_cast<float>(distance_range.max)
    );
    std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_cos_theta(-1.0f, 1.0f);

    for (unsigned int i = 0; i < config.sample_count; i++)
    {
        const float radius{ dist_radius(engine) };
        float phi{ dist_phi(engine) };
        float cos_theta{ dist_cos_theta(engine) };
        float sin_theta{ std::sqrt(1.0f - cos_theta * cos_theta) };

        std::array<float, 3> direction_vector{
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta
        };

        std::array<float, 3> sampling_position{
            reference_position.at(0) + radius * direction_vector[0],
            reference_position.at(1) + radius * direction_vector[1],
            reference_position.at(2) + radius * direction_vector[2]
        };

        out[i] = SamplingPoint{ radius, sampling_position };
    }
}
