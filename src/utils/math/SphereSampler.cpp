#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>
#include <random>
#include <string_view>
#include <stdexcept>
#include <sstream>

namespace {

void ValidateSphereSamplingInputs(
    SphereDistanceRange range,
    SphereRandomSamplingConfig random_config)
{
    if (range.min > range.max)
    {
        throw std::invalid_argument("SphereSampler: distance minimum greater than maximum");
    }
    if (range.min < 0.0 || range.max < 0.0)
    {
        throw std::invalid_argument("SphereSampler: distance range cannot be negative");
    }
    (void)random_config;
}

std::string_view GetSphereSamplingMethodName(SphereSamplingMethod method)
{
    switch (method)
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            return "RadiusUniformRandom";
        case SphereSamplingMethod::VolumeUniformRandom:
            return "VolumeUniformRandom";
    }

    return "Unknown";
}

} // namespace

SphereSamplingProfile::SphereSamplingProfile(
    SphereSamplingMethod method,
    SphereDistanceRange range,
    SphereRandomSamplingConfig random_config) :
    m_method{ method },
    m_distance_range{ range },
    m_random_config{ random_config }
{
    ValidateSphereSamplingInputs(range, random_config);
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

SphereSampler::SphereSampler() :
    m_profile{ SphereSamplingProfile::RadiusUniformRandom(SphereDistanceRange{ 0.0, 1.0 }, 10) }
{
}

void SphereSampler::Print() const
{
    const auto & config{ m_profile.GetRandomConfig() };
    const auto & distance_range{ m_profile.GetDistanceRange() };

    std::ostringstream oss;
    oss << "SphereSampler Configuration:\n"
        << " - Sampling method: " << GetSphereSamplingMethodName(m_profile.GetMethod()) << '\n'
        << " - Sample count: " << config.sample_count << '\n'
        << " - Distance range: ["
        << StringHelper::ToStringWithPrecision<double>(distance_range.min, 1) << ", "
        << StringHelper::ToStringWithPrecision<double>(distance_range.max, 1)
        << "] Angstrom.";
    Logger::Log(LogLevel::Info, oss.str());
}

void SphereSampler::SetSamplingProfile(const SphereSamplingProfile & profile)
{
    m_profile = profile;
}

std::size_t SphereSampler::GetExpectedPointCount() const
{
    return m_profile.GetRandomConfig().sample_count;
}

SamplingPointList SphereSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position) const
{
    SamplingPointList out;
    const auto & config{ m_profile.GetRandomConfig() };
    const auto & distance_range{ m_profile.GetDistanceRange() };
    switch (m_profile.GetMethod())
    {
        case SphereSamplingMethod::RadiusUniformRandom:
            GenerateRadiusUniformRandom(
                reference_position,
                distance_range,
                config,
                out);
            break;
        case SphereSamplingMethod::VolumeUniformRandom:
            GenerateVolumeUniformRandom(
                reference_position,
                distance_range,
                config,
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

        out[i] = std::make_tuple(radius, sampling_position);
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

        out[i] = std::make_tuple(radius, sampling_position);
    }
}
