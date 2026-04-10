#pragma once

#include <cstddef>
#include <variant>

#include <rhbm_gem/utils/math/SamplingTypes.hpp>

enum class SphereSamplingMethod
{
    RadiusUniformRandom,
    VolumeUniformRandom
};

struct SphereDistanceRange
{
    double min{ 0.0 };
    double max{ 1.0 };
};

struct SphereRandomSamplingConfig
{
    unsigned int sample_count{ 10 };
};

struct SphereSamplingProfile
{
    SphereSamplingMethod method{ SphereSamplingMethod::RadiusUniformRandom };
    SphereDistanceRange distance_range{};
    std::variant<SphereRandomSamplingConfig> method_config{
        SphereRandomSamplingConfig{}
    };

    static SphereSamplingProfile RadiusUniformRandom(
        SphereDistanceRange range,
        unsigned int sample_count);
    static SphereSamplingProfile VolumeUniformRandom(
        SphereDistanceRange range,
        unsigned int sample_count);
};

// SphereSampler keeps GenerateSamplingPoints() as the stable public entry point.
// Sampling behavior is selected through SphereSamplingProfile
class SphereSampler
{
    SphereSamplingProfile m_profile;

public:
    SphereSampler();
    ~SphereSampler() = default;

    void Print() const;
    void SetSamplingProfile(const SphereSamplingProfile & profile);
    const SphereSamplingProfile & GetSamplingProfile() const { return m_profile; }
    std::size_t GetExpectedPointCount() const;

    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & reference_position) const;

private:
    void GenerateVolumeUniformRandom(
        const std::array<float, 3> & reference_position,
        const SphereDistanceRange & distance_range,
        const SphereRandomSamplingConfig & config,
        SamplingPointList & out) const;
    void GenerateRadiusUniformRandom(
        const std::array<float, 3> & reference_position,
        const SphereDistanceRange & distance_range,
        const SphereRandomSamplingConfig & config,
        SamplingPointList & out) const;

};
