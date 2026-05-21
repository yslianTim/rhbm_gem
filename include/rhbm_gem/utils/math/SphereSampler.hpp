#pragma once

#include <array>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

struct SphereDistanceRange
{
    double min{ 0.0 };
    double max{ 1.0 };
};

class SphereSampler
{
    SphereSamplingMethod m_method{ SphereSamplingMethod::RadiusUniformRandom };
    SphereDistanceRange m_distance_range{ 0.0, 1.0 };
    unsigned int m_sample_count{ 10 };
    double m_radius_bin_size{ 1.0 };
    unsigned int m_samples_per_radius{ 10 };
    bool m_vary_with_radius{ false };

public:
    SphereSampler();
    ~SphereSampler() = default;

    static SphereSampler RadiusUniformRandom(
        SphereDistanceRange range,
        unsigned int sample_count);
    static SphereSampler VolumeUniformRandom(
        SphereDistanceRange range,
        unsigned int sample_count);
    static SphereSampler FibonacciDeterministic(
        SphereDistanceRange range,
        double radius_bin_size,
        unsigned int samples_per_radius,
        bool vary_with_radius = false);
    static SphereSampler AnalysisDefault(
        SphereSamplingMethod method,
        unsigned int sampling_size);

    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & center_position) const;

private:
    SphereSampler(
        SphereSamplingMethod method,
        SphereDistanceRange range,
        unsigned int sample_count,
        double radius_bin_size,
        unsigned int samples_per_radius,
        bool vary_with_radius);

    void GenerateVolumeUniformRandom(
        const std::array<float, 3> & center_position,
        SamplingPointList & out) const;
    void GenerateRadiusUniformRandom(
        const std::array<float, 3> & center_position,
        SamplingPointList & out) const;
    void GenerateFibonacciDeterministic(
        const std::array<float, 3> & center_position,
        SamplingPointList & out) const;

};
