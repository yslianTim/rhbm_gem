#pragma once

#include <array>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

class SphereSampler
{
    SphereSamplingMethod m_method;
    double m_distance_min{ 0.0 };
    double m_distance_max{ 1.0 };
    unsigned int m_sample_count{ 50 };
    double m_radius_bin_size{ 1.0 };
    bool m_vary_with_radius{ false };

public:
    ~SphereSampler() = default;

    static SphereSampler AnalysisDefault(SphereSamplingMethod method);

    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & center_position) const;

private:
    SphereSampler(
        SphereSamplingMethod method,
        double distance_min,
        double distance_max,
        unsigned int sample_count,
        double radius_bin_size,
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
