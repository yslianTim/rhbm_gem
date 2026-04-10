#pragma once

#include <rhbm_gem/utils/math/SamplingTypes.hpp>

class SphereSampler
{
    unsigned int m_sampling_size;
    double m_distance_min, m_distance_max;

public:
    SphereSampler();
    ~SphereSampler() = default;

    void Print() const;
    void SetDistanceRange(double min_value, double max_value);

    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & reference_position) const;
    unsigned int GetSampleCount() const { return m_sampling_size; }
    void SetSampleCount(unsigned int value) { m_sampling_size = value; }

private:
    void RunVolumeUniformRandomSamplingMethod(
        const std::array<float, 3> & reference_position,
        SamplingPointList & out) const;
    void RunRadiusUniformRandomSamplingMethod(
        const std::array<float, 3> & reference_position,
        SamplingPointList & out) const;

};
