#pragma once

#include <array>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

class SphereSampler
{
    unsigned int m_sampling_size;
    double m_distance_min, m_distance_max;

public:
    SphereSampler();
    ~SphereSampler() = default;

    void Print() const;

    SamplingPointList GenerateSamplingPoints(
        const std::array<float, 3> & reference_position) const;
    unsigned int GetSampleCount() const { return m_sampling_size; }
    unsigned int GetPointCount() const { return m_sampling_size; }
    void SetSampleCount(unsigned int value) { m_sampling_size = value; }

    // Transitional compatibility wrappers. Prefer GetSampleCount()/SetSampleCount().
    unsigned int GetSamplingSize() const { return GetSampleCount(); }
    void SetSamplingSize(unsigned int value) { SetSampleCount(value); }

    void SetDistanceRangeMinimum(double value) { m_distance_min = value; }
    void SetDistanceRangeMaximum(double value) { m_distance_max = value; }

private:
    void RunVolumeUniformRandomSamplingMethod(
        const std::array<float, 3> & reference_position,
        SamplingPointList & out) const;
    void RunRadiusUniformRandomSamplingMethod(
        const std::array<float, 3> & reference_position,
        SamplingPointList & out) const;

};
