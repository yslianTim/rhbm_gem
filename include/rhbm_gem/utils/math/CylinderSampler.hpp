#pragma once

#include <array>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

class CylinderSampler
{
    unsigned int m_sampling_size;
    double m_distance_min, m_distance_max;
    double m_height;

public:
    CylinderSampler();
    ~CylinderSampler() = default;

    void Print() const;
    SamplingPointList GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & axis_vector) const;
    unsigned int GetSampleCount() const { return m_sampling_size; }
    unsigned int GetPointCount() const { return m_sampling_size; }
    void SetSampleCount(unsigned int value) { m_sampling_size = value; }

    // Transitional compatibility wrappers. Prefer GetSampleCount()/SetSampleCount().
    unsigned int GetSamplingSize() const { return GetSampleCount(); }
    void SetSamplingSize(unsigned int value) { SetSampleCount(value); }

    void SetDistanceRangeMinimum(double value) { m_distance_min = value; }
    void SetDistanceRangeMaximum(double value) { m_distance_max = value; }
    void SetHeight(double value) { m_height = value; }

};
