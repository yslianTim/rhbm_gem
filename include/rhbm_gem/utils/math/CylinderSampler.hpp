#pragma once

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
    void SetDistanceRange(double min_value, double max_value);
    SamplingPointList GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & axis_vector) const;
    unsigned int GetSampleCount() const { return m_sampling_size; }
    void SetSampleCount(unsigned int value) { m_sampling_size = value; }
    void SetHeight(double value) { m_height = value; }

};
