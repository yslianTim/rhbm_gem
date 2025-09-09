#pragma once

#include <array>
#include <vector>
#include <tuple>

#include "SamplerBase.hpp"

class SphereSampler : public SamplerBase
{
    unsigned int m_sampling_size;
    double m_distance_min, m_distance_max;

public:
    SphereSampler(void);
    ~SphereSampler() = default;

    void Print(void) const override;
    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & axis_vector) const override;
    unsigned int GetSamplingSize(void) const { return m_sampling_size; }
    void SetSamplingSize(unsigned int value) { m_sampling_size = value; }
    void SetDistanceRangeMinimum(double value) { m_distance_min = value; }
    void SetDistanceRangeMaximum(double value) { m_distance_max = value; }

};
