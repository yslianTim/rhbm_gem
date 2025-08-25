#pragma once

#include <array>
#include <vector>
#include <tuple>

class SphereSampler
{
    unsigned int m_sampling_size;
    double m_distance_min, m_distance_max;

public:
    SphereSampler(void);
    ~SphereSampler() = default;

    void Print(void) const;
    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & reference_position) const;
    void SetSamplingSize(unsigned int value) { m_sampling_size = value; }
    void SetDistanceRangeMinimum(double value) { m_distance_min = value; }
    void SetDistanceRangeMaximum(double value) { m_distance_max = value; }

};
