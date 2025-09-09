#pragma once

#include <array>
#include <vector>
#include <tuple>

class CylinderSampler
{
    unsigned int m_sampling_size;
    double m_radius_min, m_radius_max;
    double m_height_min, m_height_max;

public:
    CylinderSampler(void);
    ~CylinderSampler() = default;

    void Print(void) const;
    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & axis_vector) const;
    unsigned int GetSamplingSize(void) const { return m_sampling_size; }
    void SetSamplingSize(unsigned int value) { m_sampling_size = value; }
    void SetRadiusRange(double rmin, double rmax) { m_radius_min = rmin; m_radius_max = rmax; }
    void SetHeightRange(double hmin, double hmax) { m_height_min = hmin; m_height_max = hmax; }

};
