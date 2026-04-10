#pragma once

#include <rhbm_gem/utils/math/SamplingTypes.hpp>

class GridSampler
{
    unsigned int m_sampling_size;
    float m_window_size;
    std::array<float, 3> m_reference_u_vector;

public:
    GridSampler();
    ~GridSampler() = default;

    void Print() const;
    SamplingPointList GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & plane_normal) const;
    unsigned int GetGridResolution() const { return m_sampling_size; }
    unsigned int GetPointCount() const { return m_sampling_size * m_sampling_size; }
    void SetGridResolution(unsigned int value) { m_sampling_size = value; }
    void SetWindowSize(float value) { m_window_size = value; }
    void SetReferenceUVector(const std::array<float, 3> & value) { m_reference_u_vector = value; }

};
