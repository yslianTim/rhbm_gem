#pragma once

#include <array>
#include <vector>
#include <tuple>

#include "SamplerBase.hpp"

class GridSampler : public SamplerBase
{
    unsigned int m_sampling_size;
    float m_window_size;
    std::array<float, 3> m_reference_u_vector;

public:
    GridSampler();
    ~GridSampler() = default;

    void Print() const override;
    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & axis_vector) const override;
    unsigned int GetSamplingSize() const override { return m_sampling_size * m_sampling_size; }
    void SetSamplingSize(unsigned int value) override { m_sampling_size = value; }
    void SetWindowSize(float value) { m_window_size = value; }
    void SetReferenceUVector(const std::array<float, 3> & value) { m_reference_u_vector = value; }

};
