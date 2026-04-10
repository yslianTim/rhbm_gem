#pragma once

#include <array>
#include <vector>
#include <tuple>

#include <rhbm_gem/utils/math/SamplerBase.hpp>

class SphereSampler : public SamplerBase
{
    unsigned int m_sampling_size;
    double m_distance_min, m_distance_max, m_distance_step;

public:
    SphereSampler();
    ~SphereSampler() = default;

    void Print() const override;
    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & reference_position,
        const std::array<float, 3> & axis_vector) const override;
    unsigned int GetSamplingSize() const override { return m_sampling_size; }
    void SetSamplingSize(unsigned int value) override { m_sampling_size = value; }
    void SetDistanceRangeMinimum(double value) { m_distance_min = value; }
    void SetDistanceRangeMaximum(double value) { m_distance_max = value; }

private:
    void RunVolumeUniformRandomSamplingMethod(
        const std::array<float, 3> & reference_position,
        std::vector<std::tuple<float, std::array<float, 3>>> & out) const;
    void RunRadiusUniformRandomSamplingMethod(
        const std::array<float, 3> & reference_position,
        std::vector<std::tuple<float, std::array<float, 3>>> & out) const;

};
