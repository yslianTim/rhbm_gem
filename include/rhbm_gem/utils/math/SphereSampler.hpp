#pragma once

#include <array>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

class SphereSampler
{
    SphereSamplingMethod m_method;

public:
    explicit SphereSampler(SphereSamplingMethod method);
    ~SphereSampler() = default;

    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & center_position) const;

private:
    void GenerateVolumeUniformRandom(
        const std::array<float, 3> & center_position,
        SamplingPointList & out) const;
    void GenerateRadiusUniformRandom(
        const std::array<float, 3> & center_position,
        SamplingPointList & out) const;
    void GenerateFibonacciDeterministic(
        const std::array<float, 3> & center_position,
        SamplingPointList & out) const;

};
