#pragma once

#include <array>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::sphere_sampler {

SamplingPointList GenerateSamplingPointList(
    const std::array<float, 3> & center_position,
    SphereSamplingMethod method);

SamplingPointList GenerateVolumeUniformRandom(const std::array<float, 3> & center_position);
SamplingPointList GenerateRadiusUniformRandom(const std::array<float, 3> & center_position);
SamplingPointList GenerateFibonacciDeterministic(const std::array<float, 3> & center_position);

} // namespace rhbm_gem::sphere_sampler
