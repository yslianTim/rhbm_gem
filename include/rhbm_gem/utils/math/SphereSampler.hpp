#pragma once

#include <array>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::sphere_sampler {

SamplingPointList GenerateSamplingPointList(
    const std::array<double, 3> & center_position,
    SphereSamplingMethod method);

SamplingPointList GenerateVolumeUniformRandom(const std::array<double, 3> & center_position);
SamplingPointList GenerateRadiusUniformRandom(const std::array<double, 3> & center_position);
SamplingPointList GenerateFibonacciDeterministic(const std::array<double, 3> & center_position);

} // namespace rhbm_gem::sphere_sampler
