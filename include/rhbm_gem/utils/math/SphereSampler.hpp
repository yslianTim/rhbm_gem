#pragma once

#include <array>
#include <cstddef>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::sphere_sampler {

SamplingPointList GenerateSamplingPointList(
    const std::array<double, 3> & center_position,
    SphereSamplingMethod method = SphereSamplingMethod::FibonacciDeterministic);

SamplingPointList GenerateFibonacciDeterministic(const std::array<double, 3> & center_position);

// Requires a finite positive radius and at least two points.
SamplingPointList GenerateSpiralShell(
    const std::array<double, 3> & center_position,
    double radius,
    std::size_t point_count);

} // namespace rhbm_gem::sphere_sampler
