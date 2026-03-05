#pragma once

#include <array>
#include <tuple>
#include <vector>

class SamplerBase;

namespace rhbm_gem {

class MapObject;

std::vector<std::tuple<float, float>> SampleMapValues(
    const MapObject & map_object,
    const ::SamplerBase & sampler,
    const std::array<float, 3> & position,
    const std::array<float, 3> & axis_vector = { 0.0f, 0.0f, 0.0f });

} // namespace rhbm_gem
