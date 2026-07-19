#pragma once

#include <tuple>

namespace rhbm_gem {

class MapObject;

namespace core {

// Returns (height, offset) for height * exp(-0.5 * (r / sigma)^2) + offset.
std::tuple<float, float> GetReferenceGaussianParameters(const MapObject & map_object);

} // namespace core

} // namespace rhbm_gem
