#pragma once

#include <tuple>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class MapObject;
class ModelObject;

namespace core {

// Returns (height, offset) for height * exp(-0.5 * (r / sigma)^2) + offset.
std::tuple<float, float> GetReferenceGaussianParameters(const MapObject & map_object);

// Returns accepted MapQ-style spiral-sphere points for one radial shell.
// The result may exceed num_points and is empty if all 50 attempts fail.
SamplingPointList GetRadialPointsForQScore(
    const AtomObject & atom,
    const ModelObject & model,
    double radius,
    int num_points);

} // namespace core

} // namespace rhbm_gem
