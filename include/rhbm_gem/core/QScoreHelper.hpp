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

// Calculates the MapQ-style mean-subtracted correlation for one atom.
double CalculateQScoreForAtom(
    const AtomObject & atom,
    const MapObject & map,
    const ModelObject & model,
    double sigma,
    double max_radius,
    double radial_step,
    int num_points);

// Returns the mean MapQ-style Q-score across all non-hydrogen model atoms.
double CalculateAverageQScores(const MapObject & map, const ModelObject & model);

} // namespace core

} // namespace rhbm_gem
