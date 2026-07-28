#pragma once

#include <tuple>
#include <unordered_map>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class MapObject;
class ModelObject;

namespace core {

// Returns (height, offset) for height * exp(-0.5 * (r / sigma)^2) + offset.
std::tuple<double, double> GetReferenceGaussianParameters(const MapObject & map_object);

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
    double height,
    double offset,
    double sigma=0.6,
    int num_points=8);

double CalculateQScoreForAtom(
    const LocalPotentialSampleList & sampling_entries,
    double height,
    double offset,
    double sigma=0.6);

// Returns the mean MapQ-style Q-score across all non-hydrogen model atoms.
double CalculateAverageQScores(const MapObject & map, const ModelObject & model);

// Also returns each non-hydrogen atom Q-score keyed by serial ID.
double CalculateAverageQScores(
    const MapObject & map,
    const ModelObject & model,
    std::unordered_map<int, double> & q_scores_by_serial_id);

} // namespace core

} // namespace rhbm_gem
