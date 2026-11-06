#pragma once

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <array>

namespace rhbm_gem { class AtomObject; class ModelObject; }

namespace rhbm_gem::core::detail {
// Strict QScore distance comparison; callers own neighbor discovery and radius scaling.
bool IsInsideNeighborExclusion(
    const std::array<double, 3> & sample_position,
    const std::array<double, 3> & neighbor_position,
    double exclusion_radius_squared);

SamplingPointList SampleQScoreShell(const AtomObject &, const ModelObject &, double radius, int num_points);
} // namespace rhbm_gem::core::detail
