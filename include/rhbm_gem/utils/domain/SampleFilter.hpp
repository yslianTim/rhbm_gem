#pragma once

#include <array>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::sample_filter {

// Strict QScore distance comparison; callers own neighbor discovery and radius scaling.
bool IsInsideNeighborExclusion(
    const std::array<double, 3> & sample_position,
    const std::array<double, 3> & neighbor_position,
    double exclusion_radius_squared);

void FilterSamplingPointList(
    SamplingPointList & sample_point_list,
    const std::array<double, 3> & local_position,
    const std::vector<std::array<double, 3>> & reject_position_list,
    double angle = 30.0);

LocalPotentialSampleList BuildMedianResponseSampleEntriesByRadius(
    const LocalPotentialSampleList & sample_entries);

} // namespace rhbm_gem::sample_filter
