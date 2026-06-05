#pragma once

#include <array>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::sample_filter {

void FilterSamplingPointList(
    SamplingPointList & sample_point_list,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
    double angle = 30.0);

LocalPotentialSampleList BuildMedianResponseSampleEntriesByRadius(
    const LocalPotentialSampleList & sample_entries);

LocalPotentialSampleList BuildResponseShiftedSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    double response_shift);

} // namespace rhbm_gem::sample_filter
