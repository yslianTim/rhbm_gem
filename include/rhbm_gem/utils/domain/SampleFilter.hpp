#pragma once

#include <array>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::sample_filter {

SamplingPointList FilterSamplingPointList(
    const SamplingPointList & sample_point_list,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
    double angle = 0.0);

LocalPotentialSampleList FilterLocalPotentialSampleList(LocalPotentialSampleList sample_list);

} // namespace rhbm_gem::sample_filter
