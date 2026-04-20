#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/SamplingPointAcceptanceMask.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

inline std::vector<float> BuildDefaultLocalPotentialSampleScoreList(
    std::size_t sample_count)
{
    return std::vector<float>(sample_count, 1.0f);
}

inline std::vector<float> BuildLocalPotentialSampleScoreList(
    const SamplingPointList & point_list,
    const std::vector<Eigen::VectorXd> & reject_direction_list,
    double angle = 0.0)
{
    // Callers provide score-policy inputs as plain data. This layer stays agnostic
    // to domain objects such as atoms and only owns the score calculation itself.
    return BuildSamplingPointAcceptanceMask(point_list, reject_direction_list, angle);
}

} // namespace rhbm_gem
