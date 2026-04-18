#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/SamplingPointAcceptanceMask.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

inline std::vector<float> BuildLocalPotentialSampleScoreList(
    std::size_t sample_count,
    float default_score = 1.0f)
{
    if (!std::isfinite(default_score))
    {
        throw std::invalid_argument(
            "LocalPotentialSample score must be finite.");
    }
    return std::vector<float>(sample_count, default_score);
}

inline std::vector<float> BuildLocalPotentialSampleScoreList(
    const SamplingPointList & point_list,
    const std::vector<Eigen::VectorXd> & reject_direction_list,
    double angle = 0.0)
{
    // Today the score policy is a binary acceptance mask. Keeping this decision behind a
    // LocalPotentialSample-specific entry point gives us one place to extend later.
    return BuildSamplingPointAcceptanceMask(point_list, reject_direction_list, angle);
}

} // namespace rhbm_gem
