#pragma once

#include "HaloMedianReconstruction.hpp"
#include "PostFitPeeling.hpp"
#include <rhbm_gem/data/object/PostJointQScoreResult.hpp>

namespace rhbm_gem::core::detail {

std::vector<std::pair<int, PostJointQScoreResult>> BuildPostJointQScores(
    const MapObject &, const ModelObject &, const JointProblem &, const JointAnalysisResult &,
    const HaloMedianBuildResult *, JointPostprocessingMetrics &);

} // namespace rhbm_gem::core::detail
