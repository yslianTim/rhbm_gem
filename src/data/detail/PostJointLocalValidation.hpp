#pragma once

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <string>

namespace rhbm_gem {
namespace data_internal {

struct PostJointLocalEvidenceAssessment
{
    bool valid{};
    bool group_eligible{};
    std::string reason;
};

PostJointLocalEvidenceAssessment AssessPostJointLocalResult(
    const PostJointLocalResult &, int atom_id, const LocalStageEstimate & joint_endpoint,
    const PostFitPeelingResult * peeling, const LocalPotentialSampleList & raw_samples);

} // namespace data_internal
} // namespace rhbm_gem
