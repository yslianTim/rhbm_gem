#pragma once

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <string>

namespace rhbm_gem {
struct JointAnalysisResult;
namespace data_internal {

struct PostJointLocalEvidenceAssessment
{
    bool valid{};
    bool group_eligible{};
    std::string reason;
};

std::string ValidatePostFitPeelingResult(const PostFitPeelingResult &,
    const LocalStageEstimate &, const JointAnalysisResult * = nullptr,
    const std::array<double,3> * target_position = nullptr);

PostJointLocalEvidenceAssessment AssessPostJointLocalResult(
    const PostJointLocalResult &, int atom_id, const LocalStageEstimate & joint_endpoint,
    const PostFitPeelingResult * peeling);

} // namespace data_internal
} // namespace rhbm_gem
