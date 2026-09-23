#pragma once
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <map>
#include <string>

namespace rhbm_gem {
class ModelObject;
namespace data_internal {
std::map<int, LocalStageEstimate> BuildJointStageEstimates(
    const JointAnalysisResult &, const std::string & run_id);
void ApplyJointStageEstimates(ModelObject &, const JointAnalysisResult &, const std::string & run_id);
}
}
