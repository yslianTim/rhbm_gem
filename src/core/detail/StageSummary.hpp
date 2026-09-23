#pragma once
#include <string>
#include <vector>
namespace rhbm_gem {
class ModelObject;
class AtomObject;
namespace core {
struct StageProvenance { std::string estimator, peeling_mode; };
StageProvenance CollectStageProvenance(const std::vector<const AtomObject *> &);
std::string BuildSecondStageSpotSummary(const ModelObject &);
}
}
