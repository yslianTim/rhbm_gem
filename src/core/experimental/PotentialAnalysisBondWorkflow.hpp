#pragma once

namespace rhbm_gem {

class MapObject;
class ModelObject;
struct PotentialAnalysisRequest;

} // namespace rhbm_gem

namespace rhbm_gem::experimental {

void RunPotentialAnalysisBondWorkflow(
    ModelObject & model_object,
    MapObject & map_object,
    const PotentialAnalysisRequest & options,
    int thread_size);

} // namespace rhbm_gem::experimental
