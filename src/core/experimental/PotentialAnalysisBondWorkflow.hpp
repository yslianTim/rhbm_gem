#pragma once

namespace rhbm_gem {

class MapObject;
class ModelObject;
struct PotentialAnalysisCommandOptions;

} // namespace rhbm_gem

namespace rhbm_gem::experimental {

void RunPotentialAnalysisBondWorkflow(
    ModelObject & model_object,
    MapObject & map_object,
    const PotentialAnalysisCommandOptions & options);

} // namespace rhbm_gem::experimental
