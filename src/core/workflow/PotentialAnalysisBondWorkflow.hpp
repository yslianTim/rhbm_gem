#pragma once

#include <rhbm_gem/core/PotentialAnalysisCommand.hpp>

namespace rhbm_gem {

class ModelObject;
class MapObject;

namespace detail {

struct PotentialAnalysisBondWorkflowContext
{
    ModelObject & model_object;
    MapObject & map_object;
    const PotentialAnalysisCommand::Options & options;
};

void RunPotentialAnalysisBondWorkflow(
    const PotentialAnalysisBondWorkflowContext & context);

} // namespace detail

} // namespace rhbm_gem
