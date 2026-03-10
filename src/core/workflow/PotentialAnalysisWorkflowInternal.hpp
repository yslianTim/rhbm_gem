#pragma once

#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>

namespace rhbm_gem {

class MapObject;
class ModelObject;

namespace detail {

void RunAtomSamplingWorkflow(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options);

void RunAtomGroupingWorkflow(ModelObject & model_object);

void RunLocalAtomFittingWorkflow(
    ModelObject & model_object,
    const PotentialAnalysisCommandOptions & options,
    double universal_alpha_r);

} // namespace detail

} // namespace rhbm_gem
