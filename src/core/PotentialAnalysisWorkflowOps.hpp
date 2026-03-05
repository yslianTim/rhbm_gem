#pragma once

namespace rhbm_gem {

struct PotentialAnalysisCommandOptions;
class MapObject;
class ModelObject;

void RunAtomSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options);
void RunAtomGrouping(ModelObject & model_object);
void RunLocalAtomFitting(
    ModelObject & model_object,
    const PotentialAnalysisCommandOptions & options,
    double universal_alpha_r);

void RunBondSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options);
void RunBondGrouping(ModelObject & model_object);

} // namespace rhbm_gem
