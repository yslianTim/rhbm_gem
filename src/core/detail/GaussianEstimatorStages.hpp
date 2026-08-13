#pragma once

#include <cstddef>

#include <rhbm_gem/core/GaussianEstimator.hpp>
namespace rhbm_gem::core {

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);
void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);
void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);

bool RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options);

} // namespace rhbm_gem::core
