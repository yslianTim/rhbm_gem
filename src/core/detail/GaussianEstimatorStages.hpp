#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core {

enum class LocalFittingPass
{
    FirstStage,
    ThirdStage
};

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    LocalFittingPass pass);
void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    LocalFittingPass pass);
void RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options);

} // namespace rhbm_gem::core
