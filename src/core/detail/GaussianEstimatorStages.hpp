#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core {

enum class LocalFittingPass
{
    FirstStage,
    ThirdStage
};

struct SecondStageLocalFittingInternalOptions
{
    bool enable_end_to_end_health_status{ true };
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
    const FitOptions & options,
    const SecondStageLocalFittingInternalOptions & internal_options =
        SecondStageLocalFittingInternalOptions{});

} // namespace rhbm_gem::core
