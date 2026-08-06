#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core {

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    bool use_peeling_sampling_entries);
void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    bool use_peeling_sampling_entries);
void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options,
    bool use_peeling_sampling_entries);
void RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options);

} // namespace rhbm_gem::core
