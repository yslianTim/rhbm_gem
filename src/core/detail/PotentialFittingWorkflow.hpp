#pragma once

#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <vector>

namespace rhbm_gem::core::detail {

void RunTwoStageFromPreparedSamples(
    ModelObject & model,
    const FitOptions & options);

void RunTwoStageFromPreparedSamples(
    ModelObject & model,
    const std::vector<AtomObject *> & atoms,
    const FitOptions & options);

} // namespace rhbm_gem::core::detail
