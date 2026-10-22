#include <rhbm_gem/core/GaussianEstimator.hpp>

static_assert(sizeof(rhbm_gem::core::FitOptions) > 0);
static_assert(rhbm_gem::core::FitOptions{}.estimator ==
    rhbm_gem::core::PotentialEstimator::JOINT_COMPONENTS);
