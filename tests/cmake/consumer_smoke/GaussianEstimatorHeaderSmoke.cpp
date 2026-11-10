#include <rhbm_gem/core/GaussianEstimator.hpp>

static_assert(rhbm_gem::core::FitOptions{}.first_stage_sampling_method ==
    SphereSamplingMethod::FibonacciDeterministic);

template <class T>
concept HasEstimatorField = requires(T value) { value.estimator; };

static_assert(!HasEstimatorField<rhbm_gem::core::FitOptions>);
