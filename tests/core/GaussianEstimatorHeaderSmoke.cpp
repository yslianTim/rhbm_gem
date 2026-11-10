#include "core/estimator/GaussianEstimator.hpp"

static_assert(rhbm_gem::core::FitOptions{}.first_stage_sampling_method ==
    SphereSamplingMethod::FibonacciDeterministic);

template <class T>
concept HasEstimatorField = requires(T value) { value.estimator; };

static_assert(!HasEstimatorField<rhbm_gem::core::FitOptions>);

using FirstStageRunner = void (*)(rhbm_gem::ModelObject &, const rhbm_gem::core::FitOptions &);
static_assert(static_cast<FirstStageRunner>(&rhbm_gem::core::RunFirstStageLocalAlphaTraining) != nullptr);
static_assert(static_cast<FirstStageRunner>(&rhbm_gem::core::RunFirstStageLocalFitting) != nullptr);
