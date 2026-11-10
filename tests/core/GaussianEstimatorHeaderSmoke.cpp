#include "core/estimator/GaussianEstimator.hpp"

static_assert(rhbm_gem::core::FitOptions{}.initial_stage_sampling_method ==
    SphereSamplingMethod::FibonacciDeterministic);

template <class T>
concept HasEstimatorField = requires(T value) { value.estimator; };

static_assert(!HasEstimatorField<rhbm_gem::core::FitOptions>);

using InitialStageRunner = void (*)(rhbm_gem::ModelObject &, const rhbm_gem::core::FitOptions &);
static_assert(static_cast<InitialStageRunner>(&rhbm_gem::core::RunInitialStageLocalAlphaTraining) != nullptr);
static_assert(static_cast<InitialStageRunner>(&rhbm_gem::core::RunInitialStageLocalFitting) != nullptr);
