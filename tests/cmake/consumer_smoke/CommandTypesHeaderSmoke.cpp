#include <rhbm_gem/core/CommandTypes.hpp>

template <class T>
concept HasEstimatorField = requires(T value) { value.estimator; };

static_assert(!HasEstimatorField<rhbm_gem::core::PotentialAnalysisRequest>);
