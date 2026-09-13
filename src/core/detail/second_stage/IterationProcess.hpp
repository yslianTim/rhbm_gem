#pragma once

namespace rhbm_gem {
class ModelObject;
}

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

void RunSecondStageIterations(ModelObject & model_object, const FitOptions & options);

} // namespace rhbm_gem::core::detail
