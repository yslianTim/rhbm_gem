#pragma once
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "FittingWorkset.hpp"
#include <vector>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#include <string_view>
#endif

namespace rhbm_gem { class ModelAnalysisEditor; }

namespace rhbm_gem::core::detail {
void TrainLocalAlphaForAtom(ModelAnalysisEditor &, const FitOptions &, FittingStage, AtomObject &);
void ApplyJointSeedFallback(JointInitialization &);
LocalGaussianResult FitFirstStageAtom(const AtomObject &, const FitOptions &);
void RunBatchFirstStageFromPreparedSamples(
    ModelObject &, const std::vector<AtomObject *> &, const FitOptions &);
JointInitialization RunJointFirstStageInitializationFromPreparedSamples(
    ModelObject &, const FittingWorkset &, const FitOptions &);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using FirstStageObserver = std::function<void(int, std::string_view)>;
FirstStageObserver & FirstStageObserverForTesting();
#endif
}
