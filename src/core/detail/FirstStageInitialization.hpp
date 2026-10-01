#pragma once
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "FittingWorkset.hpp"
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#include <string_view>
#endif

namespace rhbm_gem { class ModelAnalysisEditor; }

namespace rhbm_gem::core::detail {
void TrainLocalAlphaForAtom(ModelAnalysisEditor &, const FitOptions &, FittingStage, AtomObject &);
void ApplyJointSeedFallback(JointInitialization &);
enum class FirstStageMode { ExistingSamplesBatch, SampleContributorsIsolated };
LocalGaussianResult FitFirstStageAtom(const AtomObject &, const FitOptions &);
// Both modes consume prepared samples; Joint mode also collects initialization provenance.
JointInitialization RunFirstStage(ModelObject &, const FittingWorkset &, const FitOptions &, FirstStageMode);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using FirstStageObserver = std::function<void(int, std::string_view)>;
FirstStageObserver & FirstStageObserverForTesting();
#endif
}
