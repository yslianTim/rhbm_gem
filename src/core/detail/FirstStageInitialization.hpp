#pragma once
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "joint_component/JointProgress.hpp"
#include <optional>
#include <string_view>
#include <vector>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#include <string_view>
#endif

namespace rhbm_gem { class MapObject; class ModelAnalysisEditor; }

namespace rhbm_gem::core::detail {
struct FittingWorkset;
void TrainLocalAlphaForAtom(ModelAnalysisEditor &, const FitOptions &, FittingStage, AtomObject &);
void ApplyJointSeedFallback(JointInitialization &);
std::optional<JointFixedHaloParameter> MakeFixedHaloParameter(
    const JointInitializationAtom &, std::size_t contributor_index,
    std::string_view expected_atom_id, std::string & failure_reason);
RHBMBetaEstimateResult RecoverFailedFirstStageEstimate(const RHBMMemberDataset &, double alpha,
    const RHBMExecutionOptions &, const RHBMBetaEstimateResult & native);
LocalGaussianResult FitFirstStageAtom(const AtomObject &, const FitOptions &, bool enable_fixed_halo_recovery = false);
void RunBatchFirstStageFromPreparedSamples(
    ModelObject &, const std::vector<AtomObject *> &, const FitOptions &);
JointInitialization RunJointFirstStageInitializationFromPreparedSamples(
    ModelObject &, const FittingWorkset &, const FitOptions &);
JointInitialization PrepareJointWidthInitialization(
    MapObject &, ModelObject &, const FittingWorkset &, const FitOptions &);
JointFitResult FitJointComponentsInitializedWithProgress(
    const JointProblem &, const JointInitialization &,
    const joint_component::JointProgressObserver &);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using FirstStageObserver = std::function<void(int, std::string_view)>;
FirstStageObserver & FirstStageObserverForTesting();
#endif
}
