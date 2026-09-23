#pragma once
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#include <string_view>
#endif

namespace rhbm_gem::core::detail {
struct FittingWorkset
{
    std::vector<AtomObject *> contributors;
    std::vector<bool> target_mask;
};
FittingWorkset MakeJointFittingWorkset(ModelObject &, const JointProblem &);
JointInitialization RunContributorFirstStage(MapObject &, ModelObject &,
    const FittingWorkset &, const FitOptions &);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using FirstStageObserver = std::function<void(int, std::string_view)>;
FirstStageObserver & FirstStageObserverForTesting();
#endif
}
