#pragma once

#include "core/detail/ObjectiveEvaluation.hpp"
#include "core/detail/JointFitting.hpp"

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

class SecondStageObservationSession;

struct FinalDependencyPolishResult
{
    FitState state{};
    std::optional<ObjectiveBreakdown> objective{};
    bool accepted{ false };
};

FinalDependencyPolishResult RunFinalDependencyPolish(
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & objective_domain,
    const SuspiciousBlockActivity & block_activity,
    const FitState & base_state,
    BoundaryJointCorrectionWorkspaceMap & workspace_by_key,
    PerformanceCounters & performance_counters,
    SecondStageObservationSession * observation = nullptr);

} // namespace rhbm_gem::core::detail
