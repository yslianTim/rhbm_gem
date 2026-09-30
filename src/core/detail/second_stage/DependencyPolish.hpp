#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SecondStageDiagnostics.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"
#include "core/detail/second_stage/JointFitting.hpp"

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

class PerformanceCounters;

struct FinalDependencyPolishResult
{
    FitState state{};
    std::optional<ObjectiveBreakdown> objective{};
    FinalDependencyPolishDiagnostic diagnostic{};
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
    PerformanceCounters & performance_counters);

} // namespace rhbm_gem::core::detail
