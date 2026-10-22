#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"
#include "core/detail/second_stage/JointFitting.hpp"

namespace rhbm_gem::core::detail {

struct FinalDependencyPolishResult
{
    FitState state{};
    std::optional<ObjectiveBreakdown> objective{};
    bool accepted{ false };
};

FinalDependencyPolishResult RunFinalDependencyPolish(
    const SecondStageContext & context,
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & objective_domain,
    const SuspiciousBlockActivity & block_activity,
    const FitState & base_state,
    BoundaryJointCorrectionWorkspaceMap & workspace_by_key);

} // namespace rhbm_gem::core::detail
