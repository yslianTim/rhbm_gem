#pragma once

#include <stdexcept>

namespace rhbm_gem::core::detail {

enum class JointOffsetSolveStatus
{
    Converged,
    SystemBuildFailed,
    EmptySystem,
    InitialSolveFailed,
    IrlsSolveFailed,
    IrlsObjectiveDeteriorated,
    IrlsMaximumIterationsReached
};

inline bool IsJointOffsetSolveStationarityEligible(JointOffsetSolveStatus status)
{
    return status == JointOffsetSolveStatus::Converged;
}

inline bool IsJointOffsetSolveHardFailure(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return false;
    case JointOffsetSolveStatus::SystemBuildFailed:
    case JointOffsetSolveStatus::EmptySystem:
    case JointOffsetSolveStatus::InitialSolveFailed:
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return true;
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

inline const char * GetJointOffsetSolveStatusText(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
        return "converged";
    case JointOffsetSolveStatus::SystemBuildFailed:
        return "system-build-failed";
    case JointOffsetSolveStatus::EmptySystem:
        return "empty-system";
    case JointOffsetSolveStatus::InitialSolveFailed:
        return "initial-solve-failed";
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return "irls-solve-failed";
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
        return "irls-objective-deteriorated";
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return "irls-maximum-iterations-reached";
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

} // namespace rhbm_gem::core::detail
