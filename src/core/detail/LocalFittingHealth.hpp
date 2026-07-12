#pragma once

#include <stdexcept>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

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

inline bool IsJointOffsetSolveProgressEligible(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return true;
    case JointOffsetSolveStatus::SystemBuildFailed:
    case JointOffsetSolveStatus::EmptySystem:
    case JointOffsetSolveStatus::InitialSolveFailed:
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return false;
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

inline bool IsJointOffsetSolveStationarityEligible(JointOffsetSolveStatus status)
{
    return status == JointOffsetSolveStatus::Converged;
}

inline bool IsJointOffsetSolveHardFailure(JointOffsetSolveStatus status)
{
    return !IsJointOffsetSolveProgressEligible(status);
}

inline bool IsLocalGaussianRefitStatusHealthEligible(RHBMEstimationStatus status)
{
    switch (status)
    {
    case RHBMEstimationStatus::SUCCESS:
    case RHBMEstimationStatus::MAX_ITERATIONS_REACHED:
        return true;
    case RHBMEstimationStatus::SINGLE_MEMBER:
    case RHBMEstimationStatus::INSUFFICIENT_DATA:
    case RHBMEstimationStatus::NUMERICAL_FALLBACK:
        return false;
    }
    throw std::logic_error("Local Gaussian refit status is invalid.");
}

} // namespace rhbm_gem::core::detail
