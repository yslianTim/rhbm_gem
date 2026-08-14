#pragma once

#include <algorithm>
#include <stdexcept>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

#include <map>

#include "core/detail/JointOffset.hpp"
#include "core/detail/FitStateView.hpp"

namespace rhbm_gem::core::detail {

struct LocalFittingClusterHealth
{
    JointOffsetSolveStatus joint_offset_status{ JointOffsetSolveStatus::SystemBuildFailed };
    bool is_refit_stationarity_eligible{ true };

    bool IsStationarityEligible() const
    {
        return IsJointOffsetSolveStationarityEligible(joint_offset_status) &&
            is_refit_stationarity_eligible;
    }
};

using LocalFittingClusterHealthMap = std::map<ClusterKey, LocalFittingClusterHealth>;

inline bool AreLocalFittingClustersStationarityEligible(const LocalFittingClusterHealthMap & health_by_key)
{
    return std::all_of(
        health_by_key.begin(),
        health_by_key.end(),
        [](const auto & entry)
        {
            return entry.second.IsStationarityEligible();
        });
}

inline bool IsLocalGaussianRefitStatusStationarityEligible(RHBMEstimationStatus status)
{
    switch (status)
    {
    case RHBMEstimationStatus::SUCCESS:
        return true;
    case RHBMEstimationStatus::MAX_ITERATIONS_REACHED:
    case RHBMEstimationStatus::SINGLE_MEMBER:
    case RHBMEstimationStatus::INSUFFICIENT_DATA:
    case RHBMEstimationStatus::NUMERICAL_FALLBACK:
        return false;
    }
    throw std::logic_error("Local Gaussian refit status is invalid.");
}

} // namespace rhbm_gem::core::detail
