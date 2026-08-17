#pragma once

#include <algorithm>
#include <map>
#include <stdexcept>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

#include "core/detail/JointOffsetSolveStatus.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

struct ClusterHealth
{
    explicit ClusterHealth(JointOffsetSolveStatus status)
        : joint_offset_status(status)
    {
    }

    JointOffsetSolveStatus joint_offset_status;
    bool is_refit_stationarity_eligible{ true };

    bool IsStationarityEligible() const
    {
        return IsJointOffsetSolveStationarityEligible(joint_offset_status) &&
            is_refit_stationarity_eligible;
    }
};

using ClusterHealthMap = std::map<ClusterKey, ClusterHealth>;

inline bool AreClustersStationarityEligible(const ClusterHealthMap & health_by_key)
{
    return std::all_of(
        health_by_key.begin(),
        health_by_key.end(),
        [](const auto & entry)
        {
            return entry.second.IsStationarityEligible();
        });
}

inline bool IsLocalRefitStatusStationarityEligible(RHBMEstimationStatus status)
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
