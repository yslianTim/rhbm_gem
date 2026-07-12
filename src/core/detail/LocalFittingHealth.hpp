#pragma once

#include <stdexcept>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::core::detail {

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
