#pragma once

#include <cmath>

namespace rhbm_gem::core::joint_component {

struct JointConvergenceCriteria
{
    double global_ac_kkt{};
    // Threshold for the bound-projected gradient in eta = log(tau).
    double global_eta_gradient{};
};

inline constexpr JointConvergenceCriteria StrictJointConvergenceCriteria{1e-10,1e-12};
inline constexpr JointConvergenceCriteria PracticalJointConvergenceCriteria{1e-7,1e-9};

inline bool IsValidJointConvergenceCriteria(const JointConvergenceCriteria & criteria)
{
    return std::isfinite(criteria.global_ac_kkt) && criteria.global_ac_kkt>=0.0 &&
        std::isfinite(criteria.global_eta_gradient) && criteria.global_eta_gradient>=0.0;
}

} // namespace rhbm_gem::core::joint_component
