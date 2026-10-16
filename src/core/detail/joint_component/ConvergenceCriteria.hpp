#pragma once

namespace rhbm_gem::core::joint_component {

struct JointConvergenceCriteria
{
    double global_ac_kkt{};
    double global_width_gradient{};
    double confirmed_log_width_change{};
};

inline constexpr JointConvergenceCriteria StrictJointConvergenceCriteria{
    1e-10,1e-12,1e-10};
inline constexpr JointConvergenceCriteria PracticalJointConvergenceCriteria{
    1e-7,1e-9,1e-6};

bool IsValidJointConvergenceCriteria(const JointConvergenceCriteria &);

} // namespace rhbm_gem::core::joint_component
