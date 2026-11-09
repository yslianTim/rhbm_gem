#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace rhbm_gem::core::joint_component {

struct ActiveSetStep
{
    bool valid{};
    std::string status;
    Eigen::VectorXd step,gradient,multipliers;
    // Solver working bounds, not an exhaustive mask of geometrically attained bounds.
    std::vector<bool> working_bound_mask;
    double damping{},scaled_step_norm{},kkt_inf_norm{};
    std::size_t iterations{},working_set_releases{};
    std::size_t damping_solves{},bracket_solves{},bisection_solves{},damping_qr_factorizations{},
        damping_active_set_iterations{},damping_working_set_releases{};
};

// Selects the LM damping so the same constrained subproblem also satisfies
// ||D p|| <= trust_radius. The returned step is the step used by the model.
ActiveSetStep SolveBoundedTrustRegionStep(const Eigen::MatrixXd &,const Eigen::VectorXd &,
    const Eigen::VectorXd &,const Eigen::VectorXd &,double,double);

} // namespace rhbm_gem::core::joint_component
