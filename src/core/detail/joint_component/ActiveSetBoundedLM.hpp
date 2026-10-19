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
    std::vector<bool> active;
    double damping{},scaled_step_norm{},kkt_inf_norm{};
    std::size_t iterations{},releases{};
};

// Solves the damped profile least-squares subproblem with upper step bounds:
// min 0.5 ||J p + r||^2 + 0.5 lambda ||D p||^2, p <= upper_step.
// Active variables are released whenever their upper-bound multiplier has the
// wrong sign. A supplied fixed mask represents equality constraints p_i=0.
ActiveSetStep SolveBoundedDampedLeastSquares(const Eigen::MatrixXd &,const Eigen::VectorXd &,
    const Eigen::VectorXd &,const Eigen::VectorXd &,double,const std::vector<bool> & = {});

// Selects the LM damping so the same constrained subproblem also satisfies
// ||D p|| <= trust_radius. The returned step is the step used by the model.
ActiveSetStep SolveBoundedTrustRegionStep(const Eigen::MatrixXd &,const Eigen::VectorXd &,
    const Eigen::VectorXd &,const Eigen::VectorXd &,double,double,const std::vector<bool> & = {});

} // namespace rhbm_gem::core::joint_component
