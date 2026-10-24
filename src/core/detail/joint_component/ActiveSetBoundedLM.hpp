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

// Selects the LM damping so the same constrained subproblem also satisfies
// ||D p|| <= trust_radius. The returned step is the step used by the model.
ActiveSetStep SolveBoundedTrustRegionStep(const Eigen::MatrixXd &,const Eigen::VectorXd &,
    const Eigen::VectorXd &,const Eigen::VectorXd &,double,double,const std::vector<bool> & = {});

} // namespace rhbm_gem::core::joint_component
