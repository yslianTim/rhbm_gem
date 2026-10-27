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
    const Eigen::VectorXd &,const Eigen::VectorXd &,double,double);

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct BoundedLmDampingTrace
{
    double damping{};
    int phase{};
    std::string status;
    Eigen::VectorXd step,gradient,multipliers;
    std::vector<bool> active;
    std::vector<std::vector<bool>> active_set_sequence;
    std::vector<Eigen::VectorXd> step_sequence,gradient_sequence;
    std::vector<double> kkt_tolerance_sequence;
    std::vector<Eigen::Index> free_ranks;
    std::size_t iterations{},releases{};
    double scaled_step_norm{},kkt_inf_norm{};
    bool valid{},within_trust_region{},search_accepted{};
};
struct BoundedLmTestWork
{
    std::vector<BoundedLmDampingTrace> damping_trace;
    bool capture_trace{};
    std::size_t damping_solves{},bracket_solves{},bisection_solves{},active_set_iterations{},releases{};
    std::size_t dense_qr_factorizations{},workspace_allocation_sites{},workspace_peak_bytes{};
    double dense_qr_factorization_seconds{},damping_search_seconds{};
};
BoundedLmTestWork & BoundedLmWorkForTesting();
#endif

} // namespace rhbm_gem::core::joint_component
