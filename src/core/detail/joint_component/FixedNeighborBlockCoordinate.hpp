#pragma once
#include "FixedBBlockCoordinate.hpp"
#include "FixedNeighborPolicy.hpp"
#include <functional>

namespace rhbm_gem::core::joint_component {
struct FixedNeighborProfileTrial
{
    std::size_t trial_index{};
    int profile_evaluation{};
    bool accepted{};
    std::optional<int> accepted_update;
    double local_objective_before{},local_objective_after{},objective_reduction{};
    double eta_change_inf{},gradient_inf_norm{};
    double profile_seconds{},factor_seconds{},cumulative_factor_seconds{};
};
struct FixedNeighborBlockRecord
{
    std::size_t sweep{},block{},affected_rows{};
    Indices core_atoms;
    double objective_before{},objective_after{},local_objective_before{},local_objective_after{};
    double global_replay_delta{},local_global_delta_error{},objective_replay_enclosure{};
    double local_final_gradient_inf_norm{},local_final_ac_kkt{};
    double local_profile_gradient_inf_norm{},local_reference_gradient_inf_norm{},local_correction_inf_norm{};
    double local_projected_width_minimum{},local_corrected_jacobian_minimum{},local_normalized_width_minimum{};
    Eigen::Index local_projected_width_rank{},local_corrected_jacobian_rank{},local_normalized_width_rank{};
    double local_assessment_seconds{};
    Eigen::Index local_assessment_rows{},local_assessment_columns{};
    double search_seconds{},profile_factor_seconds{};
    int profile_evaluations{},accepted_updates{};
    std::vector<FixedNeighborProfileTrial> profile_trials;
    bool accepted{},local_assessment_attempted{},local_assessment_passed{},local_inner_passed{},
        local_gradient_passed{},local_correction_passed{},local_identified{},local_trust_passed{};
    std::string local_assessment_failure,local_trust_reason;
    std::string status,reason,local_search_stop_reason;
};
struct FixedNeighborBlockSweep
{
    std::size_t sweep{};
    double objective_before{},objective_after{},global_a_feasibility{},global_ac_kkt{};
    double global_width_gradient_inf_norm{},cache_replay_error{},objective_replay_error{},wall_seconds{};
    double eta_change_inf{},beta_scaled_change{};
    bool coordinate_confirmation_available{};
    double local_assessment_seconds{};
    std::size_t block_solves{},profile_evaluations{},local_assessments{},certified_local_candidates{},
        maximum_block_rows{},maximum_block_columns{},maximum_local_assessment_rows{},maximum_local_assessment_columns{},
        accepted_blocks{},unchanged_blocks{},accepted_local_updates{};
};
struct FixedNeighborPolicy
{
    std::size_t core_atoms{128},maximum_sweeps{30};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
    FixedNeighborLocalWork local_work{FixedNeighborLocalWork::Full};
    bool stop_after_stationarity{true};
    bool certify_local_candidates{};
    bool capture_local_trajectory{};
    bool stop_after_no_certified_update{};
    bool assess_final_endpoint{true};
    std::function<void(const FixedNeighborBlockSweep &)> sweep_observer;
    std::function<void(std::size_t,const BlockCoordinateState &,const FixedNeighborBlockSweep &,
        const std::vector<FixedNeighborBlockRecord> &)> state_observer;
};
struct FixedNeighborResult
{
    BlockCoordinateState state;
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{},endpoint_certified{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{};
    Assessment assessment;
    TrustEvidence endpoint_trust;
    JointFitResult fit;
};

// Search output shared by the production component route and the historical
// whole-problem wrapper. It deliberately contains no assembled JointFitResult.
struct FixedNeighborSearchResult
{
    BlockCoordinateState state;
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{},endpoint_certified{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{};
    Assessment assessment;
    TrustEvidence endpoint_trust;
};

bool IsCertifiedLocalEndpoint(const Assessment &,const TrustEvidence &);
bool IsFixedNeighborEtaChangeConfirmed(double eta_change_inf,bool has_previous_complete_sweep);
FixedNeighborSearchResult SearchFixedNeighborComponent(
    const JointProblemInput &,const JointParameterLayout &,const Domain &,VectorRef observations,
    VectorRef y,VectorRef initial_eta,const EvaluationContext &,const FixedNeighborPolicy & = {});
ComponentResult SolveFixedNeighborComponent(
    const JointProblemInput &,const JointParameterLayout &,const ComponentView &,VectorRef observations,
    VectorRef initial_b,const EvaluationContext &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr);
FixedNeighborResult SearchFixedNeighbor(const JointProblem &,VectorRef initial_eta,
    const FixedNeighborPolicy & = {});
}
