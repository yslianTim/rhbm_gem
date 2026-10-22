#pragma once
#include "Numerics.hpp"
#include "FixedNeighborPolicy.hpp"
#include "ResourceWork.hpp"
#include <functional>
#include <utility>

namespace rhbm_gem::core::joint_component {
struct ComponentProblemView;
struct BlockCoordinateState
{
    Vector eta,beta,prediction,residual;
    double objective{};
};

double BlockObjectiveReplayEnclosure(double reference);
bool WithinBlockObjectiveReplay(double error,double reference);

struct FixedNeighborProfileTrial
{
    std::size_t trial_index{};
    int profile_evaluation{};
    bool accepted{};
    std::optional<int> accepted_update;
    double local_objective_before{},local_objective_after{},objective_reduction{};
    double eta_change_inf{},gradient_inf_norm{};
    double profile_seconds{},factor_seconds{},cumulative_factor_seconds{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    double lm_actual_decrease{unavailable},lm_predicted_decrease{unavailable},lm_ratio{unavailable};
    double lm_subproblem_kkt_inf_norm{unavailable};
    bool bounded_subproblem{};
    std::string bounded_subproblem_status;
    std::size_t active_set_iterations{},active_set_releases{};
    std::vector<std::size_t> active_bound_indices;
    std::vector<double> active_bound_multipliers;
    double accepted_local_objective_improvement{unavailable};
    double weak_width_step_fraction{unavailable};
    bool lm_proposed_acceptance{};
    bool weak_width_guard_rejected{};
    std::string acceptance_reason;
#endif
};
struct FixedNeighborBlockFailureDiagnostics
{
    std::size_t component_index{},sweep_index{},block_index{},core_atom_count{},affected_row_count{},
        local_design_nonzero_entries{};
    std::vector<std::string> core_atom_ids;
    std::vector<double> local_column_norms,local_design_singular_values,initial_widths,
        width_derivative_column_norms,projected_width_singular_values,projected_width_column_norms,proposed_eta_step;
    bool initial_profile_valid{},initial_profile_trusted{},solver_valid{};
    bool width_bound_active{};
    std::string component_id,initial_profile_reason,local_search_stop_reason,solver_reason,lm_subproblem_status,
        first_invalid_atom_id,first_invalid_parameter;
    Eigen::Index local_design_rows{},local_design_columns{};
    std::optional<Eigen::Index> local_design_rank,first_invalid_coefficient;
    std::optional<int> linear_rank,linear_solve_count,active_set_releases;
    std::vector<Eigen::Index> active_amplitude_indices;
    std::optional<double> local_response_norm,local_response_min,local_response_max,parent_observation_scale,
        parent_response_norm,parent_rank_relative,release_factor,rank_relative,active_set_iteration_factor,condition_indicator,
        local_design_rank_threshold,first_invalid_coefficient_value,width_upper_bound,width_step_inf,width_gradient_inf,
        lm_actual_decrease,lm_predicted_decrease,lm_ratio,weak_width_step_fraction,projected_width_rank_threshold,
        primal_feasibility,scaled_projected_kkt,maximum_free_a_gradient_raw,maximum_free_a_gradient_scaled,
        maximum_free_c_gradient_raw,maximum_free_c_gradient_scaled,
        maximum_active_a_complementarity_violation_raw,maximum_active_a_complementarity_violation_scaled,
        replay_prediction_discrepancy,replay_gradient_discrepancy,replay_kkt_discrepancy;
    std::string failure_capture_status,replay_trial_capture_status;
};
struct FixedNeighborBlockRecord
{
    std::size_t sweep{},block{},affected_rows{};
    std::size_t recovery_attempt{};
    Indices core_atoms;
    double objective_before{},objective_after{},local_objective_before{},local_objective_after{};
    double global_replay_delta{},local_global_delta_error{},objective_replay_enclosure{};
    double local_final_gradient_inf_norm{},local_final_ac_kkt{};
    double search_seconds{},block_seconds{},profile_factor_seconds{};
    int profile_evaluations{},accepted_updates{};
    std::vector<FixedNeighborProfileTrial> profile_trials;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::size_t component_index{};
    std::string component_id;
#endif
    std::optional<FixedNeighborBlockFailureDiagnostics> failure_diagnostics;
    bool accepted{},recovery_retry{},recovery_success{};
    std::string status,reason,local_search_stop_reason;
};
struct FixedNeighborBlockSweep
{
    std::size_t sweep{};
    double objective_before{},objective_after{},global_a_feasibility{},global_ac_kkt{};
    double global_width_gradient_inf_norm{},cache_replay_error{},objective_replay_error{},wall_seconds{};
    double eta_change_inf{},beta_scaled_change{},maximum_accepted_log_width_step{},maximum_target_width{};
    bool width_bound_applied{};
    double projected_width_gradient_inf_norm{unavailable},maximum_width_bound_violation{},
        boundary_gradient_max{},boundary_multiplier_max{},width_bound_complementarity{unavailable},
        stationarity_width_gradient_inf_norm{unavailable};
    std::size_t active_width_bound_count{};
    bool coordinate_confirmation_available{};
    std::size_t block_solves{},profile_evaluations{},maximum_block_rows{},maximum_block_columns{},
        accepted_blocks{},unchanged_blocks{},accepted_local_updates{},recovery_attempts{},recovery_successes{};
};
struct FixedNeighborWork
{
    double structural_preparation_seconds{},endpoint_verification_seconds{},local_conditioning_seconds{},parameter_stability_seconds{},primary_endpoint_seconds{},
        independent_replay_seconds{},old_core_seconds{},
        effective_response_seconds{},local_search_seconds{},local_state_seconds{},
        candidate_copy_seconds{},candidate_replay_seconds{},cache_update_seconds{},sweep_replay_seconds{},
        sweep_global_state_seconds{},maximum_block_seconds{},maximum_cache_replay_error{};
    std::size_t block_visits{},block_solves{},full_candidate_replays{},candidate_state_full_copies{},
        affected_row_updates{},old_core_basis_builds{};
    std::size_t local_problem_count{},minimum_local_atoms{},maximum_local_atoms{};
    std::size_t local_conditioning_peak_memory_estimate_bytes{};
    std::size_t parameter_stability_peak_memory_estimate_bytes{},perturbation_solves{},perturbation_evaluations{};
    double total_local_atoms{};
    ProfileSearchWork local_profile_work;
};
struct FixedNeighborEndpointCertification
{
    Assessment assessment;
    TrustEvidence trust;
    MemoryEstimate verification_memory_estimate;
    bool verified() const {return trust.primary_valid && trust.passed && assessment.endpoint_verified;}
    bool qualified() const {return verified() && assessment.numerical_qualified;}
    bool resource_limited() const
    {
        return verification_memory_estimate.status!=MemoryStatus::Pass ||
            trust.reason.starts_with("resource-limit:") || assessment.failure.starts_with("resource-limit:");
    }
    std::string resource_limit_reason() const
    {
        if(verification_memory_estimate.status!=MemoryStatus::Pass)
            return std::string("resource-limit:endpoint-verification:")+verification_memory_estimate.reason;
        if(trust.reason.starts_with("resource-limit:")) return trust.reason;
        if(assessment.failure.starts_with("resource-limit:")) return assessment.failure;
        return {};
    }
};
enum class FixedNeighborStationarityState {NotStationary,CandidateStationary,ConfirmedStationary};
struct FixedNeighborDiagnosticSink
{
    std::function<void(const FixedNeighborBlockRecord &)> block;
    std::function<void(const FixedNeighborBlockSweep &)> sweep;
};
struct PreparedFixedNeighborBlock
{
    std::size_t block{};
    Indices core_atoms,affected_rows,profile_atoms,profile_rows;
    std::shared_ptr<const CompactRowMapping> row_mapping;
    Domain domain{0,{}};
    EvaluationContext context;
    Eigen::Index local_atoms{};
};
struct FixedNeighborPolicy
{
    std::size_t core_atoms{12},maximum_sweeps{30};
    JointConvergenceCriteria convergence_criteria{StrictJointConvergenceCriteria};
    bool stop_after_stationarity{true};
    bool capture_local_trajectory{};
    bool retain_diagnostic_records{true};
    bool certify_endpoint{true}; // Diagnostic facade only; production components always certify.
    bool collect_telemetry{};
    bool collect_diagnostics{};
    bool assess_local_conditioning{};
    bool assess_parameter_stability{};
    std::size_t component_index{};
    std::string component_id;
    std::function<void()> preparation_observer;
    std::function<bool(int)> claim_parameter_stability_block;
    FixedNeighborDiagnosticSink * diagnostic_sink{};
    std::function<void(const FixedNeighborBlockSweep &)> production_progress;
    std::function<void(const FixedNeighborBlockSweep &)> sweep_observer;
    std::optional<double> eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::function<bool(std::size_t,std::size_t)> force_no_trusted_profile_for_testing;
#endif
    std::function<void(std::size_t,const BlockCoordinateState &,const FixedNeighborBlockSweep &,
        const std::vector<FixedNeighborBlockRecord> &)> state_observer;
};
// Search output contains only the candidate state and search diagnostics.
struct FixedNeighborSearchResult
{
    BlockCoordinateState state;
    std::optional<FixedNeighborEndpointCertification> certification;
    std::vector<JointBlockConditioning> local_conditioning;
    std::optional<JointParameterStability> parameter_stability;
    std::vector<std::pair<double,double>> recent_sweep_changes;
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{},certification_sweep{};
    std::size_t sweep_count{},total_block_solves{},total_profile_evaluations{},total_accepted_local_updates{},
        recovery_attempts{},recovery_successes{},
        block_preparations{},domain_preparations{},mapping_preparations{};
    std::size_t prepared_mapping_entries{},prepared_mapping_bytes{};
    MemoryEstimate memory_estimate,global_state_memory_estimate;
    double search_seconds{};
    FixedNeighborWork work;
    std::optional<FixedNeighborBlockSweep> final_sweep;
};

struct FixedNeighborComponentResult
{
    FixedNeighborSearchResult search;
    bool succeeded{};
    std::string stop_reason;
};

bool IsFixedNeighborEtaChangeConfirmed(double eta_change_inf,bool has_previous_complete_sweep,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
FixedNeighborStationarityState ClassifyFixedNeighborStationarity(
    double ac_kkt,double width_gradient_inf,double eta_change_inf,bool has_previous_complete_sweep,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
bool ShouldAttemptFixedNeighborCertification(FixedNeighborStationarityState,std::size_t sweep,
    std::size_t & next_retry_sweep);
FixedNeighborSearchResult SearchFixedNeighborComponent(
    const Domain &,const JointParameterLayout &,VectorRef y,VectorRef initial_eta,const EvaluationContext &,
    const FixedNeighborPolicy & = {},
    Vector initial_beta = {});
FixedNeighborEndpointCertification CertifyFixedNeighborEndpoint(
    const Domain &,VectorRef y,const BlockCoordinateState &,const EvaluationContext &,
    ProfileSearchWork * = nullptr,std::size_t resident_override=0,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
FixedNeighborComponentResult SolveFixedNeighborComponent(
    const ComponentProblemView &,VectorRef initial_b,const EvaluationContext &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr);
std::optional<double> LastFixedNeighborObjective(const FixedNeighborSearchResult &);
std::optional<double> LastFixedNeighborGradient(const FixedNeighborSearchResult &);
}
