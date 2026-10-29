#pragma once
#include "Numerics.hpp"
#include "FixedNeighborPolicy.hpp"
#include "JointProgress.hpp"
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

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct FixedNeighborCandidateEvent
{
    std::size_t sweep{},block{};
    const BlockCoordinateState & accepted_before;
    const Indices & core_atoms,& affected_rows;
    const Vector & candidate_eta,& candidate_beta;
    const Vector & affected_prediction,& affected_residual;
    double candidate_objective{},local_objective_before{},local_objective_after{},global_replay_delta{},
        local_global_delta_error{},objective_replay_enclosure{},incremental_objective{};
    bool accepted{},unchanged{},candidate_replay_ran{};
};
enum class FixedNeighborCandidateFaultForTesting
{
    None,LocalObjectiveIncrease,ObjectiveConsistencyFailure,NaNObjective,InfiniteObjective,
    NaNParameter,InfiniteParameter,InfinitePrediction,InfiniteResidual
};
#endif

double BlockObjectiveReplayEnclosure(double reference);
bool WithinBlockObjectiveReplay(double error,double reference);

struct FixedNeighborBlockFailureDiagnostics
{
    std::size_t component_index{},sweep_index{},block_index{},core_atom_count{},affected_row_count{};
    std::vector<std::string> core_atom_ids;
    bool initial_profile_valid{},initial_profile_trusted{},solver_valid{},width_bound_active{};
    std::string component_id,initial_profile_reason,local_search_stop_reason,solver_reason,lm_subproblem_status,
        first_invalid_atom_id,first_invalid_parameter;
    std::optional<Eigen::Index> first_invalid_coefficient;
    std::optional<int> linear_rank,linear_solve_count,active_set_releases;
    std::optional<double> first_invalid_coefficient_value,width_upper_bound,width_step_inf,
        lm_actual_decrease,lm_predicted_decrease,lm_ratio,weak_width_step_fraction,
        primal_feasibility,scaled_projected_kkt,replay_prediction_discrepancy,
        replay_gradient_discrepancy,replay_kkt_discrepancy;
};
struct FixedNeighborBlockRecord
{
    std::size_t sweep{},block{},affected_rows{};
    std::size_t recovery_attempt{};
    double objective_before{},objective_after{},local_objective_before{},local_objective_after{};
    double global_replay_delta{},local_global_delta_error{},objective_replay_enclosure{};
    int accepted_updates{};
    std::optional<FixedNeighborBlockFailureDiagnostics> failure_diagnostics;
    bool accepted{},recovery_retry{},recovery_success{};
    std::string status,reason;
};
struct FixedNeighborBlockSweep
{
    std::size_t sweep{};
    double objective_before{},objective_after{},global_a_feasibility{},global_ac_kkt{};
    double global_width_gradient_inf_norm{},cache_replay_error{},objective_replay_error{},wall_seconds{};
    double eta_change_inf{},beta_scaled_change{},maximum_accepted_log_width_step{},maximum_target_width{};
    bool objective_cache_refreshed{};
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
        sweep_replay_audit_seconds{},sweep_global_state_seconds{},maximum_block_seconds{},maximum_cache_replay_error{};
    double objective_reconciliation_seconds{},maximum_objective_reconciliation_error{};
    std::size_t full_candidate_replays{},affected_candidate_replays{},candidate_replay_skipped_unchanged{},
        sweep_replays{},affected_row_updates{},affected_replay_peak_scratch_bytes{},objective_reconciliations{},
        incremental_objective_evaluations{};
    bool full_replay_memory_fallback{};
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    bool stop_after_stationarity{true};
    bool certify_endpoint{true}; // Diagnostic facade only; production components always certify.
    bool assess_local_conditioning{};
    bool assess_parameter_stability{};
    bool collect_telemetry{};
    bool collect_diagnostics{};
    std::optional<std::size_t> global_state_memory_budget_for_testing;
    std::function<void(const FixedNeighborCandidateEvent &)> candidate_observer;
    bool perturb_sweep_replay_objective_for_testing{};
    std::function<FixedNeighborCandidateFaultForTesting(std::size_t,std::size_t,
        const BlockCoordinateState &)> candidate_fault_for_testing;
#endif
    std::size_t component_index{};
    std::string component_id;
    std::function<bool(int)> claim_parameter_stability_block;
    std::function<void(const FixedNeighborBlockSweep &)> production_progress;
    std::optional<double> eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::function<bool(std::size_t,std::size_t)> force_no_trusted_profile_for_testing;
#endif
};
// Search output contains only the candidate state and search diagnostics.
struct FixedNeighborSearchResult
{
    BlockCoordinateState state;
    std::optional<FixedNeighborEndpointCertification> certification;
    std::vector<JointBlockConditioning> local_conditioning;
    std::optional<JointParameterStability> parameter_stability;
    std::size_t recent_sweep_count{};
    double maximum_eta_change{},maximum_scaled_ac_change{};
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{},certification_sweep{};
    std::size_t sweep_count{},total_block_solves{},total_profile_evaluations{},total_accepted_local_updates{},
        certification_attempts{},recovery_attempts{},recovery_successes{},
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
    const FixedNeighborPolicy & = {});
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
