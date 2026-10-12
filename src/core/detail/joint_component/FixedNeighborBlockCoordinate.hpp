#pragma once
#include "Numerics.hpp"
#include "FixedNeighborPolicy.hpp"
#include "ResourceWork.hpp"
#include <functional>

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
};
struct FixedNeighborBlockRecord
{
    std::size_t sweep{},block{},affected_rows{};
    Indices core_atoms;
    double objective_before{},objective_after{},local_objective_before{},local_objective_after{};
    double global_replay_delta{},local_global_delta_error{},objective_replay_enclosure{};
    double local_final_gradient_inf_norm{},local_final_ac_kkt{};
    double search_seconds{},block_seconds{},profile_factor_seconds{};
    int profile_evaluations{},accepted_updates{};
    std::vector<FixedNeighborProfileTrial> profile_trials;
    bool accepted{};
    std::string status,reason,local_search_stop_reason;
};
struct FixedNeighborBlockSweep
{
    std::size_t sweep{};
    double objective_before{},objective_after{},global_a_feasibility{},global_ac_kkt{};
    double global_width_gradient_inf_norm{},cache_replay_error{},objective_replay_error{},wall_seconds{};
    double eta_change_inf{},beta_scaled_change{};
    bool coordinate_confirmation_available{};
    std::size_t block_solves{},profile_evaluations{},maximum_block_rows{},maximum_block_columns{},
        accepted_blocks{},unchanged_blocks{},accepted_local_updates{};
};
struct IncrementalCandidateReplay
{
    Vector residual_delta;
    double objective{},objective_delta{},delta_condition{};
};
IncrementalCandidateReplay EvaluateCandidateDelta(VectorRef current_residual,
    const Indices & affected_rows,VectorRef old_core,VectorRef new_core,
    double current_objective,double scale);
struct FixedNeighborWork
{
    double structural_preparation_seconds{},endpoint_verification_seconds{},primary_endpoint_seconds{},
        independent_replay_seconds{},
        rank_assessment_seconds{},old_core_seconds{},
        effective_response_seconds{},local_search_seconds{},local_state_seconds{},
        candidate_copy_seconds{},candidate_replay_seconds{},cache_update_seconds{},sweep_replay_seconds{},
        sweep_global_state_seconds{},maximum_block_seconds{},incremental_candidate_replay_seconds{},
        maximum_incremental_objective_difference{},maximum_incremental_relative_objective_difference{},
        maximum_incremental_residual_difference{},maximum_incremental_delta_condition{},maximum_cache_replay_error{},
        last_incremental_candidate_objective{},last_full_candidate_objective{};
    std::size_t block_visits{},block_solves{},full_candidate_replays{},incremental_candidate_replays{},
        incremental_acceptance_mismatches{},candidate_state_full_copies{},affected_row_updates{},
        old_core_basis_builds{};
    std::size_t local_problem_count{},minimum_local_atoms{},maximum_local_atoms{};
    double total_local_atoms{};
    ProfileSearchWork local_profile_work;
};
struct FixedNeighborEndpointCertification
{
    Assessment assessment;
    TrustEvidence trust;
    MemoryEstimate verification_memory_estimate;
    bool passed() const {return trust.passed;}
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
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
    bool stop_after_stationarity{true};
    bool capture_local_trajectory{};
    bool retain_diagnostic_records{true};
    bool compare_incremental_replay{};
    bool certify_endpoint{true}; // Diagnostic facade only; production components always certify.
    bool collect_telemetry{};
    bool collect_diagnostics{};
    std::function<void()> preparation_observer;
    FixedNeighborDiagnosticSink * diagnostic_sink{};
    std::function<void(const FixedNeighborBlockSweep &)> production_progress;
    std::function<void(const FixedNeighborBlockSweep &)> sweep_observer;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    // Test-only fault injection; compiled out of production targets.
    std::function<FixedNeighborEndpointCertification(const Domain &,VectorRef,const BlockCoordinateState &,
        const EvaluationContext &,ProfileSearchWork *)> certification_for_testing;
#endif
    std::function<void(std::size_t,const BlockCoordinateState &,const FixedNeighborBlockSweep &,
        const std::vector<FixedNeighborBlockRecord> &)> state_observer;
};
// Search output contains only the candidate state and search diagnostics.
struct FixedNeighborSearchResult
{
    BlockCoordinateState state;
    std::optional<FixedNeighborEndpointCertification> certification;
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{},certification_sweep{};
    std::size_t sweep_count{},total_block_solves{},total_profile_evaluations{},total_accepted_local_updates{},
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

bool IsFixedNeighborEtaChangeConfirmed(double eta_change_inf,bool has_previous_complete_sweep);
FixedNeighborStationarityState ClassifyFixedNeighborStationarity(
    double ac_kkt,double width_gradient_inf,double eta_change_inf,bool has_previous_complete_sweep);
bool ShouldAttemptFixedNeighborCertification(FixedNeighborStationarityState,std::size_t sweep,
    std::size_t & next_retry_sweep);
FixedNeighborSearchResult SearchFixedNeighborComponent(
    const Domain &,const JointParameterLayout &,VectorRef y,VectorRef initial_eta,const EvaluationContext &,
    const FixedNeighborPolicy & = {},
    Vector initial_beta = {});
FixedNeighborEndpointCertification CertifyFixedNeighborEndpoint(
    const Domain &,VectorRef y,const BlockCoordinateState &,const EvaluationContext &,
    ProfileSearchWork * = nullptr,std::size_t resident_override=0);
FixedNeighborComponentResult SolveFixedNeighborComponent(
    const ComponentProblemView &,VectorRef initial_b,const EvaluationContext &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr);
std::optional<double> LastFixedNeighborObjective(const FixedNeighborSearchResult &);
std::optional<double> LastFixedNeighborGradient(const FixedNeighborSearchResult &);
}
