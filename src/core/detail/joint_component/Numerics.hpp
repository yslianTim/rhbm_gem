#pragma once
#include <Eigen/Dense>
#include "JointProgress.hpp"
#include "SnapshotViews.hpp"
#include "ResourceWork.hpp"
#include "ConvergenceCriteria.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <Eigen/SparseCore>
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <limits>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <chrono>
#endif

namespace rhbm_gem::core::joint_component {
using Matrix=Eigen::MatrixXd;
using Sparse=Eigen::SparseMatrix<double>;
inline constexpr double unavailable=std::numeric_limits<double>::quiet_NaN();
struct RankPolicy
{
    Eigen::Index rows{}, design_columns{}, width_columns{};
    double Relative(Eigen::Index columns) const;
    double Absolute(Eigen::Index columns,double maximum) const;
};
struct AuditPlan
{
    bool trial_details{}, expanded_if_unverified{}, precision{}, boundary{};
    bool block_precision{}, cache_precision{};
    std::vector<Eigen::Index> boundary_atoms;
    Matrix directions;
};
struct LinearPolicy
{
    double rank_relative{};
    int active_set_iteration_factor{20};
    double release_factor{128}, release_response_norm{-1};
};
struct EvaluationContext
{
    std::string snapshot_hash;
    std::shared_ptr<const VectorMap> observations;
    Identities atom_ids,row_ids;
    double scale{1};
    RankPolicy rank;
    LinearPolicy linear;
    AuditPlan audit;
    bool independent_search{};
    int profile_budget{200},update_budget{100};
    double profile_cpu_budget_seconds{};
    std::optional<double> eta_upper_bound;
};
enum class ProfileEvaluationRole {Unspecified,InitialProfile,TrialProfile,AcceptedEndpoint};
struct BoundedLmWork
{
    std::size_t subproblem_solves{},failed_subproblems{},active_set_iterations{},releases{},active_bound_occurrences{};
    double maximum_subproblem_kkt{unavailable};
    void Merge(const BoundedLmWork &);
};
struct LinearSolveWork
{
    std::size_t solve_calls{},valid_solves{},invalid_solves{},invalid_input_solves{},rank_deficient_solves{},
        nonfinite_solves{},active_set_budget_exhausted_solves{},other_invalid_solves{};
    std::size_t factorization_solves{},active_set_iterations{},active_set_releases{};
    std::size_t initial_kkt_count{},initial_kkt_le_1e10{},initial_kkt_le_1e8{},initial_kkt_le_1e7{};
    std::size_t final_kkt_count{},final_kkt_le_1e10{},final_kkt_le_1e8{},final_kkt_le_1e7{};
    std::size_t refinement_eligibility_count{},refinement_trigger_count{},refinement_iterations{},
        refinement_success_count{},refinement_failure_count{},refinement_factorizations{};
    double initial_kkt_sum{},initial_kkt_max{},final_kkt_sum{},final_kkt_max{};
    double refinement_seconds{},refinement_before_kkt_sum{},refinement_before_kkt_max{},
        refinement_after_kkt_sum{},refinement_after_kkt_max{},refinement_before_objective_sum{},
        refinement_after_objective_sum{};
};
struct ProfileRoleWork
{
    std::size_t evaluations{},derivative_preparations{},derivative_reductions{},replay_checks{};
    std::size_t local_kkt_rejections{};
    LinearSolveWork linear_solve;
    double evaluation_seconds{},basis_seconds{},linear_matrix_preparation_seconds{},linear_rhs_solve_seconds{},
        linear_certificate_seconds{},derivative_prepare_seconds{},
        derivative_reduce_seconds{},derivative_raw_assembly_seconds{},derivative_free_design_assembly_seconds{},
        derivative_factor_match_seconds{},derivative_factor_build_seconds{},derivative_factor_compact_seconds{},
        derivative_rank_seconds{},derivative_least_squares_seconds{},derivative_normal_solve_seconds{},
        derivative_cancellation_check_seconds{},derivative_cancellation_fallback_seconds{},derivative_rows_seconds{},
        derivative_jacobian_qr_seconds{},derivative_norms_seconds{},derivative_outer_overhead_seconds{},
        tiled_qr_assembly_copy_seconds{},tiled_qr_householder_seconds{},tiled_qr_rhs_transform_seconds{},
        replay_trust_seconds{};
};
struct ProfileSearchWork
{
    ProfileRoleWork total,initial_profile,trial_profile,accepted_endpoint;
    BoundedLmWork bounded_lm;
    std::size_t accepted_profile_trials{},rejected_profile_trials{},local_replay_kkt_rejections{},
        local_replay_disagreements{},invalid_gradient_stops{},width_safeguard_rejections{};
    std::size_t accepted_update_budget_stops{},profile_budget_stops{},cpu_budget_stops{},native_lm_stops{};
    std::array<std::size_t,16> lm_status_counts{};
    double total_seconds{},lm_overhead_seconds{};
    ProfileRoleWork & Role(ProfileEvaluationRole);
    void Add(ProfileEvaluationRole,const ProfileRoleWork &);
    void Merge(const ProfileSearchWork &);
    double AttributedSeconds() const;
};
EvaluationContext CreateContext(VectorRef,Eigen::Index,const std::string & = "",const AuditPlan & = {});
EvaluationContext CreateContext(std::shared_ptr<const JointProblemInput>,const std::string & = "",const AuditPlan & = {});
struct BasisValues {double gaussian{},charge{},gaussian_log_width{},charge_log_width{};};
BasisValues EvaluateKernel(double,double,double);
struct LinearResult
{
    Vector beta;
    bool valid{};
    std::string reason;
    int rank{},solves{},releases{},block_factorizations{};
    std::size_t active_set_iterations{},refinement_eligibility_count{},refinement_trigger_count{},
        refinement_iterations{},refinement_success_count{},refinement_failure_count{},refinement_factorizations{};
    double initial_projected_kkt{unavailable},final_projected_kkt{unavailable},refinement_seconds{},
        refinement_before_objective{unavailable},refinement_after_objective{unavailable};
};
struct LinearBlock {std::vector<Eigen::Index> rows,columns;};
LinearResult SolveLinear(const Sparse &,VectorRef,const Vector &,const LinearPolicy * = nullptr,
    const std::vector<LinearBlock> * = nullptr,double=0);
struct Certificate
{
    bool evaluated{},available{},feasible{},kkt_passed{};
    bool linear_solver_valid{};
    std::string linear_solver_reason;
    double projected_kkt{unavailable},rss{unavailable},objective{unavailable},residual_scale{unavailable},
        residual_rmse{unavailable},residual_max{unavailable},relative_residual{unavailable};
    std::vector<Eigen::Index> active_atoms;
    std::optional<int> linear_solves,linear_releases,free_rank,block_factorizations;
    bool width_bound_applied{},width_bound_feasible{true};
    double width_bound_projected_kkt{unavailable},maximum_width_bound_violation{},
        boundary_gradient_max{},boundary_multiplier_max{},width_bound_complementarity{unavailable};
    std::size_t active_width_bound_count{};
};
Certificate CertifyLinear(const Sparse &,VectorRef,const Vector &,double=0);
struct Endpoint
{
    Vector eta,beta,gradient;
    Certificate certificate;
    bool valid{};
    std::string reason;
};
struct Evaluation : Endpoint {Vector residual; Sparse x,derivative;};
std::optional<double> NormalizedProfileObjective(const Evaluation &,double scale);
std::optional<double> ProfileGradientInfinityNorm(const Evaluation &);
std::optional<double> ProfileGradientInfinityNorm(const Evaluation &,const EvaluationContext &);
void UpdateAcceptedProfileObjective(std::optional<double> &,const Evaluation &,double scale,bool accepted);
struct Spectrum
{
    bool available{true};
    bool threshold_sensitive{};
    std::string reason;
    Eigen::Index rank{},rows{},columns{},full_design_columns{};
    Vector singular_values,column_norms;
    double minimum{unavailable},condition{unavailable},threshold{};
};
struct WidthBoundKkt
{
    bool feasible{};
    double projected{std::numeric_limits<double>::infinity()},violation{},
        boundary_gradient{unavailable},multiplier{},complementarity{};
    std::size_t active{};
};
WidthBoundKkt CalculateWidthBoundKkt(const Vector &,const Vector &,double);
Spectrum ComputeSpectrum(const Sparse &,const RankPolicy &,Eigen::Index,bool);
Spectrum ComputeSpectrum(const Matrix &,const RankPolicy &,Eigen::Index,bool);
Evaluation EvaluateProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext *,
    const std::vector<LinearBlock> * = nullptr,
    ProfileEvaluationRole = ProfileEvaluationRole::Unspecified,ProfileSearchWork * = nullptr);
Evaluation EvaluateState(const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    ProfileEvaluationRole = ProfileEvaluationRole::Unspecified,ProfileSearchWork * = nullptr);
// Scalar/vector reductions for accepted global sweep states. Leaves x and derivative empty.
Evaluation EvaluateStateStreaming(const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    ProfileEvaluationRole = ProfileEvaluationRole::Unspecified,ProfileSearchWork * = nullptr);
MemoryEstimate EndpointVerificationMemoryEstimate(const Domain &,VectorRef,const EvaluationContext &,std::size_t resident_override=0);
struct IndependentReplayEvidence
{
    bool available{},valid{},feasible{},primary_valid{},primary_kkt_passed{},independent_kkt_passed{},
        prediction_passed{},residual_passed{},objective_passed{},kkt_agrees{},gradient_agrees{},passed{},
        assembled_prediction_available{},assembled_prediction_passed{},assembled_objective_passed{};
    std::string reason;
    double prediction_difference{unavailable},residual_difference{unavailable},objective_difference{unavailable},
        independent_objective{unavailable},
        kkt_difference{unavailable},gradient_difference{unavailable},cancellation_ratio{unavailable},
        assembled_prediction_difference{unavailable},assembled_prediction_error_ratio{unavailable},
        assembled_objective_difference{unavailable},
        primary_objective{unavailable},primary_residual_norm{unavailable},independent_residual_norm{unavailable},
        primary_ac_kkt{unavailable},independent_ac_kkt{unavailable},primary_width_gradient{unavailable},
        independent_width_gradient{unavailable};
    std::string maximum_gradient_difference_atom_id,gradient_enclosure_failure_atom_id;
    double gradient_enclosure_failure_primary{unavailable},gradient_enclosure_failure_independent{unavailable},
        gradient_enclosure_failure_difference{unavailable},gradient_enclosure_failure_bound{unavailable},
        gradient_enclosure_failure_ratio{unavailable},maximum_gradient_difference_primary{unavailable},
        maximum_gradient_difference_independent{unavailable},maximum_gradient_difference_bound{unavailable},
        maximum_gradient_difference_ratio{unavailable};
    bool width_bound_applied{},width_bound_feasible{true},primary_width_kkt_passed{true},
        independent_width_kkt_passed{true};
    double primary_width_projected_kkt{unavailable},independent_width_projected_kkt{unavailable},
        maximum_width_bound_violation{},boundary_gradient_max{},boundary_multiplier_max{},
        complementarity_residual{unavailable};
    std::size_t active_width_bounds{};
};
IndependentReplayEvidence CheckIndependentReplay(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &,
    const Vector * assembled_prediction=nullptr,const std::vector<std::size_t> * assembled_rows=nullptr,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
struct TrustEvidence
{
    std::optional<IndependentReplayEvidence> independent;
    bool primary_valid{},passed{},prediction_passed{},gradient_passed{};
    std::string reason;
    double kkt_difference{unavailable},prediction_difference{unavailable},
        gradient_difference{unavailable},cancellation_ratio{unavailable};
    bool width_bound_feasible{true};
};
struct ReplayDetailsForDiagnostics
{
    Vector prediction,residual,width_gradient;
    double objective{unavailable},scaled_kkt{unavailable};
};
TrustEvidence CheckReplay(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &,
    ReplayDetailsForDiagnostics * = nullptr);
TrustEvidence CheckTrust(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &);
TrustEvidence CheckTrust(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &,const IndependentReplayEvidence &);
struct LmTrial
{
    Vector accepted_eta,step,diagonal;
    double radius{},damping{},actual_decrease{unavailable},predicted_decrease{},ratio{unavailable};
    bool proposed_acceptance{};
    bool bounded_subproblem{};
    std::string subproblem_status;
    double subproblem_kkt_inf_norm{unavailable};
    std::size_t active_set_iterations{},active_set_releases{};
    std::vector<std::size_t> active_bound_indices;
    std::vector<double> active_bound_multipliers;
};
struct Trial
{
    Endpoint endpoint;
    int evaluation{};
    bool accepted{};
    std::optional<int> accepted_update;
    double seconds{},factor_seconds{};
    std::optional<LmTrial> lm;
    std::optional<TrustEvidence> trust;
    double weak_width_step_fraction{unavailable};
    std::string rejection_reason;
};
struct WeakWidthStepAssessment
{
    double weak_step_fraction{unavailable};
    bool available{};
};
WeakWidthStepAssessment AssessWeakWidthStep(const Matrix &,VectorRef,const RankPolicy &);
bool RejectWeakWidthTrial(double step_inf,double weak_step_fraction,double objective_before,
    double objective_after);
struct SearchResult
{
    Endpoint initial;
    Vector eta;
    std::vector<Trial> trials;
    std::optional<double> accepted_objective;
    std::optional<double> accepted_gradient_inf_norm;
    int lm_status{},evaluations{},derivatives{},accepted{};
    bool stopped{},initial_accepted{};
    std::string stop_reason;
    double seconds{};
    ProfileSearchWork profile_work;
};
SearchResult SearchProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext &,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr,ProfileSearchWork * = nullptr,
    const std::vector<bool> * fixed_width_mask=nullptr);
struct Assessment
{
    Endpoint primary;
    IndependentReplayEvidence independent_replay;
    MemoryEstimate verification_memory_estimate;
    JointConvergenceCriteria convergence_criteria{StrictJointConvergenceCriteria};
    double primary_endpoint_seconds{},independent_replay_seconds{};
    bool inner{},gradient{},endpoint_verified{},numerical_qualified{};
    std::string failure;
};
std::vector<JointCheck> AssessmentEvidence(const Assessment &,JointEvidenceScope);
JointCheckStatus MergeConvergenceStatus(JointCheckStatus,JointCheckStatus);
JointCheckStatus ConvergenceStatus(const std::vector<JointCheck> &,JointEvidenceScope,bool assembled=false);
Assessment AssessProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext &,const Vector * = nullptr,
    const Vector * assembled_prediction=nullptr,const std::vector<std::size_t> * assembled_rows=nullptr);
// Evaluations belong to this exact domain, observations and numerical policy.
Assessment AssessEvaluated(const Domain &,VectorRef,const Evaluation &,const IndependentReplayEvidence &,
    const EvaluationContext &,bool supplied=false,std::size_t resident_override=0);
// Verifies an accepted endpoint using only the streaming objective, KKT,
// width-gradient, and independent replay evidence. It never assesses global rank.
Assessment VerifyJointEndpoint(const Domain &,VectorRef,const Evaluation &,
    const IndependentReplayEvidence &,const EvaluationContext &,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct AssessmentStageWork
{
    std::string name;
    int calls{},completed_calls{};
    Eigen::Index rows{},columns{};
    double seconds{};
    bool memory_estimate_available{};
    MemoryStage memory_stage{MemoryStage::ProblemConstruction};
    MemoryStatus memory_status{MemoryStatus::Pass};
    std::size_t estimated_resident_bytes{},estimated_temporary_bytes{},estimated_peak_bytes{},budget_bytes{};
    std::string memory_preflight_reason;
    std::size_t sparse_design_bytes{},dense_factor_bytes{},factor_workspace_bytes{};
    std::size_t factor_rows{},factor_columns{},maximum_tile_rows{};
    std::optional<std::size_t> current_rss_before_bytes,current_rss_after_bytes,process_peak_rss_bytes;
};
struct AssessmentWork
{
    int assessments{};
    std::string active_stage,last_stage;
    std::vector<std::string> stage_stack;
    std::vector<AssessmentStageWork> stages;
};
AssessmentWork & AssessmentWorkForTesting();
using AssessmentStageObserverForTesting=void (*)(const AssessmentWork &,void *);
void SetAssessmentStageObserverForTesting(AssessmentStageObserverForTesting,void *);
void BeginAssessmentStageForTesting(const std::string &,Eigen::Index,Eigen::Index,const MemoryEstimate * = nullptr);
void FinishAssessmentStageForTesting(const std::string &,double,bool);
class AssessmentStageTimerForTesting
{
public:
    AssessmentStageTimerForTesting(std::string,Eigen::Index,Eigen::Index,const MemoryEstimate * = nullptr);
    ~AssessmentStageTimerForTesting();
    AssessmentStageTimerForTesting(const AssessmentStageTimerForTesting &)=delete;
    AssessmentStageTimerForTesting & operator=(const AssessmentStageTimerForTesting &)=delete;
    void Finish();
private:
    std::string name_;
    std::chrono::steady_clock::time_point started_;
    int uncaught_exceptions_{};
    bool finished_{};
};
#endif
struct ComponentView
{
    std::string id;
    std::vector<Eigen::Index> atoms,rows;
    std::shared_ptr<const PartitionMappings> mappings;
    Eigen::Index component_index{};
    Eigen::Index LocalAtom(Eigen::Index a) const {return mappings->atom_component.at(static_cast<std::size_t>(a))==component_index ? mappings->atom_to_local.at(static_cast<std::size_t>(a)) : -1;}
    Eigen::Index LocalRow(Eigen::Index r) const {return mappings->row_component.at(static_cast<std::size_t>(r))==component_index ? mappings->row_to_local.at(static_cast<std::size_t>(r)) : -1;}
    Domain domain{0,{}};
};
struct ComponentPartition
{
    std::vector<ComponentView> components;
    std::shared_ptr<const PartitionMappings> mappings;
    std::vector<Eigen::Index> constant_rows,unobserved_atoms;
};
ComponentPartition Partition(const Domain &,const Identities &);
Vector SelectValues(VectorRef,const std::vector<Eigen::Index> &);
EvaluationContext ChildContext(const EvaluationContext &,const ComponentView &,bool);
struct FixedNeighborComponentResult;
struct AssemblyResult
{
    bool available{},completed{true},mapping_valid{},consistency_passed{};
    std::string mapping_reason;
    std::vector<bool> row_mask;
    Vector eta,beta,prediction;
    double objective{unavailable};
    Assessment assessment;
    Endpoint raw;
    IndependentReplayEvidence replay;
    MemoryEstimate verification_memory_estimate;
};
AssemblyResult AssembleComponents(const Domain &,VectorRef,const ComponentPartition &,const EvaluationContext &,
    const std::vector<JointComponentResult> &,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
}
