#pragma once
#include <Eigen/Dense>
#include "SnapshotViews.hpp"
#include "ResourceWork.hpp"
#include "ConvergenceCriteria.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <Eigen/SparseCore>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <limits>

namespace rhbm_gem::core::joint_component {
using Matrix=Eigen::MatrixXd;
using Sparse=Eigen::SparseMatrix<double>;
inline constexpr double unavailable=std::numeric_limits<double>::quiet_NaN();
struct RankPolicy
{
    Eigen::Index rows{}, design_columns{}, width_columns{};
    double Relative(Eigen::Index columns) const;
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
    int profile_budget{200},update_budget{100};
    std::optional<double> eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    bool force_width_subproblem_failure_for_testing{};
#endif
};
struct BoundedLmWork
{
    std::size_t subproblem_solves{},failed_subproblems{},active_set_iterations{},releases{},active_bound_occurrences{};
    std::size_t damping_solves{},bracket_solves{},damping_qr_factorizations{};
    double maximum_subproblem_kkt{unavailable};
    void Merge(const BoundedLmWork &);
};
EvaluationContext CreateContext(std::shared_ptr<const JointProblemInput>,const std::string & = "");
struct BasisValues {double gaussian{},charge{},gaussian_log_width{},charge_log_width{};};
BasisValues EvaluateKernel(double,double,double);
struct LinearResult
{
    Vector beta;
    bool valid{};
    std::string reason;
    int rank{},solves{},releases{};
    std::size_t active_set_iterations{},refinement_eligibility_count{},refinement_trigger_count{},
        refinement_iterations{},refinement_success_count{},refinement_failure_count{},refinement_factorizations{},
        free_face_changes{},free_face_unchanged{},sparse_structure_build_count{},tile_mapping_build_count{};
};
struct LinearSolveWorkspace;
LinearResult SolveLinear(const Sparse &,VectorRef,const LinearPolicy * = nullptr,double=0,
    LinearSolveWorkspace * = nullptr);
struct Certificate
{
    bool evaluated{},available{},feasible{},kkt_passed{};
    bool linear_solver_valid{};
    std::string linear_solver_reason;
    double projected_kkt{unavailable},rss{unavailable},objective{unavailable},residual_scale{unavailable},
        residual_rmse{unavailable},residual_max{unavailable},relative_residual{unavailable};
    std::vector<Eigen::Index> active_atoms;
    std::optional<int> linear_solves,linear_releases,free_rank;
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
Evaluation EvaluateProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext *,
    LinearSolveWorkspace * = nullptr);
Evaluation EvaluateState(const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &);
// Scalar/vector reductions for accepted global sweep states. Leaves x and derivative empty.
Evaluation EvaluateStateStreaming(const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &);
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
TrustEvidence CheckReplay(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &);
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
    int lm_status{},evaluations{},derivatives{},accepted{};
    bool stopped{},initial_accepted{},width_failed{};
    std::string stop_reason;
    BoundedLmWork bounded_lm_work;
};
SearchResult SearchProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext &);
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
// Verifies an accepted endpoint using only the streaming objective, KKT,
// width-gradient, and independent replay evidence. It never assesses global rank.
Assessment VerifyJointEndpoint(const Domain &,VectorRef,const Evaluation &,
    const IndependentReplayEvidence &,const EvaluationContext &,
    const JointConvergenceCriteria & = StrictJointConvergenceCriteria);
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
