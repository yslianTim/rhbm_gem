#pragma once
#include <Eigen/Dense>
#include "JointProgress.hpp"
#include "SolverRoute.hpp"
#include "SnapshotViews.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <Eigen/SparseCore>
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
    SearchPolicy search;
    int profile_budget{200},update_budget{100};
};
enum class ProfileEvaluationRole {Unspecified,InitialProfile,TrialProfile,AcceptedEndpoint,Reference};
const char * ProfileEvaluationRoleName(ProfileEvaluationRole);
struct ProfileRoleWork
{
    std::size_t evaluations{},derivative_preparations{},derivative_reductions{},replay_checks{};
    double evaluation_seconds{},basis_seconds{},linear_matrix_preparation_seconds{},linear_symbolic_seconds{},
        linear_numeric_seconds{},linear_rhs_solve_seconds{},linear_certificate_seconds{},derivative_prepare_seconds{},
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
    ProfileRoleWork total,initial_profile,trial_profile,accepted_endpoint,reference_evaluation;
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
class LinearWorkspace;
class FreeDesignFactor;
struct LinearResult
{
    Vector beta;
    std::shared_ptr<FreeDesignFactor> factor;
    bool valid{};
    std::string reason;
    int rank{},solves{},releases{},block_factorizations{};
};
struct LinearBlock {std::vector<Eigen::Index> rows,columns;};
LinearResult SolveLinear(const Sparse &,VectorRef,const Vector &,bool=false,bool=true,
    const Sparse * = nullptr,const LinearPolicy * = nullptr,const std::vector<LinearBlock> * = nullptr,LinearWorkspace * = nullptr);
LinearResult SolveLinear(const Matrix &,VectorRef,const Vector &,bool=false,bool=false,const Sparse * = nullptr);
std::pair<Matrix,Vector> ReferenceQR(const Sparse &,const Vector &,const Vector &,VectorRef);
struct Certificate
{
    bool evaluated{},available{},feasible{},kkt_passed{};
    double projected_kkt{unavailable},rss{unavailable},objective{unavailable},residual_scale{unavailable},
        residual_rmse{unavailable},residual_max{unavailable},relative_residual{unavailable};
    std::vector<Eigen::Index> active_atoms;
    std::optional<int> linear_solves,free_rank,block_factorizations;
};
Certificate CertifyLinear(const Sparse &,VectorRef,const Vector &,double=0);
struct Endpoint
{
    Vector eta,beta,gradient;
    Certificate certificate;
    bool valid{};
    std::string reason;
};
struct Evaluation : Endpoint {Vector residual; Sparse x,derivative; std::shared_ptr<FreeDesignFactor> factor;};
std::optional<double> NormalizedProfileObjective(const Evaluation &,double scale);
std::optional<double> ProfileGradientInfinityNorm(const Evaluation &);
void UpdateAcceptedProfileObjective(std::optional<double> &,const Evaluation &,double scale,bool accepted);
struct Spectrum
{
    bool available{true};
    std::string reason;
    Eigen::Index rank{},rows{},columns{};
    Vector singular_values,column_norms;
    double minimum{unavailable},condition{unavailable},threshold{};
};
Spectrum DesignSpectrum(const Sparse &,const Vector &,const RankPolicy * = nullptr);
Spectrum ComputeSpectrum(const Sparse &,const RankPolicy &,Eigen::Index,bool);
Spectrum ComputeSpectrum(const Matrix &,const RankPolicy &,Eigen::Index,bool);
Evaluation EvaluateProfile(const Domain &,VectorRef,const Vector &,bool,const EvaluationContext *,
    const std::vector<LinearBlock> * = nullptr,LinearWorkspace * = nullptr,const void * = nullptr,
    ProfileEvaluationRole = ProfileEvaluationRole::Unspecified,ProfileSearchWork * = nullptr);
Evaluation EvaluateState(const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    ProfileEvaluationRole = ProfileEvaluationRole::Unspecified,ProfileSearchWork * = nullptr);
struct TrustEvidence
{
    std::optional<Endpoint> reference;
    bool primary_valid{},passed{},prediction_passed{},gradient_passed{};
    std::string reason;
    double coefficient_difference{unavailable},kkt_difference{unavailable},prediction_difference{unavailable},
        gradient_difference{unavailable},cancellation_ratio{unavailable};
    std::optional<Spectrum> design;
};
TrustEvidence CheckReplay(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &);
TrustEvidence CheckTrust(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &);
TrustEvidence CheckTrust(const Domain &,VectorRef,const Evaluation &,const EvaluationContext &,const Evaluation & reference);
struct LmTrial
{
    Vector accepted_eta,step,diagonal;
    double radius{},damping{},actual_decrease{unavailable},predicted_decrease{},ratio{unavailable};
    bool proposed_acceptance{};
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
};
struct SearchResult
{
    Endpoint initial;
    Vector eta;
    std::vector<Trial> trials;
    std::optional<double> accepted_objective;
    std::optional<double> accepted_gradient_inf_norm;
    int lm_status{},evaluations{},derivatives{},accepted{},references{};
    bool stopped{},initial_accepted{};
    std::string stop_reason;
    double seconds{},reference_seconds{};
    ProfileSearchWork profile_work;
};
SearchResult SearchProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext &,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr,LinearWorkspace * = nullptr,
    const void * = nullptr,ProfileSearchWork * = nullptr);
struct Assessment
{
    Endpoint primary,reference;
    std::optional<Spectrum> design,widths,normalized_widths,jacobian;
    Matrix weak_directions;
    Vector correction;
    double coefficient_difference{unavailable};
    bool inner{},gradient{},local{},identified{};
    std::string failure;
};
std::vector<JointCheck> AssessmentEvidence(const Assessment &,JointEvidenceScope);
JointCheckStatus MergeConvergenceStatus(JointCheckStatus,JointCheckStatus);
JointCheckStatus ConvergenceStatus(const std::vector<JointCheck> &,JointEvidenceScope,bool assembled=false);
Assessment AssessProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext &,const Vector * = nullptr);
// Evaluations belong to this exact domain, observations and numerical policy.
Assessment AssessEvaluated(const Domain &,VectorRef,const Evaluation &,const Evaluation &,const EvaluationContext &,bool supplied=false);
bool SameAssessmentPolicy(const EvaluationContext &,const EvaluationContext &);
// Borrowed only while the originating immutable problem and report remain alive.
struct AssessmentReuse
{
    const Domain & domain;
    VectorRef observations;
    const EvaluationContext & context;
    const Assessment & assessment;
};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct AssessmentStageWork
{
    std::string name;
    int calls{},completed_calls{};
    Eigen::Index rows{},columns{};
    double seconds{};
};
struct AssessmentWork
{
    int assessments{},reference_evaluations{};
    int compact_attempts{},compact_accepted{},compact_boundary_fallbacks{},compact_other_fallbacks{};
    std::string active_stage,last_stage;
    std::vector<std::string> stage_stack;
    std::vector<AssessmentStageWork> stages;
};
AssessmentWork & AssessmentWorkForTesting();
using AssessmentStageObserverForTesting=void (*)(const AssessmentWork &,void *);
void SetAssessmentStageObserverForTesting(AssessmentStageObserverForTesting,void *);
void BeginAssessmentStageForTesting(const std::string &,Eigen::Index,Eigen::Index);
void FinishAssessmentStageForTesting(const std::string &,double,bool);
class AssessmentStageTimerForTesting
{
public:
    AssessmentStageTimerForTesting(std::string,Eigen::Index,Eigen::Index);
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
struct ComponentResult
{
    SearchResult search;
    Assessment assessment; // Historical search-endpoint diagnostics.
    std::optional<Assessment> trusted_assessment; // Actual returned state.
    std::optional<Endpoint> trusted_state;
    std::optional<std::size_t> trusted_trial;
    double assessment_seconds{};
    std::optional<TrustEvidence> endpoint_trust;
    bool search_success{};
};
ComponentResult AssessComponentSearch(const Domain &,VectorRef,const EvaluationContext &,SearchResult,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr);
ComponentResult SolveComponent(const ComponentView &,VectorRef,const Vector &,const EvaluationContext &,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr);
struct AssemblyResult
{
    bool available{},completed{true},profile_agrees{};
    std::vector<bool> row_mask;
    Vector eta,beta,prediction;
    double objective{unavailable},profile_difference{unavailable};
    Assessment assessment;
    Endpoint raw,profile_control;
    bool profile_evaluated{};
};
AssemblyResult AssembleComponents(const Domain &,VectorRef,const ComponentPartition &,const EvaluationContext &,
    const std::vector<ComponentResult> &,const AssessmentReuse * = nullptr);
}
