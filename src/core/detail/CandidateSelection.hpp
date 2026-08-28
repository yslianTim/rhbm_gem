#pragma once

#include "core/detail/FittingModel.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointFitting.hpp"

#include <atomic>
#include <chrono>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace rhbm_gem::core::detail {

struct TrustRegionOptions
{
    double initial_radius{ 1.0 };
    double minimum_radius{ 0.0625 };
    double maximum_radius{ 4.0 };
    double shrink_factor{ 0.5 };
    double growth_factor{ 2.0 };
};

struct TrustRegionRadiusUpdate
{
    std::vector<ClusterKey> changed_key_list{};
    std::vector<ClusterKey> saturated_key_list{};
};

enum class TrustRegionRadiusAction
{
    Keep,
    Grow,
    Shrink
};

enum class TrustModelPredictionStatus
{
    Available,
    NonmaterialStep,
    ObjectiveUnavailable,
    ModelUnavailable,
    ResidualUnavailable,
    Nonfinite,
    NonpositivePrediction,
    NonmaterialPrediction
};

enum class TrustModelCandidateSource
{
    Base,
    Polish
};

enum class TrustModelTrialDisposition
{
    Accepted,
    ObjectiveRejected
};

struct TrustModelCandidateFunnel
{
    std::size_t generated_count{ 0 };
    std::size_t invalid_count{ 0 };
    std::size_t trust_skipped_count{ 0 };
    std::size_t guard_rejected_count{ 0 };
    std::size_t nonmaterial_count{ 0 };
    std::size_t objective_evaluated_count{ 0 };
    std::size_t polish_objective_evaluated_count{ 0 };
};

struct TrustModelShadowDiagnostic
{
    TrustModelPredictionStatus status{ TrustModelPredictionStatus::ModelUnavailable };
    TrustModelCandidateSource candidate_source{ TrustModelCandidateSource::Base };
    TrustModelTrialDisposition trial_disposition{
        TrustModelTrialDisposition::ObjectiveRejected };
    std::size_t search_pass{ 0 };
    std::size_t trial_number{ 0 };
    double factor{ 0.0 };
    double step_norm{ 0.0 };
    std::optional<double> actual_reduction{};
    std::optional<double> polish_reduction{};
    std::optional<double> predicted_residual_reduction{};
    std::optional<double> predicted_penalty_reduction{};
    std::optional<double> predicted_reduction{};
    std::optional<double> rho{};
    double boundary_utilization{ 0.0 };
    TrustRegionRadiusAction current_action{ TrustRegionRadiusAction::Keep };
    std::optional<TrustRegionRadiusAction> shadow_action{};
    bool objective_backtracked{ false };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    bool rejected_by_strict_polish{ false };
    bool final_local_candidate{ false };
    bool readiness_eligible{ false };
    std::size_t unselected_dependency_count{ 0 };
    double elapsed_milliseconds{ 0.0 };
};

enum class AllRejectedResolution
{
    MaximumIterations,
    BacktrackingExhausted
};

class TrustRegionStateSet
{
    TrustRegionOptions m_options{};
    std::map<ClusterKey, double> m_radius_by_key{};

public:
    explicit TrustRegionStateSet(TrustRegionOptions options = {});
    void Reconcile(const std::vector<ClusterKey> & key_list);
    double GetRadius(const ClusterKey & key) const;
    void ResetToMinimum(const std::vector<ClusterKey> & key_list);
    TrustRegionRadiusUpdate ApplyRadiusUpdates(
        const std::vector<ClusterKey> & grow_key_list,
        const std::vector<ClusterKey> & accepted_shrink_key_list,
        const std::vector<ClusterKey> & rejected_key_list,
        const std::vector<ClusterKey> & exhausted_key_list);

};


enum class SuspiciousGaussianReason
{
    None,
    InvalidModel,
    NonFiniteResponse,
    OffsetMagnitude,
    CenterSignFlip,
    RadialRebound,
    WidthGrowth,
    AmplitudeOffsetCompensation
};

using SuspiciousUpdateMask = std::vector<char>;

struct SuspiciousBlockActivity
{
    SuspiciousUpdateMask shape_fixed_atom_mask{};
    SuspiciousUpdateMask offset_fixed_atom_mask{};
    SuspiciousUpdateMask hard_failure_atom_mask{};

    SuspiciousUpdateMask BuildCombinedFixedAtomMask() const;
    bool HasActiveShape(std::size_t atom_index) const;
    bool HasActiveOffset(std::size_t atom_index) const;
};

std::size_t CountSuspiciousAtoms(const SuspiciousUpdateMask & suspicious_mask);

std::vector<double> BuildSuspiciousJointOffsetRidgeMultiplierList(const SuspiciousUpdateMask & suspicious_mask);

struct ZeroOffsetProfileDiagnostics
{
    double distance_range{ 0.0 };
    double innermost_response{ 0.0 };
    double max_abs_response{ 0.0 };
    double robust_residual_scale{ 0.0 };
    std::vector<double> radius_response_median_list{};
};

struct SuspiciousProfileAnalysis
{
    bool all_responses_finite{ true };
    std::optional<ZeroOffsetProfileDiagnostics> profile{};
};

struct SuspiciousUpdateBaseline
{
    GaussianModel3D previous_model{};
    SuspiciousProfileAnalysis previous_analysis{};
};

enum class SuspiciousUpdateMode
{
    OffsetOnly,
    PostRefit
};

struct SuspiciousGaussianAssessment
{
    SuspiciousGaussianReason reason{ SuspiciousGaussianReason::None };
    SuspiciousUpdateMode mode{ SuspiciousUpdateMode::PostRefit };
    double normalized_margin{ -std::numeric_limits<double>::infinity() };

    bool IsSuspicious() const { return reason != SuspiciousGaussianReason::None; }
};

SuspiciousGaussianAssessment AssessSuspiciousGaussianUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & candidate_model,
    const FitOptions & options,
    const SuspiciousUpdateBaseline & previous_baseline,
    SuspiciousUpdateMode mode);

SuspiciousUpdateBaseline BuildPreviousSuspiciousProfileBaseline(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const FitOptions & options);

SuspiciousUpdateMask ExpandSuspiciousSharedOffsetGroups(
    const std::vector<std::size_t> & group_id_by_position,
    const SuspiciousUpdateMask & suspicious_seed_mask);

class PerformanceCounters
{
    const bool m_quiet_mode;
    const ClusterSolverWorkspaceMap & m_solver_workspace_by_key;
    const BoundaryJointCorrectionWorkspaceMap & m_boundary_joint_correction_workspace_by_key;
    const std::chrono::steady_clock::time_point m_start_time;
    const std::size_t m_cached_sample_count;
    std::atomic<std::size_t> m_full_state_materialization_count{ 0 };
    std::atomic<std::size_t> m_gaussian_cache_hit_count{ 0 };
    std::atomic<std::size_t> m_gaussian_cache_miss_count{ 0 };
    std::atomic<std::size_t> m_objective_recomputed_sample_count{ 0 };
    std::atomic<std::size_t> m_objective_reused_sample_count{ 0 };
    std::size_t m_retired_solver_symbolic_analysis_count{ 0 };
    std::size_t m_topology_rebuild_attempt_count{ 0 };
    std::size_t m_topology_partition_change_count{ 0 };
    std::size_t m_boundary_reconciliation_attempt_count{ 0 };
    std::size_t m_boundary_reconciliation_backtracked_count{ 0 };
    std::size_t m_boundary_reconciliation_rejected_count{ 0 };
    std::size_t m_boundary_joint_correction_attempt_count{ 0 };
    std::size_t m_boundary_joint_correction_accepted_count{ 0 };
    std::size_t m_boundary_joint_correction_fallback_count{ 0 };
    std::size_t m_boundary_rescue_attempt_count{ 0 };
    std::size_t m_boundary_rescue_accepted_count{ 0 };
    std::size_t m_boundary_rescue_fallback_count{ 0 };
    std::size_t m_boundary_rescue_rejected_count{ 0 };
    std::size_t m_boundary_rescue_suspicious_exclusion_count{ 0 };
    std::size_t m_boundary_rescue_hard_failure_exclusion_count{ 0 };
    std::size_t m_boundary_rescue_invalid_proposal_exclusion_count{ 0 };
    std::size_t m_boundary_rescue_objective_unavailable_exclusion_count{ 0 };
    std::size_t m_dependency_polish_component_count{ 0 };
    std::size_t m_dependency_polish_attempt_count{ 0 };
    std::size_t m_dependency_polish_accepted_count{ 0 };
    std::size_t m_dependency_polish_fallback_count{ 0 };
    std::size_t m_dependency_polish_atom_count{ 0 };
    std::size_t m_dependency_polish_parameter_count{ 0 };
    std::size_t m_dependency_polish_round_count{ 0 };
    double m_iteration_phase_milliseconds{ 0.0 };
    double m_candidate_phase_milliseconds{ 0.0 };
    double m_topology_rebuild_milliseconds{ 0.0 };
    double m_boundary_reconciliation_milliseconds{ 0.0 };
    double m_boundary_joint_correction_milliseconds{ 0.0 };
    double m_dependency_polish_milliseconds{ 0.0 };

public:
    PerformanceCounters(
        bool quiet_mode,
        const SecondStageContext & context,
        const ClusterSolverWorkspaceMap & solver_workspace_by_key,
        const BoundaryJointCorrectionWorkspaceMap & boundary_joint_correction_workspace_by_key);

    ~PerformanceCounters();

    void RecordFullStateMaterialization();
    void RecordGaussianCacheMisses();
    void RecordGaussianCacheHits();
    void RecordObjectiveSampleEvaluation(std::size_t recomputed_sample_count, std::size_t total_sample_count);
    [[nodiscard]] std::chrono::steady_clock::time_point StartIterationPhase() const;
    void FinishIterationPhase(std::chrono::steady_clock::time_point start_time);
    [[nodiscard]] std::chrono::steady_clock::time_point StartCandidatePhase() const;
    void FinishCandidatePhase(std::chrono::steady_clock::time_point start_time);
    void RecordSolverWorkspaceReset();
    void RecordTopologyRebuild(double elapsed_milliseconds, bool partition_changed);
    void RecordBoundaryReconciliation(
        std::size_t attempt_count,
        std::size_t backtracked_count,
        std::size_t rejected_count,
        double elapsed_milliseconds);
    void RecordBoundaryJointCorrection(bool accepted, double elapsed_milliseconds);
    void RecordBoundaryRescue(bool accepted, bool used_fallback);
    void RecordBoundaryRescueExclusions(
        std::size_t suspicious_count,
        std::size_t hard_failure_count,
        std::size_t invalid_proposal_count,
        std::size_t objective_unavailable_count);
    void RecordDependencyPolish(
        std::size_t component_count,
        std::size_t attempt_count,
        std::size_t accepted_count,
        std::size_t fallback_count,
        std::size_t atom_count,
        std::size_t parameter_count,
        std::size_t round_count,
        double elapsed_milliseconds);

private:
    static double CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time);
    static std::size_t CountRawSamplingEntries(const SecondStageContext & context);
    std::size_t CountCurrentSolverSymbolicAnalyses() const;
};

constexpr double kTailValidationWeight{ 0.25 };

struct ObjectiveBreakdown
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };

    constexpr double GetTailValidationPenalty() const noexcept
    {
        return kTailValidationWeight * tail_validation_loss;
    }

    constexpr double GetTotalObjective() const noexcept
    {
        return fit_range_residual_objective + GetTailValidationPenalty() + offset_plausibility_penalty;
    }
};

struct ObjectiveTolerance
{
    double absolute_tolerance{ 0.0 };
    double relative_tolerance{ 0.0 };
};

struct AuditedState
{
    ObjectiveBreakdown objective{};
    FitState state{};
    bool uses_polish{ false };
    std::size_t source_iteration{ 0 };
};

using BestAuditState = std::optional<AuditedState>;

class CandidateEvaluationOverlay
{
    const SecondStageContext & m_context;
    const ResidualBaseline & m_baseline;
    const FitStateView & m_candidate_state;
    std::vector<char> m_changed_group_mask{};
    std::vector<std::optional<GaussianModel3D>> m_changed_group_median{};

public:
    CandidateEvaluationOverlay(
        const SecondStageContext & context,
        const ResidualBaseline & baseline,
        const FitStateView & candidate_state);

    std::optional<ResidualSample> operator()(const SampleRef & sample_ref) const;
    const FitStateView & GetState() const { return m_candidate_state; }
    const ResidualBaseline & GetBaseline() const { return m_baseline; }
};

enum class PreObjectiveFailureReason
{
    None,
    InvalidModel,
    NoCandidateWithinTrustRegion
};

enum class StabilizationTerminalReason
{
    None,
    GuardInfeasible,
    ObjectiveExhausted,
    InvalidCandidate
};

struct StabilizationTerminalDiagnostic
{
    StabilizationTerminalReason reason{ StabilizationTerminalReason::None };
    std::optional<std::size_t> guard_atom_index{};
    std::optional<SuspiciousUpdateMode> guard_mode{};
    std::optional<SuspiciousGaussianReason> guard_reason{};
};

struct ObjectiveScale
{
    double fit{ 0.0 };
    double tail{ 0.0 };
};

struct ObjectiveAttemptDiagnostic
{
    std::optional<double> accepted_factor{};
    PreObjectiveFailureReason pre_objective_failure_reason{ PreObjectiveFailureReason::None };
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<ObjectiveScale> scale{};
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
    std::optional<ObjectiveBreakdown> candidate_objective{};
    std::optional<ObjectiveBreakdown> previous_objective{};
    std::optional<ObjectiveBreakdown> best_objective{};
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    std::size_t trial_count{ 0 };
    std::size_t invalid_trial_count{ 0 };
    std::size_t trust_skipped_trial_count{ 0 };
    std::size_t guard_rejected_trial_count{ 0 };
    std::size_t objective_rejected_trial_count{ 0 };
    std::vector<StabilizationTerminalDiagnostic> terminal_diagnostic_list{};
};

struct ObjectiveDomain;

TrustRegionRadiusAction DetermineAcceptedTrustRegionRadiusAction(
    std::optional<double> first_objective_evaluated_factor,
    const ObjectiveAttemptDiagnostic & diagnostic);

TrustRegionRadiusAction DetermineTrustModelShadowAction(
    const TrustModelShadowDiagnostic & diagnostic);

TrustModelShadowDiagnostic EvaluateTrustModelShadow(
    const SecondStageContext & context,
    const ResidualBaseline & residual_baseline,
    const FitState & previous_state,
    const FitStatePatch & candidate_patch,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveDomain & objective_domain,
    const std::optional<ObjectiveBreakdown> & previous_objective,
    const std::optional<ObjectiveBreakdown> & candidate_objective,
    double trust_region_radius,
    TrustRegionRadiusAction current_action,
    TrustModelCandidateSource candidate_source,
    bool objective_backtracked);

double CalculateClusterAtomWeight(std::size_t cluster_atom_count, std::size_t active_atom_count);

std::optional<ObjectiveBreakdown> BuildObjectiveBreakdown(
    double fit_range_residual_objective,
    double tail_validation_loss,
    double offset_plausibility_penalty);

bool IsBetterAuditObjective(double candidate, double best, ObjectiveTolerance tolerance);

bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance);

struct ObjectiveClusterDomain
{
    std::vector<SampleRef> fit_sample_ref_list{};
    std::vector<SampleRef> tail_sample_ref_list{};
    std::optional<ObjectiveScale> scale{};
};

struct ObjectiveDomain
{
    std::map<ClusterKey, ObjectiveClusterDomain> cluster_by_key{};
    std::vector<ClusterKey> owner_key_by_atom_index{};
    std::vector<std::vector<char>> fit_sample_mask_by_atom{};
    std::size_t active_atom_count{ 0 };
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
};

void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset = false);

struct ClusterObjectiveState
{
    std::optional<ObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
};

using ClusterObjectiveStateMap = std::map<ClusterKey, ClusterObjectiveState>;
using ObjectiveByKey = std::map<ClusterKey, std::optional<ObjectiveBreakdown>>;

ObjectiveDomain BuildObjectiveDomain(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<ClusterKey> & cluster_key_list,
    double distance_min,
    double distance_max);


std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator);

std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const SnapshotResidualEvaluator & evaluator);

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator);

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const SnapshotResidualEvaluator & evaluator);

bool TryUpdateBestAuditState(
    const FitState & candidate_state,
    bool candidate_uses_polish,
    std::size_t source_iteration,
    const ObjectiveBreakdown & candidate_objective,
    BestAuditState & audit_state);

void ReconcileClusterObjectiveState(
    const ObjectiveByKey & previous_objective_by_key,
    ClusterObjectiveStateMap & state_by_key);

enum class BacktrackingStepStatus
{
    CandidateReady,
    InvalidCandidate,
    Exhausted
};

struct BacktrackingStep
{
    BacktrackingStepStatus status{ BacktrackingStepStatus::Exhausted };
    double factor{ 0.0 };
    std::size_t trial_number{ 0 };
};

class BacktrackingWorkspace
{
    std::size_t m_previous_state_size{ 0 };
    double m_minimum_transformed_change{ 0.0 };
    double m_next_factor{ 0.5 };
    std::size_t m_trial_number{ 1 };
    std::vector<GaussianModel3D> m_previous_model_list{};
    std::vector<GaussianModel3D> m_endpoint_model_list{};
    std::vector<double> m_previous_shared_offset_list{};
    std::vector<double> m_endpoint_shared_offset_list{};
    FitStatePatch m_candidate_patch{};

public:
    BacktrackingWorkspace(
        const SecondStageContext & context,
        const FitState & previous_state,
        const FitStatePatch & endpoint_patch,
        double minimum_transformed_change);

    BacktrackingWorkspace(const BacktrackingWorkspace &) = delete;
    BacktrackingWorkspace & operator=(const BacktrackingWorkspace &) = delete;
    BacktrackingStep BuildNextCandidate();
    FitStatePatch TakeCandidatePatch() { return std::move(m_candidate_patch); }
    const FitStatePatch & GetCandidatePatch() const { return m_candidate_patch; }

    PolishProvenance BuildCandidatePolishProvenance(
        const PolishProvenance & previous_provenance,
        const PolishProvenance & endpoint_provenance) const;

private:
    bool BuildCandidate(double factor);
    bool HasMaterialChange(std::size_t atom_position) const;
    double GetMaximumTransformedChange() const;

};

struct ClusterCandidateDiagnostic
{
    ClusterKey key{};
    ObjectiveAttemptDiagnostic attempt{};
#ifdef RHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT
    std::vector<TrustModelShadowDiagnostic> trust_model_shadow_trial_list{};
    TrustModelCandidateFunnel trust_model_candidate_funnel{};
    bool boundary_touched{ false };
#endif
    bool boundary_rescued{ false };
};

struct PolishProgress
{
    std::size_t eligible_count{ 0 };
    std::size_t accepted_count{ 0 };
    std::size_t rejected_count{ 0 };
    std::size_t skipped_count{ 0 };
};

enum class BoundaryComponentAcceptedSource
{
    None,
    Endpoint,
    JointCorrection,
    Backtracking
};

struct BoundaryComponentReconciliationDiagnostic
{
    std::vector<ClusterKey> key_list{};
    std::size_t atom_count{ 0 };
    std::size_t boundary_sample_count{ 0 };
    std::size_t trial_count{ 1 };
    std::optional<double> accepted_factor{};
    BoundaryComponentAcceptedSource accepted_source{ BoundaryComponentAcceptedSource::None };
    std::optional<BoundaryJointCorrectionStatus> joint_correction_status{};
    std::size_t interface_atom_count{ 0 };
    std::size_t shape_active_atom_count{ 0 };
    std::size_t offset_active_atom_count{ 0 };
    std::size_t offset_closure_atom_count{ 0 };
    std::size_t suspicious_candidate_atom_count{ 0 };
    std::size_t joint_parameter_count{ 0 };
    std::optional<double> joint_damping{};
    std::optional<double> maximum_normalized_trust_step{};
    std::optional<double> previous_component_objective{};
    std::optional<double> endpoint_component_objective{};
    std::optional<double> joint_reference_component_objective{};
    std::optional<double> joint_candidate_component_objective{};
    std::optional<double> candidate_component_objective{};
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rescue_candidate_cluster_count{ 0 };
    std::size_t rescued_cluster_count{ 0 };
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
    std::optional<double> component_improvement{};
    std::optional<double> global_improvement{};
    bool is_rescue_attempt{ false };
    bool accepted{ false };
    bool exhausted{ false };
};

struct CandidateSelection
{
    FitState assembled_state{};
    PolishProvenance assembled_polish_provenance{};
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    std::vector<ClusterKey> grow_trust_region_key_list{};
    std::vector<ClusterKey> shrink_trust_region_key_list{};
    std::vector<ClusterKey> exhausted_key_list{};
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<BoundaryComponentReconciliationDiagnostic> boundary_reconciliation_diagnostic_list{};
    std::optional<ObjectiveBreakdown> final_audit_objective{};
    PolishProgress polish_progress{};
};

struct CandidateSelectionInputs
{
    const SecondStageContext & context;
    const FitOptions & options;
    const ResidualBaseline & residual_baseline;
    const CouplingGraphPartition & partition;
    const ClusterHealthMap & health_by_key;
    const FitState & previous_state;
    const PolishProvenance & previous_polish_provenance;
    const FitState & operator_proposal_state;
    SuspiciousBlockActivity & block_activity;
    const std::vector<double> & ridge_multiplier_list;
    const ObjectiveDomain & objective_domain;
    const ObjectiveByKey & previous_objective_by_key;
    ClusterObjectiveStateMap & cluster_objective_state;
    const BestAuditState & best_audit_state;
    const TrustRegionStateSet & trust_region_state;
    ClusterSolverWorkspaceMap & solver_workspace_by_key;
    BoundaryJointCorrectionWorkspaceMap & boundary_joint_correction_workspace_by_key;
    int thread_size;
    PerformanceCounters & performance_counters;
};

CandidateSelection SelectClusterCandidates(const CandidateSelectionInputs & inputs);

struct FinalDependencyPolishDiagnostic
{
    std::size_t component_count{ 0 };
    std::size_t attempted_component_count{ 0 };
    std::size_t accepted_component_count{ 0 };
    std::size_t fallback_component_count{ 0 };
    std::size_t atom_count{ 0 };
    std::size_t parameter_count{ 0 };
    std::size_t round_count{ 0 };
    std::size_t suspicious_candidate_atom_count{ 0 };
    std::optional<double> objective_before{};
    std::optional<double> objective_after{};
    double elapsed_milliseconds{ 0.0 };
    struct Component
    {
        std::vector<ClusterKey> key_list{};
        std::size_t atom_count{ 0 };
        std::size_t parameter_count{ 0 };
        std::size_t round_count{ 0 };
        std::size_t suspicious_candidate_atom_count{ 0 };
        std::size_t symbolic_analysis_count{ 0 };
        std::optional<double> objective_before{};
        std::optional<double> objective_after{};
        double elapsed_milliseconds{ 0.0 };
        bool accepted{ false };
        bool fallback{ true };
    };
    std::vector<Component> component_list{};
};

struct FinalDependencyPolishResult
{
    FitState state{};
    std::optional<ObjectiveBreakdown> objective{};
    FinalDependencyPolishDiagnostic diagnostic{};
    bool accepted{ false };
};

FinalDependencyPolishResult RunFinalDependencyPolish(
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & objective_domain,
    const SuspiciousBlockActivity & block_activity,
    const TrustRegionStateSet & trust_region_state,
    const FitState & base_state,
    BoundaryJointCorrectionWorkspaceMap & workspace_by_key,
    PerformanceCounters & performance_counters);

} // namespace rhbm_gem::core::detail
