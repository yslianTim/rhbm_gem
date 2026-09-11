#pragma once

#include "core/detail/ObjectiveEvaluation.hpp"
#include "core/detail/JointFitting.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <span>
#include <vector>

namespace rhbm_gem::core::detail {

struct IterationResult;
struct ConvergenceCertificate;
struct ConvergenceDiagnostics;
struct ConvergenceAssessment;
struct SecondStageSeedSelectionRecord;
struct ClusterCandidateDiagnostic;
struct FinalDependencyPolishResult;
enum class SecondStageStopReason;
enum class FinalPolishResidualSafetyStatus;

JointCandidateObjectiveDiagnostic * BeginJointCandidateDiagnostic(
    bool quiet_mode,
    std::vector<JointCandidateObjectiveDiagnostic> & records,
    std::string_view source,
    std::optional<double> factor = std::nullopt,
    std::size_t round = 0);

void RecordJointMemberRejection(
    JointCandidateObjectiveDiagnostic * record,
    const ClusterKey & key,
    const std::optional<ObjectiveBreakdown> & previous,
    const std::optional<ObjectiveBreakdown> & best,
    const std::optional<ObjectiveBreakdown> & candidate,
    bool best_checked);

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
    void FinishIterationPhase(std::chrono::steady_clock::time_point start_time);
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

void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset = false);

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode);

using ProgressColumnWidths = std::array<std::size_t, 6>;
bool IsDebugLogLevelEnabled();
void FinishProgressLine(bool quiet_mode);
void LogSecondStageStart(bool quiet_mode);
void LogSecondStageInitializationFailure(bool quiet_mode);
void LogNoSelectedAtoms(bool quiet_mode);

void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode);
void LogFrozenBackground(const SecondStageContext & context, bool quiet_mode);

void LogRejectedClusterDiagnostics(
    bool quiet_mode,
    const std::vector<ClusterCandidateDiagnostic> & diagnostic_list);
void LogAllRejectedResolution(
    bool quiet_mode,
    const IterationResult & iteration_result);
void LogAcceptedCandidateSearchDiagnostics(
    bool quiet_mode,
    const IterationResult & iteration_result);

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count);
void LogProgressHeader(bool quiet_mode, const ProgressColumnWidths & column_widths);
void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result);

void LogUnrestrictedOperatorAssessments(
    bool quiet_mode,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    const SuspiciousBlockActivity & block_activity);
void LogConvergenceSafeguardAudit(
    bool quiet_mode,
    const IterationResult & iteration_result,
    const ConvergenceCertificate & certificate,
    const ConvergenceDiagnostics & diagnostics);

void LogAdaptiveTopologyRebuild(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    double maximum_transformed_drift,
    const GraphTopology & previous_topology,
    const GraphTopology & rebuilt_topology,
    const CouplingGraphPartition & previous_partition,
    const CouplingGraphPartition & rebuilt_partition,
    bool partition_changed);

void LogFinalDependencyPolish(
    bool quiet_mode,
    const FinalDependencyPolishResult & polish_result,
    FinalPolishResidualSafetyStatus safety_status,
    bool applied,
    const ConvergenceAssessment * candidate_certificate = nullptr);

void LogSecondStageAuditTerminal(
    bool quiet_mode,
    const SecondStageContext & context,
    SecondStageStopReason reason,
    std::size_t attempt_number,
    std::size_t accepted_iteration_count,
    const FitState & finalized_state,
    const ObjectiveDomain & comparison_objective_domain);
void LogQuarantineFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const FitState & finalized_state);
void LogConverged(
    bool quiet_mode,
    const IterationResult & iteration_result,
    const FitState & finalized_state);
void LogMaximumIterations(
    bool quiet_mode,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const BestAuditState & best_audit_state,
    const FitState & latest_state);
void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const BestAuditState & best_audit_state,
    const PolishProvenance & latest_polish_provenance,
    SecondStageStopReason stop_reason,
    bool final_uses_best_audit);

} // namespace rhbm_gem::core::detail
