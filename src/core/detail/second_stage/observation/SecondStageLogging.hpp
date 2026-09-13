#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

namespace rhbm_gem::core::detail {

class PerformanceCounters;
struct IterationObservation;
struct IterationDiagnostics;
struct FinalDependencyPolishDiagnostic;
struct SuspiciousGaussianAssessment;
struct SuspiciousBlockActivity;

void LogSecondStagePerformance(const PerformanceCounters &, std::size_t, double);

struct IterationResult;
struct FixedPointOperatorEvidence;
struct ConvergenceCertificate;
struct ConvergenceDiagnostics;
struct ConvergenceAssessment;
struct SecondStageSeedSelectionRecord;
struct ClusterCandidateDiagnostic;
struct FinalDependencyPolishResult;
enum class SecondStageStopReason;
enum class FinalPolishResidualSafetyStatus;

void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset = false);

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode);

using ProgressColumnWidths = std::array<std::size_t, 6>;
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
    const IterationObservation & iteration_result);

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count, std::size_t maximum_iterations);
void LogProgressHeader(bool quiet_mode, const ProgressColumnWidths & column_widths);
void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result, const IterationDiagnostics & diagnostics);

void LogOperatorAvailability(std::string_view phase, const FixedPointOperatorEvidence &,
    const std::vector<std::size_t> & atom_index_list);
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
    const FinalDependencyPolishDiagnostic & diagnostic,
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
