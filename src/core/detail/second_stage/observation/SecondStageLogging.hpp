#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

namespace rhbm_gem::core { struct FitOptions; }
namespace rhbm_gem::core::detail {
class SecondStageObservationSession;

class PerformanceCounters;
struct IterationDiagnostics;
struct FinalDependencyPolishDiagnostic;
struct SuspiciousGaussianAssessment;
struct SuspiciousBlockActivity;

void LogSecondStagePerformance(const PerformanceCounters &, std::size_t, double);

struct IterationResult;
struct FixedPointOperatorEvidence;
struct ConvergenceCertificate;
struct ConvergenceDiagnostics;
struct SecondStageSeedSummary;
struct FinalDependencyPolishResult;
enum class SecondStageStopReason;
std::string_view SecondStageStopReasonText(SecondStageStopReason);
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
    const SecondStageSeedSummary & selection_record_list,
    bool quiet_mode);


ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count, std::size_t maximum_iterations);
void LogProgressHeader(bool quiet_mode, const ProgressColumnWidths & column_widths);
void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result, const IterationDiagnostics & diagnostics);


void LogAdaptiveTopologyRebuild(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    double maximum_transformed_drift,
    const CouplingGraphPartition & previous_partition,
    const CouplingGraphPartition & rebuilt_partition,
    bool partition_changed);

void LogFinalDependencyPolish(
    bool quiet_mode,
    const FinalDependencyPolishResult & polish_result,
    const FinalDependencyPolishDiagnostic & diagnostic,
    FinalPolishResidualSafetyStatus safety_status,
    bool applied);

void LogQuarantineFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count);
void LogConverged(
    bool quiet_mode,
    const IterationResult & iteration_result);
void LogMaximumIterations(
    bool quiet_mode,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const BestAuditState & best_audit_state);
void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const BestAuditState & best_audit_state,
    const PolishProvenance & latest_polish_provenance,
    SecondStageStopReason stop_reason,
    bool final_uses_best_audit);

void LogDecisionAuditStart(SecondStageObservationSession &, const FitOptions &) noexcept;
void LogDecisionAuditIteration(SecondStageObservationSession &, const IterationResult &) noexcept;
void LogDecisionAuditTerminal(SecondStageObservationSession &, std::string_view reason, std::string_view source,
    const BestAuditState &, const PerformanceCounters * = nullptr) noexcept;

} // namespace rhbm_gem::core::detail
