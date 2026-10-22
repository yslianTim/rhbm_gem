#pragma once

#include "core/detail/second_stage/ObjectiveEvaluation.hpp"

#include <array>
#include <cstddef>
#include <string_view>

namespace rhbm_gem::core::detail {

bool IsDebugLogLevelEnabled();

struct IterationResult;
enum class SecondStageStopReason;
std::string_view SecondStageStopReasonText(SecondStageStopReason);
void LogGraphFallbackWarning(bool quiet_mode, bool uses_weighted_graph);

using ProgressColumnWidths = std::array<std::size_t, 2>;
void FinishProgressLine(bool quiet_mode);
void LogSecondStageStart(bool quiet_mode);
void LogSecondStageInitializationFailure(bool quiet_mode);
void LogNoSelectedAtoms(bool quiet_mode);

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count, std::size_t maximum_iterations);
void LogProgressHeader(bool quiet_mode, const ProgressColumnWidths & column_widths);
void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result);

void LogQuarantineFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count);
void LogConverged(
    bool quiet_mode,
    std::size_t accepted_iteration_count);
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

} // namespace rhbm_gem::core::detail
