#include "core/detail/second_stage/observation/SecondStageLogging.hpp"

#include "core/detail/second_stage/IterationResult.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

#include <rhbm_gem/utils/domain/Logger.hpp>

namespace rhbm_gem::core::detail {
namespace {

constexpr std::array<std::string_view, 2> kProgressHeaderList{
    "Try/Acc",
    "Cluster A/R"
};

void AppendQuarantineSummary(
    std::ostream & stream,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count)
{
    stream << "; quarantine entered/released/probation-failed/unresolved = "
        << entered_target_count << "/"
        << released_target_count << "/"
        << failed_probation_count << "/"
        << unresolved_target_count;
}

void AppendAuditSummary(std::ostringstream & stream, const AuditedState & audited_state)
{
    stream << "; best validated audit objective = "
        << std::scientific << std::setprecision(2)
        << audited_state.objective.GetTotalObjective();
}

std::string FormatProgressRow(
    const ProgressColumnWidths & column_widths,
    const std::array<std::string, 2> & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream << std::left
            << std::setw(static_cast<int>(column_widths.at(i)))
            << cell_list.at(i);
    }
    return stream.str();
}

} // namespace

std::string_view SecondStageStopReasonText(SecondStageStopReason reason)
{
    switch (reason)
    {
    case SecondStageStopReason::None:
        return "";
    case SecondStageStopReason::Quarantine:
        return "quarantine";
    case SecondStageStopReason::Converged:
        return "converged";
    case SecondStageStopReason::AuditPatience:
        return "audit-patience";
    case SecondStageStopReason::AllRejectedBacktrackingExhausted:
        return "all-rejected-backtracking-exhausted";
    case SecondStageStopReason::AllRejectedAtMaximumIterations:
    case SecondStageStopReason::RecoveryFailed:
        return "recovery-failed";
    case SecondStageStopReason::FinalCertificateFailed:
        return "final-certificate-failed";
    case SecondStageStopReason::MaximumIterations:
        return "maximum-iterations";
    }
    return "";
}

bool IsDebugLogLevelEnabled()
{
    return Logger::GetLogLevel() >= LogLevel::Debug;
}

void FinishProgressLine(bool quiet_mode)
{
    if (!quiet_mode) Logger::FinishProgressLine();
}

void LogSecondStageStart(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Info,
        "Run 2nd-stage local atom fitting with iterations...");
}

void LogSecondStageInitializationFailure(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Warning,
        "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
        "is available for every selected atom.");
}

void LogNoSelectedAtoms(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Skip 2nd-stage local atom fitting because no atoms are selected.");
}

void LogGraphFallbackWarning(bool quiet_mode, bool uses_weighted_graph)
{
    if (quiet_mode || uses_weighted_graph) return;
    Logger::Log(LogLevel::Warning,
        "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
}

ProgressColumnWidths BuildProgressColumnWidths(
    std::size_t atom_count,
    std::size_t maximum_iterations)
{
    const auto maximum_iteration_text{ std::to_string(maximum_iterations) };
    const auto maximum_atom_text{ std::to_string(atom_count) };
    const std::array<std::string, 2> maximum_cell_list{
        maximum_iteration_text + "/" + maximum_iteration_text,
        maximum_atom_text + "/" + maximum_atom_text
    };
    ProgressColumnWidths column_widths;
    for (std::size_t i = 0; i < column_widths.size(); i++)
    {
        column_widths.at(i) = std::max(
            kProgressHeaderList.at(i).size(),
            maximum_cell_list.at(i).size());
    }
    return column_widths;
}

void LogProgressHeader(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Info,
        FormatProgressRow(column_widths, { "Try/Acc", "Cluster A/R" }));
}

void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result)
{
    if (quiet_mode) return;
    const std::array<std::string, 2> cell_list{
        std::to_string(iteration_result.attempt_number) + "/" +
            std::to_string(iteration_result.accepted_iteration_count),
        std::to_string(iteration_result.accepted_key_list.size()) + "/" +
            std::to_string(iteration_result.rejected_key_list.size())
    };
    Logger::ProgressLine(FormatProgressRow(column_widths, cell_list));
}

void LogQuarantineFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after "
        << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendQuarantineSummary(
        warning_message,
        entered_target_count,
        released_target_count,
        failed_probation_count,
        unresolved_target_count);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogConverged(bool quiet_mode, std::size_t accepted_iteration_count)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Converged after " + std::to_string(accepted_iteration_count) +
        " accepted iterations.");
}

void LogMaximumIterations(
    bool quiet_mode,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const BestAuditState & best_audit_state)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendQuarantineSummary(
        warning_message,
        entered_target_count,
        released_target_count,
        failed_probation_count,
        unresolved_target_count);
    if (best_audit_state.has_value())
    {
        warning_message << "; applying best validated audit state";
        AppendAuditSummary(warning_message, *best_audit_state);
    }
    else
    {
        warning_message << "; applying latest validated state";
    }
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const BestAuditState & best_audit_state,
    const PolishProvenance & latest_polish_provenance,
    SecondStageStopReason stop_reason,
    bool final_uses_best_audit)
{
    if (quiet_mode) return;
    const auto final_uses_polish{
        final_uses_best_audit && best_audit_state.has_value() ?
            best_audit_state->uses_polish :
            std::ranges::any_of(
                latest_polish_provenance,
                [](char value) { return value != 0; })
    };
    const auto final_state_source{ final_uses_best_audit ? "best-audit" : "latest-validated" };
    std::ostringstream message;
    message << "Second-stage fitting summary: accepted_iterations="
        << accepted_iteration_count
        << ", stop_reason=" << SecondStageStopReasonText(stop_reason)
        << ", final_state=" << final_state_source
        << ", final_uses_polish=" << (final_uses_polish ? "yes" : "no") << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace rhbm_gem::core::detail
