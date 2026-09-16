#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include "core/detail/second_stage/observation/SecondStageLogging.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"

#include "core/detail/second_stage/IterationResult.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/DependencyPolish.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <set>
#include <sstream>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr std::array<std::string_view, 6> kProgressHeaderList{
    "Try/Acc",
    "Atom A/Q",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/G"
};



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
    case SecondStageStopReason::RecoveryFailed: return "recovery-failed";
    case SecondStageStopReason::FinalCertificateFailed: return "final-certificate-failed";
    case SecondStageStopReason::MaximumIterations:
        return "maximum-iterations";
    }
    return "";
}

namespace {
std::string_view GetFinalPolishResidualSafetyStatusText(
    FinalPolishResidualSafetyStatus status)
{
    switch (status)
    {
    case FinalPolishResidualSafetyStatus::NotEvaluated:
        return "not-evaluated";
    case FinalPolishResidualSafetyStatus::AbsolutePassed:
        return "absolute-passed";
    case FinalPolishResidualSafetyStatus::Failed:
        return "failed";
    case FinalPolishResidualSafetyStatus::Error:
        return "error";
    }
    throw std::invalid_argument("Unknown final polish residual safety status.");
}

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










void AppendAuditSummary(
    std::ostringstream & stream,
    const AuditedState & audited_state)
{
    const auto & objective{ audited_state.objective };
    stream << "; audit best source = ";
    if (audited_state.source_iteration != 0)
    {
        stream << "accepted iteration " << audited_state.source_iteration;
    }
    else
    {
        stream << "initial";
    }
    stream
        << std::scientific << std::setprecision(2)
        << ", fixed audit objective fit/tail-weighted/offset/total = "
        << objective.fit_range_residual_objective << "/"
        << objective.GetTailValidationPenalty() << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.GetTotalObjective()
        << ", tail raw/weight = " << objective.tail_validation_loss << "/"
        << kTailValidationWeight;
}


} // namespace

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
    Logger::Log(LogLevel::Info,
        "Second-stage local fitting summary: accepted_iterations=0, "
        "best_iteration=unavailable, stop_reason=no-valid-seed, "
        "best_audit_objective=unavailable, final_uses_polish=unavailable, "
        "final_state_source=unavailable.");
}

void LogNoSelectedAtoms(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Skip 2nd-stage local atom fitting because no atoms are selected.");
}

void LogSecondStageSeedSelections(const SecondStageSeedSummary & seeds, bool quiet_mode)
{
    if (quiet_mode || seeds.local_mdpde + seeds.global_median == 0) return;
    std::ostringstream message;
    message << "Selected second-stage initial seeds = " << seeds.local_mdpde + seeds.global_median
        << ", sources = local-mdpde:" << seeds.local_mdpde << ", global-median:" << seeds.global_median << ".";
    Logger::Log(LogLevel::Info, message.str());
}





namespace {

std::string FormatProgressMaximum(double value)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

template<typename CellList>
std::string FormatProgressRow(
    const ProgressColumnWidths & column_widths,
    const CellList & cell_list)
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

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count, std::size_t maximum_iterations)
{
    const auto maximum_iteration_text{ std::to_string(maximum_iterations) };
    const auto maximum_atom_text{ std::to_string(atom_count) };
    const auto maximum_change_text{
        FormatProgressMaximum(std::numeric_limits<double>::max())
    };
    const std::array<std::string, 6> maximum_cell_list{
        maximum_iteration_text + "/" + maximum_iteration_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text + "/" +
            maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text,
        maximum_change_text + "/" + maximum_change_text
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
    Logger::Log(
        LogLevel::Info,
        FormatProgressRow(column_widths, kProgressHeaderList));
}

void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result, const IterationDiagnostics & diagnostics)
{
    if (quiet_mode) return;
    const std::array<std::string, 6> cell_list{
        std::to_string(iteration_result.attempt_number) + "/" +
            std::to_string(iteration_result.accepted_iteration_count),
        std::to_string(iteration_result.active_atom_count) + "/" +
            std::to_string(iteration_result.quarantine_atom_count),
        std::to_string(
            iteration_result.accepted_key_list.size()) + "/" +
            std::to_string(
                iteration_result.rejected_key_list.size()),
        std::to_string(iteration_result.polish_progress.eligible_count) + "/" +
            std::to_string(iteration_result.polish_progress.accepted_count) + "/" +
            std::to_string(iteration_result.polish_progress.rejected_count) + "/" +
            std::to_string(iteration_result.polish_progress.skipped_count),
        std::to_string(iteration_result.suspicious_atom_count),
        (diagnostics.accepted_maximum_transformed_change.has_value() ?
            FormatProgressMaximum(
                *diagnostics.accepted_maximum_transformed_change) :
            std::string{ "-" }) + "/" +
            FormatProgressMaximum(
                diagnostics.proposal_maximum_transformed_change)
    };
    Logger::ProgressLine(FormatProgressRow(column_widths, cell_list));
}




void LogAdaptiveTopologyRebuild(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    double maximum_transformed_drift,
    const CouplingGraphPartition & previous_partition,
    const CouplingGraphPartition & rebuilt_partition,
    bool partition_changed)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Adaptive local-fitting topology rebuild: accepted_iteration="
        << accepted_iteration_count
        << ", trigger=drift"
        << std::scientific << std::setprecision(2)
        << ", drift=" << maximum_transformed_drift
        << ", clusters="
        << previous_partition.sample_id_list_by_key.size() << "/"
        << rebuilt_partition.sample_id_list_by_key.size()
        << ", boundary_samples="
        << previous_partition.boundary_sample_dependency_list.size() << "/"
        << rebuilt_partition.boundary_sample_dependency_list.size()
        << ", partition_changed="
        << (partition_changed ? "yes" : "no")
        << ", partition_pending=" << (partition_changed ? "yes" : "no")
        << ", objective_domain_reset=no.";
    Logger::Log(LogLevel::Info, message.str());
}

void LogFinalDependencyPolish(
    bool quiet_mode,
    const FinalDependencyPolishResult & polish_result,
    const FinalDependencyPolishDiagnostic & diagnostic,
    FinalPolishResidualSafetyStatus safety_status,
    bool applied)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    std::ostringstream message;
    message << std::scientific << std::setprecision(2)
        << "Final dependency polish: components/attempted/accepted/fallback="
        << diagnostic.component_count << "/"
        << diagnostic.attempted_component_count << "/"
        << diagnostic.accepted_component_count << "/"
        << diagnostic.component_count - diagnostic.accepted_component_count
        << ", atoms/parameters/rounds="
        << diagnostic.atom_count << "/"
        << diagnostic.parameter_count << "/"
        << diagnostic.round_count
        << ", objective before/after=";
    if (diagnostic.objective_before.has_value())
    {
        message << *diagnostic.objective_before;
    }
    else
    {
        message << "-";
    }
    message << "/";
    if (diagnostic.objective_after.has_value())
    {
        message << *diagnostic.objective_after;
    }
    else
    {
        message << "-";
    }
    message
        << ", accepted=" << (polish_result.accepted ? "yes" : "no")
        << ", residual-safety-policy=strict-fixed-point"
        << ", residual-safety="
        << GetFinalPolishResidualSafetyStatusText(safety_status)
        << ", applied=" << (applied ? "yes" : "no");
    message
        << ", elapsed_ms=" << std::fixed << std::setprecision(3)
        << diagnostic.elapsed_milliseconds << ".";
    Logger::Log(LogLevel::Info, message.str());

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

void LogConverged(
    bool quiet_mode,
    const IterationResult & iteration_result)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    const auto & transformed_change_percentile{
        iteration_result.transformed_change_percentile
    };
    std::ostringstream message;
    message
        << "Converged after " << iteration_result.accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_percentile.at(
            GaussianModel3D::LogPeakHeightCoordinateIndex())
        << ", percentile log-width change = "
        << transformed_change_percentile.at(
            GaussianModel3D::LogWidthCoordinateIndex())
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_percentile.at(
            GaussianModel3D::OffsetToPeakRatioCoordinateIndex());
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
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
    const auto * audit_state{
        best_audit_state.has_value() ? &*best_audit_state : nullptr
    };
    if (audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendAuditSummary(warning_message, *audit_state);
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
    Logger::FinishProgressLine();
    std::ostringstream message;
    message << " Second-Stage Local Fitting Summary : \n"
        << " - accepted_iterations = " << accepted_iteration_count << "\n"
        << " - best_iteration = ";
    if (!best_audit_state.has_value())
    {
        message << "unavailable\n";
    }
    else if (best_audit_state->source_iteration != 0)
    {
        message << best_audit_state->source_iteration << "\n";
    }
    else
    {
        message << "initial\n";
    }
    message << " - stop_reason = "
        << SecondStageStopReasonText(stop_reason) << "\n"
        << " - best_audit_objective = ";
    if (best_audit_state.has_value())
    {
        message << std::scientific << std::setprecision(2)
            << best_audit_state->objective.GetTotalObjective() << "\n";
    }
    else
    {
        message << "unavailable\n";
    }
    message << " - final_uses_polish = "
        << (final_uses_polish ? "yes" : "no") << "\n";
    message << " - final_state_source = "
        << (final_uses_best_audit ? "best-audit" : "latest-validated")
        << "\n";
    Logger::Log(LogLevel::Info, message.str());
}

void LogSecondStagePerformance(
    const PerformanceCounters & counters,
    std::size_t symbolic_analysis_count,
    double total_milliseconds)
{
    std::ostringstream message_info;
    message_info
        << " Second-Stage Local Fitting Performance :\n"
        << " - boundary_reconciliation_ms = "
        << std::fixed << std::setprecision(3) << counters.m_boundary_reconciliation_milliseconds << "\n"
        << " - boundary_joint_correction_ms = " << counters.m_boundary_joint_correction_milliseconds << "\n"
        << " - dependency_polish_ms = " << counters.m_dependency_polish_milliseconds << "\n"
        << " - iteration/candidate/topology/total_ms = "
        << std::fixed << std::setprecision(3)
        << counters.m_iteration_phase_milliseconds << "/"
        << counters.m_candidate_phase_milliseconds << "/"
        << counters.m_topology_rebuild_milliseconds << "/"
        << total_milliseconds << "\n";

    (void)symbolic_analysis_count;
    Logger::Log(LogLevel::Info, message_info.str());
}


void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset)
{
    if (quiet_mode) return;
    std::ostringstream message;
    message
        << (is_terminal_reset ?
            "Reset second-stage objective domain" : "Initialize second-stage objective domain")
        << ": fit/tail/offset weights = "
        << kFitRangeWeight << "/" << kTailValidationWeight << "/" << kOffsetPlausibilityPenaltyWeight
        << ", clusters = " << domain.cluster_by_key.size()
        << ", active atoms = " << domain.active_atom_count
        << ", unique fit/tail samples = " << domain.fit_sample_count << "/" << domain.tail_sample_count
        ;
    message << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message.str());
}

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode)
{
    if (quiet_mode) return;

    const auto & summary{ topology.summary };
    if (!summary.uses_weighted_graph)
    {
        Logger::Log(LogLevel::Warning,
            "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
    }
    std::ostringstream message;
    message << "Local-fitting coupling graph mode = "
        << (summary.uses_weighted_graph ? "weighted" : "binary-fallback")
        << std::scientific << std::setprecision(2)
        << ", minimum weight = " << summary.configured_minimum_weight
        << ", candidate/retained/cut edges = "
        << summary.candidate_edge_count << "/"
        << summary.retained_edge_count << "/"
        << summary.cut_edge_count
        << ", weight p50/p95/max = "
        << summary.weight_median << "/"
        << summary.weight_percentile_95 << "/"
        << summary.weight_maximum
        << ", initial components/max atoms/ratio = "
        << summary.component_count << "/"
        << summary.maximum_component_size << "/"
        << std::fixed << std::setprecision(2)
        << summary.maximum_component_ratio << ".";
    Logger::Log(LogLevel::Info, message.str());

    const auto & atom_cutoff_summary{ topology.atom_cutoff_summary };
    std::ostringstream atom_cutoff_message;
    atom_cutoff_message
        << "Local-fitting atom cutoff: atoms="
        << topology.adjacency_list.size()
        << ", limit=" << atom_cutoff_summary.maximum_atom_count_limit
        << ", clusters=" << summary.component_count
        << ", max-atoms=" << summary.maximum_component_size
        << ", cutoff-edges=" << atom_cutoff_summary.cut_edge_count << ".";
    Logger::Log(LogLevel::Info, atom_cutoff_message.str());

}

namespace {
constexpr std::array<std::string_view, static_cast<std::size_t>(AuditStage::Count)> kAuditStages{
    "proposal", "local-search", "local-polish", "boundary-endpoint", "boundary-correction", "boundary-backtracking",
    "rescue-endpoint", "rescue-correction", "rescue-backtracking", "selection-ordinary", "selection-rescue",
    "commit", "quarantine", "partition", "final-polish", "final-certification" };
constexpr std::array<std::string_view, static_cast<std::size_t>(AuditCategory::Count)> kAuditCategories{
    "none", "rejected", "unavailable", "invalid", "guard", "trust", "solver", "exhausted", "shrink",
    "enter", "retry", "release", "salvage", "rescue", "partition" };
void JsonNumber(std::ostream & out, double value)
{
    if (std::isfinite(value)) out << std::setprecision(17) << value;
    else out << "null";
}
void JsonOptional(std::ostream & out, const std::optional<double> & value)
{
    if (value) JsonNumber(out,*value); else out << "null";
}
void JsonObjective(std::ostream & out, const std::optional<ObjectiveBreakdown> & value)
{
    if (!value) { out << "{\"value\":null,\"reason\":\"unavailable\"}"; return; }
    const bool finite = std::isfinite(value->GetTotalObjective());
    out << "{\"value\":{" << "\"fit\":"; JsonNumber(out,value->fit_range_residual_objective);
    out << ",\"tail\":"; JsonNumber(out,value->tail_validation_loss);
    out << ",\"offset\":"; JsonNumber(out,value->offset_plausibility_penalty);
    out << ",\"total\":"; JsonNumber(out,value->GetTotalObjective());
    out << "},\"reason\":" << (finite ? "null" : "\"nonfinite\"") << "}";
}
template<class T> void JsonArray(std::ostream & out, const T & values)
{
    out << '['; bool first=true;
    for (const auto value : values) { if (!first) out << ','; first=false; JsonNumber(out,static_cast<double>(value)); }
    out << ']';
}
void JsonCertificate(std::ostream & out, const std::optional<ConvergenceAssessment> & assessment, std::string_view reference, bool outer = true)
{
    out << "{\"reference\":\"" << reference << "\",\"status\":\"" << (assessment ? "evaluated" : "not_evaluated") << '"';
    if (assessment)
    {
        const auto & c=assessment->certificate;
        if (outer) { out << ",\"accepted_active_p99\":"; JsonArray(out,c.accepted_active_p99); }
        out << ",\"operator_nominal_p99\":"; JsonArray(out,c.operator_nominal_p99);
        if (outer) { out << ",\"accepted_population\":"; JsonArray(out,assessment->diagnostics.accepted_active_movement.population_size_list); }
        out << ",\"operator_population\":"; JsonArray(out,assessment->diagnostics.operator_nominal_residual.population_size_list);
        out << ",\"complete\":" << c.operator_complete << ",\"qualified\":" << c.solver_qualified;
        const auto nonfinite=[](const auto & values) {
            return std::ranges::any_of(values,[](double value) { return !std::isfinite(value); });
        };
        out << ",\"residual_reason\":" << (nonfinite(c.operator_nominal_p99) || (outer && nonfinite(c.accepted_active_p99)) ?
            "\"nonfinite-residual\"" : (c.operator_complete ? "null" : "\"operator-unavailable\""));
        if (outer) out << ",\"blockers\":{\"domain\":" << c.objective_domain_changed << ",\"quarantine\":" << c.quarantine_transition
            << ",\"suspicious\":" << c.suspicious_block_fallback << ",\"rejected\":" << c.rejected_cluster << '}';
    }
    out << '}';
}
void JsonNominalSolves(std::ostream & out, const NominalSolveDiagnostics & diagnostics)
{
    out << "{\"shapes\":[";
    for (std::size_t i = 0; i < diagnostics.shapes.size(); ++i)
    {
        if (i) out << ',';
        const auto & solve{ diagnostics.shapes[i] };
        out << "{\"atom_index\":" << i << ",\"status\":";
        const auto effective_status{solve.EffectiveStatus()};
        if (effective_status) out << '"' << LocalRefitStatusText(*effective_status) << '"'; else out << "null";
        out << ",\"native_status\":";
        if (solve.status) out << '"' << LocalRefitStatusText(*solve.status) << '"'; else out << "null";
        const auto qualification{solve.Qualification()};
        out << ",\"qualification\":\"" << (qualification == RHBMSolveQualification::NativeSuccess ? "native-success" :
            qualification == RHBMSolveQualification::RefinedSuccess ? "refined-success" : "unqualified") << '"';
        if (solve.refinement)
        {
            const auto & r{*solve.refinement};
            out << ",\"refinement\":{\"accepted\":" << r.accepted << ",\"reason\":\"" << r.reason
                << "\",\"candidate_equation_evaluations\":" << r.candidate_equation_evaluations
                << ",\"reference_updates\":" << r.reference_updates << ",\"reference_stop\":\"" << r.reference_stop
                << "\",\"original_residual\":";
            JsonOptional(out,r.original_residual);
            out << ",\"candidate_residual\":"; JsonOptional(out,r.candidate_residual);
            out << ",\"reference_residual\":"; JsonOptional(out,r.reference_residual);
            out << ",\"relative_coordinate_difference\":";
            if (r.relative_coordinate_difference) JsonArray(out,*r.relative_coordinate_difference); else out << "null";
            out << ",\"weight_max_difference\":"; JsonOptional(out,r.weight_max_difference);
            out << ",\"floor_masks_equal\":";
            if (r.floor_masks_equal) out << *r.floor_masks_equal; else out << "null";
            out << '}';
        }
        out << ",\"iterations\":" << solve.diagnostics.iterations << ",\"squared_beta_change\":";
        JsonOptional(out, solve.diagnostics.squared_beta_change);
        out << ",\"relative_variance_change\":"; JsonOptional(out, solve.diagnostics.relative_variance_change);
        out << ",\"variance\":"; JsonOptional(out, solve.variance); out << '}';
    }
    out << "],\"offsets\":[";
    bool first{ true };
    for (const auto & [key, entry] : diagnostics.offsets)
    {
        if (!first) out << ','; first = false;
        out << "{\"key\":"; JsonArray(out, key);
        out << ",\"status\":\"" << JointOffsetSolveStatusText(entry.first) << "\",\"iterations\":" << entry.second.iterations;
        out << ",\"normalized_change\":"; JsonOptional(out, entry.second.normalized_change);
        out << ",\"robust_scale\":"; JsonOptional(out, entry.second.robust_scale); out << '}';
    }
    out << "]}";
}
void JsonRecovery(std::ostream & out, const RecoveryDiagnostics & diagnostic)
{
    out << "{\"attempted\":" << diagnostic.attempted << ",\"accepted\":" << diagnostic.accepted
        << ",\"reason\":\"" << diagnostic.reason << "\",\"operator_evaluations\":" << diagnostic.operator_evaluations;
    out << ",\"current_residual\":"; JsonOptional(out, diagnostic.current_residual);
    out << ",\"best_objective\":"; JsonOptional(out, diagnostic.best_objective);
    out << ",\"trials\":[";
    bool first{ true };
    for (const auto & trial : diagnostic.trials)
    {
        if (!first) out << ','; first = false;
        out << "{\"factor\":" << trial.factor << ",\"reason\":\"" << trial.reason << "\",\"residual\":";
        JsonOptional(out, trial.residual); out << ",\"objective\":"; JsonOptional(out, trial.objective); out << '}';
    }
    out << "]}";
}
void JsonBatch(std::ostream & out, const AuditBatch & batch)
{
    out << "\"stages\":{";
    for (std::size_t i=0;i<batch.stages.size();++i)
    {
        if (i) out << ','; const auto & count=batch.stages[i];
        out << '"' << kAuditStages[i] << "\":{\"total\":" << count.total << ",\"accepted\":" << count.accepted
            << ",\"rejected\":" << count.rejected << ",\"skipped\":" << count.skipped << '}';
    }
    out << "},\"anomalies\":{\"total\":" << batch.abnormal_count << ",\"shown\":" << batch.detail_count
        << ",\"omitted\":" << batch.abnormal_count-batch.detail_count << ",\"categories\":{";
    for (std::size_t i=1;i<batch.categories.size();++i)
    {
        std::size_t shown=0;
        for (std::size_t j=0;j<batch.detail_count;++j)
            if (static_cast<std::size_t>(batch.details[j].category)==i) ++shown;
        if (i!=1) out << ',';
        out << '"' << kAuditCategories[i] << "\":{\"total\":" << batch.categories[i]
            << ",\"shown\":" << shown << ",\"omitted\":" << batch.categories[i]-shown << '}';
    }
    out << "},\"details\":[";
    for (std::size_t i=0;i<batch.detail_count;++i)
    {
        const auto & e=batch.details[i]; if (i) out << ',';
        out << "{\"stage\":\"" << kAuditStages[static_cast<std::size_t>(e.stage)] << "\",\"first_atom\":" << e.first_atom
            << ",\"atom_count\":" << e.atom_count << ",\"trial\":" << e.trial << ",\"category\":\"" << kAuditCategories[static_cast<std::size_t>(e.category)]
            << "\",\"outcome\":\"" << e.outcome << "\",\"reason\":\"" << e.reason << "\",\"scope\":\"" << e.scope
            << "\",\"reference\":\"" << e.reference << "\",\"previous_checked\":" << e.previous_checked << ",\"best_checked\":" << e.best_checked
            << ",\"previous\":"; JsonObjective(out,e.previous); out << ",\"candidate\":"; JsonObjective(out,e.candidate);
        out << ",\"best\":"; JsonObjective(out,e.best); out << ",\"factor\":"; JsonOptional(out,e.factor);
        out << ",\"radius\":"; JsonOptional(out,e.radius); out << '}';
    }
    out << "]}";
}
}
void LogDecisionAuditStart(SecondStageObservationSession & session, const FitOptions & options) noexcept
{
    if (!session.Enabled()) return;
    try
    {
        std::ostringstream out; out.imbue(std::locale::classic()); out << std::boolalpha << "Second-stage audit: schema=2, payload={\"kind\":\"start\",\"version\":\"" << RHBM_GEM_AUDIT_VERSION
            << "\",\"settings\":{\"threads\":" << options.thread_size << ",\"exclude_hydrogen\":" << options.exclude_hydrogen
            << ",\"boundary_halo_depth\":" << options.second_stage_boundary_halo_depth << ",\"final_polish\":" << options.enable_second_stage_dependency_polish
            << ",\"final_polish_rounds\":" << options.second_stage_dependency_polish_max_iterations
            << ",\"failed_only_refinement\":" << options.enable_second_stage_failed_only_refinement << "}}";
        session.Write(out.str());
    }
    catch (...) { session.Disable(); }
}
void LogDecisionAuditIteration(SecondStageObservationSession & session, const IterationResult & result) noexcept
{
    const auto * data=session.Audit(); if (!data) return;
    try
    {
        std::ostringstream out; out.imbue(std::locale::classic()); out << std::boolalpha << "Second-stage audit: schema=2, payload={\"kind\":\"iteration\",\"attempt\":" << result.attempt_number
            << ",\"accepted_iterations\":" << result.accepted_iteration_count << ",\"accepted_clusters\":" << result.accepted_key_list.size()
            << ",\"rejected_clusters\":" << result.rejected_key_list.size() << ",\"objective_revision\":" << data->objective_revision
            << ",\"recovery_revision\":" << data->recovery_revision << ",\"background_revision\":" << data->background_revision
            << ",\"partition_revision\":" << data->partition_revision << ",\"objective_domain_changed\":" << result.objective_domain_changed
            << ",\"quarantine\":{\"entered\":" << data->entered << ",\"retried\":" << data->retried << ",\"released\":" << data->released
            << ",\"failed_retry\":" << data->failed_retry << "},\"score\":{\"scope\":\"global\",\"reference\":\"iteration_previous\",\"source\":\"" << data->score_source << "\",\"previous\":";
        JsonObjective(out,data->previous); out << ",\"candidate\":"; JsonObjective(out,data->candidate); out << ",\"best\":"; JsonObjective(out,data->best);
        out << "},\"selection_audit\":{";
        for (std::size_t i=0;i<2;++i)
        {
            if (i) out << ','; const auto & audit=data->selection[i];
            out << '"' << (i ? "after_rescue" : "ordinary") << "\":{\"executed\":" << audit.executed << ",\"result\":\"" << audit.result
                << "\",\"reason\":\"" << audit.reason << "\",\"evaluations\":" << audit.evaluations << ",\"removed_clusters\":" << audit.removed_clusters
                << ",\"previous\":"; JsonObjective(out,audit.previous); out << ",\"candidate\":"; JsonObjective(out,audit.candidate);
            out << ",\"best\":"; JsonObjective(out,audit.best); out << '}';
        }
        out << "},\"convergence\":"; JsonCertificate(out,data->convergence,data->convergence_reference);
        out << ",\"nominal_solves\":"; JsonNominalSolves(out, data->nominal);
        out << ",\"recovery\":"; JsonRecovery(out, data->recovery);
        out << ",\"stop_reason\":\"" << SecondStageStopReasonText(result.stop_reason) << "\","; JsonBatch(out,data->batch); out << '}';
        session.Write(out.str());
    }
    catch (...) { session.Disable(); }
}
void LogDecisionAuditTerminal(SecondStageObservationSession & session, std::string_view reason,
    std::string_view source, const BestAuditState & best, const PerformanceCounters * counters) noexcept
{
    const auto * data=session.Audit(); if (!data) return;
    try
    {
        std::ostringstream out; out.imbue(std::locale::classic()); out << std::boolalpha << "Second-stage audit: schema=2, payload={\"kind\":\"terminal\",\"stop_reason\":\"" << reason
            << "\",\"final_state_source\":\"" << source << "\",\"best_iteration\":";
        if (best) out << best->source_iteration; else out << "null";
        out << ",\"attempt\":" << data->attempt << ",\"objective_revision\":" << data->objective_revision
            << ",\"background_revision\":" << data->background_revision << ",\"partition_revision\":" << data->partition_revision
            << ",\"objective\":"; JsonObjective(out,data->final_objective);
        out << ",\"objective_scope\":\"global_frozen_domain\",\"objective_reference\":\"final_state_last_background\",\"final_polish\":{\"attempted\":" << data->polish_attempted
            << ",\"objective_accepted\":" << data->polish_accepted << ",\"operator_certified\":";
        if (data->polish_certificate) out << (data->polish_status == FinalPolishResidualSafetyStatus::AbsolutePassed); else out << "null";
        out << ",\"status\":\"" << GetFinalPolishResidualSafetyStatusText(data->polish_status) << "\",\"applied\":" << data->polish_applied << ",\"certificate\":";
        JsonCertificate(out,data->polish_certificate,"final_polish_candidate",false);
        out << "},\"final_certificate\":"; JsonCertificate(out, data->final_certificate, "persisted_state", false);
        out << ",\"final_nominal_solves\":"; JsonNominalSolves(out, data->final_nominal);
        out << ",\"recovery\":"; JsonRecovery(out, data->recovery);
        out << ",\"work_counters\":";
        if (counters) JsonArray(out,counters->AuditCounts()); else out << "null";
        out << ",\"elapsed_ms\":"; JsonNumber(out,session.ElapsedMilliseconds()); out << ','; JsonBatch(out,data->batch); out << '}';
        session.Write(out.str());
    }
    catch (...) { session.Disable(); }
}

void LogFinalStateCertificate(bool quiet, const std::optional<ConvergenceAssessment> & assessment, bool polish_applied,
    std::size_t attempts, std::size_t recovery_operators, std::size_t certificate_operators)
{
    if (quiet) return;
    std::ostringstream out;
    out.imbue(std::locale::classic()); out << std::boolalpha;
    out << "Second-stage final state: schema=1, payload={\"final_polish_applied\":" << polish_applied << ",\"certificate\":";
    JsonCertificate(out, assessment, "persisted_state", false);
    out << ",\"background_reference\":\"last_frozen_background\"";
    out << ",\"attempts\":" << attempts << ",\"recovery_operator_evaluations\":" << recovery_operators
        << ",\"certificate_operator_evaluations\":" << certificate_operators;
    out << '}';
    Logger::Log(LogLevel::Notice, out.str());
}

} // namespace rhbm_gem::core::detail
