#include <cstddef>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/PolishProvenance.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/Objective.hpp"
#include "core/detail/TerminalFailure.hpp"
#include "core/detail/IterationProcess.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/PerformanceCounters.hpp"
#include "core/detail/SecondStageInitialization.hpp"
#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/SecondStageContext.hpp"

#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rhbm_gem::core {

namespace {

using detail::FitState;

void ApplyFitState(
    ModelObject & model_object,
    const detail::SecondStageContext & context,
    const FitState & iteration_state)
{
    auto adjusted_sampling_entries_list{
        detail::BuildSecondStageAdjustedSamples(context, iteration_state)
    };

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        analysis.ApplyAtomLocalSecondStageResult(
            *context.at(i).atom,
            iteration_state.at(i),
            std::move(adjusted_sampling_entries_list.at(i)));
    }
}

void AppendOffsetSummary(std::ostringstream & stream, const FitState & state)
{
    std::size_t finite_count{ 0 };
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(state.size());
    for (const auto & result : state)
    {
        const auto offset{ result.mdpde.GetModel().GetOffset() };
        if (!std::isfinite(offset)) continue;
        finite_count++;
        absolute_offset_list.emplace_back(std::abs(offset));
    }
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
    if (!absolute_offset_list.empty())
    {
        median_absolute_offset = array_helper::ComputeMedian(absolute_offset_list);
        percentile_absolute_offset = array_helper::ComputePercentile(absolute_offset_list, 0.99);
        maximum_absolute_offset = *std::max_element(absolute_offset_list.begin(), absolute_offset_list.end());
    }
    stream << std::scientific << std::setprecision(2)
        << "; offsets finite = " << finite_count << " of " << state.size()
        << ", |C| median/p99/max = "
        << median_absolute_offset << "/"
        << percentile_absolute_offset << "/"
        << maximum_absolute_offset;
}

void AppendAuditSummary(std::ostringstream & stream, const detail::AuditedState & audited_state)
{
    const auto & objective{ audited_state.objective };
    stream << "; audit best source = ";
    if (audited_state.accepted_iteration.has_value())
    {
        stream << "accepted iteration " << *audited_state.accepted_iteration;
    }
    else
    {
        stream << "initial";
    }
    stream
        << std::scientific << std::setprecision(2)
        << ", fixed audit objective fit/tail-weighted/offset/total = "
        << objective.fit_range_residual_objective << "/"
        << objective.tail_validation_penalty << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.total_objective
        << ", tail raw/weight = "
        << objective.tail_validation_loss << "/"
        << detail::kTailValidationWeight;
}

void LogTerminalFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const detail::TerminalSummary & terminal_summary,
    const FitState & state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after " << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    detail::AppendTerminalSummary(warning_message, terminal_summary);
    AppendOffsetSummary(warning_message, state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void FinishWithNoActiveAtoms(
    ModelObject & model_object,
    const detail::SecondStageContext & context,
    const FitState & state,
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const detail::TerminalSummary & terminal_summary)
{
    ApplyFitState(model_object, context, state);
    if (terminal_summary.HasFailures())
    {
        LogTerminalFallback(quiet_mode, accepted_iteration_count, terminal_summary, state);
        return;
    }
    if (!quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Info,
            "Skip 2nd-stage local atom fitting because no atoms are selected.");
    }
}

void LogConverged(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const FitState & state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_stats.percentile_list.at(detail::kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_stats.percentile_list.at(detail::kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_stats.percentile_list.at(detail::kOffsetToPeakRatioChangeIndex);
    AppendOffsetSummary(message, state);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogMaximumIterations(
    bool quiet_mode,
    const detail::FinalStateSelection & selection,
    const detail::TerminalSummary & terminal_summary)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    detail::AppendTerminalSummary(warning_message, terminal_summary);
    if (selection.source == detail::FinalStateSource::BestAudit && selection.audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendAuditSummary(warning_message, *selection.audit_state);
    }
    else if (selection.source == detail::FinalStateSource::LatestValidated)
    {
        warning_message << "; applying latest validated state";
    }
    else
    {
        warning_message << "; no validated state is available";
    }
    AppendOffsetSummary(warning_message, selection.state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::string_view stop_reason,
    const detail::BestAuditState & best_audit_state,
    bool final_uses_polish,
    detail::FinalStateSource final_state_source)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message << "Second-stage local fitting summary: accepted_iterations="
        << accepted_iteration_count << ", best_iteration=";
    if (!best_audit_state.best.has_value())
    {
        message << "unavailable";
    }
    else if (best_audit_state.best->accepted_iteration.has_value())
    {
        message << *best_audit_state.best->accepted_iteration;
    }
    else
    {
        message << "initial";
    }
    message << ", stop_reason=" << stop_reason << ", best_audit_objective=";
    if (best_audit_state.best.has_value())
    {
        message << std::scientific << std::setprecision(8) << best_audit_state.best->objective.total_objective;
    }
    else
    {
        message << "unavailable";
    }
    message << ", final_uses_polish=";
    message << (final_uses_polish ? "yes" : "no");
    message << ", final_state_source="
        << detail::GetFinalStateSourceText(final_state_source) << ".";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace

bool RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    auto context{ detail::BuildSecondStageContext(model_object, options) };
    detail::StoreSecondStageNeighborCounts(model_object, context);
    const auto atom_size{ context.size() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto initial_state_build_result{ detail::BuildInitialFitState(context) };
    if (initial_state_build_result.failure != detail::SecondStageInitialStateBuildResult::Failure::None)
    {
        if (!options.quiet_mode)
        {
            const auto unselected_seed_failure{
                initial_state_build_result.failure == detail::SecondStageInitialStateBuildResult::Failure::UnselectedSeedUnavailable
            };
            Logger::Log(LogLevel::Warning,
                unselected_seed_failure ?
                    "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                    "is available for every unselected neighbor atom." :
                    "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                    "is available for every selected atom.");
            Logger::Log(LogLevel::Info,
                "Second-stage local fitting summary: accepted_iterations=0, "
                "best_iteration=unavailable, stop_reason=" +
                std::string(unselected_seed_failure ? "no-valid-unselected-neighbor-seed" : "no-valid-seed") +
                ", best_audit_objective=unavailable, final_uses_polish=unavailable, "
                "final_state_source=unavailable.");
        }
        return false;
    }
    detail::LogSecondStageSeedSelections(
        initial_state_build_result.selection_record_list,
        options.quiet_mode);
    detail::LogUnselectedSecondStageSeedSelections(
        context,
        initial_state_build_result.unselected_selection_record_list,
        options.quiet_mode);
    auto initial_state{ std::move(initial_state_build_result.state) };
    const auto graph_topology{
        detail::BuildSecondStageGraphTopology(context, initial_state, options.quiet_mode)
    };
    detail::LogGraphTopology(graph_topology, options.quiet_mode);
    auto iteration_state{
        detail::BuildIterationState(context, graph_topology, std::move(initial_state), options)
    };
    detail::PerformanceCounters performance_counters{
        options.quiet_mode,
        context,
        iteration_state.solver_workspace_by_key
    };
    if (iteration_state.best_audit_state.best.has_value())
    {
        performance_counters.RecordFullStateMaterialization();
    }
    detail::LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode);
    const auto progress_column_widths{ detail::BuildProgressColumnWidths(atom_size) };
    detail::LogProgressHeader(options.quiet_mode, progress_column_widths);

    for (std::size_t iter = 0; iter < detail::kMaximumIterations; iter++)
    {
        if (iteration_state.active_index_list.empty())
        {
            FinishWithNoActiveAtoms(
                model_object,
                context,
                iteration_state.previous_state,
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                iteration_state.terminal_failure_state.Summary());
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "terminal-isolation",
                iteration_state.best_audit_state,
                detail::UsesPolish(iteration_state.previous_polish_provenance),
                detail::FinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        auto iteration_result{
            detail::RunIteration(
                context,
                graph_topology,
                options,
                iter + 1,
                iteration_state,
                performance_counters)
        };
        if (iteration_result.objective_domain_changed)
        {
            detail::LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode, true);
        }
        detail::LogAcceptedBacktrackingDiagnostics(options.quiet_mode, iteration_result.diagnostics);
        if (iteration_result.outcome != detail::IterationOutcome::Accepted)
        {
            detail::LogRejectedClusterDiagnostics(
                options.quiet_mode,
                iteration_result.diagnostics.rejected_cluster_diagnostic_list);
            detail::LogIterationProgress(
                options.quiet_mode,
                progress_column_widths,
                iteration_result.progress);
            detail::LogAllRejectedResolution(
                options.quiet_mode,
                iteration_result.diagnostics.rejected_cluster_partition,
                iteration_result.diagnostics.trust_region_radius_update,
                *iteration_result.all_rejected_resolution);
            if (iteration_result.outcome == detail::IterationOutcome::Retry)
            {
                continue;
            }

            const auto final_state_selection{
                detail::SelectFinalState(
                    iteration_state.previous_state,
                    detail::UsesPolish(iteration_state.previous_polish_provenance),
                    iteration_state.best_audit_state.best)
            };
            ApplyFitState(model_object, context, final_state_selection.state);
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                detail::GetAllRejectedResolutionText(*iteration_result.all_rejected_resolution),
                iteration_state.best_audit_state,
                final_state_selection.uses_polish,
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        detail::LogRejectedClusterDiagnostics(
            options.quiet_mode,
            iteration_result.diagnostics.rejected_cluster_diagnostic_list);
        detail::LogIterationProgress(
            options.quiet_mode,
            progress_column_widths,
            iteration_result.progress);

        if (iteration_result.audit_patience_exhausted)
        {
            const auto final_state_selection{
                detail::SelectFinalState(
                    iteration_state.previous_state,
                    detail::UsesPolish(iteration_state.previous_polish_provenance),
                    iteration_state.best_audit_state.best)
            };
            ApplyFitState(model_object, context, final_state_selection.state);
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "audit-patience",
                iteration_state.best_audit_state,
                final_state_selection.uses_polish,
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iteration_result.converged)
        {
            ApplyFitState(model_object, context, iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.HasFailures())
            {
                LogTerminalFallback(
                    options.quiet_mode,
                    iteration_state.accepted_iteration_count,
                    iteration_state.terminal_failure_state.Summary(),
                    iteration_state.previous_state);
            }
            else
            {
                LogConverged(
                    options.quiet_mode,
                    iteration_state.accepted_iteration_count,
                    iteration_result.transformed_change_stats,
                    iteration_state.previous_state);
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "converged",
                iteration_state.best_audit_state,
                detail::UsesPolish(iteration_state.previous_polish_provenance),
                detail::FinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iter + 1 == detail::kMaximumIterations)
        {
            const auto final_state_selection{
                detail::SelectFinalState(
                    iteration_state.previous_state,
                    detail::UsesPolish(iteration_state.previous_polish_provenance),
                    iteration_state.best_audit_state.best)
            };
            ApplyFitState(model_object, context, final_state_selection.state);
            LogMaximumIterations(
                options.quiet_mode,
                final_state_selection,
                iteration_state.terminal_failure_state.Summary());
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "maximum-iterations",
                iteration_state.best_audit_state,
                final_state_selection.uses_polish,
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }
    }
    return false;
}
} // namespace rhbm_gem::core
