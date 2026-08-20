#include <cstddef>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/FitStateView.hpp"
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

#include <cmath>
#include <iomanip>
#include <ranges>
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
    const auto model_snapshot{
        detail::BuildSecondStageModelSnapshot(context, iteration_state)
    };

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto adjusted_sampling_entries{
            detail::BuildSecondStageAdjustedSamples(context.at(i), model_snapshot)
        };
        analysis.ApplyAtomLocalSecondStageResult(
            *context.at(i).atom,
            iteration_state.at(i),
            std::move(adjusted_sampling_entries));
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
        maximum_absolute_offset = std::ranges::max(absolute_offset_list);
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
        << detail::kTailValidationWeight;
}

void LogTerminalFallback(
    bool quiet_mode,
    const detail::IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after "
        << iteration_state.accepted_iteration_count
        << " accepted iterations with last validated states retained";
    detail::AppendTerminalSummary(
        warning_message,
        iteration_state.terminal_failure_state.Summary());
    AppendOffsetSummary(warning_message, iteration_state.previous_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogConverged(
    bool quiet_mode,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const detail::IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << iteration_state.accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_stats.percentile_list.at(detail::kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_stats.percentile_list.at(detail::kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_stats.percentile_list.at(detail::kOffsetToPeakRatioChangeIndex);
    AppendOffsetSummary(message, iteration_state.previous_state);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogMaximumIterations(bool quiet_mode, const detail::IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    detail::AppendTerminalSummary(
        warning_message,
        iteration_state.terminal_failure_state.Summary());
    const auto * audit_state{
        iteration_state.best_audit_state.has_value() ? &*iteration_state.best_audit_state : nullptr
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
    AppendOffsetSummary(
        warning_message,
        audit_state != nullptr ? audit_state->state : iteration_state.previous_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageSummary(
    bool quiet_mode,
    const detail::IterationState & iteration_state,
    std::string_view stop_reason,
    bool final_uses_best_audit)
{
    if (quiet_mode) return;

    const auto & best_audit_state{ iteration_state.best_audit_state };
    const auto final_uses_polish{
        final_uses_best_audit && best_audit_state.has_value() ?
            best_audit_state->uses_polish :
            detail::UsesPolish(iteration_state.previous_polish_provenance)
    };
    Logger::FinishProgressLine();
    std::ostringstream message;
    message << "Second-stage local fitting summary: accepted_iterations="
        << iteration_state.accepted_iteration_count << ", best_iteration=";
    if (!best_audit_state.has_value())
    {
        message << "unavailable";
    }
    else if (best_audit_state->source_iteration != 0)
    {
        message << best_audit_state->source_iteration;
    }
    else
    {
        message << "initial";
    }
    message << ", stop_reason=" << stop_reason << ", best_audit_objective=";
    if (best_audit_state.has_value())
    {
        message << std::scientific << std::setprecision(2)
            << best_audit_state->objective.GetTotalObjective();
    }
    else
    {
        message << "unavailable";
    }
    message << ", final_uses_polish=";
    message << (final_uses_polish ? "yes" : "no");
    message << ", final_state_source="
        << (final_uses_best_audit ? "best-audit" : "latest-validated") << ".";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace

bool RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    auto context{ detail::BuildSecondStageContext(model_object, options) };
    detail::StoreSecondStageNeighborCounts(model_object, context);
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    FitState initial_state;
    {
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
            initial_state_build_result.unselected_selection_record_list,
            options.quiet_mode);
        initial_state = std::move(initial_state_build_result.state);
    }
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
    if (iteration_state.best_audit_state.has_value())
    {
        performance_counters.RecordFullStateMaterialization();
    }
    detail::LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode);
    const auto progress_column_widths{ detail::BuildProgressColumnWidths(context.size()) };
    detail::LogProgressHeader(options.quiet_mode, progress_column_widths);

    std::string_view final_stop_reason;
    bool maximum_iterations_reached{ false };
    for (std::size_t iter = 0; iter < detail::kMaximumIterations; iter++)
    {
        if (iteration_state.active_index_list.empty())
        {
            ApplyFitState(model_object, context, iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.Summary().HasFailures())
            {
                LogTerminalFallback(options.quiet_mode, iteration_state);
            }
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Skip 2nd-stage local atom fitting because no atoms are selected.");
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state,
                "terminal-isolation",
                false);
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
        if (iteration_result.all_rejected_resolution.has_value())
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
                iteration_result.diagnostics.trust_region_update,
                *iteration_result.all_rejected_resolution);
            if (*iteration_result.all_rejected_resolution == detail::AllRejectedResolution::Retry)
            {
                continue;
            }

            final_stop_reason = detail::GetAllRejectedResolutionText(*iteration_result.all_rejected_resolution);
            break;
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
            final_stop_reason = "audit-patience";
            break;
        }

        if (iteration_result.converged)
        {
            ApplyFitState(model_object, context, iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.HasFailures())
            {
                LogTerminalFallback(options.quiet_mode, iteration_state);
            }
            else
            {
                LogConverged(
                    options.quiet_mode,
                    iteration_result.transformed_change_stats,
                    iteration_state);
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state,
                "converged",
                false);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iter + 1 == detail::kMaximumIterations)
        {
            final_stop_reason = "maximum-iterations";
            maximum_iterations_reached = true;
            break;
        }
    }

    if (!final_stop_reason.empty())
    {
        const auto * audit_state{
            iteration_state.best_audit_state.has_value() ? &*iteration_state.best_audit_state : nullptr
        };
        const auto & final_state{
            audit_state != nullptr ? audit_state->state : iteration_state.previous_state
        };
        ApplyFitState(model_object, context, final_state);
        if (maximum_iterations_reached)
        {
            LogMaximumIterations(options.quiet_mode, iteration_state);
        }
        LogSecondStageSummary(
            options.quiet_mode,
            iteration_state,
            final_stop_reason,
            audit_state != nullptr);
        RunGroupPotentialFitting(model_object, options, FittingStage::Second);
        return true;
    }
    return false;
}
} // namespace rhbm_gem::core
