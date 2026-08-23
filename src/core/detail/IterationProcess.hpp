#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "core/detail/JointOffset.hpp"
#include "core/detail/ClusterHealth.hpp"
#include "core/detail/ClusterSolverWorkspace.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/TerminalFailure.hpp"
#include "core/detail/Objective.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/SuspiciousUpdate.hpp"
#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/SecondStageContext.hpp"
#include "core/detail/TransformedChange.hpp"
#include "core/detail/TrustRegion.hpp"
#include "core/detail/PerformanceCounters.hpp"

namespace rhbm_gem::core::detail {

constexpr std::size_t kMaximumIterations{ 100 };
constexpr std::size_t kAuditPatience{ 3 };

struct IterationDiagnostics
{
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::size_t combined_backtracking_trial_count{ 0 };
    std::optional<double> combined_backtracking_factor{};
    bool combined_backtracking_exhausted{ false };
    TrustRegionIterationUpdate trust_region_update{};
};

struct IterationState
{
    FitState previous_state{};
    PolishProvenance previous_polish_provenance{};
    SuspiciousUpdateMask rollback_atom_mask{};
    std::vector<std::size_t> active_index_list{};
    CouplingGraphPartition graph_partition{};
    ClusterSolverWorkspaceMap solver_workspace_by_key{};
    ObjectiveDomain objective_domain{};
    BestAuditState best_audit_state{};
    TerminalFailureState terminal_failure_state{};
    ClusterObjectiveStateMap cluster_objective_state{};
    TrustRegionStateSet trust_region_state{};
    std::vector<ClusterKey> unchanged_state_exhausted_key_list{};
    std::size_t accepted_iteration_count{ 0 };
    std::size_t audit_patience_count{ 0 };
};

struct IterationProgress
{
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t terminal_atom_count{ 0 };
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rejected_cluster_count{ 0 };
    PolishProgress polish_progress{};
    std::size_t suspicious_atom_count{ 0 };
    std::optional<double> accepted_maximum_transformed_change{};
    double raw_maximum_transformed_change{ 0.0 };
};

using ProgressColumnWidths = std::array<std::size_t, 6>;

struct IterationResult
{
    IterationDiagnostics diagnostics{};
    IterationProgress progress{};
    std::optional<AllRejectedResolution> all_rejected_resolution{};
    bool objective_domain_changed{ false };
    bool converged{ false };
    bool audit_patience_exhausted{ false };
    algorithm::ParameterChangeStats transformed_change_stats{};
};

inline void AppendObjectiveBreakdown(
    std::ostringstream & stream,
    const std::optional<ObjectiveBreakdown> & breakdown)
{
    if (!breakdown.has_value())
    {
        stream << "unavailable";
        return;
    }
    stream
        << breakdown->fit_range_residual_objective << "/"
        << breakdown->GetTailValidationPenalty() << "/"
        << breakdown->offset_plausibility_penalty << "/"
        << breakdown->GetTotalObjective();
}

inline std::string_view GetPreObjectiveFailureReasonText(PreObjectiveFailureReason reason)
{
    switch (reason)
    {
    case PreObjectiveFailureReason::None:
        return "none";
    case PreObjectiveFailureReason::InvalidModel:
        return "invalid-model";
    case PreObjectiveFailureReason::PreviousSharedOffsetProjectionOutsideTrustRegion:
        return "previous-shared-offset-projection-outside-trust-region";
    case PreObjectiveFailureReason::NoCandidateWithinTrustRegion:
        return "no-candidate-within-trust-region";
    }
    return "unknown";
}

inline void LogRejectedClusterDiagnostics(
    bool quiet_mode,
    const std::vector<ClusterCandidateDiagnostic> & diagnostic_list)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug || diagnostic_list.empty())
    {
        return;
    }

    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : diagnostic_list)
    {
        std::ostringstream header;
        header
            << "Rejected local fitting cluster diagnostics: atoms = " << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/" << cluster_diagnostic.key.back()
            << ", breakdown order = fit/tail-weighted/offset/total";
        Logger::Log(LogLevel::Debug, header.str());

        const auto & diagnostic{ cluster_diagnostic.attempt };
        std::ostringstream message;
        message << std::scientific << std::setprecision(2)
            << "  fixed-point effective damping = " << diagnostic.effective_damping
            << ", trust radius/step norm = " << diagnostic.trust_region_radius << "/";

        if (diagnostic.pre_objective_failure_reason !=
            PreObjectiveFailureReason::None)
        {
            if (diagnostic.pre_objective_attempted_step_norm.has_value())
            {
                message << *diagnostic.pre_objective_attempted_step_norm;
            }
            else
            {
                message << "unavailable";
            }
            message
                << ", status = "
                << GetPreObjectiveFailureReasonText(diagnostic.pre_objective_failure_reason)
                << ", objective = not-evaluated";
            Logger::Log(LogLevel::Debug, message.str());
            continue;
        }

        message << diagnostic.trust_region_step_norm;

        if (diagnostic.is_invalid_model)
        {
            message << ", status = invalid-model";
            Logger::Log(LogLevel::Debug, message.str());
            continue;
        }

        message << ", fit/tail scales = ";
        if (diagnostic.scale.has_value())
        {
            message << diagnostic.scale->fit;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.scale.has_value() && diagnostic.tail_sample_count > 0)
        {
            message << diagnostic.scale->tail;
        }
        else
        {
            message << "empty";
        }
        message
            << ", fit/tail samples = "
            << diagnostic.fit_sample_count << "/" << diagnostic.tail_sample_count
            << ", tail raw/weight = ";
        if (diagnostic.candidate_objective.has_value())
        {
            message << diagnostic.candidate_objective->tail_validation_loss;
        }
        else
        {
            message << "unavailable";
        }
        message << "/" << kTailValidationWeight;
        message << ", candidate = ";
        AppendObjectiveBreakdown(message, diagnostic.candidate_objective);
        message << ", previous = ";
        AppendObjectiveBreakdown(message, diagnostic.previous_objective);
        message << ", best = ";
        AppendObjectiveBreakdown(message, diagnostic.best_objective);
        message << ", rejected-by = ";
        if (!diagnostic.candidate_objective.has_value())
        {
            message << "objective-unavailable";
        }
        else if (diagnostic.rejected_by_previous && diagnostic.rejected_by_best)
        {
            message << "previous+best";
        }
        else if (diagnostic.rejected_by_previous)
        {
            message << "previous";
        }
        else if (diagnostic.rejected_by_best)
        {
            message << "best";
        }
        else
        {
            message << "none";
        }
        message
            << ", backtracking trials/factor/exhausted = "
            << diagnostic.backtracking_trial_count << "/";
        if (diagnostic.accepted_backtracking_factor.has_value())
        {
            message << *diagnostic.accepted_backtracking_factor;
        }
        else
        {
            message << "-";
        }
        message << "/" << (diagnostic.backtracking_exhausted ? "yes" : "no");
        Logger::Log(LogLevel::Debug, message.str());
    }
}

inline std::string_view GetAllRejectedResolutionText(AllRejectedResolution resolution)
{
    switch (resolution)
    {
    case AllRejectedResolution::Retry:
        return "retry";
    case AllRejectedResolution::MaximumIterations:
        return "maximum-iterations";
    case AllRejectedResolution::BacktrackingExhausted:
        return "all-rejected-backtracking-exhausted";
    case AllRejectedResolution::MinimumRadius:
        return "all-rejected-minimum-radius";
    case AllRejectedResolution::NoRetryProgress:
        return "all-rejected-no-retry-progress";
    }
    return "all-rejected-no-retry-progress";
}

inline void LogAllRejectedResolution(
    bool quiet_mode,
    const TrustRegionIterationUpdate & trust_region_update,
    AllRejectedResolution resolution)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "All-rejected local fitting resolution: outcome = "
        << GetAllRejectedResolutionText(resolution)
        << ", exhausted/retryable/radius-changed/radius-saturated = "
        << trust_region_update.rejected_cluster_partition.exhausted_key_list.size() << "/"
        << trust_region_update.rejected_cluster_partition.retryable_key_list.size() << "/"
        << trust_region_update.radius_update.changed_key_list.size() << "/"
        << trust_region_update.radius_update.saturated_key_list.size() << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

inline void LogAcceptedBacktrackingDiagnostics(
    bool quiet_mode,
    const IterationDiagnostics & diagnostics)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }
    const auto has_local_backtracking{
        std::any_of(
            diagnostics.accepted_cluster_diagnostic_list.begin(),
            diagnostics.accepted_cluster_diagnostic_list.end(),
            [](const ClusterCandidateDiagnostic & diagnostic)
            {
                return diagnostic.attempt.backtracking_trial_count > 1;
            })
    };
    if (!has_local_backtracking && diagnostics.combined_backtracking_trial_count <= 1)
    {
        return;
    }
    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : diagnostics.accepted_cluster_diagnostic_list)
    {
        const auto & diagnostic{ cluster_diagnostic.attempt };
        if (diagnostic.backtracking_trial_count <= 1) continue;
        std::ostringstream message;
        message
            << "Accepted local fitting objective backtracking: atoms = " << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/" << cluster_diagnostic.key.back()
            << ", trials/factor = " << diagnostic.backtracking_trial_count << "/";
        if (diagnostic.accepted_backtracking_factor.has_value())
        {
            message << *diagnostic.accepted_backtracking_factor;
        }
        else
        {
            message << "-";
        }
        message << ", fixed fit/tail scales = ";
        if (diagnostic.scale.has_value())
        {
            message << diagnostic.scale->fit;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.scale.has_value() && diagnostic.tail_sample_count > 0)
        {
            message << diagnostic.scale->tail;
        }
        else
        {
            message << "empty";
        }
        message << ".";
        Logger::Log(LogLevel::Debug, message.str());
    }
    if (diagnostics.combined_backtracking_trial_count <= 1) return;
    std::ostringstream message;
    message << "Combined-objective backtracking: trials/factor/exhausted = "
        << diagnostics.combined_backtracking_trial_count << "/";
    if (diagnostics.combined_backtracking_factor.has_value())
    {
        message << *diagnostics.combined_backtracking_factor;
    }
    else
    {
        message << "-";
    }
    message << "/" << (diagnostics.combined_backtracking_exhausted ? "yes" : "no") << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

inline std::string FormatProgressMaximum(double value)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

inline constexpr std::array<std::string_view, 6> kProgressHeaderList
{
    "Try/Acc",
    "Atom A/T",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/R"
};

inline std::string FormatProgressRow(
    const ProgressColumnWidths & column_widths,
    const auto & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream << std::left << std::setw(static_cast<int>(column_widths.at(i))) << cell_list.at(i);
    }
    return stream.str();
}

inline ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_size)
{
    const auto maximum_iteration_text{ std::to_string(kMaximumIterations) };
    const auto maximum_atom_text{ std::to_string(atom_size) };
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

inline void LogProgressHeader(bool quiet_mode, const ProgressColumnWidths & column_widths)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Info, FormatProgressRow(column_widths, kProgressHeaderList));
}

inline void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationProgress & progress)
{
    if (quiet_mode) return;

    const std::array<std::string, 6> cell_list{
        std::to_string(progress.attempt_number) + "/" +
            std::to_string(progress.accepted_iteration_count),
        std::to_string(progress.active_atom_count) + "/" +
            std::to_string(progress.terminal_atom_count),
        std::to_string(progress.accepted_cluster_count) + "/" +
            std::to_string(progress.rejected_cluster_count),
        std::to_string(progress.polish_progress.eligible_count) + "/" +
            std::to_string(progress.polish_progress.accepted_count) + "/" +
            std::to_string(progress.polish_progress.rejected_count) + "/" +
        std::to_string(progress.polish_progress.skipped_count),
        std::to_string(progress.suspicious_atom_count),
        (progress.accepted_maximum_transformed_change.has_value() ?
            FormatProgressMaximum(*progress.accepted_maximum_transformed_change) :
            std::string{ "-" }) + "/" +
            FormatProgressMaximum(progress.raw_maximum_transformed_change)
    };
    Logger::ProgressLine(FormatProgressRow(column_widths, cell_list));
}

struct RawIterationResult
{
    FitState state{};
    SuspiciousUpdateMask rollback_atom_mask{};
    ClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_stationarity_eligible{ false };
};

inline std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const AtomContext & atom_context,
    const LocalGaussianResult & previous_result,
    const GaussianModel3D & offset_model,
    const std::vector<double> & adjusted_response_list,
    const FitOptions & options)
{
    auto adjusted_sampling_entries{
        BuildSecondStageAdjustedSamples(atom_context, adjusted_response_list)
    };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto previous_baseline{
        BuildPreviousSuspiciousProfileBaseline(adjusted_sampling_entries, previous_model, options)
    };
    const auto is_candidate_acceptable = [&](
        const GaussianModel3D & model,
        SuspiciousUpdateMode mode)
    {
        const auto reason{ EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                model,
                options,
                previous_baseline,
                mode) };
        return reason == SuspiciousGaussianReason::None;
    };
    try
    {
        auto candidate_result{
            EstimateLocalGaussianPrepared(
                atom_context.refit_design_template,
                adjusted_response_list,
                atom_context.alpha_r,
                options,
                offset_model)
        };
        const auto & candidate_model{ candidate_result.mdpde.GetModel() };
        if (is_candidate_acceptable(candidate_model, SuspiciousUpdateMode::PostRefit))
        {
            const auto is_stationarity_eligible{
                candidate_result.fit_result.has_value() &&
                IsLocalRefitStatusStationarityEligible(candidate_result.fit_result->status)
            };
            return LocalAtomRefitResult{
                std::move(candidate_result),
                is_stationarity_eligible
            };
        }
    }
    catch (const std::exception &)
    {
    }

    auto result{ previous_result };
    result.ols = GaussianModel3DWithUncertainty{
        result.ols.GetModel().WithOffset(offset_model.GetOffset()),
        result.ols.GetStandardDeviationModel()
    };
    result.mdpde = GaussianModel3DWithUncertainty{
        result.mdpde.GetModel().WithOffset(offset_model.GetOffset()),
        result.mdpde.GetStandardDeviationModel()
    };
    if (!is_candidate_acceptable(result.mdpde.GetModel(), SuspiciousUpdateMode::OffsetOnly))
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{ std::move(result), false };
}

inline RawIterationResult RunRawIteration(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list,
    ClusterSolverWorkspaceMap & solver_workspace_by_key)
{
    auto current_model_snapshot{
        BuildSecondStageModelSnapshot(context, previous_state)
    };
    const auto is_debug_logging_enabled{ Logger::GetLogLevel() >= LogLevel::Debug };
    const auto log_debug_diagnostics{ !options.quiet_mode && is_debug_logging_enabled };
    std::vector<JointOffsetSolveResult> joint_offset_result_list(cluster_key_list.size());
    std::vector<std::exception_ptr> joint_offset_exception_list(cluster_key_list.size());
    const auto solve_joint_offset = [&](std::size_t cluster_position)
    {
        try
        {
            joint_offset_result_list.at(cluster_position) = EstimateJointOffsets(
                context,
                cluster_key_list.at(cluster_position),
                current_model_snapshot,
                ridge_multiplier_list,
                solver_workspace_by_key.at(
                    cluster_key_list.at(cluster_position)).joint_offset,
                log_debug_diagnostics);
        }
        catch (...)
        {
            joint_offset_exception_list.at(cluster_position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    const bool parallel_joint_offsets{
        !is_debug_logging_enabled &&
        options.thread_size > 1 &&
        cluster_key_list.size() > 1
    };
    if (parallel_joint_offsets)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
        for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
        {
            solve_joint_offset(cluster_position);
        }
    }
    else
#endif
    {
        for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
        {
            solve_joint_offset(cluster_position);
        }
    }
    for (const auto & exception : joint_offset_exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    ClusterHealthMap health_by_key;
    for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        const auto & result{ joint_offset_result_list.at(cluster_position) };
        for (std::size_t i = 0; i < key.size(); i++)
        {
            const auto atom_index{ key.at(i) };
            current_model_snapshot.selected.at(atom_index) =
                GetFitModel(current_model_snapshot.selected, atom_index)
                    .WithOffset(result.offset(static_cast<Eigen::Index>(i)));
        }
        health_by_key.emplace(key, ClusterHealth{ result.status });
    }

    auto iteration_state{ previous_state };
    SuspiciousUpdateMask rollback_atom_mask(context.size(), 0);
    std::vector<std::size_t> group_id_by_atom_index;
    group_id_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_id_by_atom_index.emplace_back(atom_context.group_id);
    }
    for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        std::vector<std::size_t> group_id_by_position;
        SuspiciousUpdateMask suspicious_seed_mask(key.size(), 0);
        group_id_by_position.reserve(key.size());
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            group_id_by_position.emplace_back(group_id_by_atom_index.at(atom_index));
            if (EvaluateSuspiciousOffsetUpdate(
                    context.at(atom_index).raw_sampling_entries,
                    previous_state.at(atom_index).mdpde.GetModel(),
                    GetFitModel(current_model_snapshot.selected, atom_index),
                    options) != SuspiciousGaussianReason::None)
            {
                suspicious_seed_mask.at(position) = 1;
            }
        }
        const auto cluster_rollback_mask{
            ExpandSuspiciousSharedOffsetGroups(group_id_by_position, suspicious_seed_mask)
        };
        for (std::size_t position = 0; position < key.size(); position++)
        {
            if (cluster_rollback_mask.at(position) == 0) continue;
            rollback_atom_mask.at(key.at(position)) = 1;
        }
    }
    for (std::size_t atom_index = 0; atom_index < rollback_atom_mask.size(); atom_index++)
    {
        if (rollback_atom_mask.at(atom_index) == 0) continue;
        current_model_snapshot.selected.at(atom_index) = previous_state.at(atom_index).mdpde.GetModel();
    }

    FittedGaussianSnapshot refit_model_snapshot{
        BuildGroupMedianModelList(
            group_id_by_atom_index,
            current_model_snapshot.selected)
    };
    const auto refit_model_bundle{
        BuildSecondStageModelSnapshot(context, std::move(refit_model_snapshot))
    };
    const auto refit_response_cache{
        BuildSecondStageAdjustedResponseCache(context, refit_model_bundle)
    };
    std::vector<std::size_t> refit_atom_index_list;
    for (const auto & key : cluster_key_list)
    {
        for (const auto atom_index : key)
        {
            if (rollback_atom_mask.at(atom_index) == 0)
            {
                refit_atom_index_list.emplace_back(atom_index);
            }
        }
    }
    std::vector<std::optional<LocalAtomRefitResult>> refit_result_list(refit_atom_index_list.size());
    std::vector<std::exception_ptr> refit_exception_list(refit_atom_index_list.size());
#ifdef USE_OPENMP
    const bool parallel_refits{
        !is_debug_logging_enabled &&
        options.thread_size > 1 &&
        refit_atom_index_list.size() > 1
    };
#else
    const bool parallel_refits{ false };
#endif
    FitOptions refit_options{ options };
    if (parallel_refits)
    {
        refit_options.thread_size = 1;
    }
    const auto run_refit = [&](std::size_t refit_position)
    {
        const auto atom_index{ refit_atom_index_list.at(refit_position) };
        try
        {
            refit_result_list.at(refit_position) =
                FitAtomWithJointOffsetFallback(
                    context.at(atom_index),
                    previous_state.at(atom_index),
                    GetFitModel(refit_model_bundle.selected, atom_index),
                    refit_response_cache.at(atom_index),
                    refit_options);
        }
        catch (...)
        {
            refit_exception_list.at(refit_position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (parallel_refits)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
        for (std::size_t refit_position = 0; refit_position < refit_atom_index_list.size(); refit_position++)
        {
            run_refit(refit_position);
        }
    }
    else
#endif
    {
        for (std::size_t refit_position = 0; refit_position < refit_atom_index_list.size(); refit_position++)
        {
            run_refit(refit_position);
        }
    }
    for (const auto & exception : refit_exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    std::size_t refit_position{ 0 };
    for (const auto & key : cluster_key_list)
    {
        auto & health{ health_by_key.at(key) };
        bool has_post_refit_suspicious_atom{ false };
        for (const auto atom_index : key)
        {
            if (rollback_atom_mask.at(atom_index) != 0) continue;

            auto refit_result{ std::move(refit_result_list.at(refit_position++)) };
            if (!refit_result.has_value())
            {
                health.is_refit_stationarity_eligible = false;
                has_post_refit_suspicious_atom = true;
                continue;
            }
            if (!refit_result->is_stationarity_eligible)
            {
                health.is_refit_stationarity_eligible = false;
            }
            iteration_state.at(atom_index) = std::move(refit_result->result);
        }
        if (has_post_refit_suspicious_atom)
        {
            for (const auto atom_index : key)
            {
                rollback_atom_mask.at(atom_index) = 1;
            }
        }
    }
    for (std::size_t atom_index = 0; atom_index < rollback_atom_mask.size(); atom_index++)
    {
        if (rollback_atom_mask.at(atom_index) == 0) continue;
        iteration_state.at(atom_index) = previous_state.at(atom_index);
    }

    return RawIterationResult{
        std::move(iteration_state),
        std::move(rollback_atom_mask),
        std::move(health_by_key)
    };
}

inline IterationState BuildIterationState(
    const SecondStageContext & context,
    const GraphTopology & graph_topology,
    FitState initial_state,
    const FitOptions & options)
{
    IterationState iteration_state;
    iteration_state.previous_state = std::move(initial_state);
    iteration_state.previous_polish_provenance.assign(context.size(), 0);
    iteration_state.rollback_atom_mask.assign(context.size(), 0);
    iteration_state.terminal_failure_state = TerminalFailureState(context.size());
    iteration_state.active_index_list = iteration_state.terminal_failure_state.BuildEligibleActiveIndexList();
    iteration_state.graph_partition = BuildGraphPartition(graph_topology, iteration_state.active_index_list);
    const auto cluster_key_list{ BuildGraphClusterKeyList(iteration_state.graph_partition) };
    ResetClusterSolverWorkspace(
        cluster_key_list,
        iteration_state.solver_workspace_by_key);
    const auto initial_model_snapshot{
        BuildSecondStageModelSnapshot(context, iteration_state.previous_state)
    };
    iteration_state.objective_domain = BuildObjectiveDomain(
        context,
        initial_model_snapshot,
        cluster_key_list,
        options.distance_min,
        options.distance_max);
    const auto initial_audit_objective{
        EvaluateAuditObjective(
            iteration_state.objective_domain,
            SnapshotResidualEvaluator{ context, initial_model_snapshot })
    };
    if (initial_audit_objective.has_value())
    {
        TryUpdateBestAuditState(
            iteration_state.previous_state,
            UsesPolish(iteration_state.previous_polish_provenance),
            0,
            *initial_audit_objective,
            iteration_state.best_audit_state);
    }
    return iteration_state;
}

inline IterationResult RunIteration(
    const SecondStageContext & context,
    const GraphTopology & graph_topology,
    const FitOptions & options,
    std::size_t attempt_number,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    const auto & previous_state{ iteration_state.previous_state };
    const auto & active_index_list{ iteration_state.active_index_list };
    const auto & graph_partition{ iteration_state.graph_partition };
    const auto cluster_key_list{ BuildGraphClusterKeyList(graph_partition) };
    const auto & objective_domain{ iteration_state.objective_domain };

    const auto residual_baseline{
        BuildResidualBaseline(context, previous_state)
    };
    performance_counters.RecordGaussianCacheMisses();

    const auto previous_objective_by_key{
        BuildObjectiveByKey(graph_partition, objective_domain, residual_baseline)
    };
    ReconcileClusterObjectiveState(previous_objective_by_key, iteration_state.cluster_objective_state);
    iteration_state.trust_region_state.Reconcile(cluster_key_list);

    const auto joint_offset_ridge_multiplier_list{
        BuildSuspiciousJointOffsetRidgeMultiplierList(iteration_state.rollback_atom_mask)
    };

    const auto iteration_phase_start{
        performance_counters.StartIterationPhase()
    };
    auto raw_iteration_result{
        RunRawIteration(
            context,
            cluster_key_list,
            previous_state,
            options,
            joint_offset_ridge_multiplier_list,
            iteration_state.solver_workspace_by_key)
    };
    performance_counters.FinishIterationPhase(iteration_phase_start);
    performance_counters.RecordGaussianCacheHits();

    const auto iteration_suspicious_atom_count{
        CountSuspiciousAtoms(raw_iteration_result.rollback_atom_mask)
    };
    const auto has_suspicious_offset_fallback{ iteration_suspicious_atom_count > 0 };
    const auto is_stationarity_eligible{
        AreClustersStationarityEligible(raw_iteration_result.health_by_key)
    };
    const auto & raw_state{ raw_iteration_result.state };
    const auto raw_fixed_point_change_summary{
        SummarizeTransformedChanges(raw_state, previous_state, active_index_list)
    };
    auto working_cluster_objective_state{
        iteration_state.cluster_objective_state
    };
    const CandidateSelectionInputs candidate_inputs{
        .context = context,
        .residual_baseline = residual_baseline,
        .partition = graph_partition,
        .health_by_key = raw_iteration_result.health_by_key,
        .previous_state = previous_state,
        .previous_polish_provenance = iteration_state.previous_polish_provenance,
        .raw_state = raw_state,
        .rollback_atom_mask = raw_iteration_result.rollback_atom_mask,
        .ridge_multiplier_list = joint_offset_ridge_multiplier_list,
        .unchanged_state_exhausted_key_list =
            std::span<const ClusterKey>{ iteration_state.unchanged_state_exhausted_key_list },
        .objective_domain = objective_domain,
        .previous_objective_by_key = previous_objective_by_key,
        .cluster_objective_state = working_cluster_objective_state,
        .trust_region_state = iteration_state.trust_region_state,
        .solver_workspace_by_key = iteration_state.solver_workspace_by_key,
        .thread_size = options.thread_size,
        .performance_counters = performance_counters
    };
    const auto candidate_phase_start{ performance_counters.StartCandidatePhase() };
    auto selection{ SelectClusterCandidates(candidate_inputs) };
    performance_counters.FinishCandidatePhase(candidate_phase_start);
    performance_counters.RecordFullStateMaterialization();

    const auto * best_audit_objective{
        iteration_state.best_audit_state.has_value() ? &iteration_state.best_audit_state->objective : nullptr
    };
    const auto combined_check{
        EvaluateCombinedCandidateObjective(
            context,
            residual_baseline,
            graph_partition,
            previous_state,
            selection.assembled_state,
            selection.accepted_key_list,
            objective_domain,
            best_audit_objective,
            performance_counters)
    };
    selection.combined_backtracking_objective = combined_check.candidate_objective;
    auto combined_objective_accepted{ combined_check.accepted };
    if (!combined_objective_accepted)
    {
        combined_objective_accepted =
            TryBacktrackCombinedCandidate(
                candidate_inputs,
                combined_check.previous_objective.has_value() ?
                    &*combined_check.previous_objective : nullptr,
                best_audit_objective,
                iteration_state.cluster_objective_state,
                selection);
    }
    if (!combined_objective_accepted)
    {
        if (selection.combined_backtracking_exhausted)
        {
            for (const auto & key : selection.accepted_key_list)
            {
                if (std::ranges::find(
                        selection.backtracking_exhausted_key_list,
                        key) == selection.backtracking_exhausted_key_list.end())
                {
                    selection.backtracking_exhausted_key_list.emplace_back(key);
                }
            }
        }
        RejectCombinedCandidate(
            previous_state,
            iteration_state.previous_polish_provenance,
            selection);
    }
    else
    {
        iteration_state.cluster_objective_state = std::move(working_cluster_objective_state);
    }

    auto assembled_state{ std::move(selection.assembled_state) };
    auto assembled_polish_provenance{ std::move(selection.assembled_polish_provenance) };
    const auto has_new_terminal_failures{
        iteration_state.terminal_failure_state.IsolatePersistentFailures(
            selection.accepted_key_list,
            raw_iteration_result.rollback_atom_mask,
            raw_iteration_result.health_by_key,
            assembled_state,
            previous_state,
            iteration_state.previous_polish_provenance,
            assembled_polish_provenance)
    };
    const auto assembled_uses_polish{
        UsesPolish(assembled_polish_provenance)
    };
    bool objective_domain_changed{ false };
    if (has_new_terminal_failures)
    {
        auto remaining_active_index_list{
            iteration_state.terminal_failure_state.BuildEligibleActiveIndexList()
        };
        if (!remaining_active_index_list.empty())
        {
            auto remaining_graph_partition{
                BuildGraphPartition(graph_topology, remaining_active_index_list)
            };
            auto remaining_cluster_key_list{
                BuildGraphClusterKeyList(remaining_graph_partition)
            };
            const auto assembled_model_snapshot{
                BuildSecondStageModelSnapshot(context, assembled_state)
            };
            iteration_state.objective_domain = BuildObjectiveDomain(
                context,
                assembled_model_snapshot,
                remaining_cluster_key_list,
                options.distance_min,
                options.distance_max);
            iteration_state.cluster_objective_state.clear();
            const auto remaining_objective_by_key{
                BuildObjectiveByKey(
                    remaining_graph_partition,
                    iteration_state.objective_domain,
                    SnapshotResidualEvaluator{ context, assembled_model_snapshot })
            };
            ReconcileClusterObjectiveState(
                remaining_objective_by_key,
                iteration_state.cluster_objective_state);
            const auto reset_audit_objective{
                EvaluateAuditObjective(
                    iteration_state.objective_domain,
                    SnapshotResidualEvaluator{ context, assembled_model_snapshot })
            };
            iteration_state.best_audit_state.reset();
            if (reset_audit_objective.has_value())
            {
                TryUpdateBestAuditState(
                    assembled_state,
                    assembled_uses_polish,
                    iteration_state.accepted_iteration_count + 1,
                    *reset_audit_objective,
                    iteration_state.best_audit_state);
            }
            iteration_state.active_index_list = std::move(remaining_active_index_list);
            performance_counters.RecordSolverWorkspaceReset();
            ResetClusterSolverWorkspace(
                remaining_cluster_key_list,
                iteration_state.solver_workspace_by_key);
            iteration_state.graph_partition = std::move(remaining_graph_partition);
            objective_domain_changed = true;
        }
        else
        {
            iteration_state.active_index_list.clear();
        }
    }

    auto trust_region_iteration_update{
        iteration_state.trust_region_state.UpdateAfterIteration(
            selection.grow_trust_region_key_list,
            selection.rejected_key_list,
            selection.backtracking_exhausted_key_list)
    };

    IterationResult result;
    result.objective_domain_changed = objective_domain_changed;
    if (!selection.accepted_key_list.empty())
    {
        result.diagnostics.accepted_cluster_diagnostic_list = std::move(selection.accepted_cluster_diagnostic_list);
    }
    result.diagnostics.rejected_cluster_diagnostic_list = std::move(selection.rejected_cluster_diagnostic_list);
    result.diagnostics.combined_backtracking_trial_count = selection.combined_backtracking_trial_count;
    result.diagnostics.combined_backtracking_factor = selection.combined_backtracking_factor;
    result.diagnostics.combined_backtracking_exhausted = selection.combined_backtracking_exhausted;
    result.diagnostics.trust_region_update = std::move(trust_region_iteration_update);
    iteration_state.rollback_atom_mask = std::move(raw_iteration_result.rollback_atom_mask);
    result.progress = IterationProgress{
        attempt_number,
        iteration_state.accepted_iteration_count,
        context.size() - iteration_state.terminal_failure_state.AtomCount(),
        iteration_state.terminal_failure_state.AtomCount(),
        selection.accepted_key_list.size(),
        selection.rejected_key_list.size(),
        selection.polish_progress,
        iteration_suspicious_atom_count,
        std::nullopt,
        GetMaximumTransformedChange(raw_fixed_point_change_summary)
    };

    if (selection.accepted_key_list.empty())
    {
        result.all_rejected_resolution = ResolveAllRejected(
            attempt_number >= kMaximumIterations,
            result.diagnostics.trust_region_update.rejected_cluster_partition,
            result.diagnostics.trust_region_update.radius_update);
        if (*result.all_rejected_resolution == AllRejectedResolution::Retry)
        {
            for (const auto & key :
                result.diagnostics.trust_region_update.rejected_cluster_partition.exhausted_key_list)
            {
                if (std::ranges::find(
                        iteration_state.unchanged_state_exhausted_key_list,
                        key) == iteration_state.unchanged_state_exhausted_key_list.end())
                {
                    iteration_state.unchanged_state_exhausted_key_list.emplace_back(key);
                }
            }
        }
        return result;
    }

    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            assembled_state,
            previous_state,
            iteration_state.active_index_list)
    };
    iteration_state.accepted_iteration_count++;
    bool improved_best_audit{ false };
    if (!objective_domain_changed)
    {
        auto candidate_audit_objective{ selection.combined_backtracking_objective };
        if (!candidate_audit_objective.has_value())
        {
            const auto candidate_model_snapshot{
                BuildSecondStageModelSnapshot(context, assembled_state)
            };
            candidate_audit_objective = EvaluateAuditObjective(
                iteration_state.objective_domain,
                SnapshotResidualEvaluator{ context, candidate_model_snapshot });
        }
        if (candidate_audit_objective.has_value())
        {
            improved_best_audit = TryUpdateBestAuditState(
                assembled_state,
                assembled_uses_polish,
                iteration_state.accepted_iteration_count,
                *candidate_audit_objective,
                iteration_state.best_audit_state);
        }
    }
    if (improved_best_audit)
    {
        performance_counters.RecordFullStateMaterialization();
    }
    const auto changed_rejected_trust_radius{
        !selection.rejected_key_list.empty() &&
        !result.diagnostics.trust_region_update.radius_update.changed_key_list.empty()
    };
    if (objective_domain_changed || improved_best_audit || changed_rejected_trust_radius)
    {
        iteration_state.audit_patience_count = 0;
    }
    else
    {
        iteration_state.audit_patience_count++;
    }

    result.progress.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.progress.accepted_maximum_transformed_change =
        GetMaximumTransformedChange(transformed_change_summary);
    result.transformed_change_stats = transformed_change_summary.percentile_stats;
    result.audit_patience_exhausted = iteration_state.audit_patience_count >= kAuditPatience;
    result.converged =
        is_stationarity_eligible &&
        !has_suspicious_offset_fallback &&
        selection.rejected_key_list.empty() &&
        IsTransformedChangeConverged(transformed_change_summary) &&
        IsTransformedChangeConverged(raw_fixed_point_change_summary);

    iteration_state.previous_state = std::move(assembled_state);
    iteration_state.previous_polish_provenance = std::move(assembled_polish_provenance);
    iteration_state.unchanged_state_exhausted_key_list.clear();
    return result;
}

} // namespace rhbm_gem::core::detail
