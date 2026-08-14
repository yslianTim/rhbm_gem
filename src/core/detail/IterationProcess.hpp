#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/Logger.hpp>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "core/detail/JointOffset.hpp"
#include "core/detail/ClusterHealth.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/LocalFittingTerminalFailure.hpp"
#include "core/detail/Objective.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/SuspiciousUpdate.hpp"
#include "core/detail/ReusableWeightedRidgeSolver.hpp"
#include "core/detail/ScopedEigenThreadCount.hpp"
#include "core/detail/SecondStageContext.hpp"
#include "core/detail/TransformedChange.hpp"
#include "core/detail/TrustRegion.hpp"
#include "core/detail/LocalFittingPerformanceCounters.hpp"

namespace rhbm_gem::core::detail {

constexpr std::size_t kMaximumIterations{ 100 };
constexpr std::size_t kAuditPatience{ 3 };

enum class IterationOutcome
{
    Accepted,
    Retry,
    AllRejected
};

struct IterationDiagnostics
{
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    std::vector<RejectedClusterDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<RejectedClusterDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<ClusterKey> backtracking_exhausted_key_list{};
    std::size_t combined_backtracking_trial_count{ 0 };
    std::optional<double> combined_backtracking_factor{};
    std::optional<ObjectiveBreakdown> combined_backtracking_objective{};
    bool combined_backtracking_exhausted{ false };
    PolishProgress polish_progress{};
    RejectedClusterPartition rejected_cluster_partition{};
    TrustRegionRadiusUpdate trust_region_radius_update{};
};

struct IterationState
{
    FitState previous_state{};
    PolishProvenance previous_polish_provenance{};
    SuspiciousUpdateMask rollback_atom_mask{};
    std::vector<std::size_t> active_index_list{};
    CouplingGraphPartition graph_partition{};
    std::vector<ClusterKey> cluster_key_list{};
    ClusterSolverWorkspaceMap solver_workspace_by_key{};
    ObjectiveDomain objective_domain{};
    BestAuditState best_audit_state{};
    LocalFittingTerminalFailureState terminal_failure_state{};
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
    std::optional<double> raw_maximum_transformed_change{};
};

using ProgressColumnWidths = std::array<std::size_t, 6>;

struct IterationResult
{
    IterationOutcome outcome{ IterationOutcome::Accepted };
    IterationDiagnostics diagnostics{};
    IterationProgress progress{};
    std::optional<AllRejectedResolution> all_rejected_resolution{};
    bool objective_domain_changed{ false };
    bool converged{ false };
    bool audit_patience_exhausted{ false };
    std::size_t suspicious_atom_count{ 0 };
    std::optional<double> accepted_maximum_transformed_change{};
    std::optional<double> raw_maximum_transformed_change{};
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
        << breakdown->tail_validation_penalty << "/"
        << breakdown->offset_plausibility_penalty << "/"
        << breakdown->total_objective;
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
    const FitOptions & options,
    const std::vector<RejectedClusterDiagnostic> & diagnostic_list)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug || diagnostic_list.empty())
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
        if (diagnostic.fit_scale.has_value())
        {
            message << *diagnostic.fit_scale;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.tail_scale.has_value())
        {
            message << *diagnostic.tail_scale;
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
    const FitOptions & options,
    const RejectedClusterPartition & partition,
    const TrustRegionRadiusUpdate & radius_update,
    AllRejectedResolution resolution)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "All-rejected local fitting resolution: outcome = "
        << GetAllRejectedResolutionText(resolution)
        << ", exhausted/retryable/radius-changed/radius-saturated = "
        << partition.exhausted_key_list.size() << "/"
        << partition.retryable_key_list.size() << "/"
        << radius_update.changed_key_list.size() << "/"
        << radius_update.saturated_key_list.size() << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

inline void LogAcceptedBacktrackingDiagnostics(
    const FitOptions & options,
    const IterationDiagnostics & selection)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }
    const auto has_local_backtracking{
        std::any_of(
            selection.accepted_cluster_diagnostic_list.begin(),
            selection.accepted_cluster_diagnostic_list.end(),
            [&](const RejectedClusterDiagnostic & diagnostic)
            {
                return diagnostic.attempt.backtracking_trial_count > 1 &&
                    std::find(
                        selection.accepted_key_list.begin(),
                        selection.accepted_key_list.end(),
                        diagnostic.key) != selection.accepted_key_list.end();
            })
    };
    if (!has_local_backtracking && selection.combined_backtracking_trial_count <= 1)
    {
        return;
    }
    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : selection.accepted_cluster_diagnostic_list)
    {
        const auto & diagnostic{ cluster_diagnostic.attempt };
        if (diagnostic.backtracking_trial_count <= 1 ||
            std::find(
                selection.accepted_key_list.begin(),
                selection.accepted_key_list.end(),
                cluster_diagnostic.key) == selection.accepted_key_list.end())
        {
            continue;
        }
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
        if (diagnostic.fit_scale.has_value())
        {
            message << *diagnostic.fit_scale;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.tail_scale.has_value())
        {
            message << *diagnostic.tail_scale;
        }
        else
        {
            message << "empty";
        }
        message << ".";
        Logger::Log(LogLevel::Debug, message.str());
    }
    if (selection.combined_backtracking_trial_count <= 1) return;
    std::ostringstream message;
    message
        << "Combined-objective backtracking: trials/factor/exhausted = "
        << selection.combined_backtracking_trial_count << "/";
    if (selection.combined_backtracking_factor.has_value())
    {
        message << *selection.combined_backtracking_factor;
    }
    else
    {
        message << "-";
    }
    message << "/"
        << (selection.combined_backtracking_exhausted ? "yes" : "no")
        << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

inline std::string FormatProgressMaximum(const std::optional<double> & value)
{
    if (!value.has_value()) return "-";
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << *value;
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
    const std::array<std::string, 6> & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream << std::left << std::setw(static_cast<int>(column_widths.at(i)))
            << cell_list.at(i);
    }
    return stream.str();
}

inline ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_size)
{
    const auto maximum_iteration_text{
        std::to_string(kMaximumIterations)
    };
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

inline void LogProgressHeader(const FitOptions & options, const ProgressColumnWidths & column_widths)
{
    if (options.quiet_mode) return;
    std::array<std::string, 6> header_list;
    for (std::size_t i = 0; i < header_list.size(); i++)
    {
        header_list.at(i) = kProgressHeaderList.at(i);
    }
    Logger::Log(LogLevel::Info, FormatProgressRow(column_widths, header_list));
}

inline void LogIterationProgress(
    const FitOptions & options,
    const ProgressColumnWidths & column_widths,
    const IterationProgress & progress)
{
    if (options.quiet_mode) return;

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
        FormatProgressMaximum(progress.accepted_maximum_transformed_change) + "/" +
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
    const SecondStageContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const GaussianModel3D & offset_model,
    const std::vector<double> & adjusted_response_list,
    const FitOptions & options)
{
    auto adjusted_sampling_entries{
        BuildSecondStageAdjustedSamples(context, atom_index, adjusted_response_list)
    };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto previous_baseline{
        local_fitting_suspicious_internal::BuildPreviousSuspiciousProfileBaseline(
            adjusted_sampling_entries,
            previous_model,
            options)
    };
    const auto is_post_refit_candidate_acceptable = [&](const GaussianModel3D & model)
    {
        const auto reason{ local_fitting_suspicious_internal::EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                previous_model,
                model,
                options,
                previous_baseline,
                true) };
        return reason == SuspiciousGaussianReason::None;
    };
    const auto is_offset_only_fallback_acceptable =
        [&](const GaussianModel3D & model)
    {
        const auto reason{ local_fitting_suspicious_internal::EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                previous_model,
                model,
                options,
                previous_baseline,
                false) };
        return reason == SuspiciousGaussianReason::None;
    };
    try
    {
        auto candidate_result{
            EstimateLocalGaussianPrepared(
                context.at(atom_index).refit_design_template,
                adjusted_response_list,
                context.at(atom_index).alpha_r,
                options,
                offset_model)
        };
        const auto candidate_model{ candidate_result.mdpde.GetModel() };
        if (is_post_refit_candidate_acceptable(candidate_model))
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
    if (!is_offset_only_fallback_acceptable(result.mdpde.GetModel()))
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
    const auto selected_atom_size{ context.size() };
    auto current_model_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(previous_state))
    };
    const auto is_debug_logging_enabled{
        Logger::GetLogLevel() >= LogLevel::Debug
    };
    const auto log_debug_diagnostics{
        !options.quiet_mode && is_debug_logging_enabled
    };
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
        ScopedEigenThreadCount eigen_thread_guard{ 1 };
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
                current_model_snapshot.selected.at(atom_index)
                    .WithOffset(result.offset(static_cast<Eigen::Index>(i)));
        }
        health_by_key.emplace(key, ClusterHealth{ result.status });
    }

    auto iteration_state{ previous_state };
    SuspiciousUpdateMask rollback_atom_mask(selected_atom_size, 0);
    std::vector<GroupKey> group_key_by_atom_index;
    group_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_key_by_atom_index.emplace_back(atom_context.group_key);
    }
    for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        std::vector<GroupKey> group_key_by_position;
        SuspiciousUpdateMask suspicious_seed_mask(key.size(), 0);
        group_key_by_position.reserve(key.size());
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            group_key_by_position.emplace_back(group_key_by_atom_index.at(atom_index));
            if (EvaluateSuspiciousOffsetUpdate(
                    context.at(atom_index).raw_sampling_entries,
                    previous_state.at(atom_index).mdpde.GetModel(),
                    current_model_snapshot.selected.at(atom_index),
                    options) != SuspiciousGaussianReason::None)
            {
                suspicious_seed_mask.at(position) = 1;
            }
        }
        const auto cluster_rollback_mask{
            ExpandSuspiciousSharedOffsetGroups(group_key_by_position, suspicious_seed_mask)
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

    auto refit_model_snapshot{
        BuildLocalFittingGroupMedianModelList(
            group_key_by_atom_index,
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
                    context,
                    atom_index,
                    previous_state.at(atom_index),
                    refit_model_bundle.selected.at(atom_index),
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
        ScopedEigenThreadCount eigen_thread_guard{ 1 };
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

    std::vector<std::size_t> post_refit_suspicious_seed_atom_index_list;
    std::size_t refit_position{ 0 };
    for (const auto & key : cluster_key_list)
    {
        auto & health{ health_by_key.at(key) };
        for (const auto atom_index : key)
        {
            if (rollback_atom_mask.at(atom_index) != 0) continue;

            auto refit_result{ std::move(refit_result_list.at(refit_position++)) };
            if (!refit_result.has_value())
            {
                health.is_refit_stationarity_eligible = false;
                post_refit_suspicious_seed_atom_index_list.emplace_back(atom_index);
                continue;
            }
            if (!refit_result->is_stationarity_eligible)
            {
                health.is_refit_stationarity_eligible = false;
            }
            iteration_state.at(atom_index) = std::move(refit_result->result);
        }
    }
    ExpandPostRefitSuspiciousClusters(
        cluster_key_list,
        post_refit_suspicious_seed_atom_index_list,
        rollback_atom_mask);
    for (std::size_t atom_index = 0; atom_index < rollback_atom_mask.size(); atom_index++)
    {
        if (rollback_atom_mask.at(atom_index) == 0) continue;
        iteration_state.at(atom_index) = previous_state.at(atom_index);
    }

    RawIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
    iteration_result.rollback_atom_mask = std::move(rollback_atom_mask);
    iteration_result.health_by_key = std::move(health_by_key);
    return iteration_result;
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
    iteration_state.terminal_failure_state =
        LocalFittingTerminalFailureState(context.size());
    iteration_state.active_index_list =
        iteration_state.terminal_failure_state.BuildEligibleActiveIndexList();
    iteration_state.graph_partition = BuildGraphPartition(
        graph_topology,
        iteration_state.active_index_list);
    iteration_state.cluster_key_list = BuildGraphClusterKeyList(iteration_state.graph_partition);
    ResetClusterSolverWorkspace(
        iteration_state.cluster_key_list,
        iteration_state.solver_workspace_by_key);
    iteration_state.objective_domain = BuildObjectiveDomain(
        context,
        iteration_state.previous_state,
        iteration_state.graph_partition,
        options);
    iteration_state.best_audit_state = BuildInitialBestAuditState(
        context,
        iteration_state.previous_state,
        iteration_state.previous_polish_provenance,
        std::nullopt,
        iteration_state.objective_domain);
    return iteration_state;
}

inline IterationResult RunIteration(
    const SecondStageContext & context,
    const GraphTopology & graph_topology,
    const FitOptions & options,
    std::size_t attempt_number,
    std::size_t cached_sample_count,
    IterationState & iteration_state,
    LocalFittingPerformanceCounters & performance_counters)
{
    const auto & previous_state{ iteration_state.previous_state };
    const auto & active_index_list{ iteration_state.active_index_list };
    const auto & graph_partition{ iteration_state.graph_partition };
    const auto & cluster_key_list{ iteration_state.cluster_key_list };
    const auto & objective_domain{ iteration_state.objective_domain };

    const auto previous_model_snapshot{
        BuildSecondStageModelSnapshot(context, previous_state)
    };
    const auto residual_baseline{
        BuildResidualBaseline(context, previous_state, previous_model_snapshot)
    };
    performance_counters.gaussian_cache_miss_count += cached_sample_count;

    const auto previous_objective_by_key{
        BuildObjectiveByKey(
            context,
            previous_state,
            graph_partition,
            objective_domain,
            residual_baseline)
    };
    ReconcileClusterObjectiveState(
        graph_partition,
        previous_objective_by_key,
        iteration_state.cluster_objective_state);
    iteration_state.trust_region_state.Reconcile(cluster_key_list);

    const auto joint_offset_ridge_multiplier_list{
        BuildSuspiciousJointOffsetRidgeMultiplierList(
            context.size(),
            iteration_state.rollback_atom_mask)
    };

    const auto iteration_phase_start{ std::chrono::steady_clock::now() };
    auto raw_iteration_result{
        RunRawIteration(
            context,
            cluster_key_list,
            previous_state,
            options,
            joint_offset_ridge_multiplier_list,
            iteration_state.solver_workspace_by_key)
    };
    performance_counters.iteration_phase_milliseconds +=
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - iteration_phase_start).count();
    performance_counters.gaussian_cache_hit_count += cached_sample_count;

    const auto iteration_suspicious_atom_count{
        CountSuspiciousAtoms(raw_iteration_result.rollback_atom_mask)
    };
    const auto has_suspicious_offset_fallback{
        iteration_suspicious_atom_count > 0
    };
    const auto is_stationarity_eligible{
        AreClustersStationarityEligible(raw_iteration_result.health_by_key)
    };
    const auto & raw_state{ raw_iteration_result.state };
    const auto raw_fixed_point_change_summary{
        SummarizeTransformedChanges(raw_state, previous_state, active_index_list)
    };
    const auto previous_transformed_estimation_list{
        BuildTransformedEstimationList(previous_state)
    };
    const auto raw_transformed_estimation_list{
        BuildTransformedEstimationList(raw_state)
    };

    auto working_cluster_objective_state{
        iteration_state.cluster_objective_state
    };
    const auto candidate_phase_start{ std::chrono::steady_clock::now() };
    auto selection{
        SelectClusterCandidates(
            context,
            previous_model_snapshot,
            residual_baseline,
            graph_partition,
            raw_iteration_result.health_by_key,
            previous_state,
            iteration_state.previous_polish_provenance,
            raw_state,
            previous_transformed_estimation_list,
            raw_transformed_estimation_list,
            raw_iteration_result.rollback_atom_mask,
            joint_offset_ridge_multiplier_list,
            iteration_state.unchanged_state_exhausted_key_list,
            objective_domain,
            previous_objective_by_key,
            working_cluster_objective_state,
            iteration_state.trust_region_state,
            iteration_state.solver_workspace_by_key,
            options.thread_size,
            performance_counters)
    };
    performance_counters.candidate_phase_milliseconds +=
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - candidate_phase_start).count();
    performance_counters.full_state_materialization_count++;

    const auto combined_changed_key_list{ selection.accepted_key_list };
    const auto combined_check{
        EvaluateCombinedCandidateObjective(
            context,
            previous_model_snapshot,
            residual_baseline,
            graph_partition,
            previous_state,
            selection.assembled_state,
            selection.accepted_key_list,
            objective_domain,
            iteration_state.best_audit_state,
            &performance_counters)
    };
    selection.combined_backtracking_objective = combined_check.candidate_objective;
    auto combined_objective_accepted{
        !combined_check.guard_required || combined_check.accepted
    };
    if (!combined_objective_accepted)
    {
        combined_objective_accepted =
            TryBacktrackCombinedCandidate(
                context,
                previous_model_snapshot,
                residual_baseline,
                graph_partition,
                previous_state,
                iteration_state.previous_polish_provenance,
                objective_domain,
                previous_objective_by_key,
                combined_check.previous_objective,
                iteration_state.best_audit_state,
                iteration_state.cluster_objective_state,
                working_cluster_objective_state,
                selection,
                performance_counters);
    }
    if (!combined_objective_accepted)
    {
        RejectCombinedCandidate(
            previous_state,
            iteration_state.previous_polish_provenance,
            cluster_key_list,
            selection);
        if (selection.combined_backtracking_exhausted)
        {
            for (const auto & key : combined_changed_key_list)
            {
                if (std::find(
                        selection.backtracking_exhausted_key_list.begin(),
                        selection.backtracking_exhausted_key_list.end(),
                        key) == selection.backtracking_exhausted_key_list.end())
                {
                    selection.backtracking_exhausted_key_list.emplace_back(key);
                }
            }
        }
    }
    else
    {
        iteration_state.cluster_objective_state = std::move(working_cluster_objective_state);
    }

    auto assembled_state{ std::move(selection.assembled_state) };
    auto assembled_polish_provenance{ std::move(selection.assembled_polish_provenance) };
    const auto terminal_key_list{
        iteration_state.terminal_failure_state.IsolatePersistentFailures(
            selection.accepted_key_list,
            raw_iteration_result.rollback_atom_mask,
            raw_iteration_result.health_by_key,
            assembled_state,
            previous_state,
            iteration_state.previous_polish_provenance,
            assembled_polish_provenance)
    };
    bool objective_domain_changed{ false };
    if (!terminal_key_list.empty())
    {
        auto remaining_active_index_list{
            iteration_state.terminal_failure_state.BuildEligibleActiveIndexList()
        };
        if (!remaining_active_index_list.empty())
        {
            auto remaining_graph_partition{
                BuildGraphPartition(graph_topology, remaining_active_index_list)
            };
            iteration_state.objective_domain = BuildObjectiveDomain(
                context,
                assembled_state,
                remaining_graph_partition,
                options);
            iteration_state.cluster_objective_state.clear();
            const auto assembled_model_snapshot{
                BuildSecondStageModelSnapshot(context, assembled_state)
            };
            const auto remaining_objective_by_key{
                BuildObjectiveByKey(
                    context,
                    assembled_state,
                    remaining_graph_partition,
                    iteration_state.objective_domain,
                    assembled_model_snapshot)
            };
            ReconcileClusterObjectiveState(
                remaining_graph_partition,
                remaining_objective_by_key,
                iteration_state.cluster_objective_state);
            ResetBestAuditAfterObjectiveDomainChange(
                context,
                assembled_state,
                assembled_polish_provenance,
                iteration_state.accepted_iteration_count + 1,
                iteration_state.objective_domain,
                iteration_state.best_audit_state);
            iteration_state.active_index_list = std::move(remaining_active_index_list);
            iteration_state.cluster_key_list = BuildGraphClusterKeyList(remaining_graph_partition);
            ResetClusterSolverWorkspace(
                iteration_state.cluster_key_list,
                iteration_state.solver_workspace_by_key);
            iteration_state.graph_partition = std::move(remaining_graph_partition);
            objective_domain_changed = true;
        }
        else
        {
            iteration_state.active_index_list.clear();
        }
    }

    const auto trust_region_iteration_update{
        iteration_state.trust_region_state.UpdateAfterIteration(
            selection.grow_trust_region_key_list,
            selection.rejected_key_list,
            selection.backtracking_exhausted_key_list)
    };
    const auto & rejected_cluster_partition{
        trust_region_iteration_update.rejected_cluster_partition
    };
    const auto & trust_region_radius_update{
        trust_region_iteration_update.radius_update
    };

    IterationResult result;
    result.objective_domain_changed = objective_domain_changed;
    result.suspicious_atom_count = iteration_suspicious_atom_count;
    result.raw_maximum_transformed_change = GetMaximumTransformedChange(raw_fixed_point_change_summary);
    result.diagnostics.accepted_key_list = std::move(selection.accepted_key_list);
    result.diagnostics.rejected_key_list = std::move(selection.rejected_key_list);
    result.diagnostics.accepted_cluster_diagnostic_list = std::move(selection.accepted_cluster_diagnostic_list);
    result.diagnostics.rejected_cluster_diagnostic_list = std::move(selection.rejected_cluster_diagnostic_list);
    result.diagnostics.backtracking_exhausted_key_list = std::move(selection.backtracking_exhausted_key_list);
    result.diagnostics.combined_backtracking_trial_count = selection.combined_backtracking_trial_count;
    result.diagnostics.combined_backtracking_factor = selection.combined_backtracking_factor;
    result.diagnostics.combined_backtracking_objective = selection.combined_backtracking_objective;
    result.diagnostics.combined_backtracking_exhausted = selection.combined_backtracking_exhausted;
    result.diagnostics.polish_progress = selection.polish_progress;
    result.diagnostics.rejected_cluster_partition = rejected_cluster_partition;
    result.diagnostics.trust_region_radius_update = trust_region_radius_update;
    iteration_state.rollback_atom_mask = std::move(raw_iteration_result.rollback_atom_mask);
    result.progress = IterationProgress{
        attempt_number,
        iteration_state.accepted_iteration_count,
        context.size() - iteration_state.terminal_failure_state.AtomCount(),
        iteration_state.terminal_failure_state.AtomCount(),
        result.diagnostics.accepted_key_list.size(),
        result.diagnostics.rejected_key_list.size(),
        result.diagnostics.polish_progress,
        result.suspicious_atom_count,
        std::nullopt,
        result.raw_maximum_transformed_change
    };

    if (result.diagnostics.accepted_key_list.empty())
    {
        result.outcome = IterationOutcome::AllRejected;
        result.all_rejected_resolution = ResolveAllRejected(
            attempt_number >= kMaximumIterations,
            result.diagnostics.rejected_cluster_partition,
            result.diagnostics.trust_region_radius_update);
        if (*result.all_rejected_resolution == AllRejectedResolution::Retry)
        {
            for (const auto & key :
                result.diagnostics.rejected_cluster_partition.exhausted_key_list)
            {
                if (std::find(
                        iteration_state.unchanged_state_exhausted_key_list.begin(),
                        iteration_state.unchanged_state_exhausted_key_list.end(),
                        key) == iteration_state.unchanged_state_exhausted_key_list.end())
                {
                    iteration_state.unchanged_state_exhausted_key_list.emplace_back(key);
                }
            }
            result.outcome = IterationOutcome::Retry;
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
    const auto improved_best_audit{
        TryUpdateBestAuditState(
            context,
            assembled_state,
            assembled_polish_provenance,
            iteration_state.accepted_iteration_count,
            iteration_state.objective_domain,
            iteration_state.best_audit_state,
            objective_domain_changed ?
                std::nullopt : result.diagnostics.combined_backtracking_objective)
    };
    if (improved_best_audit)
    {
        performance_counters.full_state_materialization_count++;
    }
    const auto changed_rejected_trust_radius{
        !result.diagnostics.rejected_key_list.empty() &&
        !result.diagnostics.trust_region_radius_update.changed_key_list.empty()
    };
    if (objective_domain_changed || improved_best_audit || changed_rejected_trust_radius)
    {
        iteration_state.audit_patience_count = 0;
    }
    else
    {
        iteration_state.audit_patience_count++;
    }

    result.accepted_maximum_transformed_change =
        GetMaximumTransformedChange(transformed_change_summary);
    result.progress.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.progress.accepted_maximum_transformed_change = result.accepted_maximum_transformed_change;
    result.transformed_change_stats = transformed_change_summary.percentile_stats;
    result.audit_patience_exhausted = iteration_state.audit_patience_count >= kAuditPatience;
    result.converged =
        is_stationarity_eligible &&
        !has_suspicious_offset_fallback &&
        result.diagnostics.rejected_key_list.empty() &&
        IsTransformedChangeConverged(transformed_change_summary) &&
        IsTransformedChangeConverged(raw_fixed_point_change_summary);

    iteration_state.previous_state = std::move(assembled_state);
    iteration_state.previous_polish_provenance = std::move(assembled_polish_provenance);
    iteration_state.unchanged_state_exhausted_key_list.clear();
    result.outcome = IterationOutcome::Accepted;
    return result;
}

} // namespace rhbm_gem::core::detail
