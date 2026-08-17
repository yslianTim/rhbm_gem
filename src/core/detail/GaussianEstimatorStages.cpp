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
#include "core/detail/TransformedGaussianEvaluation.hpp"
#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/SecondStageContext.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core {

namespace {

using detail::BuildSecondStageModelSnapshot;
using detail::SampleRef;
using detail::PerformanceCounters;
using detail::PolishProvenance;
using detail::FitState;
using detail::SecondStageContext;
using detail::SecondStageInitialStateBuildResult;

constexpr double kOffsetSummaryPercentile{ 0.99 };

using detail::kMaximumIterations;
using detail::BuildProgressColumnWidths;
using detail::LogProgressHeader;
using detail::LogIterationProgress;
using detail::LogRejectedClusterDiagnostics;
using detail::LogAcceptedBacktrackingDiagnostics;
using detail::LogAllRejectedResolution;
using detail::GetAllRejectedResolutionText;
using detail::BuildSecondStageAdjustedSamples;
using detail::ObjectiveDomain;
using detail::AuditedState;
using detail::BestAuditState;
using detail::TerminalSummary;
using detail::AppendTerminalSummary;
using detail::UsesPolish;
using detail::kFitRangeWeight;
using detail::kTailValidationWeight;
using detail::kOffsetPlausibilityPenaltyWeight;
using detail::BuildSecondStageContext;
using detail::StoreSecondStageNeighborCounts;
using detail::BuildInitialFitState;
using detail::GetSecondStageSeedSourceText;
using detail::SecondStageSeedSelectionRecord;
using detail::BuildIterationState;
using detail::IterationOutcome;
using detail::RunIteration;

enum class FinalStateSource
{
    BestAudit,
    LatestValidated,
    Unavailable
};

struct FinalStateSelection
{
    const FitState * state{ nullptr };
    const PolishProvenance * polish_provenance{ nullptr };
    const AuditedState * audit_state{ nullptr };
    FinalStateSource source{ FinalStateSource::Unavailable };
};

FinalStateSelection SelectFinalState(
    const FitState & latest_validated_state,
    const PolishProvenance & latest_validated_polish_provenance,
    const std::optional<AuditedState> & audited_state)
{
    if (audited_state.has_value())
    {
        return FinalStateSelection{
            &audited_state->state,
            &audited_state->polish_provenance,
            &*audited_state,
            FinalStateSource::BestAudit
        };
    }
    return FinalStateSelection{
        &latest_validated_state,
        &latest_validated_polish_provenance,
        nullptr,
        FinalStateSource::LatestValidated
    };
}

std::string_view GetFinalStateSourceText(FinalStateSource source)
{
    switch (source)
    {
    case FinalStateSource::BestAudit:
        return "best-audit";
    case FinalStateSource::LatestValidated:
        return "latest-validated";
    case FinalStateSource::Unavailable:
        return "unavailable";
    }
    return "unavailable";
}

struct OffsetStat
{
    std::size_t atom_count{ 0 };
    std::size_t finite_count{ 0 };
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
};

detail::GraphTopology BuildGraphTopology(
    const SecondStageContext & context,
    const FitState & initial_state,
    bool quiet_mode)
{
    std::size_t total_sample_count{ 0 };
    for (const auto & atom_context : context)
    {
        total_sample_count += atom_context.raw_sampling_entries.size();
    }
    const std::size_t total_work{ total_sample_count + 1 };
    std::size_t completed_work{ 0 };
    const std::string progress_message{ " Build local-fitting coupling topology" };
    int last_progress_percent{ -1 };
    const auto update_progress = [&]()
    {
        if (quiet_mode) return;
        const auto progress_percent{ static_cast<int>(
            100.0 * static_cast<double>(completed_work) / static_cast<double>(total_work)) };
        if (progress_percent == last_progress_percent) return;
        last_progress_percent = progress_percent;
        Logger::ProgressPercent(
            completed_work,
            total_work,
            50,
            progress_message);
    };
    update_progress();

    detail::CouplingGraphBuilder builder{ context.size() };
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, initial_state)
    };
    std::vector<std::optional<detail::TransformedModelInvariants>> selected_model_invariants;
    selected_model_invariants.reserve(model_snapshot.selected.model_list.size());
    for (const auto & model : model_snapshot.selected.model_list)
    {
        selected_model_invariants.emplace_back(detail::BuildTransformedModelInvariants(model));
    }
    std::vector<std::optional<detail::TransformedModelInvariants>> unselected_model_invariants;
    unselected_model_invariants.reserve(model_snapshot.unselected.model_list.size());
    for (const auto & model : model_snapshot.unselected.model_list)
    {
        unselected_model_invariants.emplace_back(detail::BuildTransformedModelInvariants(model));
    }
    const auto evaluate_model = [](
        const std::optional<detail::TransformedModelInvariants> & invariants,
        double distance)
    {
        if (!invariants.has_value())
        {
            return std::optional<detail::TransformedResponse>{};
        }
        return detail::EvaluateTransformedResponse(*invariants, distance);
    };
    const auto invalid_jacobian{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    std::vector<detail::GraphParticipant> participant_list;
    participant_list.reserve(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        const auto & atom_context{ context.at(atom_index) };
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            const auto target_evaluation{ evaluate_model(
                selected_model_invariants.at(atom_index),
                static_cast<double>(sample.point.distance)) };
            participant_list.clear();
            participant_list.emplace_back(
                detail::GraphParticipant{
                    atom_index,
                    target_evaluation.has_value() ?
                        target_evaluation->jacobian : invalid_jacobian
                });
            for (auto neighbor_iter = atom_context.NeighborBegin(sample_index);
                neighbor_iter != atom_context.NeighborEnd(sample_index);
                ++neighbor_iter)
            {
                const auto & neighbor_atom_sample{ *neighbor_iter };
                const auto neighbor_evaluation{ neighbor_atom_sample.is_selected ?
                    evaluate_model(
                        selected_model_invariants.at(neighbor_atom_sample.atom_index),
                        neighbor_atom_sample.distance) :
                    evaluate_model(
                        unselected_model_invariants.at(neighbor_atom_sample.atom_index),
                        neighbor_atom_sample.distance) };
                const auto jacobian{
                    neighbor_evaluation.has_value() ? neighbor_evaluation->jacobian : invalid_jacobian
                };
                if (neighbor_atom_sample.is_selected)
                {
                    participant_list.emplace_back(
                        detail::GraphParticipant{
                            neighbor_atom_sample.atom_index,
                            jacobian
                        });
                    continue;
                }
                const auto selected_group_id{
                    context.unselected_atom_list.at(neighbor_atom_sample.atom_index).selected_group_id
                };
                if (!selected_group_id.has_value()) continue;
                for (const auto selected_index :
                    context.selected_atom_index_list_by_group.at(*selected_group_id))
                {
                    participant_list.emplace_back(
                        detail::GraphParticipant{
                            selected_index,
                            jacobian
                        });
                }
            }
            builder.AddSample(
                SampleRef{ atom_index, sample_index },
                participant_list);
            completed_work++;
            update_progress();
        }
    }

    std::vector<detail::ResidueKey> residue_key_by_atom_index;
    residue_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        residue_key_by_atom_index.emplace_back(atom_context.residue_key);
    }
    const auto topology{ builder.BuildTopology(std::move(residue_key_by_atom_index)) };
    completed_work++;
    update_progress();
    return topology;
}

void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    constexpr std::array<detail::SecondStageSeedSource, 4> source_list{
        detail::SecondStageSeedSource::GroupPosterior,
        detail::SecondStageSeedSource::GroupPrior,
        detail::SecondStageSeedSource::GroupMedian,
        detail::SecondStageSeedSource::GlobalMedian
    };
    std::array<std::size_t, source_list.size()> source_count{};
    for (const auto & record : selection_record_list)
    {
        source_count.at(static_cast<std::size_t>(record.source))++;
    }

    std::ostringstream summary;
    summary << "Selected second-stage initial seeds = "
        << selection_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < source_list.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedSourceText(source_list.at(i))
            << ":" << source_count.at(i);
    }
    summary << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        std::ostringstream detail_message;
        detail_message << "Second-stage seed selection: atom index = "
            << record.atom_index
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", original MDPDE A/B/C = "
            << record.original_model.GetAmplitude() << "/"
            << record.original_model.GetWidth() << "/"
            << record.original_model.GetOffset()
            << ", selected A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

void LogUnselectedSecondStageSeedSelections(
    const SecondStageContext & context,
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    std::size_t group_median_count{ 0 };
    std::size_t global_median_count{ 0 };
    for (const auto & record : selection_record_list)
    {
        if (record.source == detail::SecondStageSeedSource::GroupMedian)
        {
            group_median_count++;
        }
        else if (record.source == detail::SecondStageSeedSource::GlobalMedian)
        {
            global_median_count++;
        }
    }

    std::ostringstream summary;
    summary << "Unselected second-stage neighbor seeds = " << selection_record_list.size()
        << ", sources = group-median:" << group_median_count
        << ", global-median:" << global_median_count << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        const auto & unselected_atom_contributor{
            context.unselected_atom_list.at(record.atom_index)
        };
        std::ostringstream detail_message;
        detail_message
            << "Unselected second-stage neighbor seed selection: serial ID = "
            << unselected_atom_contributor.atom->GetSerialID()
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", seed A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

void LogGraphTopology(const detail::GraphTopology & topology, bool quiet_mode)
{
    if (quiet_mode) return;

    const auto & summary{ topology.summary };
    if (!summary.uses_weighted_graph)
    {
        Logger::Log(LogLevel::Warning,
            "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
        if (Logger::GetLogLevel() >= LogLevel::Debug)
        {
            Logger::Log(LogLevel::Debug,
                "Local-fitting weighted threshold sensitivity is unavailable in binary fallback mode.");
        }
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

    const auto & residue_cutoff_summary{ topology.residue_cutoff_summary };
    std::ostringstream residue_cutoff_message;
    residue_cutoff_message
        << "Local-fitting residue cutoff: residues="
        << residue_cutoff_summary.residue_count
        << ", limit=" << residue_cutoff_summary.maximum_residue_count_limit
        << ", clusters=" << residue_cutoff_summary.cluster_count
        << ", max-residues=" << residue_cutoff_summary.maximum_residue_count
        << ", cutoff-edges=" << residue_cutoff_summary.cut_edge_count << ".";
    Logger::Log(LogLevel::Info, residue_cutoff_message.str());

    for (const auto & sensitivity : summary.threshold_sensitivity_list)
    {
        std::ostringstream sensitivity_message;
        sensitivity_message
            << std::scientific << std::setprecision(2)
            << "Coupling sensitivity: threshold=" << sensitivity.minimum_weight
            << ", retained/cut="
            << sensitivity.retained_edge_count << "/"
            << sensitivity.cut_edge_count
            << ", components/max-atoms/ratio="
            << sensitivity.component_count << "/"
            << sensitivity.maximum_component_size << "/"
            << std::fixed << std::setprecision(2)
            << sensitivity.maximum_component_ratio << ".";
        Logger::Log(LogLevel::Info, sensitivity_message.str());
    }
}

void LogObjectiveDomain(const ObjectiveDomain & domain, bool quiet_mode, bool is_terminal_reset = false)
{
    if (quiet_mode) return;
    std::vector<double> fit_scale_list;
    std::vector<double> tail_scale_list;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        static_cast<void>(key);
        if (!cluster_domain.scale.has_value()) continue;
        fit_scale_list.emplace_back(cluster_domain.scale->fit);
        if (cluster_domain.scale->tail.has_value())
        {
            tail_scale_list.emplace_back(*cluster_domain.scale->tail);
        }
    }
    const auto append_scale_summary = [](
        std::ostringstream & message,
        const std::vector<double> & scale_list)
    {
        if (scale_list.empty())
        {
            message << "unavailable";
            return;
        }
        message
            << array_helper::ComputePercentile(scale_list, 0.5) << "/"
            << array_helper::ComputePercentile(scale_list, 0.99) << "/"
            << *std::max_element(scale_list.begin(), scale_list.end());
    };

    std::ostringstream message;
    message
        << (is_terminal_reset ?
            "Reset second-stage objective domain" : "Initialize second-stage objective domain")
        << ": fit/tail/offset weights = "
        << kFitRangeWeight << "/"
        << kTailValidationWeight << "/"
        << kOffsetPlausibilityPenaltyWeight
        << ", clusters = " << domain.cluster_by_key.size()
        << ", active atoms = " << domain.active_atom_count
        << ", unique fit/tail samples = "
        << domain.fit_sample_count << "/"
        << domain.tail_sample_count
        << ", fixed fit scale median/p99/max = ";
    append_scale_summary(message, fit_scale_list);
    message << ", fixed tail scale median/p99/max = ";
    append_scale_summary(message, tail_scale_list);
    message << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message.str());
}

void ApplyFitState(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitState & iteration_state)
{
    auto adjusted_sampling_entries_list{
        BuildSecondStageAdjustedSamples(context, iteration_state)
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

OffsetStat SummarizeOffsets(const FitState & state)
{
    OffsetStat stats;
    stats.atom_count = state.size();
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(state.size());
    for (const auto & result : state)
    {
        const auto offset{ result.mdpde.GetModel().GetOffset() };
        if (!std::isfinite(offset)) continue;
        stats.finite_count++;
        absolute_offset_list.emplace_back(std::abs(offset));
    }
    if (absolute_offset_list.empty()) return stats;

    stats.median_absolute_offset = array_helper::ComputeMedian(absolute_offset_list);
    stats.percentile_absolute_offset = array_helper::ComputePercentile(
        absolute_offset_list,
        kOffsetSummaryPercentile);
    stats.maximum_absolute_offset = *std::max_element(absolute_offset_list.begin(), absolute_offset_list.end());
    return stats;
}

void AppendOffsetSummary(std::ostringstream & stream, const OffsetStat & stats)
{
    stream
        << std::scientific << std::setprecision(2)
        << "; offsets finite = " << stats.finite_count
        << " of " << stats.atom_count
        << ", |C| median/p99/max = "
        << stats.median_absolute_offset << "/"
        << stats.percentile_absolute_offset << "/"
        << stats.maximum_absolute_offset;
}

void AppendAuditSummary(std::ostringstream & stream, const AuditedState & audited_state)
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
        << kTailValidationWeight;
}

void LogTerminalFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const TerminalSummary & terminal_summary,
    const OffsetStat & offset_stats)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Completed local fitting after " << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendTerminalSummary(warning_message, terminal_summary);
    AppendOffsetSummary(warning_message, offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void FinishWithNoActiveAtoms(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitState & state,
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const TerminalSummary & terminal_summary)
{
    ApplyFitState(model_object, context, state);
    const auto offset_stats{ SummarizeOffsets(state) };
    if (terminal_summary.HasFailures())
    {
        LogTerminalFallback(quiet_mode, accepted_iteration_count, terminal_summary, offset_stats);
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
    const OffsetStat & offset_stats)
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
    AppendOffsetSummary(message, offset_stats);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogMaximumIterations(
    bool quiet_mode,
    const FinalStateSelection & final_state_selection,
    const TerminalSummary & terminal_summary,
    const OffsetStat & applied_offset_stats)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendTerminalSummary(warning_message, terminal_summary);
    if (final_state_selection.source == FinalStateSource::BestAudit &&
        final_state_selection.audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendAuditSummary(warning_message, *final_state_selection.audit_state);
    }
    else if (final_state_selection.source == FinalStateSource::LatestValidated)
    {
        warning_message << "; applying latest validated state";
    }
    else
    {
        warning_message << "; no validated state is available";
    }
    AppendOffsetSummary(warning_message, applied_offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::string_view stop_reason,
    const BestAuditState & best_audit_state,
    std::optional<bool> final_uses_polish,
    FinalStateSource final_state_source)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Second-stage local fitting summary: accepted_iterations="
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
        message << std::scientific << std::setprecision(8)
            << best_audit_state.best->objective.total_objective;
    }
    else
    {
        message << "unavailable";
    }
    message << ", final_uses_polish=";
    if (!final_uses_polish.has_value())
    {
        message << "unavailable";
    }
    else
    {
        message << (*final_uses_polish ? "yes" : "no");
    }
    message
        << ", final_state_source="
        << GetFinalStateSourceText(final_state_source)
        << ".";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace

bool RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    auto context{ BuildSecondStageContext(model_object, options) };
    StoreSecondStageNeighborCounts(model_object, context);
    const auto atom_size{ context.size() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto initial_state_build_result{ BuildInitialFitState(context) };
    if (initial_state_build_result.failure !=
        SecondStageInitialStateBuildResult::Failure::None)
    {
        if (!options.quiet_mode)
        {
            const auto unselected_seed_failure{
                initial_state_build_result.failure ==
                SecondStageInitialStateBuildResult::Failure::UnselectedSeedUnavailable
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
                    std::string(unselected_seed_failure ?
                        "no-valid-unselected-neighbor-seed" :
                        "no-valid-seed") +
                    ", "
                "best_audit_objective=unavailable, final_uses_polish=unavailable, "
                "final_state_source=unavailable.");
        }
        return false;
    }
    LogSecondStageSeedSelections(initial_state_build_result.selection_record_list, options.quiet_mode);
    LogUnselectedSecondStageSeedSelections(
        context,
        initial_state_build_result.unselected_selection_record_list,
        options.quiet_mode);
    auto initial_state{ std::move(initial_state_build_result.state) };
    const auto graph_topology{
        BuildGraphTopology(context, initial_state, options.quiet_mode)
    };
    LogGraphTopology(graph_topology, options.quiet_mode);
    auto iteration_state{
        BuildIterationState(context, graph_topology, std::move(initial_state), options)
    };
    PerformanceCounters performance_counters{
        options.quiet_mode,
        context,
        iteration_state.solver_workspace_by_key
    };
    if (iteration_state.best_audit_state.best.has_value())
    {
        performance_counters.RecordFullStateMaterialization();
    }
    LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode);
    const auto progress_column_widths{ BuildProgressColumnWidths(atom_size) };
    LogProgressHeader(options.quiet_mode, progress_column_widths);

    for (std::size_t iter = 0; iter < kMaximumIterations; iter++)
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
                UsesPolish(iteration_state.previous_polish_provenance),
                FinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        auto iteration_result{
            RunIteration(
                context,
                graph_topology,
                options,
                iter + 1,
                iteration_state,
                performance_counters)
        };
        if (iteration_result.objective_domain_changed)
        {
            LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode, true);
        }
        LogAcceptedBacktrackingDiagnostics(options.quiet_mode, iteration_result.diagnostics);
        if (iteration_result.outcome != IterationOutcome::Accepted)
        {
            LogRejectedClusterDiagnostics(
                options.quiet_mode,
                iteration_result.diagnostics.rejected_cluster_diagnostic_list);
            LogIterationProgress(
                options.quiet_mode,
                progress_column_widths,
                iteration_result.progress);
            LogAllRejectedResolution(
                options.quiet_mode,
                iteration_result.diagnostics.rejected_cluster_partition,
                iteration_result.diagnostics.trust_region_radius_update,
                *iteration_result.all_rejected_resolution);
            if (iteration_result.outcome == IterationOutcome::Retry)
            {
                continue;
            }

            const auto final_state_selection{
                SelectFinalState(
                    iteration_state.previous_state,
                    iteration_state.previous_polish_provenance,
                    iteration_state.best_audit_state.best)
            };
            ApplyFitState(model_object, context, *final_state_selection.state);
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                GetAllRejectedResolutionText(*iteration_result.all_rejected_resolution),
                iteration_state.best_audit_state,
                UsesPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        LogRejectedClusterDiagnostics(
            options.quiet_mode,
            iteration_result.diagnostics.rejected_cluster_diagnostic_list);
        LogIterationProgress(
            options.quiet_mode,
            progress_column_widths,
            iteration_result.progress);

        if (iteration_result.audit_patience_exhausted)
        {
            const auto final_state_selection{
                SelectFinalState(
                    iteration_state.previous_state,
                    iteration_state.previous_polish_provenance,
                    iteration_state.best_audit_state.best)
            };
            ApplyFitState(
                model_object,
                context,
                *final_state_selection.state);
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "audit-patience",
                iteration_state.best_audit_state,
                UsesPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iteration_result.converged)
        {
            const auto accepted_offset_stats{
                SummarizeOffsets(iteration_state.previous_state)
            };
            ApplyFitState(
                model_object,
                context,
                iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.HasFailures())
            {
                LogTerminalFallback(
                    options.quiet_mode,
                    iteration_state.accepted_iteration_count,
                    iteration_state.terminal_failure_state.Summary(),
                    accepted_offset_stats);
            }
            else
            {
                LogConverged(
                    options.quiet_mode,
                    iteration_state.accepted_iteration_count,
                    iteration_result.transformed_change_stats,
                    accepted_offset_stats);
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "converged",
                iteration_state.best_audit_state,
                UsesPolish(iteration_state.previous_polish_provenance),
                FinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iter + 1 == kMaximumIterations)
        {
            const auto final_state_selection{
                SelectFinalState(
                    iteration_state.previous_state,
                    iteration_state.previous_polish_provenance,
                    iteration_state.best_audit_state.best)
            };
            ApplyFitState(
                model_object,
                context,
                *final_state_selection.state);
            LogMaximumIterations(
                options.quiet_mode,
                final_state_selection,
                iteration_state.terminal_failure_state.Summary(),
                SummarizeOffsets(*final_state_selection.state));
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                "maximum-iterations",
                iteration_state.best_audit_state,
                UsesPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }
    }
    return false;
}
} // namespace rhbm_gem::core
