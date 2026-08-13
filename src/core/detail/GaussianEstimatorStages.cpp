#include <cstddef>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingCandidateEvaluationOverlay.hpp"
#include "core/detail/LocalFittingResidualEvaluation.hpp"
#include "core/detail/LocalFittingObjective.hpp"
#include "core/detail/LocalFittingTerminalFailure.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "core/detail/LocalFittingIteration.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointOffset.hpp"
#include "core/detail/JointPolish.hpp"
#include "core/detail/LocalFittingObjectiveAttemptDiagnostic.hpp"
#include "core/detail/LocalFittingPerformanceCounters.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/SecondStageLocalFittingInitialization.hpp"
#include "core/detail/LocalFittingStateView.hpp"
#include "core/detail/TrustRegion.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/ReusableWeightedRidgeSolver.hpp"
#include "core/detail/ScopedEigenThreadCount.hpp"
#include "core/detail/SecondStageLocalFittingContext.hpp"
#include "data/detail/AtomClassifier.hpp"
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
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {

namespace {

using detail::FittedGaussianSnapshot;
using detail::BuildFittedGaussianSnapshot;
using detail::BuildSecondStageModelSnapshot;
using detail::LocalFittingCandidateEvaluationOverlay;
using detail::LocalFittingClusterKey;
using detail::LocalFittingClusterSolverWorkspace;
using detail::LocalFittingClusterSolverWorkspaceMap;
using detail::LocalFittingObjectiveAttemptDiagnostic;
using detail::LocalFittingObjectiveSampleRef;
using detail::LocalFittingPerformanceCounters;
using detail::LocalFittingPolishProvenance;
using detail::LocalFittingPreObjectiveFailureReason;
using detail::LocalFittingResidualBaseline;
using detail::LocalFittingState;
using detail::LocalFittingStatePatch;
using detail::LocalFittingStateView;
using detail::SecondStageAtomContext;
using detail::SecondStageLocalFittingContext;
using detail::SecondStageModelSnapshot;
using detail::ReusableWeightedRidgeSolver;
using detail::ScopedEigenThreadCount;

constexpr std::size_t kLocalFittingMaximumIterations{ 100 };
constexpr std::size_t kLocalFittingAuditPatience{ 3 };
constexpr bool kApplyLocalFittingBestIteration{ true };

using detail::JointOffsetSolveStatus;
using detail::IsJointOffsetSolveStationarityEligible;
using detail::GetJointOffsetSolveStatusText;
using detail::ResetLocalFittingClusterSolverWorkspace;
using detail::BuildSecondStageAdjustedResponseCache;
using detail::BuildSecondStageAdjustedSamples;
using detail::BuildLocalFittingResidualBaseline;
using detail::SummarizeLocalFittingTransformedChanges;
using detail::IsLocalFittingTransformedChangeConverged;
using detail::BuildLocalFittingTransformedEstimationList;
using detail::kLocalFittingChangePercentile;
using detail::LocalFittingObjectiveDomain;
using detail::LocalFittingClusterObjectiveState;
using detail::LocalFittingClusterObjectiveStateMap;
using detail::LocalFittingObjectiveByKey;
using detail::LocalFittingCombinedObjectiveCheck;
using detail::LocalFittingAuditedState;
using detail::LocalFittingBestAuditState;
using detail::BuildLocalFittingObjectiveDomain;
using detail::BuildLocalFittingObjectiveByKey;
using detail::EvaluateLocalFittingAuditObjective;
using detail::EvaluateLocalFittingCombinedObjective;
using detail::BuildInitialLocalFittingBestAuditState;
using detail::ResetLocalFittingBestAuditAfterObjectiveDomainChange;
using detail::TryUpdateLocalFittingBestAuditState;
using detail::ReconcileLocalFittingClusterObjectiveState;
using detail::PersistentTerminalFailureState;
using detail::PersistentTerminalFailureStateMap;
using detail::LocalFittingTerminalSummary;
using detail::UpdatePersistentTerminalFailureState;
using detail::ApplyTerminalFallbackClusters;
using detail::BuildEligibleLocalFittingActiveIndexList;
using detail::RejectedClusterDiagnostic;
using detail::PolishProgress;
using detail::CandidateSelection;
using detail::RejectCombinedCandidate;
using detail::SelectClusterCandidates;
using detail::TryBacktrackCombinedCandidate;
using detail::BuildStatePatch;
using detail::kLocalFittingFitRangeWeight;
using detail::kLocalFittingTailValidationWeight;
using detail::kLocalFittingOffsetPlausibilityPenaltyWeight;
using detail::BuildSecondStageLocalFittingContext;
using detail::StoreSecondStageNeighborCounts;
using detail::BuildInitialLocalFittingState;
using detail::GetSecondStageSeedSourceText;
using detail::SecondStageSeedSelectionRecord;
using detail::RunLocalFittingIteration;

struct LocalFittingFinalStateSelection
{
    const LocalFittingState * state{ nullptr };
    const LocalFittingPolishProvenance * polish_provenance{ nullptr };
    const LocalFittingAuditedState * audit_state{ nullptr };
    detail::LocalFittingFinalStateSource source{
        detail::LocalFittingFinalStateSource::Unavailable
    };
};

LocalFittingFinalStateSelection SelectLocalFittingFinalState(
    const LocalFittingState & latest_validated_state,
    const LocalFittingPolishProvenance & latest_validated_polish_provenance,
    const std::optional<LocalFittingAuditedState> & audited_state)
{
    const auto source{ detail::SelectLocalFittingFinalStateSource(
        kApplyLocalFittingBestIteration,
        true,
        audited_state.has_value()) };
    if (source == detail::LocalFittingFinalStateSource::BestAudit)
    {
        return LocalFittingFinalStateSelection{
            &audited_state->state,
            &audited_state->polish_provenance,
            &*audited_state,
            source
        };
    }
    return LocalFittingFinalStateSelection{
        &latest_validated_state,
        &latest_validated_polish_provenance,
        nullptr,
        source
    };
}

std::string_view GetLocalFittingFinalStateSourceText(detail::LocalFittingFinalStateSource source)
{
    switch (source)
    {
    case detail::LocalFittingFinalStateSource::BestAudit:
        return "best-audit";
    case detail::LocalFittingFinalStateSource::LatestValidated:
        return "latest-validated";
    case detail::LocalFittingFinalStateSource::Unavailable:
        return "unavailable";
    }
    return "unavailable";
}

bool UsesLocalFittingPolish(const LocalFittingPolishProvenance & provenance)
{
    return std::any_of(
        provenance.begin(),
        provenance.end(),
        [](char is_polished)
        {
            return is_polished != 0;
        });
}

struct LocalFittingOffsetStats
{
    std::size_t atom_count{ 0 };
    std::size_t finite_count{ 0 };
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
};

struct LocalFittingIterationProgress
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

using LocalFittingProgressColumnWidths = std::array<std::size_t, 6>;

detail::GraphTopology BuildGraphTopology(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state,
    const FitOptions & options)
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
        if (options.quiet_mode) return;
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
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(initial_state))
    };
    std::vector<std::optional<detail::TransformedModelInvariants>>
        selected_model_invariants;
    selected_model_invariants.reserve(model_snapshot.selected.size());
    for (const auto & model : model_snapshot.selected)
    {
        selected_model_invariants.emplace_back(
            detail::BuildTransformedModelInvariants(model));
    }
    std::vector<std::optional<detail::TransformedModelInvariants>>
        unselected_model_invariants;
    unselected_model_invariants.reserve(model_snapshot.unselected.size());
    for (const auto & model : model_snapshot.unselected)
    {
        unselected_model_invariants.emplace_back(
            detail::BuildTransformedModelInvariants(model));
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
                const auto & neighbor_sample{ *neighbor_iter };
                const auto neighbor_evaluation{ neighbor_sample.is_selected ?
                    evaluate_model(
                        selected_model_invariants.at(neighbor_sample.atom_index),
                        neighbor_sample.distance) :
                    evaluate_model(
                        unselected_model_invariants.at(neighbor_sample.atom_index),
                        neighbor_sample.distance) };
                const auto jacobian{
                    neighbor_evaluation.has_value() ?
                        neighbor_evaluation->jacobian : invalid_jacobian
                };
                if (neighbor_sample.is_selected)
                {
                    participant_list.emplace_back(
                        detail::GraphParticipant{
                            neighbor_sample.atom_index,
                            jacobian
                        });
                    continue;
                }
                const auto selected_group_id{
                    context.unselected_atom_list.at(
                        neighbor_sample.atom_index).selected_group_id
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
                LocalFittingObjectiveSampleRef{ atom_index, sample_index },
                participant_list);
            completed_work++;
            update_progress();
        }
    }

    std::vector<detail::GraphResidueKey>
        residue_key_by_atom_index;
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
    const FitOptions & options)
{
    if (options.quiet_mode || selection_record_list.empty()) return;

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
    const SecondStageLocalFittingContext & context,
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    const FitOptions & options)
{
    if (options.quiet_mode || selection_record_list.empty()) return;

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
    summary << "Unselected second-stage neighbor seeds = "
        << selection_record_list.size()
        << ", sources = group-median:" << group_median_count
        << ", global-median:" << global_median_count << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        const auto & contributor{
            context.unselected_atom_list.at(record.atom_index)
        };
        std::ostringstream detail_message;
        detail_message
            << "Unselected second-stage neighbor seed selection: serial ID = "
            << contributor.atom->GetSerialID()
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", seed A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

void LogGraphTopology(
    const detail::GraphTopology & topology,
    const FitOptions & options)
{
    if (options.quiet_mode) return;

    const auto & summary{ topology.summary };
    if (!summary.uses_weighted_graph)
    {
        Logger::Log(
            LogLevel::Warning,
            "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
        if (Logger::GetLogLevel() >= LogLevel::Debug)
        {
            Logger::Log(
                LogLevel::Debug,
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

void LogLocalFittingObjectiveDomain(
    const LocalFittingObjectiveDomain & domain,
    const FitOptions & options,
    bool is_terminal_reset = false)
{
    if (options.quiet_mode) return;
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
            "Reset second-stage objective domain" :
            "Initialize second-stage objective domain")
        << ": fit/tail/offset weights = "
        << kLocalFittingFitRangeWeight << "/"
        << kLocalFittingTailValidationWeight << "/"
        << kLocalFittingOffsetPlausibilityPenaltyWeight
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

void ApplyLocalFittingState(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & iteration_state)
{
    if (context.size() != iteration_state.size())
    {
        throw std::invalid_argument(
            "Local fitting context and state sizes are inconsistent.");
    }

    const auto model_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(iteration_state))
    };
    const auto adjusted_response_cache{
        BuildSecondStageAdjustedResponseCache(context, model_snapshot)
    };
    std::vector<LocalPotentialSampleList> adjusted_sampling_entries_list;
    adjusted_sampling_entries_list.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        adjusted_sampling_entries_list.emplace_back(
            BuildSecondStageAdjustedSamples(
                context,
                i,
                adjusted_response_cache.at(i)));
    }

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        analysis.ApplyAtomLocalSecondStageResult(
            *context.at(i).atom,
            iteration_state.at(i),
            std::move(adjusted_sampling_entries_list.at(i)));
    }
}

LocalFittingOffsetStats SummarizeLocalFittingOffsets(const LocalFittingState & state)
{
    LocalFittingOffsetStats stats;
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
    stats.percentile_absolute_offset = array_helper::ComputePercentile(absolute_offset_list, kLocalFittingChangePercentile);
    stats.maximum_absolute_offset = *std::max_element(absolute_offset_list.begin(), absolute_offset_list.end());
    return stats;
}

void AppendLocalFittingOffsetSummary(std::ostringstream & stream, const LocalFittingOffsetStats & stats)
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

void AppendLocalFittingAuditSummary(std::ostringstream & stream, const LocalFittingAuditedState & audited_state)
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
        << kLocalFittingTailValidationWeight;
}

void AppendLocalFittingTerminalSummary(std::ostringstream & stream, const LocalFittingTerminalSummary & summary)
{
    if (summary.suspicious_atom_count > 0)
    {
        stream << "; terminal suspicious rollback fallback clusters/atoms = "
            << summary.suspicious_cluster_count
            << "/" << summary.suspicious_atom_count;
    }
    if (summary.joint_offset_failure_atom_count > 0)
    {
        stream << "; terminal joint-offset failure fallback clusters/atoms = "
            << summary.joint_offset_failure_cluster_count
            << "/" << summary.joint_offset_failure_atom_count;
        if (!summary.joint_offset_failure_status_count.empty())
        {
            stream << ", statuses = ";
            bool is_first_status{ true };
            for (const auto & [status, count] :
                summary.joint_offset_failure_status_count)
            {
                if (!is_first_status) stream << ",";
                stream << GetJointOffsetSolveStatusText(status) << ":" << count;
                is_first_status = false;
            }
        }
    }
}

void LogLocalFittingTerminalFallback(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Completed local fitting after " << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    AppendLocalFittingOffsetSummary(warning_message, offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void FinishLocalFittingWithNoActiveAtoms(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingTerminalSummary & terminal_summary)
{
    ApplyLocalFittingState(model_object, context, state);
    const auto offset_stats{ SummarizeLocalFittingOffsets(state) };
    if (terminal_summary.AtomCount() > 0)
    {
        LogLocalFittingTerminalFallback(
            options,
            accepted_iteration_count,
            terminal_summary,
            offset_stats);
        return;
    }
    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(
            LogLevel::Info,
            "Skip 2nd-stage local atom fitting because no atoms are selected.");
    }
}

void AppendLocalFittingObjectiveBreakdown(
    std::ostringstream & stream,
    const std::optional<detail::LocalFittingObjectiveBreakdown> & breakdown)
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

std::string_view GetLocalFittingPreObjectiveFailureReasonText(
    LocalFittingPreObjectiveFailureReason reason)
{
    switch (reason)
    {
    case LocalFittingPreObjectiveFailureReason::None:
        return "none";
    case LocalFittingPreObjectiveFailureReason::InvalidModel:
        return "invalid-model";
    case LocalFittingPreObjectiveFailureReason::PreviousSharedOffsetProjectionOutsideTrustRegion:
        return "previous-shared-offset-projection-outside-trust-region";
    case LocalFittingPreObjectiveFailureReason::NoCandidateWithinTrustRegion:
        return "no-candidate-within-trust-region";
    }
    return "unknown";
}

void LogRejectedLocalFittingClusterDiagnostics(
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
            << "Rejected local fitting cluster diagnostics: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/" << cluster_diagnostic.key.back()
            << ", breakdown order = fit/tail-weighted/offset/total";
        Logger::Log(LogLevel::Debug, header.str());

        const auto & diagnostic{ cluster_diagnostic.attempt };
        std::ostringstream message;
        message
            << std::scientific << std::setprecision(2)
            << "  fixed-point effective damping = "
            << diagnostic.effective_damping
            << ", trust radius/step norm = "
            << diagnostic.trust_region_radius << "/";

        if (diagnostic.pre_objective_failure_reason !=
            LocalFittingPreObjectiveFailureReason::None)
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
                << GetLocalFittingPreObjectiveFailureReasonText(
                    diagnostic.pre_objective_failure_reason)
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
            << diagnostic.fit_sample_count << "/"
            << diagnostic.tail_sample_count
            << ", tail raw/weight = ";
        if (diagnostic.candidate_objective.has_value())
        {
            message << diagnostic.candidate_objective->tail_validation_loss;
        }
        else
        {
            message << "unavailable";
        }
        message << "/" << kLocalFittingTailValidationWeight;
        message << ", candidate = ";
        AppendLocalFittingObjectiveBreakdown(message, diagnostic.candidate_objective);
        message << ", previous = ";
        AppendLocalFittingObjectiveBreakdown(message, diagnostic.previous_objective);
        message << ", best = ";
        AppendLocalFittingObjectiveBreakdown(message, diagnostic.best_objective);
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

std::string_view GetAllRejectedResolutionText(
    detail::AllRejectedResolution resolution)
{
    switch (resolution)
    {
    case detail::AllRejectedResolution::Retry:
        return "retry";
    case detail::AllRejectedResolution::MaximumIterations:
        return "maximum-iterations";
    case detail::AllRejectedResolution::BacktrackingExhausted:
        return "all-rejected-backtracking-exhausted";
    case detail::AllRejectedResolution::MinimumRadius:
        return "all-rejected-minimum-radius";
    case detail::AllRejectedResolution::NoRetryProgress:
        return "all-rejected-no-retry-progress";
    }
    return "all-rejected-no-retry-progress";
}

void LogAllRejectedResolution(
    const FitOptions & options,
    const detail::RejectedClusterPartition & partition,
    const detail::TrustRegionRadiusUpdate & radius_update,
    detail::AllRejectedResolution resolution)
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

void LogAcceptedLocalFittingBacktrackingDiagnostics(
    const FitOptions & options,
    const CandidateSelection & selection)
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
    if (!has_local_backtracking &&
        selection.combined_backtracking_trial_count <= 1)
    {
        return;
    }
    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic :
        selection.accepted_cluster_diagnostic_list)
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
            << "Accepted local fitting objective backtracking: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/"
            << cluster_diagnostic.key.back()
            << ", trials/factor = "
            << diagnostic.backtracking_trial_count << "/";
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

std::string FormatLocalFittingProgressMaximum(
    const std::optional<double> & value)
{
    if (!value.has_value()) return "-";
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << *value;
    return stream.str();
}

std::optional<double> SummarizeLocalFittingProgressMaximum(
    const std::vector<double> & maximum_list)
{
    if (maximum_list.empty()) return std::nullopt;
    return *std::max_element(maximum_list.begin(), maximum_list.end());
}

constexpr std::array<std::string_view, 6> kLocalFittingProgressHeaderList{
    "Try/Acc",
    "Atom A/T",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/R"
};

std::string FormatLocalFittingProgressRow(
    const LocalFittingProgressColumnWidths & column_widths,
    const std::array<std::string, 6> & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream
            << std::left
            << std::setw(static_cast<int>(column_widths.at(i)))
            << cell_list.at(i);
    }
    return stream.str();
}

LocalFittingProgressColumnWidths BuildLocalFittingProgressColumnWidths(
    std::size_t atom_size)
{
    const auto maximum_iteration_text{ std::to_string(kLocalFittingMaximumIterations) };
    const auto maximum_atom_text{ std::to_string(atom_size) };
    const auto maximum_change_text{
        FormatLocalFittingProgressMaximum(std::numeric_limits<double>::max())
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

    LocalFittingProgressColumnWidths column_widths;
    for (std::size_t i = 0; i < column_widths.size(); i++)
    {
        column_widths.at(i) = std::max(
            kLocalFittingProgressHeaderList.at(i).size(),
            maximum_cell_list.at(i).size());
    }
    return column_widths;
}

void LogLocalFittingProgressHeader(
    const FitOptions & options,
    const LocalFittingProgressColumnWidths & column_widths)
{
    if (options.quiet_mode) return;
    std::array<std::string, 6> header_list;
    for (std::size_t i = 0; i < header_list.size(); i++)
    {
        header_list.at(i) = kLocalFittingProgressHeaderList.at(i);
    }
    Logger::Log(
        LogLevel::Info,
        FormatLocalFittingProgressRow(column_widths, header_list));
}

void LogLocalFittingIterationProgress(
    const FitOptions & options,
    const LocalFittingProgressColumnWidths & column_widths,
    const LocalFittingIterationProgress & progress)
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
        FormatLocalFittingProgressMaximum(progress.accepted_maximum_transformed_change) + "/" +
            FormatLocalFittingProgressMaximum(progress.raw_maximum_transformed_change)
    };
    Logger::ProgressLine(FormatLocalFittingProgressRow(column_widths, cell_list));
}

void LogLocalFittingConverged(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

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
    AppendLocalFittingOffsetSummary(message, offset_stats);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogLocalFittingMaximumIterations(
    const FitOptions & options,
    detail::LocalFittingFinalStateSource final_state_source,
    const LocalFittingAuditedState * applied_audit_state,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & applied_offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    if (final_state_source == detail::LocalFittingFinalStateSource::BestAudit &&
        applied_audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendLocalFittingAuditSummary(warning_message, *applied_audit_state);
    }
    else if (final_state_source == detail::LocalFittingFinalStateSource::LatestValidated)
    {
        warning_message << "; applying latest validated state";
    }
    else
    {
        warning_message << "; no validated state is available";
    }
    AppendLocalFittingOffsetSummary(warning_message, applied_offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageLocalFittingSummary(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    std::string_view stop_reason,
    const LocalFittingBestAuditState & best_audit_state,
    std::optional<bool> final_uses_polish,
    detail::LocalFittingFinalStateSource final_state_source)
{
    if (options.quiet_mode) return;

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
        << GetLocalFittingFinalStateSourceText(final_state_source)
        << ".";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace

bool RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options)
{
    auto context{ BuildSecondStageLocalFittingContext(model_object, options) };
    StoreSecondStageNeighborCounts(model_object, context);
    const auto atom_size{ context.size() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
        Logger::Log(
            LogLevel::Info,
            kApplyLocalFittingBestIteration ?
                "Second-stage best-iteration application: enabled." :
                "Second-stage best-iteration application: disabled.");
    }

    bool unselected_seed_failure{ false };
    auto initial_state_build_result{
        BuildInitialLocalFittingState(
            context,
            unselected_seed_failure)
    };
    if (!initial_state_build_result.has_value())
    {
        if (!options.quiet_mode)
        {
            Logger::Log(
                LogLevel::Warning,
                unselected_seed_failure ?
                    "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                    "is available for every unselected neighbor atom." :
                    "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                    "is available for every selected atom.");
            Logger::Log(
                LogLevel::Info,
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
    LogSecondStageSeedSelections(initial_state_build_result->selection_record_list, options);
    LogUnselectedSecondStageSeedSelections(
        context,
        initial_state_build_result->unselected_selection_record_list,
        options);
    auto previous_state{ std::move(initial_state_build_result->state) };
    LocalFittingPolishProvenance previous_polish_provenance(atom_size, 0);
    const auto graph_topology{
        BuildGraphTopology(context, previous_state, options)
    };
    LogGraphTopology(graph_topology, options);
    std::vector<char> terminal_fallback_atom_mask(atom_size, 0);
    auto active_index_list{
        BuildEligibleLocalFittingActiveIndexList(terminal_fallback_atom_mask)
    };
    auto graph_partition{
        detail::BuildGraphPartition(graph_topology, active_index_list)
    };
    auto cluster_key_list{
        detail::BuildGraphClusterKeyList(graph_partition)
    };
    LocalFittingClusterSolverWorkspaceMap solver_workspace_by_key;
    ResetLocalFittingClusterSolverWorkspace(
        cluster_key_list,
        solver_workspace_by_key);
    LocalFittingPerformanceCounters performance_counters{
        options,
        solver_workspace_by_key
    };
    const auto cached_sample_count{
        std::accumulate(
            context.begin(),
            context.end(),
            std::size_t{ 0 },
            [](std::size_t count, const SecondStageAtomContext & atom_context)
            {
                return count + atom_context.raw_sampling_entries.size();
            })
    };
    auto objective_domain{
        BuildLocalFittingObjectiveDomain(
            context,
            previous_state,
            graph_partition,
            options)
    };
    LogLocalFittingObjectiveDomain(objective_domain, options);
    auto best_audit_state{
        BuildInitialLocalFittingBestAuditState(
            context,
            previous_state,
            previous_polish_provenance,
            std::nullopt,
            objective_domain)
    };
    if (best_audit_state.best.has_value())
    {
        performance_counters.full_state_materialization_count++;
    }

    std::vector<char> rollback_atom_mask(atom_size, 0);
    PersistentTerminalFailureStateMap persistent_terminal_failure_state_by_key;
    LocalFittingTerminalSummary terminal_summary;
    LocalFittingClusterObjectiveStateMap cluster_objective_state;
    detail::TrustRegionStateSet trust_region_state;
    const auto progress_column_widths{ BuildLocalFittingProgressColumnWidths(atom_size) };
    LogLocalFittingProgressHeader(options, progress_column_widths);

    std::size_t accepted_iteration_count{ 0 };
    std::size_t audit_patience_count{ 0 };
    std::vector<LocalFittingClusterKey> unchanged_state_exhausted_key_list;
    for (std::size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        if (active_index_list.empty())
        {
            FinishLocalFittingWithNoActiveAtoms(
                model_object,
                context,
                previous_state,
                options,
                accepted_iteration_count,
                terminal_summary);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "terminal-isolation",
                best_audit_state,
                UsesLocalFittingPolish(previous_polish_provenance),
                detail::LocalFittingFinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        const auto previous_model_snapshot{
            BuildSecondStageModelSnapshot(
                context,
                BuildFittedGaussianSnapshot(previous_state))
        };
        const auto residual_baseline{
            BuildLocalFittingResidualBaseline(
                context,
                previous_state,
                previous_model_snapshot)
        };
        performance_counters.gaussian_cache_miss_count += cached_sample_count;

        const auto previous_objective_by_key{
            BuildLocalFittingObjectiveByKey(
                context,
                previous_state,
                graph_partition,
                objective_domain,
                residual_baseline)
        };
        ReconcileLocalFittingClusterObjectiveState(
            graph_partition,
            previous_objective_by_key,
            cluster_objective_state);
        trust_region_state.Reconcile(cluster_key_list);

        std::vector<double> joint_offset_ridge_multiplier_list(atom_size, 1.0);
        for (std::size_t atom_index = 0; atom_index < atom_size; atom_index++)
        {
            if (rollback_atom_mask.at(atom_index) == 0) continue;
            joint_offset_ridge_multiplier_list.at(atom_index) =
                detail::kSuspiciousJointOffsetRidgeMultiplier;
        }

        const auto iteration_phase_start{ std::chrono::steady_clock::now() };
        auto iteration_result{
            RunLocalFittingIteration(
                context,
                cluster_key_list,
                previous_state,
                options,
                joint_offset_ridge_multiplier_list,
                solver_workspace_by_key)
        };
        performance_counters.iteration_phase_milliseconds +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() -
                iteration_phase_start).count();
        performance_counters.gaussian_cache_hit_count += cached_sample_count;
        const auto current_health_by_key{ std::move(iteration_result.health_by_key) };
        rollback_atom_mask = std::move(iteration_result.rollback_atom_mask);
        const auto iteration_suspicious_atom_count{
            static_cast<std::size_t>(std::count_if(
                rollback_atom_mask.begin(),
                rollback_atom_mask.end(),
                [](char is_suspicious)
                {
                    return is_suspicious != 0;
                }))
        };
        const auto has_suspicious_offset_fallback{ iteration_suspicious_atom_count > 0 };

        std::size_t stationarity_ineligible_cluster_count{ 0 };
        std::vector<LocalFittingClusterKey> polish_eligible_key_list;
        for (const auto & [key, health] : current_health_by_key)
        {
            if (!IsJointOffsetSolveStationarityEligible(
                    health.joint_offset_status) ||
                !health.is_refit_stationarity_eligible)
            {
                stationarity_ineligible_cluster_count++;
                continue;
            }
            const auto contains_suspicious_atom{
                std::any_of(
                    key.begin(),
                    key.end(),
                    [&](std::size_t atom_index)
                    {
                        return rollback_atom_mask.at(atom_index) != 0;
                    })
            };
            if (!contains_suspicious_atom)
            {
                polish_eligible_key_list.emplace_back(key);
            }
        }
        const auto raw_state{ std::move(iteration_result.state) };
        const auto raw_fixed_point_change_summary{
            SummarizeLocalFittingTransformedChanges(raw_state, previous_state, active_index_list)
        };
        const auto previous_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(previous_state)
        };
        const auto raw_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(raw_state)
        };

        auto working_cluster_objective_state{ cluster_objective_state };
        const auto candidate_phase_start{ std::chrono::steady_clock::now() };
        auto selection{
            SelectClusterCandidates(
                context,
                previous_model_snapshot,
                residual_baseline,
                graph_partition,
                polish_eligible_key_list,
                previous_state,
                previous_polish_provenance,
                raw_state,
                previous_transformed_estimation_list,
                raw_transformed_estimation_list,
                rollback_atom_mask,
                joint_offset_ridge_multiplier_list,
                unchanged_state_exhausted_key_list,
                objective_domain,
                previous_objective_by_key,
                working_cluster_objective_state,
                trust_region_state,
                solver_workspace_by_key,
                options.thread_size,
                performance_counters)
        };
        performance_counters.candidate_phase_milliseconds +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() -
                candidate_phase_start).count();
        performance_counters.full_state_materialization_count++;

        const auto needs_combined_objective_guard{
            graph_partition.boundary_sample_count > 0 && !selection.accepted_key_list.empty()
        };
        const auto combined_changed_key_list{ selection.accepted_key_list };
        std::optional<detail::LocalFittingObjectiveBreakdown>
            previous_audit_objective;
        LocalFittingCombinedObjectiveCheck combined_check;
        if (needs_combined_objective_guard)
        {
            previous_audit_objective =
                EvaluateLocalFittingAuditObjective(
                    context,
                    previous_state,
                    objective_domain,
                    residual_baseline);
            LocalFittingClusterKey changed_atom_index_list;
            for (const auto & key : combined_changed_key_list)
            {
                changed_atom_index_list.insert(
                    changed_atom_index_list.end(),
                    key.begin(),
                    key.end());
            }
            const auto combined_patch{
                BuildStatePatch(
                    selection.assembled_state,
                    std::move(changed_atom_index_list))
            };
            const LocalFittingStateView combined_state_view{
                previous_state,
                combined_patch
            };
            const LocalFittingCandidateEvaluationOverlay combined_overlay{
                context,
                previous_model_snapshot,
                residual_baseline,
                combined_state_view,
                combined_patch
            };
            const auto affected_sample_ref_list{
                detail::BuildGraphAffectedSampleUnion(
                    graph_partition,
                    combined_changed_key_list)
            };
            combined_check = EvaluateLocalFittingCombinedObjective(
                context,
                previous_state,
                combined_state_view,
                residual_baseline,
                combined_overlay,
                combined_patch.atom_index_list,
                affected_sample_ref_list,
                objective_domain,
                best_audit_state,
                previous_audit_objective,
                &performance_counters);
            selection.combined_backtracking_objective =
                combined_check.candidate_objective;
        }
        auto combined_objective_accepted{
            !needs_combined_objective_guard || combined_check.accepted
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
                    previous_polish_provenance,
                    objective_domain,
                    previous_objective_by_key,
                    previous_audit_objective,
                    best_audit_state,
                    cluster_objective_state,
                    working_cluster_objective_state,
                    selection,
                    performance_counters);
        }
        if (!combined_objective_accepted)
        {
            RejectCombinedCandidate(
                previous_state,
                previous_polish_provenance,
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
            cluster_objective_state = std::move(working_cluster_objective_state);
        }
        LogAcceptedLocalFittingBacktrackingDiagnostics(options, selection);

        auto assembled_state{ std::move(selection.assembled_state) };
        auto assembled_polish_provenance{ std::move(selection.assembled_polish_provenance) };
        const auto terminal_failure_by_key{
            UpdatePersistentTerminalFailureState(
                selection.accepted_key_list,
                rollback_atom_mask,
                current_health_by_key,
                assembled_state,
                previous_state,
                persistent_terminal_failure_state_by_key)
        };

        std::vector<LocalFittingClusterKey> terminal_key_list;
        for (const auto & [key, reason] : terminal_failure_by_key)
        {
            terminal_key_list.emplace_back(key);
            if (std::holds_alternative<detail::PersistentSuspiciousRollbackReason>(reason))
            {
                terminal_summary.suspicious_cluster_count++;
                terminal_summary.suspicious_atom_count += key.size();
                continue;
            }
            const auto status{ std::get<JointOffsetSolveStatus>(reason) };
            terminal_summary.joint_offset_failure_cluster_count++;
            terminal_summary.joint_offset_failure_atom_count += key.size();
            terminal_summary.joint_offset_failure_status_count[status]++;
        }
        auto objective_domain_changed{ false };
        ApplyTerminalFallbackClusters(
            terminal_key_list,
            previous_state,
            previous_polish_provenance,
            terminal_fallback_atom_mask,
            assembled_state,
            assembled_polish_provenance);
        if (!terminal_key_list.empty())
        {
            for (const auto & key : terminal_key_list)
            {
                for (const auto atom_index : key)
                {
                    rollback_atom_mask.at(atom_index) = 0;
                }
            }
            auto remaining_active_index_list{
                BuildEligibleLocalFittingActiveIndexList(terminal_fallback_atom_mask)
            };
            if (!remaining_active_index_list.empty())
            {
                auto remaining_graph_partition{
                    detail::BuildGraphPartition(
                        graph_topology,
                        remaining_active_index_list)
                };
                objective_domain = BuildLocalFittingObjectiveDomain(
                    context,
                    assembled_state,
                    remaining_graph_partition,
                    options);
                LogLocalFittingObjectiveDomain(
                    objective_domain,
                    options,
                    true);
                cluster_objective_state.clear();
                const auto assembled_model_snapshot{
                    BuildSecondStageModelSnapshot(
                        context,
                        BuildFittedGaussianSnapshot(assembled_state))
                };
                const auto remaining_objective_by_key{
                    BuildLocalFittingObjectiveByKey(
                        context,
                        assembled_state,
                        remaining_graph_partition,
                        objective_domain,
                        assembled_model_snapshot)
                };
                ReconcileLocalFittingClusterObjectiveState(
                    remaining_graph_partition,
                    remaining_objective_by_key,
                    cluster_objective_state);
                ResetLocalFittingBestAuditAfterObjectiveDomainChange(
                    context,
                    assembled_state,
                    assembled_polish_provenance,
                    accepted_iteration_count + 1,
                    objective_domain,
                    best_audit_state);
                active_index_list = std::move(remaining_active_index_list);
                cluster_key_list =
                    detail::BuildGraphClusterKeyList(
                        remaining_graph_partition);
                ResetLocalFittingClusterSolverWorkspace(
                    cluster_key_list,
                    solver_workspace_by_key);
                graph_partition = std::move(remaining_graph_partition);
                objective_domain_changed = true;
            }
            else
            {
                active_index_list.clear();
            }
        }

        const auto trust_region_iteration_update{
            trust_region_state.UpdateAfterIteration(
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
        const auto terminal_atom_count{ terminal_summary.AtomCount() };
        LocalFittingIterationProgress progress{
            iter + 1,
            accepted_iteration_count,
            atom_size - terminal_atom_count,
            terminal_atom_count,
            selection.accepted_key_list.size(),
            selection.rejected_key_list.size(),
            selection.polish_progress,
            iteration_suspicious_atom_count,
            std::nullopt,
            SummarizeLocalFittingProgressMaximum(raw_fixed_point_change_summary.maximum_list)
        };
        if (selection.accepted_key_list.empty())
        {
            const auto all_rejected_resolution{
                detail::ResolveAllRejected(
                    iter + 1 >= kLocalFittingMaximumIterations,
                    rejected_cluster_partition,
                    trust_region_radius_update)
            };
            LogRejectedLocalFittingClusterDiagnostics(options, selection.rejected_cluster_diagnostic_list);
            LogLocalFittingIterationProgress(options, progress_column_widths, progress);
            LogAllRejectedResolution(
                options,
                rejected_cluster_partition,
                trust_region_radius_update,
                all_rejected_resolution);
            if (all_rejected_resolution == detail::AllRejectedResolution::Retry)
            {
                for (const auto & key : rejected_cluster_partition.exhausted_key_list)
                {
                    if (std::find(
                            unchanged_state_exhausted_key_list.begin(),
                            unchanged_state_exhausted_key_list.end(),
                            key) == unchanged_state_exhausted_key_list.end())
                    {
                        unchanged_state_exhausted_key_list.emplace_back(key);
                    }
                }
                continue;
            }

            const auto final_state_selection{
                SelectLocalFittingFinalState(
                    previous_state,
                    previous_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(
                model_object,
                context,
                *final_state_selection.state);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                GetAllRejectedResolutionText(all_rejected_resolution),
                best_audit_state,
                UsesLocalFittingPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(
                model_object,
                options,
                FittingStage::Second);
            return true;
        }

        const auto transformed_change_summary{
            SummarizeLocalFittingTransformedChanges(assembled_state, previous_state, active_index_list)
        };

        accepted_iteration_count++;
        const auto improved_best_audit{
            TryUpdateLocalFittingBestAuditState(
                context,
                assembled_state,
                assembled_polish_provenance,
                accepted_iteration_count,
                objective_domain,
                best_audit_state,
                objective_domain_changed ?
                    std::nullopt : selection.combined_backtracking_objective)
        };
        if (improved_best_audit)
        {
            performance_counters.full_state_materialization_count++;
        }
        audit_patience_count = objective_domain_changed ? 0 :
            detail::AdvanceLocalFittingAuditPatience(
                audit_patience_count,
                improved_best_audit,
                !selection.rejected_key_list.empty() &&
                    !trust_region_radius_update.changed_key_list.empty());
        progress.accepted_iteration_count = accepted_iteration_count;
        progress.accepted_maximum_transformed_change = SummarizeLocalFittingProgressMaximum(transformed_change_summary.maximum_list);
        LogRejectedLocalFittingClusterDiagnostics(options, selection.rejected_cluster_diagnostic_list);
        LogLocalFittingIterationProgress(options, progress_column_widths, progress);

        if (audit_patience_count >= kLocalFittingAuditPatience)
        {
            const auto final_state_selection{
                SelectLocalFittingFinalState(
                    assembled_state,
                    assembled_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(
                model_object,
                context,
                *final_state_selection.state);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "audit-patience",
                best_audit_state,
                UsesLocalFittingPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(
                model_object,
                options,
                FittingStage::Second);
            return true;
        }

        const auto converged{
            stationarity_ineligible_cluster_count == 0 &&
            !has_suspicious_offset_fallback &&
            selection.rejected_key_list.empty() &&
            IsLocalFittingTransformedChangeConverged(
                transformed_change_summary.percentile_stats,
                transformed_change_summary.maximum_list) &&
            IsLocalFittingTransformedChangeConverged(
                raw_fixed_point_change_summary.percentile_stats,
                raw_fixed_point_change_summary.maximum_list)
        };
        if (converged)
        {
            const auto accepted_offset_stats{
                SummarizeLocalFittingOffsets(assembled_state)
            };
            ApplyLocalFittingState(model_object, context, assembled_state);
            if (terminal_summary.AtomCount() > 0)
            {
                LogLocalFittingTerminalFallback(
                    options,
                    accepted_iteration_count,
                    terminal_summary,
                    accepted_offset_stats);
            }
            else
            {
                LogLocalFittingConverged(
                    options,
                    accepted_iteration_count,
                    transformed_change_summary.percentile_stats,
                    accepted_offset_stats);
            }
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "converged",
                best_audit_state,
                UsesLocalFittingPolish(assembled_polish_provenance),
                detail::LocalFittingFinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            const auto final_state_selection{
                SelectLocalFittingFinalState(
                    assembled_state,
                    assembled_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(
                model_object,
                context,
                *final_state_selection.state);
            LogLocalFittingMaximumIterations(
                options,
                final_state_selection.source,
                final_state_selection.audit_state,
                terminal_summary,
                SummarizeLocalFittingOffsets(*final_state_selection.state));
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "maximum-iterations",
                best_audit_state,
                UsesLocalFittingPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }
        unchanged_state_exhausted_key_list.clear();
        previous_state = std::move(assembled_state);
        previous_polish_provenance = std::move(assembled_polish_provenance);
    }
    return false;
}
} // namespace rhbm_gem::core
