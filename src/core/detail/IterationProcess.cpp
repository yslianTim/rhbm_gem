#include "core/detail/IterationProcess.hpp"

#include "core/detail/CandidateSelection.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ostream>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>

#include "data/detail/AtomClassifier.hpp"

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core::detail {

namespace {

constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{
    2.0 * kNeighborContributionDistanceMax
};

struct SecondStageSeedSelectionRecord
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

struct UnselectedSecondStageSeedSelectionRecord
{
    int atom_serial_id{ 0 };
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3D selected_model{};
};

struct SecondStageInitialStateBuildResult
{
    FitState state{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
    std::vector<UnselectedSecondStageSeedSelectionRecord>
        unselected_selection_record_list{};
    enum class Failure
    {
        None,
        SelectedSeedUnavailable,
        UnselectedSeedUnavailable
    } failure{ Failure::None };
};

struct QuarantineSummary
{
    std::size_t entered_target_count{ 0 };
    std::size_t released_target_count{ 0 };
    std::size_t failed_probation_count{ 0 };
};

struct QuarantineState
{
    QuarantineFailureStateMap state_by_target{};
    std::vector<QuarantineTarget> probation_target_list{};
    QuarantineSummary summary{};
    bool force_probation{ false };
    QuarantineState() = default;

    explicit QuarantineState(std::size_t atom_count)
        : m_atom_count(atom_count)
    {
    }

    SuspiciousBlockActivity BeginIteration(std::size_t accepted_iteration_count);
    SuspiciousBlockActivity BuildFinalActivity() const;
    QuarantineStateTransition UpdateAfterIteration(
        const SecondStageContext & context,
        const std::vector<ClusterKey> & accepted_key_list,
        const SuspiciousBlockActivity & block_activity,
        std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
        const ClusterHealthMap & health_by_key,
        FitState & assembled_state,
        const FitState & previous_state,
        const PolishProvenance & previous_polish_provenance,
        PolishProvenance & assembled_polish_provenance,
        std::size_t accepted_iteration_count);
    void NotifyTopologyChanged() { force_probation = true; }
    std::size_t AtomCount() const;
    std::size_t TargetCount() const
    {
        return static_cast<std::size_t>(std::ranges::count_if(
            state_by_target,
            [](const auto & entry) { return entry.second.quarantined; }));
    }
    bool HasFailures() const { return TargetCount() != 0; }

private:
    std::size_t m_atom_count{ 0 };
};

constexpr std::size_t kMaximumIterations{ 100 };
constexpr std::size_t kAuditPatience{ 3 };

struct IterationDiagnostics
{
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<BoundaryComponentReconciliationDiagnostic> boundary_reconciliation_diagnostic_list{};
    TrustRegionIterationUpdate trust_region_update{};
};

struct IterationState
{
    FitState previous_state{};
    FitState topology_reference_state{};
    PolishProvenance previous_polish_provenance{};
    SuspiciousUpdateMask rollback_atom_mask{};
    std::vector<std::size_t> active_index_list{};
    CouplingGraphPartition graph_partition{};
    ClusterSolverWorkspaceMap solver_workspace_by_key{};
    BoundaryJointCorrectionWorkspaceMap boundary_joint_correction_workspace_by_key{};
    ObjectiveDomain objective_domain{};
    BestAuditState best_audit_state{};
    QuarantineState quarantine_state{};
    ClusterObjectiveStateMap cluster_objective_state{};
    TrustRegionStateSet trust_region_state{};
    std::vector<ClusterKey> unchanged_state_exhausted_key_list{};
    std::size_t accepted_iteration_count{ 0 };
    std::size_t accepted_iterations_since_topology_rebuild{ 0 };
    std::size_t audit_patience_count{ 0 };
};

struct IterationProgress
{
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t quarantine_atom_count{ 0 };
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

struct RawIterationResult
{
    FitState state{};
    SuspiciousUpdateMask rollback_atom_mask{};
    SuspiciousBlockActivity block_activity{};
    SuspiciousBlockActivity failure_block_activity{};
    std::vector<SuspiciousGaussianAssessment> assessment_by_atom{};
    std::vector<std::optional<RHBMEstimationStatus>> local_refit_status_by_atom{};
    ClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_stationarity_eligible{ false };
    bool is_boundary_correction_eligible{ false };
    bool shape_fixed{ false };
    SuspiciousGaussianAssessment assessment{};
    std::optional<RHBMEstimationStatus> attempted_refit_status{};
};

static std::optional<SecondStageSeedSelection>
SelectValidSecondStageSeedCandidate(
    SecondStageSeedSource source,
    const std::optional<GaussianModel3DWithUncertainty> & candidate)
{
    if (!candidate.has_value() || !IsValidSecondStageGaussianModel(candidate->GetModel()))
    {
        return std::nullopt;
    }
    return SecondStageSeedSelection{ source, *candidate };
}

} // namespace

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(const SecondStageSeedCandidates & candidates)
{
    if (const auto selected{
            SelectValidSecondStageSeedCandidate(
                SecondStageSeedSource::GroupPosterior,
                candidates.group_posterior) })
    {
        return selected;
    }
    if (const auto selected{
            SelectValidSecondStageSeedCandidate(
                SecondStageSeedSource::GroupPrior,
                candidates.group_prior) })
    {
        return selected;
    }
    if (const auto selected{
            SelectValidSecondStageSeedCandidate(
                SecondStageSeedSource::GroupMedian,
                candidates.group_median) })
    {
        return selected;
    }
    return SelectValidSecondStageSeedCandidate(
        SecondStageSeedSource::GlobalMedian,
        candidates.global_median);
}

AdaptiveTopologyRebuildDecision EvaluateAdaptiveTopologyRebuildTrigger(
    const FitState & accepted_state,
    const FitState & topology_reference_state,
    const std::vector<std::size_t> & active_index_list,
    std::size_t accepted_iterations_since_rebuild)
{
    const auto drift_summary{
        SummarizeTransformedChanges(
            accepted_state,
            topology_reference_state,
            active_index_list)
    };
    const auto maximum_transformed_drift{
        GetMaximumTransformedChange(drift_summary)
    };
    if (maximum_transformed_drift >= kAdaptiveTopologyRebuildDriftThreshold)
    {
        return AdaptiveTopologyRebuildDecision{
            AdaptiveTopologyRebuildTrigger::Drift,
            maximum_transformed_drift
        };
    }
    if (accepted_iterations_since_rebuild >= kAdaptiveTopologyRebuildAcceptedIterationInterval)
    {
        return AdaptiveTopologyRebuildDecision{
            AdaptiveTopologyRebuildTrigger::Interval,
            maximum_transformed_drift
        };
    }
    return AdaptiveTopologyRebuildDecision{
        AdaptiveTopologyRebuildTrigger::None,
        maximum_transformed_drift
    };
}

ConvergenceSafeguardPredicates EvaluateConvergenceSafeguardPredicates(
    bool stationarity_eligible,
    const TransformedChangeSummary & accepted_change,
    const TransformedChangeSummary & raw_change)
{
    return ConvergenceSafeguardPredicates{
        stationarity_eligible,
        ConvergenceChangePredicates{
            IsTransformedPercentileConverged(accepted_change),
            IsTransformedMaximumConverged(accepted_change)
        },
        ConvergenceChangePredicates{
            IsTransformedPercentileConverged(raw_change),
            IsTransformedMaximumConverged(raw_change)
        }
    };
}

TransformedChangeIndexListByParameter BuildActiveBlockChangeIndexLists(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousBlockActivity & block_activity)
{
    TransformedChangeIndexListByParameter result;
    for (const auto atom_index : atom_index_list)
    {
        if (block_activity.HasActiveShape(atom_index))
        {
            result.at(kLogPeakHeightChangeIndex).emplace_back(atom_index);
            result.at(kLogWidthChangeIndex).emplace_back(atom_index);
        }
        if (block_activity.HasActiveOffset(atom_index))
        {
            result.at(kOffsetToPeakRatioChangeIndex).emplace_back(atom_index);
        }
    }
    return result;
}

ConvergenceStationarityAudit EvaluateConvergenceStationarityAudit(
    const ClusterHealthMap & health_by_key)
{
    ConvergenceStationarityAudit result;
    result.active_block_eligible = std::ranges::all_of(
        health_by_key | std::views::values,
        &ClusterHealth::is_active_block_stationarity_eligible);
    result.full_cluster_eligible = AreClustersStationarityEligible(health_by_key);
    for (const auto & health : health_by_key | std::views::values)
    {
        if (!health.is_active_block_stationarity_eligible)
        {
            result.active_block_ineligible_cluster_count++;
        }
        if (!health.is_refit_stationarity_eligible)
        {
            result.refit_ineligible_cluster_count++;
        }
        if (health.joint_offset_status != JointOffsetSolveStatus::Converged)
        {
            if (IsJointOffsetSolveHardFailure(health.joint_offset_status))
            {
                result.hard_joint_failure_cluster_count++;
            }
            else
            {
                result.soft_joint_nonconverged_cluster_count++;
            }
        }
    }
    return result;
}

namespace {

static SecondStageContext BuildSecondStageContext(
    const ModelObject & model_object,
    const FitOptions & options)
{
    SecondStageContext context;
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.selected_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.selected_atom_list.emplace_back(AtomContext{ atom });
    }
    std::unordered_map<GroupKey, std::size_t> selected_group_id_by_key;
    selected_group_id_by_key.reserve(context.size());
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        atom_index_map.emplace(context.at(i).atom, i);
    }
    std::unordered_map<const AtomObject *, std::size_t> unselected_atom_index_map;
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto group_key{ data_internal::GetGroupKey(atom) };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        atom_context.raw_sampling_entries = local_view.GetRawSamplingEntries(false);
        atom_context.initial_result = local_view.GetGaussianResult(FittingStage::Second);
        atom_context.group_prior = analysis_view.FindAtomGroupPriorWithUncertainty(FittingStage::Second, *atom);
        atom_context.alpha_r = local_view.GetAlphaR(FittingStage::Second);
        atom_context.refit_design_template = BuildLocalGaussianDesignTemplate(
            atom_context.raw_sampling_entries,
            options.distance_min,
            options.distance_max);
        auto [group_iter, inserted]{
            selected_group_id_by_key.emplace(group_key, context.selected_atom_index_list_by_group.size())
        };
        if (inserted)
        {
            context.selected_atom_index_list_by_group.emplace_back();
        }
        atom_context.group_id = group_iter->second;
        context.selected_atom_index_list_by_group.at(atom_context.group_id).emplace_back(atom_index);
    }

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };
        std::unordered_set<const AtomObject *> neighbor_atom_set;

        atom_context.neighbor_atom_sample_offset_list.reserve(atom_context.raw_sampling_entries.size() + 1);
        atom_context.neighbor_atom_sample_offset_list.emplace_back(0);
        atom_context.neighbor_atom_sample_list.reserve(
            atom_context.raw_sampling_entries.size() * neighbor_atom_list.size());
        for (std::size_t sample_index = 0; sample_index < atom_context.raw_sampling_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            for (auto * neighbor_atom : neighbor_atom_list)
            {
                if (options.exclude_hydrogen && neighbor_atom->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(sample.point.position, neighbor_atom->GetPositionRef()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;
                neighbor_atom_set.emplace(neighbor_atom);

                const auto selected_iter{ atom_index_map.find(neighbor_atom) };
                if (selected_iter != atom_index_map.end())
                {
                    atom_context.neighbor_atom_sample_list.emplace_back(
                        NeighborAtomSample{
                            true,
                            selected_iter->second,
                            distance
                        });
                    continue;
                }

                auto unselected_atom_contributor_iter{
                    unselected_atom_index_map.find(neighbor_atom)
                };
                if (unselected_atom_contributor_iter == unselected_atom_index_map.end())
                {
                    const auto unselected_atom_contributor_index{
                        context.unselected_atom_list.size()
                    };
                    const auto group_key{
                        data_internal::GetGroupKey(neighbor_atom)
                    };
                    const auto selected_group_iter{
                        selected_group_id_by_key.find(group_key)
                    };
                    context.unselected_atom_list.emplace_back(
                        UnselectedAtomContributor{
                            neighbor_atom->GetSerialID(),
                            selected_group_iter ==
                                selected_group_id_by_key.end() ?
                                std::nullopt :
                                std::optional<std::size_t>{ selected_group_iter->second }
                        });
                    unselected_atom_contributor_iter = unselected_atom_index_map.emplace(
                        neighbor_atom,
                        unselected_atom_contributor_index).first;
                }
                atom_context.neighbor_atom_sample_list.emplace_back(
                    NeighborAtomSample{
                        false,
                        unselected_atom_contributor_iter->second,
                        distance
                    });
            }
            atom_context.neighbor_atom_sample_offset_list.emplace_back(
                atom_context.neighbor_atom_sample_list.size());
        }
        atom_context.neighbor_count_for_peeling = static_cast<int>(neighbor_atom_set.size());
    }

    return context;
}

static void StoreSecondStageNeighborCounts(
    ModelObject & model_object,
    const SecondStageContext & context)
{
    auto analysis{ model_object.EditAnalysis() };
    for (const auto & atom_context : context)
    {
        analysis.SetAtomLocalNeighborCountForPeeling(
            *atom_context.atom,
            atom_context.neighbor_count_for_peeling);
    }
}

static std::optional<GaussianModel3DWithUncertainty>
BuildValidGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list)
{
    const auto median_model{ BuildGaussianParameterMedian(model_list) };
    if (!median_model.has_value()) return std::nullopt;
    return GaussianModel3DWithUncertainty{
        *median_model,
        GaussianModel3DUncertainty{}
    };
}

static SecondStageInitialStateBuildResult BuildInitialFitState(
    SecondStageContext & context)
{
    SecondStageInitialStateBuildResult build_result;
    auto & state{ build_result.state };
    state.resize(context.size());
    std::vector<std::vector<GaussianModel3D>> models_by_group(context.selected_atom_index_list_by_group.size());
    std::vector<GaussianModel3D> global_models;
    global_models.reserve(context.size());

    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        state.at(i) = atom_context.initial_result;

        const auto & result{ state.at(i) };
        const auto direct_selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    result.posterior,
                    atom_context.group_prior,
                    std::nullopt,
                    std::nullopt
                })
        };
        if (!direct_selection.has_value()) continue;

        models_by_group.at(atom_context.group_id).emplace_back(direct_selection->model.GetModel());
        global_models.emplace_back(direct_selection->model.GetModel());
    }

    std::vector<std::optional<GaussianModel3DWithUncertainty>> median_by_group(models_by_group.size());
    for (std::size_t group_id = 0; group_id < models_by_group.size(); group_id++)
    {
        median_by_group.at(group_id) =
            BuildValidGaussianParameterMedian(models_by_group.at(group_id));
    }
    const auto global_median{ BuildValidGaussianParameterMedian(global_models) };

    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto & result{ state.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        const auto & atom_context{ context.at(i) };
        const auto & group_median{ median_by_group.at(atom_context.group_id) };
        const auto selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    result.posterior,
                    atom_context.group_prior,
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value())
        {
            build_result.failure = SecondStageInitialStateBuildResult::Failure::SelectedSeedUnavailable;
            return build_result;
        }

        result.mdpde = selection->model;
        build_result.selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                selection->source,
                original_model,
                selection->model.GetModel()
            });
    }

    for (std::size_t i = 0; i < context.unselected_atom_list.size(); i++)
    {
        auto & unselected_atom_contributor{
            context.unselected_atom_list.at(i)
        };
        std::optional<GaussianModel3DWithUncertainty> group_median;
        if (unselected_atom_contributor.selected_group_id.has_value() &&
            *unselected_atom_contributor.selected_group_id < median_by_group.size() &&
            median_by_group.at(*unselected_atom_contributor.selected_group_id).has_value())
        {
            group_median = *median_by_group.at(*unselected_atom_contributor.selected_group_id);
        }
        const auto selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    std::nullopt,
                    std::nullopt,
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value())
        {
            build_result.failure = SecondStageInitialStateBuildResult::Failure::UnselectedSeedUnavailable;
            return build_result;
        }

        unselected_atom_contributor.initial_seed = selection->model;
        build_result.unselected_selection_record_list.emplace_back(
            UnselectedSecondStageSeedSelectionRecord{
                unselected_atom_contributor.atom_serial_id,
                selection->source,
                selection->model.GetModel()
            });
    }
    return build_result;
}

static const char * GetSecondStageSeedSourceText(
    SecondStageSeedSource source)
{
    switch (source)
    {
    case SecondStageSeedSource::GroupPosterior:
        return "group-posterior";
    case SecondStageSeedSource::GroupPrior:
        return "group-prior";
    case SecondStageSeedSource::GroupMedian:
        return "group-median";
    case SecondStageSeedSource::GlobalMedian:
        return "global-median";
    }
    throw std::logic_error("Unknown second-stage seed source.");
}

static void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    constexpr std::array<SecondStageSeedSource, 4> source_list{
        SecondStageSeedSource::GroupPosterior,
        SecondStageSeedSource::GroupPrior,
        SecondStageSeedSource::GroupMedian,
        SecondStageSeedSource::GlobalMedian
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

    for (std::size_t atom_index = 0; atom_index < selection_record_list.size(); atom_index++)
    {
        const auto & record{ selection_record_list.at(atom_index) };
        std::ostringstream detail_message;
        detail_message << "Second-stage seed selection: atom index = "
            << atom_index
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

static void LogUnselectedSecondStageSeedSelections(
    const std::vector<UnselectedSecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    std::size_t group_median_count{ 0 };
    std::size_t global_median_count{ 0 };
    for (const auto & record : selection_record_list)
    {
        if (record.source == SecondStageSeedSource::GroupMedian)
        {
            group_median_count++;
        }
        else if (record.source == SecondStageSeedSource::GlobalMedian)
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
        std::ostringstream detail_message;
        detail_message
            << "Unselected second-stage neighbor seed selection: serial ID = "
            << record.atom_serial_id
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", seed A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

static void AppendQuarantineSummary(
    std::ostream & stream,
    const QuarantineSummary & summary,
    std::size_t unresolved_target_count)
{
    stream << "; quarantine entered/released/probation-failed/unresolved = "
        << summary.entered_target_count << "/"
        << summary.released_target_count << "/"
        << summary.failed_probation_count << "/"
        << unresolved_target_count;
}

} // namespace

QuarantineStateTransition UpdateQuarantineFailureState(
    const std::vector<QuarantineFailureObservation> & observation_list,
    const std::vector<QuarantineTarget> & successful_probation_target_list,
    std::size_t accepted_iteration_count,
    QuarantineFailureStateMap & state_by_target)
{
    QuarantineStateTransition transition;
    const std::set<QuarantineTarget> successful_target_set{
        successful_probation_target_list.begin(),
        successful_probation_target_list.end()
    };
    for (auto iter = state_by_target.begin(); iter != state_by_target.end();)
    {
        auto & [target, state]{ *iter };
        if (!state.probation_active)
        {
            ++iter;
            continue;
        }
        if (successful_target_set.contains(target))
        {
            transition.released_target_list.emplace_back(target);
            iter = state_by_target.erase(iter);
            continue;
        }
        state.probation_active = false;
        state.probation_count++;
        state.probation_exhausted =
            state.probation_count >= kQuarantineMaximumProbationCount;
        state.next_probation_iteration =
            accepted_iteration_count + kQuarantineProbationCooldown;
        transition.failed_probation_target_list.emplace_back(target);
        ++iter;
    }

    std::map<QuarantineTarget, QuarantineFailureReason> observation_by_target;
    for (const auto & observation : observation_list)
    {
        observation_by_target.try_emplace(observation.target, observation.reason);
    }
    for (auto iter = state_by_target.begin(); iter != state_by_target.end();)
    {
        if (!iter->second.quarantined && !observation_by_target.contains(iter->first))
        {
            iter = state_by_target.erase(iter);
            continue;
        }
        ++iter;
    }

    for (const auto & [target, reason] : observation_by_target)
    {
        auto [iter, inserted]{
            state_by_target.try_emplace(
                target,
                QuarantineFailureState{
                    reason,
                    0,
                    0,
                    0,
                    false,
                    false,
                    false
                })
        };
        auto & state{ iter->second };
        if (state.quarantined) continue;
        if (!inserted && state.reason != reason)
        {
            state.reason = reason;
            state.stable_iteration_count = 0;
        }
        state.stable_iteration_count++;
        if (state.stable_iteration_count >= kPersistentQuarantineFailureIterationLimit)
        {
            state.quarantined = true;
            state.next_probation_iteration =
                accepted_iteration_count + kQuarantineProbationCooldown;
            transition.entered_target_list.emplace_back(target);
        }
    }
    return transition;
}

namespace {

static void ApplyQuarantineFallbackTargets(
    const std::vector<QuarantineTarget> & target_list,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    FitState & assembled_state,
    PolishProvenance & assembled_polish_provenance)
{
    for (const auto & target : target_list)
    {
        for (const auto atom_index : target.atom_index_list)
        {
            if (target.kind == QuarantineTargetKind::HardFailureCluster)
            {
                assembled_state.at(atom_index) = previous_state.at(atom_index);
                assembled_polish_provenance.at(atom_index) =
                    previous_polish_provenance.at(atom_index);
                continue;
            }
            if (target.kind == QuarantineTargetKind::OffsetGroup)
            {
                const auto previous_offset{
                    previous_state.at(atom_index).mdpde.GetModel().GetOffset()
                };
                assembled_state.at(atom_index).ols = GaussianModel3DWithUncertainty{
                    assembled_state.at(atom_index).ols.GetModel().WithOffset(previous_offset),
                    assembled_state.at(atom_index).ols.GetStandardDeviationModel()
                };
                assembled_state.at(atom_index).mdpde = GaussianModel3DWithUncertainty{
                    assembled_state.at(atom_index).mdpde.GetModel().WithOffset(previous_offset),
                    assembled_state.at(atom_index).mdpde.GetStandardDeviationModel()
                };
                continue;
            }
            const auto assembled_offset{
                assembled_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            assembled_state.at(atom_index).ols = GaussianModel3DWithUncertainty{
                previous_state.at(atom_index).ols.GetModel().WithOffset(assembled_offset),
                previous_state.at(atom_index).ols.GetStandardDeviationModel()
            };
            assembled_state.at(atom_index).mdpde = GaussianModel3DWithUncertainty{
                previous_state.at(atom_index).mdpde.GetModel().WithOffset(assembled_offset),
                previous_state.at(atom_index).mdpde.GetStandardDeviationModel()
            };
            assembled_polish_provenance.at(atom_index) =
                previous_polish_provenance.at(atom_index);
        }
    }
}

SuspiciousBlockActivity QuarantineState::BeginIteration(
    std::size_t accepted_iteration_count)
{
    SuspiciousBlockActivity activity{
        SuspiciousUpdateMask(m_atom_count, 0),
        SuspiciousUpdateMask(m_atom_count, 0),
        SuspiciousUpdateMask(m_atom_count, 0)
    };
    probation_target_list.clear();
    std::vector<QuarantineTarget> due_target_list;
    for (auto & [target, state] : state_by_target)
    {
        if (!state.quarantined) continue;
        const auto probation_due{
            !state.probation_exhausted &&
            (force_probation || accepted_iteration_count >= state.next_probation_iteration)
        };
        if (probation_due)
        {
            due_target_list.emplace_back(target);
        }
    }
    std::ranges::sort(
        due_target_list,
        [](const auto & lhs, const auto & rhs)
        {
            if (lhs.kind != rhs.kind) return lhs.kind > rhs.kind;
            return lhs.atom_index_list < rhs.atom_index_list;
        });
    std::set<std::size_t> selected_probation_atom_index_set;
    for (const auto & target : due_target_list)
    {
        const auto overlaps_selected{
            std::ranges::any_of(
                target.atom_index_list,
                [&](const auto atom_index)
                {
                    return selected_probation_atom_index_set.contains(atom_index);
                })
        };
        if (overlaps_selected) continue;
        state_by_target.at(target).probation_active = true;
        probation_target_list.emplace_back(target);
        selected_probation_atom_index_set.insert(
            target.atom_index_list.begin(),
            target.atom_index_list.end());
    }
    for (auto & [target, state] : state_by_target)
    {
        if (!state.quarantined || state.probation_active) continue;
        const auto shadowed_by_broader_probation{
            std::ranges::any_of(
                probation_target_list,
                [&](const auto & probation_target)
                {
                    return probation_target.kind > target.kind &&
                        std::ranges::any_of(
                            target.atom_index_list,
                            [&](const auto atom_index)
                            {
                                return std::ranges::binary_search(
                                    probation_target.atom_index_list,
                                    atom_index);
                            });
                })
        };
        if (shadowed_by_broader_probation) continue;
        for (const auto atom_index : target.atom_index_list)
        {
            if (target.kind != QuarantineTargetKind::OffsetGroup)
            {
                activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            if (target.kind != QuarantineTargetKind::ShapeAtom)
            {
                activity.offset_fixed_atom_mask.at(atom_index) = 1;
            }
            if (target.kind == QuarantineTargetKind::HardFailureCluster)
            {
                activity.hard_failure_atom_mask.at(atom_index) = 1;
            }
        }
    }
    force_probation = false;
    return activity;
}

SuspiciousBlockActivity QuarantineState::BuildFinalActivity() const
{
    SuspiciousBlockActivity activity{
        SuspiciousUpdateMask(m_atom_count, 0),
        SuspiciousUpdateMask(m_atom_count, 0),
        SuspiciousUpdateMask(m_atom_count, 0)
    };
    for (const auto & [target, state] : state_by_target)
    {
        if (!state.quarantined) continue;
        for (const auto atom_index : target.atom_index_list)
        {
            if (target.kind != QuarantineTargetKind::OffsetGroup)
            {
                activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            if (target.kind != QuarantineTargetKind::ShapeAtom)
            {
                activity.offset_fixed_atom_mask.at(atom_index) = 1;
            }
            if (target.kind == QuarantineTargetKind::HardFailureCluster)
            {
                activity.hard_failure_atom_mask.at(atom_index) = 1;
            }
        }
    }
    return activity;
}

static bool ContainsAcceptedTargetAtoms(
    const std::vector<ClusterKey> & accepted_key_list,
    const QuarantineTarget & target)
{
    for (const auto atom_index : target.atom_index_list)
    {
        if (std::ranges::none_of(
                accepted_key_list,
                [&](const auto & key)
                {
                    return std::ranges::binary_search(key, atom_index);
                }))
        {
            return false;
        }
    }
    return true;
}

static bool DoesFailureObservationAffectTarget(
    const QuarantineFailureObservation & observation,
    const QuarantineTarget & target)
{
    const auto overlaps{
        std::ranges::any_of(
            target.atom_index_list,
            [&](const auto atom_index)
            {
                return std::ranges::binary_search(
                    observation.target.atom_index_list,
                    atom_index);
            })
    };
    if (!overlaps) return false;
    if (target.kind == QuarantineTargetKind::HardFailureCluster) return true;
    if (target.kind == QuarantineTargetKind::OffsetGroup)
    {
        return observation.target.kind != QuarantineTargetKind::ShapeAtom;
    }
    return observation.target.kind != QuarantineTargetKind::OffsetGroup;
}

static bool IsGuardSafeNonMaterialStationarity(
    const QuarantineTarget & target,
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    const ClusterHealthMap & health_by_key)
{
    for (const auto atom_index : target.atom_index_list)
    {
        const auto health_iter{
            std::ranges::find_if(
                health_by_key,
                [&](const auto & entry)
                {
                    return std::ranges::binary_search(entry.first, atom_index);
                })
        };
        if (health_iter == health_by_key.end() ||
            !health_iter->second.IsStationarityEligible() ||
            atom_index >= assessment_by_atom.size() ||
            assessment_by_atom[atom_index].IsSuspicious() ||
            assessment_by_atom[atom_index].damping_factor != 1.0)
        {
            return false;
        }
        if (target.kind != QuarantineTargetKind::OffsetGroup &&
            block_activity.shape_fixed_atom_mask.at(atom_index) == 0)
        {
            return false;
        }
        if (target.kind != QuarantineTargetKind::ShapeAtom &&
            block_activity.offset_fixed_atom_mask.at(atom_index) == 0)
        {
            return false;
        }
    }
    return true;
}

QuarantineStateTransition QuarantineState::UpdateAfterIteration(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & accepted_key_list,
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    const ClusterHealthMap & health_by_key,
    FitState & assembled_state,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    PolishProvenance & assembled_polish_provenance,
    std::size_t accepted_iteration_count)
{
    std::vector<QuarantineFailureObservation> observation_list;
    std::set<std::size_t> observed_offset_group_id_set;
    for (const auto & key : accepted_key_list)
    {
        const auto transformed_change_summary{
            SummarizeTransformedChanges(assembled_state, previous_state, key)
        };
        if (!IsTransformedPercentileConverged(transformed_change_summary)) continue;
        const auto & health{ health_by_key.at(key) };
        if (IsJointOffsetSolveHardFailure(health.joint_offset_status))
        {
            observation_list.emplace_back(QuarantineFailureObservation{
                QuarantineTarget{ QuarantineTargetKind::HardFailureCluster, key },
                health.joint_offset_status
            });
            continue;
        }
        for (const auto atom_index : key)
        {
            if (block_activity.shape_fixed_atom_mask.at(atom_index) != 0)
            {
                if (atom_index < assessment_by_atom.size() &&
                    assessment_by_atom[atom_index].reason != SuspiciousGaussianReason::None)
                {
                    observation_list.emplace_back(QuarantineFailureObservation{
                        QuarantineTarget{
                            QuarantineTargetKind::ShapeAtom,
                            { atom_index }
                        },
                        assessment_by_atom[atom_index].reason
                    });
                }
            }
            if (block_activity.offset_fixed_atom_mask.at(atom_index) == 0) continue;
            const auto group_id{ context.at(atom_index).group_id };
            if (!observed_offset_group_id_set.emplace(group_id).second) continue;
            std::vector<std::size_t> group_atom_index_list;
            for (std::size_t group_atom_index = 0;
                group_atom_index < context.size();
                group_atom_index++)
            {
                if (context.at(group_atom_index).group_id == group_id)
                {
                    group_atom_index_list.emplace_back(group_atom_index);
                }
            }
            SuspiciousGaussianReason reason{ SuspiciousGaussianReason::None };
            for (const auto group_atom_index : group_atom_index_list)
            {
                if (group_atom_index < assessment_by_atom.size() &&
                    assessment_by_atom[group_atom_index].reason != SuspiciousGaussianReason::None)
                {
                    reason = assessment_by_atom[group_atom_index].reason;
                    break;
                }
            }
            if (reason == SuspiciousGaussianReason::None) continue;
            observation_list.emplace_back(QuarantineFailureObservation{
                QuarantineTarget{
                    QuarantineTargetKind::OffsetGroup,
                    std::move(group_atom_index_list)
                },
                reason
            });
        }
    }

    std::vector<QuarantineTarget> successful_probation_target_list;
    for (const auto & target : probation_target_list)
    {
        const auto has_affecting_observation{
            std::ranges::any_of(
                observation_list,
                [&](const auto & observation)
                {
                    return DoesFailureObservationAffectTarget(observation, target);
                })
        };
        const auto accepted_material_proposal{
            ContainsAcceptedTargetAtoms(accepted_key_list, target)
        };
        if (!has_affecting_observation &&
            (accepted_material_proposal ||
                IsGuardSafeNonMaterialStationarity(
                    target,
                    block_activity,
                    assessment_by_atom,
                    health_by_key)))
        {
            successful_probation_target_list.emplace_back(target);
        }
    }
    auto transition{
        UpdateQuarantineFailureState(
            observation_list,
            successful_probation_target_list,
            accepted_iteration_count,
            state_by_target)
    };
    std::vector<QuarantineTarget> fallback_target_list{
        transition.entered_target_list
    };
    fallback_target_list.insert(
        fallback_target_list.end(),
        transition.failed_probation_target_list.begin(),
        transition.failed_probation_target_list.end());
    ApplyQuarantineFallbackTargets(
        fallback_target_list,
        previous_state,
        previous_polish_provenance,
        assembled_state,
        assembled_polish_provenance);
    summary.entered_target_count += transition.entered_target_list.size();
    summary.released_target_count += transition.released_target_list.size();
    summary.failed_probation_count += transition.failed_probation_target_list.size();
    return transition;
}

std::size_t QuarantineState::AtomCount() const
{
    std::set<std::size_t> atom_index_set;
    for (const auto & [target, state] : state_by_target)
    {
        if (!state.quarantined) continue;
        atom_index_set.insert(
            target.atom_index_list.begin(),
            target.atom_index_list.end());
    }
    return atom_index_set.size();
}

static void AppendObjectiveBreakdown(
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

static std::string_view GetPreObjectiveFailureReasonText(
    PreObjectiveFailureReason reason)
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

static std::string_view GetSuspiciousGaussianReasonText(
    SuspiciousGaussianReason reason)
{
    switch (reason)
    {
    case SuspiciousGaussianReason::None: return "none";
    case SuspiciousGaussianReason::InvalidModel: return "invalid-model";
    case SuspiciousGaussianReason::NonFiniteResponse: return "non-finite-response";
    case SuspiciousGaussianReason::OffsetMagnitude: return "offset-magnitude";
    case SuspiciousGaussianReason::CenterSignFlip: return "center-sign-flip";
    case SuspiciousGaussianReason::RadialRebound: return "radial-rebound";
    case SuspiciousGaussianReason::WidthGrowth: return "width-growth";
    case SuspiciousGaussianReason::AmplitudeOffsetCompensation:
        return "amplitude-offset-compensation";
    }
    return "none";
}

static void LogRejectedClusterDiagnostics(
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

static std::string_view GetAllRejectedResolutionText(
    AllRejectedResolution resolution)
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

static void LogAllRejectedResolution(
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

static void LogAcceptedBacktrackingDiagnostics(
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
    const auto has_boundary_reconciliation_diagnostic{
        !diagnostics.boundary_reconciliation_diagnostic_list.empty()
    };
    if (!has_local_backtracking && !has_boundary_reconciliation_diagnostic) return;
    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : diagnostics.accepted_cluster_diagnostic_list)
    {
        const auto & diagnostic{ cluster_diagnostic.attempt };
        if (cluster_diagnostic.boundary_rescued)
        {
            std::ostringstream rescue_message;
            rescue_message
                << "Accepted local fitting cluster after boundary rescue: atoms = "
                << cluster_diagnostic.key.size()
                << ", key first/last = "
                << cluster_diagnostic.key.front() << "/"
                << cluster_diagnostic.key.back() << ".";
            Logger::Log(LogLevel::Debug, rescue_message.str());
        }
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
    for (const auto & diagnostic : diagnostics.boundary_reconciliation_diagnostic_list)
    {
        const auto accepted_source_text = [&]() -> std::string_view
        {
            switch (diagnostic.accepted_source)
            {
            case BoundaryComponentAcceptedSource::None:
                return "none";
            case BoundaryComponentAcceptedSource::Endpoint:
                return "endpoint";
            case BoundaryComponentAcceptedSource::JointCorrection:
                return "joint-correction";
            case BoundaryComponentAcceptedSource::Backtracking:
                return "backtracking";
            }
            return "none";
        };
        std::ostringstream message;
        message
            << "Boundary-component reconciliation: clusters/atoms/boundary-samples = "
            << diagnostic.key_list.size() << "/"
            << diagnostic.atom_count << "/"
            << diagnostic.boundary_sample_count
            << ", accepted/rescue-candidates/rescued = "
            << diagnostic.accepted_cluster_count << "/"
            << diagnostic.rescue_candidate_cluster_count << "/"
            << diagnostic.rescued_cluster_count
            << ", locally-deteriorated/max-delta="
            << diagnostic.locally_deteriorated_member_count << "/"
            << diagnostic.maximum_local_deterioration
            << ", component/global-improvement=";
        if (diagnostic.component_improvement.has_value())
        {
            message << *diagnostic.component_improvement;
        }
        else
        {
            message << "-";
        }
        message << "/";
        if (diagnostic.global_improvement.has_value())
        {
            message << *diagnostic.global_improvement;
        }
        else
        {
            message << "-";
        }
        message
            << ", mode=" << (diagnostic.is_rescue_attempt ? "rescue" : "accepted-only")
            << ", trials/factor/accepted/exhausted = "
            << diagnostic.trial_count << "/";
        if (diagnostic.accepted_factor.has_value())
        {
            message << *diagnostic.accepted_factor;
        }
        else
        {
            message << "-";
        }
        message
            << "/" << (diagnostic.accepted ? "yes" : "no")
            << "/" << (diagnostic.exhausted ? "yes" : "no")
            << ", accepted_source=" << accepted_source_text()
            << ", objectives previous/endpoint/final=";
        const auto append_objective = [&](const std::optional<double> & objective)
        {
            if (objective.has_value())
            {
                message << *objective;
            }
            else
            {
                message << "-";
            }
        };
        append_objective(diagnostic.previous_component_objective);
        message << "/";
        append_objective(diagnostic.endpoint_component_objective);
        message << "/";
        append_objective(diagnostic.candidate_component_objective);
        message << ".";
        Logger::Log(LogLevel::Debug, message.str());
        if (!diagnostic.joint_correction_status.has_value()) continue;
        std::ostringstream correction_message;
        correction_message << std::scientific << std::setprecision(2)
            << "Boundary-interface joint correction: direct-interface/shape-active/offset-active/closure/parameters = "
            << diagnostic.interface_atom_count << "/"
            << diagnostic.shape_active_atom_count << "/"
            << diagnostic.offset_active_atom_count << "/"
            << diagnostic.offset_closure_atom_count << "/"
            << diagnostic.joint_parameter_count
            << ", status="
            << GetBoundaryJointCorrectionStatusText(*diagnostic.joint_correction_status)
            << ", damping/trust=";
        if (diagnostic.joint_damping.has_value())
        {
            correction_message << *diagnostic.joint_damping;
        }
        else
        {
            correction_message << "-";
        }
        correction_message << "/";
        if (diagnostic.maximum_normalized_trust_step.has_value())
        {
            correction_message << *diagnostic.maximum_normalized_trust_step;
        }
        else
        {
            correction_message << "-";
        }
        correction_message << ", reference/candidate=";
        if (diagnostic.joint_reference_component_objective.has_value())
        {
            correction_message << *diagnostic.joint_reference_component_objective;
        }
        else
        {
            correction_message << "-";
        }
        correction_message << "/";
        if (diagnostic.joint_candidate_component_objective.has_value())
        {
            correction_message << *diagnostic.joint_candidate_component_objective;
        }
        else
        {
            correction_message << "-";
        }
        const auto correction_accepted{
            diagnostic.accepted_source == BoundaryComponentAcceptedSource::JointCorrection
        };
        std::string_view correction_outcome{ "failed" };
        if (correction_accepted)
        {
            correction_outcome = diagnostic.endpoint_component_objective.has_value() ?
                "accepted-over-endpoint" : "accepted-over-previous";
        }
        else if (diagnostic.accepted_source == BoundaryComponentAcceptedSource::Endpoint)
        {
            correction_outcome = "fallback-endpoint";
        }
        correction_message << ", suspicious="
            << diagnostic.suspicious_candidate_atom_count
            << ", accepted=" << (correction_accepted ? "yes" : "no")
            << ", outcome=" << correction_outcome << ".";
        Logger::Log(LogLevel::Debug, correction_message.str());
    }
}

static std::string FormatProgressMaximum(double value)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

constexpr std::array<std::string_view, 6> kProgressHeaderList
{
    "Try/Acc",
    "Atom A/Q",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/R"
};

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

static ProgressColumnWidths BuildProgressColumnWidths(
    std::size_t atom_size)
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

static void LogProgressHeader(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Info, FormatProgressRow(column_widths, kProgressHeaderList));
}

static void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationProgress & progress)
{
    if (quiet_mode) return;

    const std::array<std::string, 6> cell_list{
        std::to_string(progress.attempt_number) + "/" +
            std::to_string(progress.accepted_iteration_count),
        std::to_string(progress.active_atom_count) + "/" +
            std::to_string(progress.quarantine_atom_count),
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

static void AppendAuditValues(
    std::ostringstream & message,
    const std::vector<double> & value_list)
{
    for (std::size_t i = 0; i < value_list.size(); i++)
    {
        if (i != 0) message << "/";
        message << value_list.at(i);
    }
}

static void AppendAuditPopulation(
    std::ostringstream & message,
    const std::array<std::size_t, kTransformedChangeSize> & population_size_list)
{
    for (std::size_t i = 0; i < population_size_list.size(); i++)
    {
        if (i != 0) message << "/";
        message << population_size_list.at(i);
    }
}

static void LogConvergenceSafeguardAudit(
    bool quiet_mode,
    std::size_t attempt_number,
    const IterationProgress & progress,
    const ConvergenceSafeguardPredicates & predicates,
    const ConvergenceSafeguardPredicates & shadow_predicates,
    const TransformedChangeSummary & accepted_change,
    const TransformedChangeSummary & raw_change,
    const TransformedChangeSummary & accepted_shadow_change,
    const TransformedChangeSummary & raw_shadow_change,
    const ConvergenceStationarityAudit & stationarity,
    const IterationDiagnostics & diagnostics,
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    bool accepted_equals_raw,
    bool assembled_uses_polish,
    bool has_quarantine_transition,
    bool objective_domain_changed)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;

    const auto trust_limited_count{ static_cast<std::size_t>(std::ranges::count_if(
        diagnostics.accepted_cluster_diagnostic_list,
        [](const auto & diagnostic)
        {
            return diagnostic.attempt.effective_damping < 1.0;
        })) };
    const auto cluster_backtracked_count{ static_cast<std::size_t>(std::ranges::count_if(
        diagnostics.accepted_cluster_diagnostic_list,
        [](const auto & diagnostic)
        {
            return diagnostic.attempt.accepted_backtracking_factor.has_value() &&
                *diagnostic.attempt.accepted_backtracking_factor < 1.0;
        })) };
    const auto boundary_backtracked_count{ static_cast<std::size_t>(std::ranges::count_if(
        diagnostics.boundary_reconciliation_diagnostic_list,
        [](const auto & diagnostic)
        {
            return diagnostic.accepted &&
                diagnostic.accepted_source == BoundaryComponentAcceptedSource::Backtracking;
        })) };
    const auto accepted_boundary_count{ static_cast<std::size_t>(std::ranges::count_if(
        diagnostics.boundary_reconciliation_diagnostic_list,
        &BoundaryComponentReconciliationDiagnostic::accepted)) };
    const auto rescued_boundary_count{ static_cast<std::size_t>(std::ranges::count_if(
        diagnostics.boundary_reconciliation_diagnostic_list,
        [](const auto & diagnostic)
        {
            return diagnostic.accepted && diagnostic.is_rescue_attempt;
        })) };
    const auto damped_atom_count{ static_cast<std::size_t>(std::ranges::count_if(
        assessment_by_atom,
        [](const auto & assessment)
        {
            return assessment.damping_factor != 1.0;
        })) };
    const auto shape_fixed_count{ static_cast<std::size_t>(
        std::ranges::count(block_activity.shape_fixed_atom_mask, 1)) };
    const auto offset_fixed_count{ static_cast<std::size_t>(
        std::ranges::count(block_activity.offset_fixed_atom_mask, 1)) };
    const auto hard_fixed_count{ static_cast<std::size_t>(
        std::ranges::count(block_activity.hard_failure_atom_mask, 1)) };
    std::array<std::size_t, 5> local_refit_status_count{};
    std::size_t unavailable_local_refit_status_count{ 0 };
    for (const auto status : local_refit_status_by_atom)
    {
        if (!status.has_value())
        {
            unavailable_local_refit_status_count++;
            continue;
        }
        switch (*status)
        {
        case RHBMEstimationStatus::SUCCESS:
            local_refit_status_count.at(0)++;
            break;
        case RHBMEstimationStatus::MAX_ITERATIONS_REACHED:
            local_refit_status_count.at(1)++;
            break;
        case RHBMEstimationStatus::SINGLE_MEMBER:
            local_refit_status_count.at(2)++;
            break;
        case RHBMEstimationStatus::INSUFFICIENT_DATA:
            local_refit_status_count.at(3)++;
            break;
        case RHBMEstimationStatus::NUMERICAL_FALLBACK:
            local_refit_status_count.at(4)++;
            break;
        }
    }

    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "Convergence safeguard audit: try=" << attempt_number
        << ", acc=" << progress.accepted_iteration_count
        << ", atoms=" << progress.active_atom_count + progress.quarantine_atom_count
        << ", quarantine=" << progress.quarantine_atom_count
        << ", population=";
    AppendAuditPopulation(message, accepted_change.population_size_list);
    message << ", shadow-population=";
    AppendAuditPopulation(message, accepted_shadow_change.population_size_list);
    message
        << ", predicates[s/a99/amax/r99/rmax]="
        << predicates.stationarity_eligible << "/"
        << predicates.accepted.percentile_converged << "/"
        << predicates.accepted.maximum_converged << "/"
        << predicates.raw.percentile_converged << "/"
        << predicates.raw.maximum_converged
        << ", shadow-predicates[s/a99/amax/r99/rmax]="
        << shadow_predicates.stationarity_eligible << "/"
        << shadow_predicates.accepted.percentile_converged << "/"
        << shadow_predicates.accepted.maximum_converged << "/"
        << shadow_predicates.raw.percentile_converged << "/"
        << shadow_predicates.raw.maximum_converged
        << ", accepted-p99=";
    AppendAuditValues(message, accepted_change.percentile_stats.percentile_list);
    message << ", accepted-max=";
    AppendAuditValues(message, accepted_change.maximum_list);
    message << ", raw-p99=";
    AppendAuditValues(message, raw_change.percentile_stats.percentile_list);
    message << ", raw-max=";
    AppendAuditValues(message, raw_change.maximum_list);
    message << ", shadow-accepted-p99=";
    AppendAuditValues(message, accepted_shadow_change.percentile_stats.percentile_list);
    message << ", shadow-accepted-max=";
    AppendAuditValues(message, accepted_shadow_change.maximum_list);
    message << ", shadow-raw-p99=";
    AppendAuditValues(message, raw_shadow_change.percentile_stats.percentile_list);
    message << ", shadow-raw-max=";
    AppendAuditValues(message, raw_shadow_change.maximum_list);
    message
        << ", accepted-equals-raw=" << accepted_equals_raw
        << ", path[trust/backtrack/polish/boundary/rescue]="
        << trust_limited_count << "/"
        << cluster_backtracked_count + boundary_backtracked_count << "/"
        << assembled_uses_polish << "/"
        << accepted_boundary_count << "/"
        << rescued_boundary_count
        << ", stationarity[current/full/active-ineligible/refit-ineligible/soft-joint/hard-joint]="
        << stationarity.active_block_eligible << "/"
        << stationarity.full_cluster_eligible << "/"
        << stationarity.active_block_ineligible_cluster_count << "/"
        << stationarity.refit_ineligible_cluster_count << "/"
        << stationarity.soft_joint_nonconverged_cluster_count << "/"
        << stationarity.hard_joint_failure_cluster_count
        << ", local-status[success/max-iter/single/insufficient/numerical/unavailable]="
        << local_refit_status_count.at(0) << "/"
        << local_refit_status_count.at(1) << "/"
        << local_refit_status_count.at(2) << "/"
        << local_refit_status_count.at(3) << "/"
        << local_refit_status_count.at(4) << "/"
        << unavailable_local_refit_status_count
        << ", fixed[shape/offset/hard]="
        << shape_fixed_count << "/" << offset_fixed_count << "/" << hard_fixed_count
        << ", damped-atoms=" << damped_atom_count
        << ", blockers[suspicious/rejected/quarantine-transition/domain-change]="
        << progress.suspicious_atom_count << "/"
        << progress.rejected_cluster_count << "/"
        << has_quarantine_transition << "/"
        << objective_domain_changed << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Debug, message.str());
}

constexpr std::size_t kGuardAwareDampingMaximumHalvings{ 8 };

struct GuardAwareDampingResult
{
    std::optional<GaussianModel3D> model{};
    SuspiciousGaussianAssessment assessment{};
};

static std::optional<GaussianModel3D> InterpolateGaussianShape(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & endpoint_model,
    double fixed_offset,
    double factor)
{
    const auto previous_coordinates{ EncodeTransformedCoordinates(previous_model) };
    const auto endpoint_coordinates{ EncodeTransformedCoordinates(endpoint_model) };
    if (!previous_coordinates.has_value() || !endpoint_coordinates.has_value() ||
        !std::isfinite(factor) || factor < 0.0 || factor > 1.0)
    {
        return std::nullopt;
    }
    Eigen::Vector3d shape_coordinates{
        (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) +
            factor * (
                (*endpoint_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) -
                (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex))),
        (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) +
            factor * (
                (*endpoint_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) -
                (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex))),
        0.0
    };
    const auto shape_model{ DecodeTransformedCoordinates(shape_coordinates) };
    if (!shape_model.has_value()) return std::nullopt;
    const auto candidate_model{ shape_model->WithOffset(fixed_offset) };
    return IsValidSecondStageGaussianModel(candidate_model) ?
        std::optional<GaussianModel3D>{ candidate_model } : std::nullopt;
}

static bool HasMaterialShapeChange(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model)
{
    const auto previous_coordinates{ EncodeTransformedCoordinates(previous_model) };
    const auto candidate_coordinates{ EncodeTransformedCoordinates(candidate_model) };
    if (!previous_coordinates.has_value() || !candidate_coordinates.has_value()) return false;
    return std::abs(
               (*candidate_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) -
               (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex))) >=
            kTransformedChangeTolerance ||
        std::abs(
               (*candidate_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) -
               (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex))) >=
            kTransformedChangeTolerance;
}

static GuardAwareDampingResult FindGuardSafeShapeCandidate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & endpoint_model,
    const FitOptions & options,
    const SuspiciousUpdateBaseline & previous_baseline)
{
    GuardAwareDampingResult result;
    double factor{ 1.0 };
    for (std::size_t halving_count = 0;
        halving_count <= kGuardAwareDampingMaximumHalvings;
        halving_count++, factor *= 0.5)
    {
        const auto candidate_model{
            InterpolateGaussianShape(
                previous_model,
                endpoint_model,
                endpoint_model.GetOffset(),
                factor)
        };
        if (!candidate_model.has_value()) continue;
        auto assessment{
            AssessSuspiciousGaussianUpdate(
                sample_entries,
                *candidate_model,
                options,
                previous_baseline,
                SuspiciousUpdateMode::PostRefit)
        };
        assessment.damping_trial_count = halving_count + 1;
        assessment.damping_factor = factor;
        result.assessment = assessment;
        if (assessment.IsSuspicious()) continue;
        if (!HasMaterialShapeChange(previous_model, *candidate_model)) return result;
        result.model = *candidate_model;
        return result;
    }
    return result;
}

static std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
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
    std::optional<SuspiciousGaussianAssessment> failed_shape_assessment;
    std::optional<RHBMEstimationStatus> attempted_refit_status;
    bool has_guard_safe_nonmaterial_stationarity{ false };
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
        if (candidate_result.fit_result.has_value())
        {
            attempted_refit_status = candidate_result.fit_result->status;
        }
        const auto damped_candidate{
            FindGuardSafeShapeCandidate(
                adjusted_sampling_entries,
                previous_model,
                candidate_result.mdpde.GetModel(),
                options,
                previous_baseline)
        };
        if (damped_candidate.model.has_value())
        {
            const auto damping_factor{ damped_candidate.assessment.damping_factor };
            candidate_result.mdpde = GaussianModel3DWithUncertainty{
                *damped_candidate.model,
                candidate_result.mdpde.GetStandardDeviationModel()
            };
            if (const auto damped_ols{
                    InterpolateGaussianShape(
                        previous_model,
                        candidate_result.ols.GetModel(),
                        offset_model.GetOffset(),
                        damping_factor) })
            {
                candidate_result.ols = GaussianModel3DWithUncertainty{
                    *damped_ols,
                    candidate_result.ols.GetStandardDeviationModel()
                };
            }
            const auto is_stationarity_eligible{
                damping_factor == 1.0 &&
                candidate_result.fit_result.has_value() &&
                IsLocalRefitStatusStationarityEligible(candidate_result.fit_result->status)
            };
            return LocalAtomRefitResult{
                std::move(candidate_result),
                is_stationarity_eligible,
                true,
                false,
                damped_candidate.assessment,
                attempted_refit_status
            };
        }
        has_guard_safe_nonmaterial_stationarity =
            !damped_candidate.assessment.IsSuspicious() &&
            damped_candidate.assessment.damping_factor == 1.0 &&
            candidate_result.fit_result.has_value() &&
            IsLocalRefitStatusStationarityEligible(candidate_result.fit_result->status);
        failed_shape_assessment = damped_candidate.assessment;
    }
    catch (const std::exception &)
    {
        failed_shape_assessment = SuspiciousGaussianAssessment{
            SuspiciousGaussianReason::InvalidModel,
            SuspiciousUpdateMode::PostRefit,
            std::numeric_limits<double>::infinity(),
            1,
            1.0
        };
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
    auto fallback_assessment{
        AssessSuspiciousGaussianUpdate(
            adjusted_sampling_entries,
            result.mdpde.GetModel(),
            options,
            previous_baseline,
            SuspiciousUpdateMode::OffsetOnly)
    };
    if (fallback_assessment.IsSuspicious())
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{
        std::move(result),
        has_guard_safe_nonmaterial_stationarity,
        true,
        true,
        failed_shape_assessment.value_or(fallback_assessment),
        attempted_refit_status
    };
}

static RawIterationResult RunRawIteration(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list,
    const SuspiciousBlockActivity & quarantine_activity,
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

    auto iteration_state{ previous_state };
    SuspiciousBlockActivity block_activity{
        SuspiciousUpdateMask(context.size(), 0),
        SuspiciousUpdateMask(context.size(), 0),
        SuspiciousUpdateMask(context.size(), 0)
    };
    std::vector<SuspiciousGaussianAssessment> assessment_by_atom(context.size());
    std::vector<std::optional<RHBMEstimationStatus>>
        local_refit_status_by_atom(context.size());
    std::vector<std::size_t> group_id_by_atom_index;
    group_id_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_id_by_atom_index.emplace_back(atom_context.group_id);
    }
    ClusterHealthMap health_by_key;
    for (std::size_t cluster_position = 0;
        cluster_position < cluster_key_list.size();
        cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        const auto & offset_result{ joint_offset_result_list.at(cluster_position) };
        auto [health_iter, inserted]{
            health_by_key.emplace(key, ClusterHealth{ offset_result.status })
        };
        static_cast<void>(inserted);
        auto & health{ health_iter->second };
        std::map<std::size_t, std::vector<std::size_t>> position_list_by_group;
        for (std::size_t position = 0; position < key.size(); position++)
        {
            position_list_by_group[group_id_by_atom_index.at(key.at(position))].emplace_back(position);
        }
        if (IsJointOffsetSolveHardFailure(offset_result.status))
        {
            health.is_refit_stationarity_eligible = false;
            health.is_active_block_stationarity_eligible =
                std::ranges::none_of(
                    key,
                    [&](const auto atom_index)
                    {
                        return quarantine_activity.HasActiveShape(atom_index) ||
                            quarantine_activity.HasActiveOffset(atom_index);
                    });
            for (const auto atom_index : key)
            {
                block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                block_activity.hard_failure_atom_mask.at(atom_index) = 1;
                current_model_snapshot.selected.at(atom_index) =
                    previous_state.at(atom_index).mdpde.GetModel();
            }
            continue;
        }

        for (const auto & [group_id, position_list] : position_list_by_group)
        {
            static_cast<void>(group_id);
            const auto is_quarantined_group{
                std::ranges::any_of(
                    position_list,
                    [&](const auto position)
                    {
                        const auto atom_index{ key.at(position) };
                        return !quarantine_activity.HasActiveOffset(atom_index);
                    })
            };
            if (is_quarantined_group)
            {
                health.is_refit_stationarity_eligible = false;
                for (const auto position : position_list)
                {
                    const auto atom_index{ key.at(position) };
                    block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                    current_model_snapshot.selected.at(atom_index) =
                        previous_state.at(atom_index).mdpde.GetModel();
                }
                continue;
            }
            bool accepted{ false };
            std::vector<GaussianModel3D> accepted_model_list;
            std::vector<SuspiciousGaussianAssessment> accepted_assessment_list;
            std::vector<SuspiciousGaussianAssessment> last_assessment_list;
            bool has_guard_safe_nonmaterial_stationarity{ false };
            double factor{ 1.0 };
            for (std::size_t halving_count = 0;
                halving_count <= kGuardAwareDampingMaximumHalvings;
                halving_count++, factor *= 0.5)
            {
                bool all_safe{ true };
                std::vector<GaussianModel3D> candidate_model_list;
                std::vector<SuspiciousGaussianAssessment> candidate_assessment_list;
                candidate_model_list.reserve(position_list.size());
                candidate_assessment_list.reserve(position_list.size());
                for (const auto position : position_list)
                {
                    const auto atom_index{ key.at(position) };
                    const auto & previous_model{
                        previous_state.at(atom_index).mdpde.GetModel()
                    };
                    const auto proposed_offset{
                        offset_result.offset(static_cast<Eigen::Index>(position))
                    };
                    const auto candidate_model{
                        previous_model.WithOffset(
                            previous_model.GetOffset() +
                            factor * (proposed_offset - previous_model.GetOffset()))
                    };
                    auto assessment{
                        AssessSuspiciousGaussianUpdate(
                            context.at(atom_index).raw_sampling_entries,
                            candidate_model,
                            options,
                            BuildPreviousSuspiciousProfileBaseline(
                                context.at(atom_index).raw_sampling_entries,
                                previous_model,
                                options),
                            SuspiciousUpdateMode::OffsetOnly)
                    };
                    assessment.damping_trial_count = halving_count + 1;
                    assessment.damping_factor = factor;
                    all_safe = all_safe && !assessment.IsSuspicious();
                    candidate_model_list.emplace_back(candidate_model);
                    candidate_assessment_list.emplace_back(assessment);
                }
                last_assessment_list = candidate_assessment_list;
                if (!all_safe) continue;
                bool has_material_group_change{ false };
                for (std::size_t member_position = 0;
                    member_position < candidate_model_list.size();
                    member_position++)
                {
                    const auto atom_index{
                        key.at(position_list.at(member_position))
                    };
                    has_material_group_change = has_material_group_change ||
                        IsTransformedChangeMaterial(
                            CalculateTransformedChange(
                                candidate_model_list.at(member_position),
                                previous_state.at(atom_index).mdpde.GetModel()),
                            kTransformedChangeTolerance);
                }
                if (!has_material_group_change)
                {
                    has_guard_safe_nonmaterial_stationarity = halving_count == 0;
                    break;
                }
                accepted = true;
                accepted_model_list = std::move(candidate_model_list);
                accepted_assessment_list = std::move(candidate_assessment_list);
                if (halving_count != 0)
                {
                    health.is_refit_stationarity_eligible = false;
                    health.is_active_block_stationarity_eligible = false;
                }
                break;
            }
            for (std::size_t member_position = 0;
                member_position < position_list.size();
                member_position++)
            {
                const auto atom_index{ key.at(position_list.at(member_position)) };
                if (accepted)
                {
                    current_model_snapshot.selected.at(atom_index) =
                        accepted_model_list.at(member_position);
                    assessment_by_atom.at(atom_index) =
                        accepted_assessment_list.at(member_position);
                }
                else
                {
                    block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                    current_model_snapshot.selected.at(atom_index) =
                        previous_state.at(atom_index).mdpde.GetModel();
                    if (member_position < last_assessment_list.size())
                    {
                        assessment_by_atom.at(atom_index) =
                            last_assessment_list.at(member_position);
                    }
                }
            }
            if (!accepted && !has_guard_safe_nonmaterial_stationarity)
            {
                health.is_refit_stationarity_eligible = false;
                health.is_active_block_stationarity_eligible = false;
            }
        }
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
            refit_atom_index_list.emplace_back(atom_index);
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
        for (const auto atom_index : key)
        {
            auto refit_result{ std::move(refit_result_list.at(refit_position++)) };
            if (!refit_result.has_value())
            {
                health.is_refit_stationarity_eligible = false;
                if (quarantine_activity.HasActiveShape(atom_index) ||
                    quarantine_activity.HasActiveOffset(atom_index))
                {
                    health.is_active_block_stationarity_eligible = false;
                }
                block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
                const auto failed_group_id{ group_id_by_atom_index.at(atom_index) };
                for (std::size_t group_atom_index = 0;
                    group_atom_index < group_id_by_atom_index.size();
                    group_atom_index++)
                {
                    if (group_id_by_atom_index.at(group_atom_index) == failed_group_id)
                    {
                        block_activity.offset_fixed_atom_mask.at(group_atom_index) = 1;
                    }
                }
                assessment_by_atom.at(atom_index) = SuspiciousGaussianAssessment{
                    SuspiciousGaussianReason::InvalidModel,
                    SuspiciousUpdateMode::PostRefit,
                    std::numeric_limits<double>::infinity(),
                    1,
                    1.0
                };
                continue;
            }
            local_refit_status_by_atom.at(atom_index) =
                refit_result->attempted_refit_status;
            if (!refit_result->is_stationarity_eligible)
            {
                health.is_refit_stationarity_eligible = false;
                if (quarantine_activity.HasActiveShape(atom_index))
                {
                    health.is_active_block_stationarity_eligible = false;
                }
            }
            if (!refit_result->is_boundary_correction_eligible)
            {
                health.is_boundary_correction_eligible = false;
            }
            if (refit_result->shape_fixed)
            {
                block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            assessment_by_atom.at(atom_index) = refit_result->assessment;
            iteration_state.at(atom_index) = std::move(refit_result->result);
        }
    }
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        if (block_activity.offset_fixed_atom_mask.at(atom_index) == 0 &&
            block_activity.hard_failure_atom_mask.at(atom_index) == 0)
        {
            continue;
        }
        const auto previous_offset{
            previous_state.at(atom_index).mdpde.GetModel().GetOffset()
        };
        iteration_state.at(atom_index).ols = GaussianModel3DWithUncertainty{
            iteration_state.at(atom_index).ols.GetModel().WithOffset(previous_offset),
            iteration_state.at(atom_index).ols.GetStandardDeviationModel()
        };
        iteration_state.at(atom_index).mdpde = GaussianModel3DWithUncertainty{
            iteration_state.at(atom_index).mdpde.GetModel().WithOffset(previous_offset),
            iteration_state.at(atom_index).mdpde.GetStandardDeviationModel()
        };
    }

    const auto failure_block_activity{ block_activity };
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        if (!quarantine_activity.HasActiveShape(atom_index))
        {
            const auto accepted_offset{
                iteration_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            iteration_state.at(atom_index).ols = GaussianModel3DWithUncertainty{
                previous_state.at(atom_index).ols.GetModel().WithOffset(accepted_offset),
                previous_state.at(atom_index).ols.GetStandardDeviationModel()
            };
            iteration_state.at(atom_index).mdpde = GaussianModel3DWithUncertainty{
                previous_state.at(atom_index).mdpde.GetModel().WithOffset(accepted_offset),
                previous_state.at(atom_index).mdpde.GetStandardDeviationModel()
            };
            block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
        }
        if (!quarantine_activity.HasActiveOffset(atom_index))
        {
            const auto previous_offset{
                previous_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            iteration_state.at(atom_index).ols = GaussianModel3DWithUncertainty{
                iteration_state.at(atom_index).ols.GetModel().WithOffset(previous_offset),
                iteration_state.at(atom_index).ols.GetStandardDeviationModel()
            };
            iteration_state.at(atom_index).mdpde = GaussianModel3DWithUncertainty{
                iteration_state.at(atom_index).mdpde.GetModel().WithOffset(previous_offset),
                iteration_state.at(atom_index).mdpde.GetStandardDeviationModel()
            };
            block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
        }
        if (quarantine_activity.hard_failure_atom_mask.at(atom_index) != 0)
        {
            iteration_state.at(atom_index) = previous_state.at(atom_index);
            block_activity.hard_failure_atom_mask.at(atom_index) = 1;
        }
    }
    auto rollback_atom_mask{ block_activity.BuildCombinedFixedAtomMask() };

    return RawIterationResult{
        std::move(iteration_state),
        std::move(rollback_atom_mask),
        std::move(block_activity),
        failure_block_activity,
        std::move(assessment_by_atom),
        std::move(local_refit_status_by_atom),
        std::move(health_by_key)
    };
}

using GraphEdgeSet = std::set<std::pair<std::size_t, std::size_t>>;

static GraphEdgeSet BuildGraphEdgeSet(const GraphTopology & topology)
{
    GraphEdgeSet edge_set;
    for (std::size_t atom_index = 0; atom_index < topology.adjacency_list.size(); atom_index++)
    {
        for (const auto neighbor_index : topology.adjacency_list.at(atom_index))
        {
            if (atom_index < neighbor_index)
            {
                edge_set.emplace(atom_index, neighbor_index);
            }
        }
    }
    return edge_set;
}

static std::size_t CountGraphEdgeDifference(const GraphEdgeSet & source, const GraphEdgeSet & destination)
{
    return static_cast<std::size_t>(
        std::ranges::count_if(
            source,
            [&](const auto & edge)
            {
                return !destination.contains(edge);
            }));
}

static bool AreGraphPartitionsEqual(const CouplingGraphPartition & lhs, const CouplingGraphPartition & rhs)
{
    return lhs.boundary_sample_count == rhs.boundary_sample_count &&
        lhs.boundary_sample_dependency_list == rhs.boundary_sample_dependency_list &&
        lhs.sample_id_list_by_key == rhs.sample_id_list_by_key;
}

static std::string_view GetAdaptiveTopologyTriggerText(
    AdaptiveTopologyRebuildTrigger trigger)
{
    switch (trigger)
    {
        case AdaptiveTopologyRebuildTrigger::Drift: return "drift";
        case AdaptiveTopologyRebuildTrigger::Interval: return "interval";
        case AdaptiveTopologyRebuildTrigger::None: return "none";
    }
    return "unknown";
}

static void ResetIterationStateForPartition(
    const SecondStageContext & context,
    const FitOptions & options,
    const FitState & accepted_state,
    bool accepted_uses_polish,
    std::size_t source_iteration,
    CouplingGraphPartition partition,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    const auto cluster_key_list{ BuildGraphClusterKeyList(partition) };
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, accepted_state)
    };
    iteration_state.objective_domain = BuildObjectiveDomain(
        context,
        model_snapshot,
        cluster_key_list,
        options.distance_min,
        options.distance_max);
    iteration_state.cluster_objective_state.clear();
    const auto objective_by_key{
        BuildObjectiveByKey(
            partition,
            iteration_state.objective_domain,
            SnapshotResidualEvaluator{ context, model_snapshot })
    };
    ReconcileClusterObjectiveState(
        objective_by_key,
        iteration_state.cluster_objective_state);
    const auto audit_objective{
        EvaluateAuditObjective(
            iteration_state.objective_domain,
            SnapshotResidualEvaluator{ context, model_snapshot })
    };
    iteration_state.best_audit_state.reset();
    if (audit_objective.has_value())
    {
        TryUpdateBestAuditState(
            accepted_state,
            accepted_uses_polish,
            source_iteration,
            *audit_objective,
            iteration_state.best_audit_state);
    }
    iteration_state.trust_region_state.Reconcile(cluster_key_list);
    performance_counters.RecordSolverWorkspaceReset();
    ResetClusterSolverWorkspace(
        cluster_key_list,
        iteration_state.solver_workspace_by_key);
    iteration_state.boundary_joint_correction_workspace_by_key.clear();
    iteration_state.graph_partition = std::move(partition);
    iteration_state.unchanged_state_exhausted_key_list.clear();
    iteration_state.audit_patience_count = 0;
}

static bool TryRebuildAdaptiveTopology(
    const SecondStageContext & context,
    const FitOptions & options,
    const FitState & accepted_state,
    bool accepted_uses_polish,
    GraphTopology & graph_topology,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    const auto decision{
        EvaluateAdaptiveTopologyRebuildTrigger(
            accepted_state,
            iteration_state.topology_reference_state,
            iteration_state.active_index_list,
            iteration_state.accepted_iterations_since_topology_rebuild)
    };
    if (decision.trigger == AdaptiveTopologyRebuildTrigger::None) return false;

    if (!options.quiet_mode) Logger::FinishProgressLine();
    const auto rebuild_start{ std::chrono::steady_clock::now() };
    auto rebuilt_topology{
        BuildAdaptiveSecondStageGraphTopology(
            context,
            accepted_state,
            graph_topology,
            options.quiet_mode)
    };
    auto rebuilt_partition{
        BuildGraphPartition(
            rebuilt_topology,
            iteration_state.active_index_list)
    };
    const auto partition_changed{
        !AreGraphPartitionsEqual(
            iteration_state.graph_partition,
            rebuilt_partition)
    };
    const auto previous_edge_set{ BuildGraphEdgeSet(graph_topology) };
    const auto rebuilt_edge_set{ BuildGraphEdgeSet(rebuilt_topology) };
    const auto removed_edge_count{
        CountGraphEdgeDifference(previous_edge_set, rebuilt_edge_set)
    };
    const auto added_edge_count{
        CountGraphEdgeDifference(rebuilt_edge_set, previous_edge_set)
    };
    const auto elapsed_milliseconds{
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - rebuild_start).count()
    };
    performance_counters.RecordTopologyRebuild(elapsed_milliseconds, partition_changed);

    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        std::ostringstream message;
        message
            << "Adaptive local-fitting topology rebuild: accepted_iteration="
            << iteration_state.accepted_iteration_count
            << ", trigger=" << GetAdaptiveTopologyTriggerText(decision.trigger)
            << std::scientific << std::setprecision(2)
            << ", drift=" << decision.maximum_transformed_drift
            << ", clusters="
            << iteration_state.graph_partition.sample_id_list_by_key.size()
            << "/" << rebuilt_partition.sample_id_list_by_key.size()
            << ", boundary_samples="
            << iteration_state.graph_partition.boundary_sample_count
            << "/" << rebuilt_partition.boundary_sample_count
            << ", edges_added/removed="
            << added_edge_count << "/" << removed_edge_count
            << ", partition_changed="
            << (partition_changed ? "yes" : "no")
            << ", objective_domain_reset="
            << (partition_changed ? "yes" : "no") << ".";
        Logger::Log(LogLevel::Info, message.str());
    }

    graph_topology = std::move(rebuilt_topology);
    iteration_state.topology_reference_state = accepted_state;
    iteration_state.accepted_iterations_since_topology_rebuild = 0;
    performance_counters.RecordFullStateMaterialization();
    if (partition_changed)
    {
        ResetIterationStateForPartition(
            context,
            options,
            accepted_state,
            accepted_uses_polish,
            iteration_state.accepted_iteration_count,
            std::move(rebuilt_partition),
            iteration_state,
            performance_counters);
    }
    return partition_changed;
}

static IterationState BuildIterationState(
    const SecondStageContext & context,
    const GraphTopology & graph_topology,
    FitState initial_state,
    const FitOptions & options)
{
    IterationState iteration_state;
    iteration_state.previous_state = std::move(initial_state);
    iteration_state.topology_reference_state = iteration_state.previous_state;
    iteration_state.previous_polish_provenance.assign(context.size(), 0);
    iteration_state.rollback_atom_mask.assign(context.size(), 0);
    iteration_state.quarantine_state = QuarantineState(context.size());
    iteration_state.active_index_list.resize(context.size());
    std::iota(
        iteration_state.active_index_list.begin(),
        iteration_state.active_index_list.end(),
        0);
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

static IterationResult RunIteration(
    const SecondStageContext & context,
    GraphTopology & graph_topology,
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

    const auto quarantine_activity{
        iteration_state.quarantine_state.BeginIteration(
            iteration_state.accepted_iteration_count)
    };
    std::set<std::size_t> probation_atom_index_set;
    for (const auto & target : iteration_state.quarantine_state.probation_target_list)
    {
        probation_atom_index_set.insert(
            target.atom_index_list.begin(),
            target.atom_index_list.end());
    }
    std::vector<ClusterKey> probation_key_list;
    for (const auto & key : cluster_key_list)
    {
        if (std::ranges::any_of(
                key,
                [&](const auto atom_index)
                {
                    return probation_atom_index_set.contains(atom_index);
                }))
        {
            probation_key_list.emplace_back(key);
        }
    }
    iteration_state.trust_region_state.ResetToMinimum(probation_key_list);
    auto ridge_atom_mask{ iteration_state.rollback_atom_mask };
    const auto quarantine_atom_mask{
        quarantine_activity.BuildCombinedFixedAtomMask()
    };
    for (std::size_t atom_index = 0; atom_index < ridge_atom_mask.size(); atom_index++)
    {
        if (quarantine_atom_mask.at(atom_index) != 0)
        {
            ridge_atom_mask.at(atom_index) = 1;
        }
        if (probation_atom_index_set.contains(atom_index))
        {
            ridge_atom_mask.at(atom_index) = 1;
        }
    }
    const auto joint_offset_ridge_multiplier_list{
        BuildSuspiciousJointOffsetRidgeMultiplierList(ridge_atom_mask)
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
            quarantine_activity,
            iteration_state.solver_workspace_by_key)
    };
    performance_counters.FinishIterationPhase(iteration_phase_start);
    performance_counters.RecordGaussianCacheHits();

    if (!options.quiet_mode && Logger::GetLogLevel() >= LogLevel::Debug)
    {
        for (std::size_t atom_index = 0;
            atom_index < raw_iteration_result.assessment_by_atom.size();
            atom_index++)
        {
            const auto & assessment{
                raw_iteration_result.assessment_by_atom.at(atom_index)
            };
            const auto has_fixed_block{
                raw_iteration_result.block_activity.shape_fixed_atom_mask.at(atom_index) != 0 ||
                raw_iteration_result.block_activity.offset_fixed_atom_mask.at(atom_index) != 0 ||
                raw_iteration_result.block_activity.hard_failure_atom_mask.at(atom_index) != 0
            };
            if (!assessment.IsSuspicious() &&
                assessment.damping_factor == 1.0 &&
                !has_fixed_block)
            {
                continue;
            }
            std::ostringstream message;
            message << std::scientific << std::setprecision(2)
                << "Guard-aware local fitting: atom=" << atom_index
                << ", reason=" << GetSuspiciousGaussianReasonText(assessment.reason)
                << ", margin=" << assessment.normalized_margin
                << ", trials/factor=" << assessment.damping_trial_count
                << "/" << assessment.damping_factor
                << ", shape-fixed/offset-fixed/hard-fixed="
                << (raw_iteration_result.block_activity.shape_fixed_atom_mask.at(atom_index) != 0 ? "yes" : "no")
                << "/"
                << (raw_iteration_result.block_activity.offset_fixed_atom_mask.at(atom_index) != 0 ? "yes" : "no")
                << "/"
                << (raw_iteration_result.block_activity.hard_failure_atom_mask.at(atom_index) != 0 ? "yes" : "no")
                << ".";
            Logger::Log(LogLevel::Debug, message.str());
        }
    }

    const auto iteration_failure_atom_mask{
        raw_iteration_result.failure_block_activity.BuildCombinedFixedAtomMask()
    };
    const auto iteration_suspicious_atom_count{
        CountSuspiciousAtoms(iteration_failure_atom_mask)
    };
    const auto has_suspicious_offset_fallback{ iteration_suspicious_atom_count > 0 };
    const auto is_stationarity_eligible{
        std::ranges::all_of(
            raw_iteration_result.health_by_key | std::views::values,
            &ClusterHealth::is_active_block_stationarity_eligible)
    };
    const auto & raw_state{ raw_iteration_result.state };
    const auto raw_fixed_point_change_summary{
        SummarizeTransformedChanges(raw_state, previous_state, active_index_list)
    };
    const CandidateSelectionInputs candidate_inputs{
        .context = context,
        .options = options,
        .residual_baseline = residual_baseline,
        .partition = graph_partition,
        .health_by_key = raw_iteration_result.health_by_key,
        .previous_state = previous_state,
        .previous_polish_provenance = iteration_state.previous_polish_provenance,
        .raw_state = raw_state,
        .rollback_atom_mask = raw_iteration_result.rollback_atom_mask,
        .block_activity = raw_iteration_result.block_activity,
        .assessment_by_atom = raw_iteration_result.assessment_by_atom,
        .ridge_multiplier_list = joint_offset_ridge_multiplier_list,
        .unchanged_state_exhausted_key_list = std::span<const ClusterKey>{ iteration_state.unchanged_state_exhausted_key_list },
        .objective_domain = objective_domain,
        .previous_objective_by_key = previous_objective_by_key,
        .cluster_objective_state = iteration_state.cluster_objective_state,
        .best_audit_state = iteration_state.best_audit_state,
        .trust_region_state = iteration_state.trust_region_state,
        .solver_workspace_by_key = iteration_state.solver_workspace_by_key,
        .boundary_joint_correction_workspace_by_key = iteration_state.boundary_joint_correction_workspace_by_key,
        .thread_size = options.thread_size,
        .performance_counters = performance_counters
    };
    auto selection{ SelectClusterCandidates(candidate_inputs) };

    auto assembled_state{ std::move(selection.assembled_state) };
    auto assembled_polish_provenance{ std::move(selection.assembled_polish_provenance) };
    auto trust_region_iteration_update{
        iteration_state.trust_region_state.UpdateAfterIteration(
            selection.grow_trust_region_key_list,
            selection.rejected_key_list,
            selection.backtracking_exhausted_key_list)
    };
    const auto quarantine_transition{
        iteration_state.quarantine_state.UpdateAfterIteration(
            context,
            selection.accepted_key_list,
            raw_iteration_result.failure_block_activity,
            raw_iteration_result.assessment_by_atom,
            raw_iteration_result.health_by_key,
            assembled_state,
            previous_state,
            iteration_state.previous_polish_provenance,
            assembled_polish_provenance,
            iteration_state.accepted_iteration_count + 1)
    };
    const auto has_quarantine_transition{
        !quarantine_transition.entered_target_list.empty() ||
        !quarantine_transition.released_target_list.empty() ||
        !quarantine_transition.failed_probation_target_list.empty()
    };
    if (has_quarantine_transition) selection.final_audit_objective.reset();
    const auto assembled_uses_polish{
        UsesPolish(assembled_polish_provenance)
    };
    bool objective_domain_changed{ false };

    IterationResult result;
    result.objective_domain_changed = objective_domain_changed;
    if (!selection.accepted_key_list.empty())
    {
        result.diagnostics.accepted_cluster_diagnostic_list = std::move(selection.accepted_cluster_diagnostic_list);
    }
    result.diagnostics.rejected_cluster_diagnostic_list = std::move(selection.rejected_cluster_diagnostic_list);
    result.diagnostics.boundary_reconciliation_diagnostic_list = selection.boundary_reconciliation_diagnostic_list;
    result.diagnostics.trust_region_update = std::move(trust_region_iteration_update);
    iteration_state.rollback_atom_mask = std::move(raw_iteration_result.rollback_atom_mask);
    result.progress = IterationProgress{
        attempt_number,
        iteration_state.accepted_iteration_count,
        context.size() - iteration_state.quarantine_state.AtomCount(),
        iteration_state.quarantine_state.AtomCount(),
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
    const auto safeguard_predicates{
        EvaluateConvergenceSafeguardPredicates(
            is_stationarity_eligible,
            transformed_change_summary,
            raw_fixed_point_change_summary)
    };
    iteration_state.accepted_iteration_count++;
    iteration_state.accepted_iterations_since_topology_rebuild++;
    if (TryRebuildAdaptiveTopology(
            context,
            options,
            assembled_state,
            assembled_uses_polish,
            graph_topology,
            iteration_state,
            performance_counters))
    {
        objective_domain_changed = true;
        iteration_state.quarantine_state.NotifyTopologyChanged();
    }
    result.objective_domain_changed = objective_domain_changed;
    bool improved_best_audit{ false };
    if (!objective_domain_changed)
    {
        auto candidate_audit_objective{ selection.final_audit_objective };
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
        !objective_domain_changed &&
        !has_quarantine_transition &&
        !has_suspicious_offset_fallback &&
        selection.rejected_key_list.empty() &&
        safeguard_predicates.Converged();

    if (!options.quiet_mode && Logger::GetLogLevel() >= LogLevel::Debug)
    {
        const auto active_block_change_index_list{
            BuildActiveBlockChangeIndexLists(
                iteration_state.active_index_list,
                raw_iteration_result.block_activity)
        };
        const auto accepted_shadow_change_summary{
            SummarizeTransformedChangesByParameter(
                assembled_state,
                previous_state,
                active_block_change_index_list)
        };
        const auto raw_shadow_change_summary{
            SummarizeTransformedChangesByParameter(
                raw_state,
                previous_state,
                active_block_change_index_list)
        };
        const auto shadow_safeguard_predicates{
            EvaluateConvergenceSafeguardPredicates(
                is_stationarity_eligible,
                accepted_shadow_change_summary,
                raw_shadow_change_summary)
        };
        const auto accepted_raw_change_summary{
            SummarizeTransformedChanges(
                assembled_state,
                raw_state,
                iteration_state.active_index_list)
        };
        LogConvergenceSafeguardAudit(
            options.quiet_mode,
            attempt_number,
            result.progress,
            safeguard_predicates,
            shadow_safeguard_predicates,
            transformed_change_summary,
            raw_fixed_point_change_summary,
            accepted_shadow_change_summary,
            raw_shadow_change_summary,
            EvaluateConvergenceStationarityAudit(
                raw_iteration_result.health_by_key),
            result.diagnostics,
            raw_iteration_result.block_activity,
            raw_iteration_result.assessment_by_atom,
            raw_iteration_result.local_refit_status_by_atom,
            GetMaximumTransformedChange(accepted_raw_change_summary) == 0.0,
            assembled_uses_polish,
            has_quarantine_transition,
            objective_domain_changed);
    }

    iteration_state.previous_state = std::move(assembled_state);
    iteration_state.previous_polish_provenance = std::move(assembled_polish_provenance);
    iteration_state.unchanged_state_exhausted_key_list.clear();
    return result;
}

void ApplyFitState(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitState & iteration_state)
{
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, iteration_state)
    };

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto adjusted_sampling_entries{
            BuildSecondStageAdjustedSamples(context.at(i), model_snapshot)
        };
        analysis.ApplyAtomLocalSecondStageResult(
            *context.at(i).atom,
            iteration_state.at(i),
            std::move(adjusted_sampling_entries));
    }
}

static void LogFinalDependencyPolish(
    bool quiet_mode,
    const FinalDependencyPolishResult & polish_result)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    const auto & diagnostic{ polish_result.diagnostic };
    std::ostringstream message;
    message << std::scientific << std::setprecision(2)
        << "Final dependency polish: components/attempted/accepted/fallback="
        << diagnostic.component_count << "/"
        << diagnostic.attempted_component_count << "/"
        << diagnostic.accepted_component_count << "/"
        << diagnostic.fallback_component_count
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
    message << ", accepted=" << (polish_result.accepted ? "yes" : "no")
        << ", elapsed_ms=" << std::fixed << std::setprecision(3)
        << diagnostic.elapsed_milliseconds << ".";
    Logger::Log(LogLevel::Info, message.str());

    if (Logger::GetLogLevel() < LogLevel::Debug) return;
    for (std::size_t position = 0;
        position < diagnostic.component_list.size();
        position++)
    {
        const auto & component{ diagnostic.component_list.at(position) };
        std::ostringstream component_message;
        component_message << std::scientific << std::setprecision(2)
            << "Final dependency polish component " << position + 1
            << ": clusters/atoms/parameters/rounds="
            << component.key_list.size() << "/"
            << component.atom_count << "/"
            << component.parameter_count << "/"
            << component.round_count
            << ", suspicious/symbolic="
            << component.suspicious_candidate_atom_count << "/"
            << component.symbolic_analysis_count
            << ", objective before/after=";
        if (component.objective_before.has_value())
        {
            component_message << *component.objective_before;
        }
        else
        {
            component_message << "-";
        }
        component_message << "/";
        if (component.objective_after.has_value())
        {
            component_message << *component.objective_after;
        }
        else
        {
            component_message << "-";
        }
        component_message
            << ", accepted/fallback="
            << (component.accepted ? "yes" : "no") << "/"
            << (component.fallback ? "yes" : "no")
            << ", elapsed_ms=" << std::fixed << std::setprecision(3)
            << component.elapsed_milliseconds << ".";
        Logger::Log(LogLevel::Debug, component_message.str());
    }
}

static void FinalizeSecondStageState(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & graph_topology,
    IterationState & iteration_state,
    bool use_best_audit_state,
    PerformanceCounters & performance_counters)
{
    const auto & base_state{
        use_best_audit_state && iteration_state.best_audit_state.has_value() ?
            iteration_state.best_audit_state->state :
            iteration_state.previous_state
    };
    const auto final_block_activity{
        iteration_state.quarantine_state.BuildFinalActivity()
    };
    auto polish_result{
        RunFinalDependencyPolish(
            context,
            options,
            graph_topology,
            iteration_state.graph_partition,
            iteration_state.objective_domain,
            final_block_activity,
            iteration_state.trust_region_state,
            base_state,
            iteration_state.boundary_joint_correction_workspace_by_key,
            performance_counters)
    };
    LogFinalDependencyPolish(options.quiet_mode, polish_result);
    if (polish_result.accepted && polish_result.objective.has_value())
    {
        if (use_best_audit_state &&
            iteration_state.best_audit_state.has_value())
        {
            iteration_state.best_audit_state->state = polish_result.state;
            iteration_state.best_audit_state->objective = *polish_result.objective;
            iteration_state.best_audit_state->uses_polish = true;
        }
        else
        {
            const auto previous_state{ iteration_state.previous_state };
            iteration_state.previous_state = polish_result.state;
            if (iteration_state.previous_polish_provenance.size() != context.size())
            {
                iteration_state.previous_polish_provenance.resize(context.size(), 0);
            }
            for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
            {
                if (IsTransformedChangeMaterial(
                        CalculateTransformedChange(
                            iteration_state.previous_state.at(atom_index).mdpde.GetModel(),
                            previous_state.at(atom_index).mdpde.GetModel()),
                        kTransformedChangeTolerance))
                {
                    iteration_state.previous_polish_provenance.at(atom_index) = 1;
                }
            }
            TryUpdateBestAuditState(
                iteration_state.previous_state,
                true,
                iteration_state.accepted_iteration_count,
                *polish_result.objective,
                iteration_state.best_audit_state);
        }
    }
    const auto & final_state{
        use_best_audit_state && iteration_state.best_audit_state.has_value() ?
            iteration_state.best_audit_state->state :
            iteration_state.previous_state
    };
    ApplyFitState(model_object, context, final_state);
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

void AppendAuditSummary(std::ostringstream & stream, const AuditedState & audited_state)
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

void LogQuarantineFallback(
    bool quiet_mode,
    const IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after "
        << iteration_state.accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendQuarantineSummary(
        warning_message,
        iteration_state.quarantine_state.summary,
        iteration_state.quarantine_state.TargetCount());
    AppendOffsetSummary(warning_message, iteration_state.previous_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogConverged(
    bool quiet_mode,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << iteration_state.accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_stats.percentile_list.at(kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_stats.percentile_list.at(kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_stats.percentile_list.at(kOffsetToPeakRatioChangeIndex);
    AppendOffsetSummary(message, iteration_state.previous_state);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogMaximumIterations(bool quiet_mode, const IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendQuarantineSummary(
        warning_message,
        iteration_state.quarantine_state.summary,
        iteration_state.quarantine_state.TargetCount());
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
    const IterationState & iteration_state,
    std::string_view stop_reason,
    bool final_uses_best_audit)
{
    if (quiet_mode) return;

    const auto & best_audit_state{ iteration_state.best_audit_state };
    const auto final_uses_polish{
        final_uses_best_audit && best_audit_state.has_value() ?
            best_audit_state->uses_polish :
            UsesPolish(iteration_state.previous_polish_provenance)
    };
    Logger::FinishProgressLine();
    std::ostringstream message;
    message << " Second-Stage Local Fitting Summary : \n"
        << " - accepted_iterations = " << iteration_state.accepted_iteration_count << "\n"
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
    message << " - stop_reason = " << stop_reason << "\n"
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
    message << " - final_uses_polish = " << (final_uses_polish ? "yes" : "no") << "\n";
    message << " - final_state_source = " << (final_uses_best_audit ? "best-audit" : "latest-validated") << "\n";
    Logger::Log(LogLevel::Info, message.str());
}

static bool RunSecondStageIterations(ModelObject & model_object, const FitOptions & options)
{
    if (options.enable_second_stage_dependency_polish &&
        options.second_stage_dependency_polish_max_iterations == 0)
    {
        throw std::invalid_argument(
            "Second-stage dependency polish maximum iterations must be positive when enabled.");
    }
    auto context{ BuildSecondStageContext(model_object, options) };
    StoreSecondStageNeighborCounts(model_object, context);
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    FitState initial_state;
    {
        auto initial_state_build_result{ BuildInitialFitState(context) };
        if (initial_state_build_result.failure != SecondStageInitialStateBuildResult::Failure::None)
        {
            if (!options.quiet_mode)
            {
                const auto unselected_seed_failure{
                    initial_state_build_result.failure == SecondStageInitialStateBuildResult::Failure::UnselectedSeedUnavailable
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
        LogSecondStageSeedSelections(
            initial_state_build_result.selection_record_list,
            options.quiet_mode);
        LogUnselectedSecondStageSeedSelections(
            initial_state_build_result.unselected_selection_record_list,
            options.quiet_mode);
        initial_state = std::move(initial_state_build_result.state);
    }
    auto graph_topology{
        BuildSecondStageGraphTopology(context, initial_state, options.quiet_mode)
    };
    LogGraphTopology(graph_topology, options.quiet_mode);
    auto iteration_state{
        BuildIterationState(context, graph_topology, std::move(initial_state), options)
    };
    PerformanceCounters performance_counters{
        options.quiet_mode,
        context,
        iteration_state.solver_workspace_by_key,
        iteration_state.boundary_joint_correction_workspace_by_key
    };
    if (iteration_state.best_audit_state.has_value())
    {
        performance_counters.RecordFullStateMaterialization();
    }
    LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode);
    const auto progress_column_widths{ BuildProgressColumnWidths(context.size()) };
    LogProgressHeader(options.quiet_mode, progress_column_widths);

    std::string_view final_stop_reason;
    bool maximum_iterations_reached{ false };
    for (std::size_t iter = 0; iter < kMaximumIterations; iter++)
    {
        if (iteration_state.active_index_list.empty())
        {
            FinalizeSecondStageState(
                model_object,
                context,
                options,
                graph_topology,
                iteration_state,
                false,
                performance_counters);
            if (iteration_state.quarantine_state.HasFailures())
            {
                LogQuarantineFallback(options.quiet_mode, iteration_state);
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
                    "quarantine",
                false);
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
        if (iteration_result.all_rejected_resolution.has_value())
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
                iteration_result.diagnostics.trust_region_update,
                *iteration_result.all_rejected_resolution);
            if (*iteration_result.all_rejected_resolution == AllRejectedResolution::Retry)
            {
                continue;
            }

            final_stop_reason = GetAllRejectedResolutionText(*iteration_result.all_rejected_resolution);
            break;
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
            final_stop_reason = "audit-patience";
            break;
        }

        if (iteration_result.converged)
        {
            FinalizeSecondStageState(
                model_object,
                context,
                options,
                graph_topology,
                iteration_state,
                false,
                performance_counters);
            if (iteration_state.quarantine_state.HasFailures())
            {
                LogQuarantineFallback(options.quiet_mode, iteration_state);
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
            return true;
        }

        if (iter + 1 == kMaximumIterations)
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
        FinalizeSecondStageState(
            model_object,
            context,
            options,
            graph_topology,
            iteration_state,
            audit_state != nullptr,
            performance_counters);
        if (maximum_iterations_reached)
        {
            LogMaximumIterations(options.quiet_mode, iteration_state);
        }
        LogSecondStageSummary(
            options.quiet_mode,
            iteration_state,
            final_stop_reason,
            audit_state != nullptr);
        return true;
    }
    return false;
}

} // namespace

} // namespace rhbm_gem::core::detail

namespace rhbm_gem::core {

bool RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto fitted{ detail::RunSecondStageIterations(model_object, options) };
    if (fitted)
    {
        RunGroupPotentialFitting(model_object, options, FittingStage::Second);
    }
    return fitted;
}

} // namespace rhbm_gem::core
