#include "core/detail/IterationProcess.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <ostream>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

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

struct TerminalSummary
{
    std::size_t suspicious_cluster_count{ 0 };
    std::size_t suspicious_atom_count{ 0 };
    std::size_t joint_offset_failure_cluster_count{ 0 };
    std::size_t joint_offset_failure_atom_count{ 0 };
    std::map<JointOffsetSolveStatus, std::size_t>
        joint_offset_failure_status_count{};

    std::size_t AtomCount() const
    {
        return suspicious_atom_count + joint_offset_failure_atom_count;
    }

    bool HasFailures() const
    {
        return AtomCount() > 0;
    }
};

static std::vector<ClusterKey> AccumulateTerminalFailureSummary(
    const TerminalPersistentFailureMap & terminal_failure_by_key,
    TerminalSummary & terminal_summary);

static void ApplyTerminalFallbackClusters(
    const std::vector<ClusterKey> & terminal_key_list,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    std::vector<char> & terminal_atom_mask,
    FitState & assembled_state,
    PolishProvenance & assembled_polish_provenance);

struct TerminalFailureState
{
    PersistentTerminalFailureStateMap persistent_state_by_key{};
    std::vector<char> terminal_atom_mask{};
    TerminalSummary terminal_summary{};
    TerminalFailureState() = default;

    explicit TerminalFailureState(std::size_t atom_count)
        : terminal_atom_mask(atom_count, 0)
    {
    }

    bool HasFailures() const { return terminal_summary.HasFailures(); }
    std::size_t AtomCount() const { return terminal_summary.AtomCount(); }

    std::vector<std::size_t> BuildEligibleActiveIndexList() const;

    bool IsolatePersistentFailures(
        const std::vector<ClusterKey> & accepted_key_list,
        SuspiciousUpdateMask & suspicious_atom_mask,
        const ClusterHealthMap & health_by_key,
        FitState & assembled_state,
        const FitState & previous_state,
        const PolishProvenance & previous_polish_provenance,
        PolishProvenance & assembled_polish_provenance);
};

constexpr std::size_t kMaximumIterations{ 100 };
constexpr std::size_t kAuditPatience{ 3 };

struct IterationDiagnostics
{
    std::vector<ClusterCandidateDiagnostic>
        accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic>
        rejected_cluster_diagnostic_list{};
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

std::vector<std::size_t> TerminalFailureState::BuildEligibleActiveIndexList()
    const
{
    std::vector<std::size_t> active_index_list;
    active_index_list.reserve(terminal_atom_mask.size());
    for (std::size_t atom_index = 0;
        atom_index < terminal_atom_mask.size();
        atom_index++)
    {
        if (terminal_atom_mask.at(atom_index) == 0)
        {
            active_index_list.emplace_back(atom_index);
        }
    }
    return active_index_list;
}

bool TerminalFailureState::IsolatePersistentFailures(
    const std::vector<ClusterKey> & accepted_key_list,
    SuspiciousUpdateMask & suspicious_atom_mask,
    const ClusterHealthMap & health_by_key,
    FitState & assembled_state,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    PolishProvenance & assembled_polish_provenance)
{
    const auto terminal_failure_by_key{
        UpdatePersistentTerminalFailureState(
            accepted_key_list,
            suspicious_atom_mask,
            health_by_key,
            assembled_state,
            previous_state,
            persistent_state_by_key)
    };
    const auto terminal_key_list{
        AccumulateTerminalFailureSummary(
            terminal_failure_by_key,
            terminal_summary)
    };
    ApplyTerminalFallbackClusters(
        terminal_key_list,
        previous_state,
        previous_polish_provenance,
        terminal_atom_mask,
        assembled_state,
        assembled_polish_provenance);
    if (!terminal_key_list.empty())
    {
        ClearSuspiciousUpdateMaskForClusters(
            terminal_key_list,
            suspicious_atom_mask);
    }
    return !terminal_key_list.empty();
}

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

static void AppendTerminalSummary(
    std::ostream & stream,
    const TerminalSummary & summary)
{
    if (summary.suspicious_atom_count > 0)
    {
        stream << "; terminal suspicious rollback fallback clusters/atoms = "
            << summary.suspicious_cluster_count << "/" << summary.suspicious_atom_count;
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

static std::vector<ClusterKey> AccumulateTerminalFailureSummary(
    const TerminalPersistentFailureMap & terminal_failure_by_key,
    TerminalSummary & terminal_summary)
{
    std::vector<ClusterKey> terminal_key_list;
    terminal_key_list.reserve(terminal_failure_by_key.size());
    for (const auto & [key, reason] : terminal_failure_by_key)
    {
        terminal_key_list.emplace_back(key);
        if (std::holds_alternative<PersistentSuspiciousRollbackReason>(reason))
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
    return terminal_key_list;
}

TerminalPersistentFailureMap UpdatePersistentTerminalFailureState(
    const std::vector<ClusterKey> & accepted_key_list,
    const SuspiciousUpdateMask & suspicious_atom_mask,
    const ClusterHealthMap & health_by_key,
    const FitState & assembled_state,
    const FitState & previous_state,
    PersistentTerminalFailureStateMap & state_by_key)
{
    PersistentTerminalFailureStateMap next_state_by_key;
    TerminalPersistentFailureMap terminal_failure_by_key;
    for (const auto & [key, health] : health_by_key)
    {
        if (std::ranges::find(accepted_key_list, key) ==
            accepted_key_list.end())
        {
            continue;
        }

        auto cluster_suspicious_atom_index_list{
            CollectSuspiciousAtomIndices(key, suspicious_atom_mask)
        };
        PersistentTerminalFailureReason reason;
        if (!cluster_suspicious_atom_index_list.empty())
        {
            reason = std::move(cluster_suspicious_atom_index_list);
        }
        else
        {
            const auto status{ health.joint_offset_status };
            if (!IsJointOffsetSolveHardFailure(status)) continue;
            reason = status;
        }

        const auto transformed_change_summary{
            SummarizeTransformedChanges(assembled_state, previous_state, key)
        };
        if (!IsTransformedPercentileConverged(transformed_change_summary))
        {
            continue;
        }

        PersistentTerminalFailureState next_state{ std::move(reason), 1 };
        const auto previous_iter{ state_by_key.find(key) };
        if (previous_iter != state_by_key.end() && previous_iter->second.reason == next_state.reason)
        {
            next_state.stable_iteration_count = previous_iter->second.stable_iteration_count + 1;
        }

        if (next_state.stable_iteration_count >= kPersistentTerminalFailureIterationLimit)
        {
            terminal_failure_by_key.emplace(key, std::move(next_state.reason));
            continue;
        }
        next_state_by_key.emplace(key, std::move(next_state));
    }
    state_by_key = std::move(next_state_by_key);
    return terminal_failure_by_key;
}

static void ApplyTerminalFallbackClusters(
    const std::vector<ClusterKey> & terminal_key_list,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    std::vector<char> & terminal_atom_mask,
    FitState & assembled_state,
    PolishProvenance & assembled_polish_provenance)
{
    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            terminal_atom_mask.at(atom_index) = 1;
            assembled_state.at(atom_index) = previous_state.at(atom_index);
            assembled_polish_provenance.at(atom_index) = previous_polish_provenance.at(atom_index);
        }
    }
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

static std::string FormatProgressMaximum(double value)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

namespace {

constexpr std::array<std::string_view, 6> kProgressHeaderList
{
    "Try/Acc",
    "Atom A/T",
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

} // namespace

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

static RawIterationResult RunRawIteration(
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

static IterationState BuildIterationState(
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

static IterationResult RunIteration(
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

namespace {

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

void LogTerminalFallback(
    bool quiet_mode,
    const IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after "
        << iteration_state.accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendTerminalSummary(
        warning_message,
        iteration_state.terminal_failure_state.terminal_summary);
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
    AppendTerminalSummary(
        warning_message,
        iteration_state.terminal_failure_state.terminal_summary);
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

static bool RunSecondStageIterations(ModelObject & model_object, const FitOptions & options)
{
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
    const auto graph_topology{
        BuildSecondStageGraphTopology(context, initial_state, options.quiet_mode)
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
            ApplyFitState(model_object, context, iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.HasFailures())
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
        return true;
    }
    return false;
}

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
