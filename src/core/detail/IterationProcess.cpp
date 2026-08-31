#include "core/detail/IterationProcess.hpp"

#include "core/detail/PreparedLocalGaussianFit.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "data/detail/AtomClassifier.hpp"

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

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core::detail {

namespace {

constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{
    2.0 * kNeighborContributionDistanceMax
};
constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
constexpr std::size_t kMaximumIterations{ 100 };
constexpr std::size_t kAuditPatience{ 3 };
constexpr double kLegacyMaximumTransformedChangeTolerance{ 1.0e-3 };
constexpr double kConvergencePercentile{ 0.99 };
constexpr std::array<SecondStageSeedSource, 4> kSecondStageSeedSourceList{
    SecondStageSeedSource::GroupPosterior,
    SecondStageSeedSource::GroupPrior,
    SecondStageSeedSource::GroupMedian,
    SecondStageSeedSource::GlobalMedian
};
constexpr std::array<std::string_view, 6> kProgressHeaderList{
    "Try/Acc",
    "Atom A/Q",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/G"
};

using ProgressColumnWidths = std::array<std::size_t, 6>;
using GraphEdgeSet = std::set<std::pair<std::size_t, std::size_t>>;

namespace iteration_internal {

bool UsesPolish(const PolishProvenance & provenance)
{
    return std::ranges::any_of(
        provenance,
        [](char value) { return value != 0; });
}

std::vector<double> BuildSuspiciousJointOffsetRidgeMultiplierList(
    const SuspiciousUpdateMask & suspicious_mask)
{
    std::vector<double> ridge_multiplier_list(suspicious_mask.size(), 1.0);
    for (std::size_t atom_index = 0; atom_index < suspicious_mask.size(); atom_index++)
    {
        if (suspicious_mask.at(atom_index) != 0)
        {
            ridge_multiplier_list.at(atom_index) =
                kSuspiciousJointOffsetRidgeMultiplier;
        }
    }
    return ridge_multiplier_list;
}

void ResetClusterSolverWorkspace(
    const std::vector<ClusterKey> & cluster_key_list,
    ClusterSolverWorkspaceMap & workspace_by_key)
{
    workspace_by_key.clear();
    for (const auto & key : cluster_key_list)
    {
        workspace_by_key.try_emplace(key);
    }
}

} // namespace iteration_internal

void SetLocalResultOffset(LocalGaussianResult & result, double offset)
{
    result.ols = WithPreservedUncertaintyOffset(result.ols, offset);
    result.mdpde = WithPreservedUncertaintyOffset(result.mdpde, offset);
}

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

struct SecondStageInitializationResult
{
    SecondStageContext context{};
    FitState state{};
    std::vector<int> neighbor_count_list{};
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
        std::span<const ClusterCandidateDiagnostic> accepted_diagnostic_list,
        std::span<const ClusterCandidateDiagnostic> rejected_diagnostic_list,
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

struct IterationDiagnostics
{
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<BoundaryComponentReconciliationDiagnostic> boundary_reconciliation_diagnostic_list{};
    TrustRegionRadiusUpdate trust_region_update{};
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
    double operator_maximum_transformed_change{ 0.0 };
};

struct IterationResult
{
    IterationDiagnostics diagnostics{};
    IterationProgress progress{};
    std::optional<AllRejectedResolution> all_rejected_resolution{};
    bool objective_domain_changed{ false };
    bool converged{ false };
    bool audit_patience_exhausted{ false };
    TransformedChange transformed_change_percentile{};
};

enum class OperatorEndpointStatus : char
{
    NotRequested,
    Available,
    OffsetSolverFailure,
    InvalidOffset,
    ShapeRefitFailure
};

struct FixedPointOperatorEvidence
{
    FitState state{};
    std::vector<OperatorEndpointStatus> shape_status_by_atom{};
    std::vector<OperatorEndpointStatus> offset_status_by_atom{};
    bool shadow_shape_refit_performed{ false };
};

struct IterationProposalResult
{
    FitState operator_proposal_state{};
    FixedPointOperatorEvidence fixed_point_operator{};
    SuspiciousBlockActivity block_activity{};
    std::vector<SuspiciousGaussianAssessment> assessment_by_atom{};
    std::vector<std::optional<RHBMEstimationStatus>> local_refit_status_by_atom{};
    ClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    std::optional<GaussianModel3DWithUncertainty> unrestricted_mdpde{};
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

struct OffsetGroupEntry
{
    ClusterKey cluster_key{};
    ClusterKey atom_index_list{};
};

static void ValidateActiveCoordinateInputs(
    std::size_t atom_count,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity)
{
    if (group_id_by_atom_index.size() != atom_count ||
        block_activity.shape_fixed_atom_mask.size() != atom_count ||
        block_activity.offset_fixed_atom_mask.size() != atom_count ||
        block_activity.hard_failure_atom_mask.size() != atom_count ||
        quarantine_activity.shape_fixed_atom_mask.size() != atom_count ||
        quarantine_activity.offset_fixed_atom_mask.size() != atom_count ||
        quarantine_activity.hard_failure_atom_mask.size() != atom_count)
    {
        throw std::invalid_argument(
            "Active-coordinate convergence inputs are inconsistent.");
    }
}

static std::vector<OffsetGroupEntry> BuildOffsetGroupEntries(
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index)
{
    std::vector<OffsetGroupEntry> result;
    for (const auto & cluster_key : cluster_key_list)
    {
        std::map<std::size_t, ClusterKey> atom_index_list_by_group;
        for (const auto atom_index : cluster_key)
        {
            atom_index_list_by_group[group_id_by_atom_index.at(atom_index)]
                .emplace_back(atom_index);
        }
        for (auto & [group_id, atom_index_list] : atom_index_list_by_group)
        {
            static_cast<void>(group_id);
            result.emplace_back(OffsetGroupEntry{
                cluster_key,
                std::move(atom_index_list)
            });
        }
    }
    return result;
}

static ConvergenceCertificate SummarizeFixedPointOperator(
    const FixedPointOperatorEvidence & evidence,
    const FitState & previous_state,
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index)
{
    SuspiciousBlockActivity nominal_activity{
        SuspiciousUpdateMask(previous_state.size(), 0),
        SuspiciousUpdateMask(previous_state.size(), 0),
        SuspiciousUpdateMask(previous_state.size(), 0)
    };
    ConvergenceCertificate result;
    result.operator_nominal_population = BuildActiveCoordinatePopulation(
        atom_index_list,
        cluster_key_list,
        group_id_by_atom_index,
        nominal_activity,
        nominal_activity);
    result.operator_shadow_shape_refit_performed =
        evidence.shadow_shape_refit_performed;

    std::vector<TransformedChange> change_list;
    change_list.reserve(previous_state.size());
    for (std::size_t atom_index = 0;
        atom_index < previous_state.size(); atom_index++)
    {
        auto change{ CalculateTransformedChange(
            GetFitModel(evidence.state, atom_index),
            GetFitModel(previous_state, atom_index)) };
        if (evidence.shape_status_by_atom.at(atom_index) !=
            OperatorEndpointStatus::Available)
        {
            change.at(kLogPeakHeightChangeIndex) =
                std::numeric_limits<double>::infinity();
            change.at(kLogWidthChangeIndex) =
                std::numeric_limits<double>::infinity();
        }
        if (evidence.offset_status_by_atom.at(atom_index) !=
            OperatorEndpointStatus::Available)
        {
            change.at(kOffsetToPeakRatioChangeIndex) =
                std::numeric_limits<double>::infinity();
        }
        change_list.emplace_back(std::move(change));
    }
    result.operator_nominal_residual = SummarizeActiveDofChanges(
        change_list, result.operator_nominal_population);

    for (const auto atom_index : atom_index_list)
    {
        const auto shape_available{
            evidence.shape_status_by_atom.at(atom_index) ==
            OperatorEndpointStatus::Available
        };
        for (const auto coordinate : {
                kLogPeakHeightChangeIndex, kLogWidthChangeIndex })
        {
            if (!shape_available)
            {
                result.operator_unavailable_count.at(coordinate)++;
            }
            else if (std::abs(change_list.at(atom_index).at(coordinate)) >=
                kLegacyMaximumTransformedChangeTolerance)
            {
                result.operator_tail_count.at(coordinate)++;
            }
        }
        if (!shape_available)
        {
            const auto status{ evidence.shape_status_by_atom.at(atom_index) };
            if (status == OperatorEndpointStatus::OffsetSolverFailure)
            {
                result.operator_unavailable_reason_count.at(0)++;
            }
            else if (status == OperatorEndpointStatus::InvalidOffset)
            {
                result.operator_unavailable_reason_count.at(1)++;
            }
            else
            {
                result.operator_unavailable_reason_count.at(2)++;
            }
        }
    }
    for (const auto & group :
        result.operator_nominal_population.active_offset_group_atom_index_list)
    {
        bool available{ true };
        double maximum_change{ 0.0 };
        for (const auto atom_index : group)
        {
            available = available &&
                evidence.offset_status_by_atom.at(atom_index) ==
                    OperatorEndpointStatus::Available;
            maximum_change = std::max(
                maximum_change,
                std::abs(change_list.at(atom_index).at(
                    kOffsetToPeakRatioChangeIndex)));
        }
        if (!available)
        {
            result.operator_unavailable_count.at(
                kOffsetToPeakRatioChangeIndex)++;
            const auto failed_atom{ *std::ranges::find_if(
                group,
                [&](std::size_t atom_index)
                {
                    return evidence.offset_status_by_atom.at(atom_index) !=
                        OperatorEndpointStatus::Available;
                }) };
            const auto status{
                evidence.offset_status_by_atom.at(failed_atom)
            };
            if (status == OperatorEndpointStatus::OffsetSolverFailure)
            {
                result.operator_unavailable_reason_count.at(0)++;
            }
            else
            {
                result.operator_unavailable_reason_count.at(1)++;
            }
        }
        else if (maximum_change >= kLegacyMaximumTransformedChangeTolerance)
        {
            result.operator_tail_count.at(kOffsetToPeakRatioChangeIndex)++;
        }
    }
    return result;
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

static SecondStageInitializationResult BuildSecondStageInitialization(
    const ModelObject & model_object,
    const FitOptions & options)
{
    SecondStageInitializationResult build_result;
    auto & context{ build_result.context };
    auto & state{ build_result.state };

    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.selected_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.selected_atom_list.emplace_back(AtomContext{ atom });
    }
    state.resize(context.size());
    std::vector<std::optional<GaussianModel3DWithUncertainty>>
        group_prior_list(context.size());

    std::unordered_map<GroupKey, std::size_t> selected_group_id_by_key;
    selected_group_id_by_key.reserve(context.size());
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        atom_index_map.emplace(context.at(i).atom, i);
    }
    for (std::size_t atom_index = 0;
        atom_index < context.size();
        atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto group_key{ data_internal::GetGroupKey(atom) };
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        atom_context.raw_sampling_entries =
            local_view.GetRawSamplingEntries(false);
        state.at(atom_index) =
            local_view.GetGaussianResult(FittingStage::Second);
        group_prior_list.at(atom_index) =
            analysis_view.FindAtomGroupPriorWithUncertainty(
                FittingStage::Second,
                *atom);
        atom_context.alpha_r =
            local_view.GetAlphaR(FittingStage::Second);
        atom_context.refit_design = PreparedLocalGaussianDesign{
            atom_context.raw_sampling_entries,
            options.distance_min,
            options.distance_max
        };
        auto [group_iter, inserted]{
            selected_group_id_by_key.emplace(
                group_key,
                context.selected_atom_index_list_by_group.size())
        };
        if (inserted)
        {
            context.selected_atom_index_list_by_group.emplace_back();
        }
        atom_context.group_id = group_iter->second;
        context.selected_atom_index_list_by_group.at(
            atom_context.group_id).emplace_back(atom_index);
    }

    std::unordered_map<const AtomObject *, std::size_t>
        unselected_atom_index_map;
    std::vector<int> unselected_atom_serial_id_list;
    build_result.neighbor_count_list.reserve(context.size());
    for (std::size_t atom_index = 0;
        atom_index < context.size();
        atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{
            atom->FindNeighborAtoms(kNeighborAtomSearchRange)
        };
        std::unordered_set<const AtomObject *> neighbor_atom_set;

        atom_context.neighbor_atom_sample_offset_list.reserve(
            atom_context.raw_sampling_entries.size() + 1);
        atom_context.neighbor_atom_sample_offset_list.emplace_back(0);
        atom_context.neighbor_atom_sample_list.reserve(
            atom_context.raw_sampling_entries.size() *
            neighbor_atom_list.size());
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{
                atom_context.raw_sampling_entries.at(sample_index)
            };
            for (auto * neighbor_atom : neighbor_atom_list)
            {
                if (options.exclude_hydrogen &&
                    neighbor_atom->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto distance{
                    array_helper::ComputeNorm(
                        sample.point.position,
                        neighbor_atom->GetPositionRef())
                };
                if (distance > kNeighborContributionDistanceMax) continue;
                neighbor_atom_set.emplace(neighbor_atom);

                const auto selected_iter{
                    atom_index_map.find(neighbor_atom)
                };
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

                auto contributor_iter{
                    unselected_atom_index_map.find(neighbor_atom)
                };
                if (contributor_iter == unselected_atom_index_map.end())
                {
                    const auto contributor_index{
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
                            selected_group_iter ==
                                selected_group_id_by_key.end() ?
                                std::nullopt :
                                std::optional<std::size_t>{
                                    selected_group_iter->second
                                }
                        });
                    unselected_atom_serial_id_list.emplace_back(
                        neighbor_atom->GetSerialID());
                    contributor_iter = unselected_atom_index_map.emplace(
                        neighbor_atom,
                        contributor_index).first;
                }
                atom_context.neighbor_atom_sample_list.emplace_back(
                    NeighborAtomSample{
                        false,
                        contributor_iter->second,
                        distance
                    });
            }
            atom_context.neighbor_atom_sample_offset_list.emplace_back(
                atom_context.neighbor_atom_sample_list.size());
        }
        build_result.neighbor_count_list.emplace_back(
            static_cast<int>(neighbor_atom_set.size()));
    }

    std::vector<std::vector<GaussianModel3D>> models_by_group(
        context.selected_atom_index_list_by_group.size());
    std::vector<GaussianModel3D> global_models;
    global_models.reserve(context.size());

    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        const auto & result{ state.at(i) };
        const auto direct_selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    result.posterior,
                    group_prior_list.at(i),
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
                    group_prior_list.at(i),
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value())
        {
            build_result.failure =
                SecondStageInitializationResult::Failure::
                    SelectedSeedUnavailable;
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
            build_result.failure =
                SecondStageInitializationResult::Failure::
                    UnselectedSeedUnavailable;
            return build_result;
        }

        unselected_atom_contributor.initial_seed = selection->model.GetModel();
        build_result.unselected_selection_record_list.emplace_back(
            UnselectedSecondStageSeedSelectionRecord{
                unselected_atom_serial_id_list.at(i),
                selection->source,
                selection->model.GetModel()
            });
    }
    return build_result;
}

static void StoreSecondStageNeighborCounts(
    ModelObject & model_object,
    const SecondStageContext & context,
    const std::vector<int> & neighbor_count_list)
{
    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t atom_index = 0;
        atom_index < context.size();
        atom_index++)
    {
        analysis.SetAtomLocalNeighborCountForPeeling(
            *context.at(atom_index).atom,
            neighbor_count_list.at(atom_index));
    }
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

    std::array<std::size_t, kSecondStageSeedSourceList.size()> source_count{};
    for (const auto & record : selection_record_list)
    {
        source_count.at(static_cast<std::size_t>(record.source))++;
    }

    std::ostringstream summary;
    summary << "Selected second-stage initial seeds = "
        << selection_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < kSecondStageSeedSourceList.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedSourceText(kSecondStageSeedSourceList.at(i))
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
                SetLocalResultOffset(assembled_state.at(atom_index), previous_offset);
                continue;
            }
            const auto assembled_offset{
                assembled_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            assembled_state.at(atom_index).ols = WithPreservedUncertaintyOffset(
                previous_state.at(atom_index).ols,
                assembled_offset);
            assembled_state.at(atom_index).mdpde = WithPreservedUncertaintyOffset(
                previous_state.at(atom_index).mdpde,
                assembled_offset);
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

static bool HasAcceptedMaterialTargetChange(
    const FitState & assembled_state,
    const FitState & previous_state,
    const QuarantineTarget & target)
{
    for (const auto atom_index : target.atom_index_list)
    {
        const auto change{ CalculateTransformedChange(
            assembled_state.at(atom_index).mdpde.GetModel(),
            previous_state.at(atom_index).mdpde.GetModel()) };
        if (target.kind == QuarantineTargetKind::ShapeAtom)
        {
            if (std::max(change.at(0), change.at(1)) >=
                kTransformedChangeTolerance)
            {
                return true;
            }
        }
        else if (target.kind == QuarantineTargetKind::OffsetGroup)
        {
            if (change.at(2) >= kTransformedChangeTolerance) return true;
        }
        else if (IsTransformedChangeMaterial(change, kTransformedChangeTolerance))
        {
            return true;
        }
    }
    return false;
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

static bool IsGuardSafeNonMaterialSolverQualifiedEndpoint(
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
            !health_iter->second.IsSolverQualified() ||
            atom_index >= assessment_by_atom.size() ||
            assessment_by_atom[atom_index].IsSuspicious())
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
    std::span<const ClusterCandidateDiagnostic> accepted_diagnostic_list,
    std::span<const ClusterCandidateDiagnostic> rejected_diagnostic_list,
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
    for (const auto & [key, health] : health_by_key)
    {
        if (!IsJointOffsetSolveHardFailure(health.joint_offset_status)) continue;
        observation_list.emplace_back(QuarantineFailureObservation{
            QuarantineTarget{ QuarantineTargetKind::HardFailureCluster, key },
            health.joint_offset_status
        });
    }
    const auto append_terminal_observations = [&](const auto & diagnostic_list)
    {
        for (const auto & diagnostic : diagnostic_list)
        {
            for (const auto & terminal :
                diagnostic.attempt.terminal_diagnostic_list)
            {
                const auto reason{ terminal.reason };
                if (reason == StabilizationTerminalReason::None) continue;
                const StabilizationTerminalFailure failure{
                    reason,
                    terminal.guard_reason
                };
                if (reason == StabilizationTerminalReason::GuardInfeasible &&
                    terminal.guard_atom_index.has_value() &&
                    terminal.guard_mode.has_value())
                {
                    const auto atom_index{ *terminal.guard_atom_index };
                    if (*terminal.guard_mode == SuspiciousUpdateMode::OffsetOnly)
                    {
                        std::vector<std::size_t> group_atom_index_list;
                        const auto group_id{ context.at(atom_index).group_id };
                        for (std::size_t member_index = 0;
                            member_index < context.size(); member_index++)
                        {
                            if (context.at(member_index).group_id == group_id)
                            {
                                group_atom_index_list.emplace_back(member_index);
                            }
                        }
                        observation_list.emplace_back(QuarantineFailureObservation{
                            QuarantineTarget{
                                QuarantineTargetKind::OffsetGroup,
                                std::move(group_atom_index_list)
                            },
                            failure
                        });
                    }
                    else
                    {
                        observation_list.emplace_back(QuarantineFailureObservation{
                            QuarantineTarget{
                                QuarantineTargetKind::ShapeAtom,
                                { atom_index }
                            },
                            failure
                        });
                    }
                    continue;
                }
                observation_list.emplace_back(QuarantineFailureObservation{
                    QuarantineTarget{
                        QuarantineTargetKind::HardFailureCluster,
                        diagnostic.key
                    },
                    failure
                });
            }
        }
    };
    append_terminal_observations(accepted_diagnostic_list);
    append_terminal_observations(rejected_diagnostic_list);

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
            HasAcceptedMaterialTargetChange(
                assembled_state,
                previous_state,
                target)
        };
        if (!has_affecting_observation &&
            (accepted_material_proposal ||
                IsGuardSafeNonMaterialSolverQualifiedEndpoint(
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

static std::string_view GetStabilizationTerminalReasonText(
    StabilizationTerminalReason reason)
{
    switch (reason)
    {
    case StabilizationTerminalReason::None: return "none";
    case StabilizationTerminalReason::GuardInfeasible: return "guard-infeasible";
    case StabilizationTerminalReason::ObjectiveExhausted: return "objective-exhausted";
    case StabilizationTerminalReason::InvalidCandidate: return "invalid-candidate";
    }
    return "unknown";
}

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
static std::string_view GetTrustModelPredictionStatusText(
    TrustModelPredictionStatus status)
{
    switch (status)
    {
    case TrustModelPredictionStatus::Available: return "available";
    case TrustModelPredictionStatus::NonmaterialStep: return "nonmaterial-step";
    case TrustModelPredictionStatus::ObjectiveUnavailable: return "objective-unavailable";
    case TrustModelPredictionStatus::ModelUnavailable: return "model-unavailable";
    case TrustModelPredictionStatus::ResidualUnavailable: return "residual-unavailable";
    case TrustModelPredictionStatus::Nonfinite: return "nonfinite";
    case TrustModelPredictionStatus::NonpositivePrediction: return "nonpositive-prediction";
    case TrustModelPredictionStatus::NonmaterialPrediction: return "nonmaterial-prediction";
    }
    return "model-unavailable";
}

static std::string_view GetTrustModelCandidateSourceText(
    TrustModelCandidateSource source)
{
    return source == TrustModelCandidateSource::Polish ? "polish" : "base";
}

static std::string_view GetTrustModelTrialDispositionText(
    TrustModelTrialDisposition disposition)
{
    return disposition == TrustModelTrialDisposition::Accepted ?
        "accepted" : "objective-rejected";
}

static std::string_view GetTrustRegionRadiusActionText(
    TrustRegionRadiusAction action)
{
    switch (action)
    {
    case TrustRegionRadiusAction::Keep: return "keep";
    case TrustRegionRadiusAction::Grow: return "grow";
    case TrustRegionRadiusAction::Shrink: return "shrink";
    }
    return "keep";
}

static void AppendTrustModelOptionalValue(
    std::ostringstream & stream,
    const std::optional<double> & value)
{
    if (value.has_value()) stream << *value;
    else stream << "-";
}

static void LogTrustModelShadowDiagnostics(
    bool quiet_mode,
    const IterationDiagnostics & diagnostics,
    const IterationProgress & progress)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    const auto log_records = [&](const auto & diagnostic_list, std::string_view disposition)
    {
        for (const auto & cluster_diagnostic : diagnostic_list)
        {
            const auto & funnel{ cluster_diagnostic.trust_model_candidate_funnel };
            std::ostringstream funnel_message;
            funnel_message
                << "Trust-model funnel: schema=1"
                << ", try=" << progress.attempt_number
                << ", acc=" << progress.accepted_iteration_count
                << ", atoms=" << cluster_diagnostic.key.size()
                << ", key-first=" << cluster_diagnostic.key.front()
                << ", key-last=" << cluster_diagnostic.key.back()
                << ", disposition=" << disposition
                << ", generated=" << funnel.generated_count
                << ", invalid=" << funnel.invalid_count
                << ", trust-skipped=" << funnel.trust_skipped_count
                << ", guard-rejected=" << funnel.guard_rejected_count
                << ", nonmaterial=" << funnel.nonmaterial_count
                << ", objective-evaluated=" << funnel.objective_evaluated_count
                << ", polish-objective-evaluated="
                << funnel.polish_objective_evaluated_count;
            Logger::Log(LogLevel::Debug, funnel_message.str());

            for (const auto & diagnostic :
                cluster_diagnostic.trust_model_shadow_trial_list)
            {
                std::ostringstream message;
                message << std::scientific << std::setprecision(17)
                    << "Trust-model shadow: schema=2"
                    << ", try=" << progress.attempt_number
                    << ", acc=" << progress.accepted_iteration_count
                    << ", atoms=" << cluster_diagnostic.key.size()
                    << ", key-first=" << cluster_diagnostic.key.front()
                    << ", key-last=" << cluster_diagnostic.key.back()
                    << ", disposition=" << disposition
                    << ", boundary-touched=" << cluster_diagnostic.boundary_touched
                    << ", boundary-rescued=" << cluster_diagnostic.boundary_rescued
                    << ", readiness-eligible=" << diagnostic.readiness_eligible
                    << ", final-local-candidate=" << diagnostic.final_local_candidate
                    << ", status=" << GetTrustModelPredictionStatusText(diagnostic.status)
                    << ", source=" << GetTrustModelCandidateSourceText(
                        diagnostic.candidate_source)
                    << ", search-pass=" << diagnostic.search_pass
                    << ", trial=" << diagnostic.trial_number
                    << ", factor=" << diagnostic.factor
                    << ", trial-disposition=" << GetTrustModelTrialDispositionText(
                        diagnostic.trial_disposition)
                    << ", rejected-by-previous=" << diagnostic.rejected_by_previous
                    << ", rejected-by-best=" << diagnostic.rejected_by_best
                    << ", rejected-by-strict-polish="
                    << diagnostic.rejected_by_strict_polish
                    << ", step-norm=" << diagnostic.step_norm
                    << ", actual-reduction=";
                AppendTrustModelOptionalValue(message, diagnostic.actual_reduction);
                message << ", polish-reduction=";
                AppendTrustModelOptionalValue(message, diagnostic.polish_reduction);
                message << ", predicted-residual-reduction=";
                AppendTrustModelOptionalValue(
                    message, diagnostic.predicted_residual_reduction);
                message << ", predicted-penalty-reduction=";
                AppendTrustModelOptionalValue(
                    message, diagnostic.predicted_penalty_reduction);
                message << ", predicted-reduction=";
                AppendTrustModelOptionalValue(message, diagnostic.predicted_reduction);
                message << ", rho=";
                AppendTrustModelOptionalValue(message, diagnostic.rho);
                message
                    << ", boundary-utilization=" << diagnostic.boundary_utilization
                    << ", current-action="
                    << GetTrustRegionRadiusActionText(diagnostic.current_action)
                    << ", shadow-action=";
                if (diagnostic.shadow_action.has_value())
                {
                    message << GetTrustRegionRadiusActionText(
                        *diagnostic.shadow_action);
                }
                else
                {
                    message << "suppressed";
                }
                message
                    << ", objective-backtracked=" << diagnostic.objective_backtracked
                    << ", unselected-dependencies="
                    << diagnostic.unselected_dependency_count
                    << ", elapsed-ms=" << diagnostic.elapsed_milliseconds;
                Logger::Log(LogLevel::Debug, message.str());
            }
        }
    };
    Logger::FinishProgressLine();
    log_records(diagnostics.accepted_cluster_diagnostic_list, "accepted");
    log_records(diagnostics.rejected_cluster_diagnostic_list, "rejected");
}
#endif

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
            << "  unified accepted factor = ";
        if (diagnostic.accepted_factor.has_value())
        {
            message << *diagnostic.accepted_factor;
        }
        else
        {
            message << "-";
        }
        message
            << ", trials[total/invalid/trust-skipped/guard-rejected/objective-rejected] = "
            << diagnostic.trial_count << "/"
            << diagnostic.invalid_trial_count << "/"
            << diagnostic.trust_skipped_trial_count << "/"
            << diagnostic.guard_rejected_trial_count << "/"
            << diagnostic.objective_rejected_trial_count
            << ", terminal = "
            << (diagnostic.terminal_diagnostic_list.empty() ? "none" :
                GetStabilizationTerminalReasonText(
                    diagnostic.terminal_diagnostic_list.back().reason))
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
        Logger::Log(LogLevel::Debug, message.str());
    }
}

static std::string_view GetAllRejectedResolutionText(
    AllRejectedResolution resolution)
{
    switch (resolution)
    {
    case AllRejectedResolution::MaximumIterations:
        return "maximum-iterations";
    case AllRejectedResolution::BacktrackingExhausted:
        return "all-rejected-backtracking-exhausted";
    }
    return "all-rejected-backtracking-exhausted";
}

static void LogAllRejectedResolution(
    bool quiet_mode,
    const TrustRegionRadiusUpdate & trust_region_update,
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
        << ", radius-changed/radius-saturated = "
        << trust_region_update.changed_key_list.size() << "/"
        << trust_region_update.saturated_key_list.size() << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

static void LogAcceptedCandidateSearchDiagnostics(
    bool quiet_mode,
    const IterationDiagnostics & diagnostics)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }
    const auto has_local_search{
        std::any_of(
            diagnostics.accepted_cluster_diagnostic_list.begin(),
            diagnostics.accepted_cluster_diagnostic_list.end(),
            [](const ClusterCandidateDiagnostic & diagnostic)
            {
                return diagnostic.attempt.trial_count > 1;
            })
    };
    const auto has_boundary_reconciliation_diagnostic{
        !diagnostics.boundary_reconciliation_diagnostic_list.empty()
    };
    if (!has_local_search && !has_boundary_reconciliation_diagnostic) return;
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
        if (diagnostic.trial_count <= 1) continue;
        std::ostringstream message;
        message
            << "Accepted local fitting candidate search: atoms = " << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/" << cluster_diagnostic.key.back()
            << ", trials/factor = " << diagnostic.trial_count << "/";
        if (diagnostic.accepted_factor.has_value())
        {
            message << *diagnostic.accepted_factor;
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
            FormatProgressMaximum(progress.operator_maximum_transformed_change)
    };
    Logger::ProgressLine(FormatProgressRow(column_widths, cell_list));
}

static void AppendAuditValues(
    std::ostringstream & message,
    const TransformedChange & value_list)
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
    const ConvergenceCertificate & certificate,
    const TransformedChangeSummary & accepted_production_change,
    const ActiveCoordinatePopulation & active_population,
    const IterationDiagnostics & diagnostics,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    bool accepted_equals_operator,
    bool assembled_uses_polish,
    bool has_quarantine_transition,
    bool objective_domain_changed)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;

    const auto & solver_qualification{ certificate.solver_qualification };

    const auto accepted_limited_count{ static_cast<std::size_t>(std::ranges::count_if(
        diagnostics.accepted_cluster_diagnostic_list,
        [](const auto & diagnostic)
        {
            return diagnostic.attempt.accepted_factor.has_value() &&
                *diagnostic.attempt.accepted_factor < 1.0;
        })) };
    const auto sum_trial_diagnostic = [&](const auto projection)
    {
        std::size_t count{ 0 };
        for (const auto & diagnostic : diagnostics.accepted_cluster_diagnostic_list)
        {
            count += projection(diagnostic.attempt);
        }
        for (const auto & diagnostic : diagnostics.rejected_cluster_diagnostic_list)
        {
            count += projection(diagnostic.attempt);
        }
        return count;
    };
    const auto trial_count{ sum_trial_diagnostic(
        [](const auto & diagnostic) { return diagnostic.trial_count; }) };
    const auto invalid_trial_count{ sum_trial_diagnostic(
        [](const auto & diagnostic) { return diagnostic.invalid_trial_count; }) };
    const auto trust_skipped_trial_count{ sum_trial_diagnostic(
        [](const auto & diagnostic) { return diagnostic.trust_skipped_trial_count; }) };
    const auto guard_rejected_trial_count{ sum_trial_diagnostic(
        [](const auto & diagnostic) { return diagnostic.guard_rejected_trial_count; }) };
    const auto objective_rejected_trial_count{ sum_trial_diagnostic(
        [](const auto & diagnostic) { return diagnostic.objective_rejected_trial_count; }) };
    const auto terminal_count{ sum_trial_diagnostic(
        [](const auto & diagnostic)
        {
            return diagnostic.terminal_diagnostic_list.size();
        }) };
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

    const auto selected_atom_count{
        progress.active_atom_count + progress.quarantine_atom_count
    };
    const auto ratio = [selected_atom_count](std::size_t count)
    {
        return selected_atom_count == 0 ? 0.0 :
            static_cast<double>(count) / static_cast<double>(selected_atom_count);
    };
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "Convergence safeguard audit: schema=9"
        << ", certificate-definition=1, try=" << attempt_number
        << ", acc=" << progress.accepted_iteration_count
        << ", atoms=" << progress.active_atom_count + progress.quarantine_atom_count
        << ", quarantine=" << progress.quarantine_atom_count
        << ", accepted-active-population=";
    AppendAuditPopulation(message, accepted_production_change.population_size_list);
    message << ", operator-nominal-population=";
    AppendAuditPopulation(
        message, certificate.operator_nominal_residual.population_size_list);
    message
        << ", certificate[solver/accepted-p99/operator-complete/operator-p99/invariants/orthogonal/production]="
        << solver_qualification.solver_qualified << "/"
        << certificate.AcceptedPercentilePassed() << "/"
        << certificate.OperatorComplete() << "/"
        << certificate.OperatorPercentilePassed() << "/"
        << certificate.InvariantsClear() << "/"
        << certificate.blockers.Clear() << "/"
        << certificate.ProductionConverged()
        << ", accepted-active-p99=";
    AppendAuditValues(
        message,
        accepted_production_change.percentile_list);
    message << ", accepted-active-max=";
    AppendAuditValues(message, accepted_production_change.maximum_list);
    message << ", operator-nominal-residual-p99=";
    AppendAuditValues(
        message,
        certificate.operator_nominal_residual.percentile_list);
    message << ", operator-nominal-residual-max=";
    AppendAuditValues(message, certificate.operator_nominal_residual.maximum_list);
    message << ", operator-nominal-unavailable[height/width/offset]=";
    AppendAuditPopulation(message, certificate.operator_unavailable_count);
    message << ", operator-nominal-unavailable-reasons[offset-solver/invalid-offset/shape-refit]=";
    AppendAuditPopulation(message, certificate.operator_unavailable_reason_count);
    message << ", operator-nominal-tail[height/width/offset]=";
    AppendAuditPopulation(message, certificate.operator_tail_count);
    message
        << ", residual-state=" << GetFixedPointResidualInterpretationText(
            EvaluateFixedPointResidualInterpretation(
                certificate.OperatorComplete(),
                solver_qualification.solver_qualified,
                accepted_production_change,
                certificate.operator_nominal_residual))
        << ", operator-shadow-refit="
        << certificate.operator_shadow_shape_refit_performed
        << ", accepted-equals-operator=" << accepted_equals_operator
        << ", unified-search[trials/invalid/trust-skipped/guard-rejected/objective-rejected/accepted-limited/terminal]="
        << trial_count << "/"
        << invalid_trial_count << "/"
        << trust_skipped_trial_count << "/"
        << guard_rejected_trial_count << "/"
        << objective_rejected_trial_count << "/"
        << accepted_limited_count << "/"
        << terminal_count
        << ", path[limited/polish/boundary/rescue]="
        << accepted_limited_count << "/"
        << assembled_uses_polish << "/"
        << accepted_boundary_count << "/"
        << rescued_boundary_count
        << ", joint-status[converged/system-build/empty/initial-solve/irls-solve/objective-deteriorated/max-iter]="
        << solver_qualification.joint_offset_status_count.at(0) << "/"
        << solver_qualification.joint_offset_status_count.at(1) << "/"
        << solver_qualification.joint_offset_status_count.at(2) << "/"
        << solver_qualification.joint_offset_status_count.at(3) << "/"
        << solver_qualification.joint_offset_status_count.at(4) << "/"
        << solver_qualification.joint_offset_status_count.at(5) << "/"
        << solver_qualification.joint_offset_status_count.at(6)
        << ", qualification[production/solver/restricted/all-fixed/active-shape/solver-shape/soft-shape/hard-shape/fixed-shape/quarantine-shape/active-offset/solver-offset/soft-offset/hard-offset/fixed-offset/quarantine-offset/mixed-offset]="
        << solver_qualification.production_qualified << "/"
        << solver_qualification.solver_qualified << "/"
        << solver_qualification.restricted_active_set << "/"
        << solver_qualification.all_fixed << "/"
        << solver_qualification.active_shape_count << "/"
        << solver_qualification.qualified_shape_count << "/"
        << solver_qualification.soft_unqualified_shape_count << "/"
        << solver_qualification.hard_failure_shape_count << "/"
        << solver_qualification.fixed_shape_count << "/"
        << solver_qualification.quarantined_shape_count << "/"
        << solver_qualification.active_offset_group_count << "/"
        << solver_qualification.qualified_offset_group_count << "/"
        << solver_qualification.soft_unqualified_offset_group_count << "/"
        << solver_qualification.hard_failure_offset_group_count << "/"
        << solver_qualification.fixed_offset_group_count << "/"
        << solver_qualification.quarantined_offset_group_count << "/"
        << solver_qualification.mixed_offset_group_count
        << ", local-status[success/max-iter/single/insufficient/numerical/unavailable]="
        << local_refit_status_count.at(0) << "/"
        << local_refit_status_count.at(1) << "/"
        << local_refit_status_count.at(2) << "/"
        << local_refit_status_count.at(3) << "/"
        << local_refit_status_count.at(4) << "/"
        << unavailable_local_refit_status_count
        << ", offset-groups[total/active/fixed/quarantine/mixed]="
        << active_population.total_offset_group_count << "/"
        << solver_qualification.active_offset_group_count << "/"
        << active_population.fixed_offset_group_count << "/"
        << active_population.quarantined_offset_group_count << "/"
        << active_population.mixed_offset_group_count
        << ", ratios[shape-active/offset-active/quarantine]="
        << ratio(active_population.active_atom_index_list_by_parameter.at(
            kLogPeakHeightChangeIndex).size()) << "/"
        << ratio(active_population.active_atom_index_list_by_parameter.at(
            kOffsetToPeakRatioChangeIndex).size()) << "/"
        << ratio(progress.quarantine_atom_count)
        << ", certificate-blockers[objective-domain/quarantine-transition/suspicious-offset/rejected-cluster]="
        << certificate.blockers.objective_domain_changed << "/"
        << certificate.blockers.quarantine_transition << "/"
        << certificate.blockers.suspicious_offset_fallback << "/"
        << certificate.blockers.rejected_cluster
        << ", limiters[guard/fixed/quarantine/trust/objective/reject/polish/boundary/rescue]="
        << guard_rejected_trial_count << "/"
        << shape_fixed_count + offset_fixed_count + hard_fixed_count << "/"
        << progress.quarantine_atom_count << "/"
        << trust_skipped_trial_count << "/"
        << objective_rejected_trial_count + boundary_backtracked_count << "/"
        << progress.rejected_cluster_count << "/"
        << assembled_uses_polish << "/"
        << accepted_boundary_count << "/"
        << rescued_boundary_count
        << ", fixed[shape/offset/hard]="
        << shape_fixed_count << "/" << offset_fixed_count << "/" << hard_fixed_count
        << ", blockers[suspicious/rejected/quarantine-transition/domain-change]="
        << progress.suspicious_atom_count << "/"
        << progress.rejected_cluster_count << "/"
        << has_quarantine_transition << "/"
        << objective_domain_changed << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Debug, message.str());
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
    std::optional<GaussianModel3DWithUncertainty> unrestricted_mdpde;
    try
    {
        auto candidate_result{
            atom_context.refit_design.Estimate(
                adjusted_response_list,
                atom_context.alpha_r,
                options.thread_size,
                offset_model)
        };
        if (candidate_result.fit_result.has_value())
        {
            attempted_refit_status = candidate_result.fit_result->status;
        }
        if (IsValidSecondStageGaussianModel(candidate_result.mdpde.GetModel()))
        {
            unrestricted_mdpde = candidate_result.mdpde;
        }
        const auto assessment{
            AssessSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                candidate_result.mdpde.GetModel(),
                options,
                previous_baseline,
                SuspiciousUpdateMode::PostRefit)
        };
        if (unrestricted_mdpde.has_value())
        {
            return LocalAtomRefitResult{
                std::move(candidate_result),
                std::move(unrestricted_mdpde),
                assessment,
                attempted_refit_status
            };
        }
        failed_shape_assessment = assessment;
    }
    catch (const std::exception &)
    {
        failed_shape_assessment = SuspiciousGaussianAssessment{
            SuspiciousGaussianReason::InvalidModel,
            SuspiciousUpdateMode::PostRefit,
            std::numeric_limits<double>::infinity()
        };
    }

    auto result{ previous_result };
    SetLocalResultOffset(result, offset_model.GetOffset());
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
        std::move(unrestricted_mdpde),
        failed_shape_assessment.value_or(fallback_assessment),
        attempted_refit_status
    };
}

static std::vector<std::optional<GaussianModel3DWithUncertainty>>
RunUnrestrictedShapeRefits(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & operator_offset_state,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const FitOptions & options)
{
    const auto operator_refit_snapshot{
        BuildGroupMedianModelList(
            group_id_by_atom_index,
            BuildFittedGaussianSnapshot(operator_offset_state))
    };
    const auto operator_model_bundle{
        BuildSecondStageModelSnapshot(context, operator_refit_snapshot)
    };
    const auto adjusted_response_cache{
        BuildSecondStageAdjustedResponseCache(context, operator_model_bundle)
    };
    std::vector<std::size_t> atom_index_list;
    for (const auto & key : cluster_key_list)
    {
        atom_index_list.insert(atom_index_list.end(), key.begin(), key.end());
    }
    std::vector<std::optional<GaussianModel3DWithUncertainty>> result(
        context.size());
    int refit_thread_size{ options.thread_size };
#ifdef USE_OPENMP
    const bool parallel_refits{
        Logger::GetLogLevel() < LogLevel::Debug &&
        options.thread_size > 1 && atom_index_list.size() > 1
    };
    if (parallel_refits) refit_thread_size = 1;
#pragma omp parallel for schedule(dynamic) if(parallel_refits) num_threads(options.thread_size)
#endif
    for (std::size_t position = 0; position < atom_index_list.size(); position++)
    {
        const auto atom_index{ atom_index_list.at(position) };
        try
        {
            const auto candidate{
                context.at(atom_index).refit_design.Estimate(
                    adjusted_response_cache.at(atom_index),
                    context.at(atom_index).alpha_r,
                    refit_thread_size,
                    GetFitModel(operator_model_bundle.selected, atom_index))
            };
            if (IsValidSecondStageGaussianModel(candidate.mdpde.GetModel()))
            {
                result.at(atom_index) = candidate.mdpde;
            }
        }
        catch (const std::exception &)
        {
        }
    }
    return result;
}

static IterationProposalResult RunProposalIteration(
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

    FixedPointOperatorEvidence fixed_point_operator;
    fixed_point_operator.state = previous_state;
    fixed_point_operator.shape_status_by_atom.assign(
        context.size(), OperatorEndpointStatus::NotRequested);
    fixed_point_operator.offset_status_by_atom.assign(
        context.size(), OperatorEndpointStatus::NotRequested);
    for (std::size_t cluster_position = 0;
        cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        const auto & offset_result{
            joint_offset_result_list.at(cluster_position)
        };
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            if (IsJointOffsetSolveHardFailure(offset_result.status))
            {
                fixed_point_operator.offset_status_by_atom.at(atom_index) =
                    OperatorEndpointStatus::OffsetSolverFailure;
                continue;
            }
            const auto proposed_offset{
                offset_result.offset(static_cast<Eigen::Index>(position))
            };
            const auto operator_model{
                previous_state.at(atom_index).mdpde.GetModel().WithOffset(
                    proposed_offset)
            };
            if (!IsValidSecondStageGaussianModel(operator_model))
            {
                fixed_point_operator.offset_status_by_atom.at(atom_index) =
                    OperatorEndpointStatus::InvalidOffset;
                continue;
            }
            auto & operator_result{ fixed_point_operator.state.at(atom_index) };
            SetLocalResultOffset(operator_result, proposed_offset);
            fixed_point_operator.offset_status_by_atom.at(atom_index) =
                OperatorEndpointStatus::Available;
        }
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
            health.all_local_refits_solver_qualified = false;
            health.production_convergence_qualified =
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
                health.all_local_refits_solver_qualified = false;
                for (const auto position : position_list)
                {
                    const auto atom_index{ key.at(position) };
                    block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                    current_model_snapshot.selected.at(atom_index) =
                        previous_state.at(atom_index).mdpde.GetModel();
                }
                continue;
            }
            bool accepted{ true };
            std::vector<GaussianModel3D> accepted_model_list;
            std::vector<SuspiciousGaussianAssessment> accepted_assessment_list;
            accepted_model_list.reserve(position_list.size());
            accepted_assessment_list.reserve(position_list.size());
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
                    previous_model.WithOffset(proposed_offset)
                };
                if (!IsValidSecondStageGaussianModel(candidate_model))
                {
                    accepted = false;
                    break;
                }
                accepted_model_list.emplace_back(candidate_model);
                accepted_assessment_list.emplace_back(
                    AssessSuspiciousGaussianUpdate(
                        context.at(atom_index).raw_sampling_entries,
                        candidate_model,
                        options,
                        BuildPreviousSuspiciousProfileBaseline(
                            context.at(atom_index).raw_sampling_entries,
                            previous_model,
                            options),
                        SuspiciousUpdateMode::OffsetOnly));
            }
            for (std::size_t member_position = 0;
                member_position < position_list.size(); member_position++)
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
                }
            }
            if (!accepted)
            {
                health.all_local_refits_solver_qualified = false;
                health.production_convergence_qualified = false;
            }
        }
    }

    bool operator_offsets_complete{ true };
    bool operator_offsets_match_proposal{ true };
    for (const auto & key : cluster_key_list)
    {
        for (const auto atom_index : key)
        {
            if (fixed_point_operator.offset_status_by_atom.at(atom_index) !=
                OperatorEndpointStatus::Available)
            {
                operator_offsets_complete = false;
                operator_offsets_match_proposal = false;
                continue;
            }
            operator_offsets_match_proposal =
                operator_offsets_match_proposal &&
                fixed_point_operator.state.at(atom_index).mdpde.GetModel()
                        .GetOffset() ==
                    current_model_snapshot.selected.at(atom_index).GetOffset();
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
                if (operator_offsets_complete && operator_offsets_match_proposal)
                {
                    fixed_point_operator.shape_status_by_atom.at(atom_index) =
                        OperatorEndpointStatus::ShapeRefitFailure;
                }
                health.all_local_refits_solver_qualified = false;
                if (quarantine_activity.HasActiveShape(atom_index) ||
                    quarantine_activity.HasActiveOffset(atom_index))
                {
                    health.production_convergence_qualified = false;
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
                    std::numeric_limits<double>::infinity()
                };
                continue;
            }
            if (operator_offsets_complete && operator_offsets_match_proposal)
            {
                if (refit_result->unrestricted_mdpde.has_value())
                {
                    fixed_point_operator.state.at(atom_index).mdpde =
                        *refit_result->unrestricted_mdpde;
                    fixed_point_operator.shape_status_by_atom.at(atom_index) =
                        OperatorEndpointStatus::Available;
                }
                else
                {
                    fixed_point_operator.shape_status_by_atom.at(atom_index) =
                        OperatorEndpointStatus::ShapeRefitFailure;
                }
            }
            local_refit_status_by_atom.at(atom_index) =
                refit_result->attempted_refit_status;
            const auto shape_solver_qualified{
                refit_result->attempted_refit_status.has_value() &&
                IsLocalRefitStatusSolverQualified(
                    *refit_result->attempted_refit_status)
            };
            if (!shape_solver_qualified)
            {
                health.all_local_refits_solver_qualified = false;
                if (quarantine_activity.HasActiveShape(atom_index))
                {
                    health.production_convergence_qualified = false;
                }
            }
            if (!refit_result->unrestricted_mdpde.has_value())
            {
                block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            assessment_by_atom.at(atom_index) = refit_result->assessment;
            iteration_state.at(atom_index) = std::move(refit_result->result);
        }
    }
    if (!operator_offsets_complete)
    {
        const auto unavailable_reason{
            std::ranges::any_of(
                fixed_point_operator.offset_status_by_atom,
                [](OperatorEndpointStatus status)
                {
                    return status == OperatorEndpointStatus::OffsetSolverFailure;
                }) ?
                OperatorEndpointStatus::OffsetSolverFailure :
                OperatorEndpointStatus::InvalidOffset
        };
        for (const auto & key : cluster_key_list)
        {
            for (const auto atom_index : key)
            {
                fixed_point_operator.shape_status_by_atom.at(atom_index) =
                    unavailable_reason;
            }
        }
    }
    else if (!operator_offsets_match_proposal)
    {
        fixed_point_operator.shadow_shape_refit_performed = true;
        const auto unrestricted_shape_list{
            RunUnrestrictedShapeRefits(
                context,
                cluster_key_list,
                fixed_point_operator.state,
                group_id_by_atom_index,
                options)
        };
        for (const auto & key : cluster_key_list)
        {
            for (const auto atom_index : key)
            {
                if (unrestricted_shape_list.at(atom_index).has_value())
                {
                    fixed_point_operator.state.at(atom_index).mdpde =
                        *unrestricted_shape_list.at(atom_index);
                    fixed_point_operator.shape_status_by_atom.at(atom_index) =
                        OperatorEndpointStatus::Available;
                }
                else
                {
                    fixed_point_operator.shape_status_by_atom.at(atom_index) =
                        OperatorEndpointStatus::ShapeRefitFailure;
                }
            }
        }
    }
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        if (block_activity.hard_failure_atom_mask.at(atom_index) != 0)
        {
            iteration_state.at(atom_index) = previous_state.at(atom_index);
            continue;
        }
        if (block_activity.offset_fixed_atom_mask.at(atom_index) == 0)
        {
            continue;
        }
        const auto previous_offset{
            previous_state.at(atom_index).mdpde.GetModel().GetOffset()
        };
        SetLocalResultOffset(iteration_state.at(atom_index), previous_offset);
    }

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        if (!quarantine_activity.HasActiveShape(atom_index))
        {
            const auto accepted_offset{
                iteration_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            iteration_state.at(atom_index).ols = WithPreservedUncertaintyOffset(
                previous_state.at(atom_index).ols,
                accepted_offset);
            iteration_state.at(atom_index).mdpde = WithPreservedUncertaintyOffset(
                previous_state.at(atom_index).mdpde,
                accepted_offset);
            block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
        }
        if (!quarantine_activity.HasActiveOffset(atom_index))
        {
            const auto previous_offset{
                previous_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            SetLocalResultOffset(iteration_state.at(atom_index), previous_offset);
            block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
        }
        if (quarantine_activity.hard_failure_atom_mask.at(atom_index) != 0)
        {
            iteration_state.at(atom_index) = previous_state.at(atom_index);
            block_activity.hard_failure_atom_mask.at(atom_index) = 1;
        }
    }
    return IterationProposalResult{
        std::move(iteration_state),
        std::move(fixed_point_operator),
        std::move(block_activity),
        std::move(assessment_by_atom),
        std::move(local_refit_status_by_atom),
        std::move(health_by_key)
    };
}

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
    return lhs.boundary_sample_dependency_list.size() ==
            rhs.boundary_sample_dependency_list.size() &&
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
    iteration_internal::ResetClusterSolverWorkspace(
        cluster_key_list,
        iteration_state.solver_workspace_by_key);
    iteration_state.boundary_joint_correction_workspace_by_key.clear();
    iteration_state.graph_partition = std::move(partition);
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
            << iteration_state.graph_partition.boundary_sample_dependency_list.size()
            << "/" << rebuilt_partition.boundary_sample_dependency_list.size()
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
    iteration_internal::ResetClusterSolverWorkspace(
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
            iteration_internal::UsesPolish(
                iteration_state.previous_polish_provenance),
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
        iteration_internal::BuildSuspiciousJointOffsetRidgeMultiplierList(
            ridge_atom_mask)
    };
    const auto debug_convergence_audit_enabled{
        !options.quiet_mode && Logger::GetLogLevel() >= LogLevel::Debug };
    const auto iteration_phase_start{
        std::chrono::steady_clock::now()
    };
    auto raw_iteration_result{
        RunProposalIteration(
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

    if (debug_convergence_audit_enabled)
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
            if (!assessment.IsSuspicious() && !has_fixed_block)
            {
                continue;
            }
            std::ostringstream message;
            message << std::scientific << std::setprecision(2)
                << "Unrestricted operator assessment: atom=" << atom_index
                << ", reason=" << GetSuspiciousGaussianReasonText(assessment.reason)
                << ", margin=" << assessment.normalized_margin
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

    const auto & operator_proposal_state{
        raw_iteration_result.operator_proposal_state
    };
    const auto operator_proposal_change_summary{
        SummarizeTransformedChanges(
            operator_proposal_state, previous_state, active_index_list)
    };
    const CandidateSelectionInputs candidate_inputs{
        .context = context,
        .options = options,
        .residual_baseline = residual_baseline,
        .partition = graph_partition,
        .health_by_key = raw_iteration_result.health_by_key,
        .previous_state = previous_state,
        .previous_polish_provenance = iteration_state.previous_polish_provenance,
        .operator_proposal_state = operator_proposal_state,
        .block_activity = raw_iteration_result.block_activity,
        .ridge_multiplier_list = joint_offset_ridge_multiplier_list,
        .objective_domain = objective_domain,
        .previous_objective_by_key = previous_objective_by_key,
        .cluster_objective_state = iteration_state.cluster_objective_state,
        .best_audit_state = iteration_state.best_audit_state,
        .trust_region_state = iteration_state.trust_region_state,
        .solver_workspace_by_key = iteration_state.solver_workspace_by_key,
        .boundary_joint_correction_workspace_by_key = iteration_state.boundary_joint_correction_workspace_by_key,
        .performance_counters = performance_counters
    };
    auto selection{ SelectClusterCandidates(candidate_inputs) };
    const auto iteration_failure_atom_mask{
        BuildSuspiciousFailureAtomMask(
            raw_iteration_result.block_activity,
            raw_iteration_result.assessment_by_atom)
    };
    const auto iteration_suspicious_atom_count{
        static_cast<std::size_t>(std::ranges::count_if(
            iteration_failure_atom_mask,
            [](char value) { return value != 0; }))
    };
    const auto has_suspicious_offset_fallback{ iteration_suspicious_atom_count > 0 };

    auto assembled_state{ std::move(selection.assembled_state) };
    auto assembled_polish_provenance{ std::move(selection.assembled_polish_provenance) };
    auto trust_region_update{
        iteration_state.trust_region_state.ApplyRadiusUpdates(
            selection.grow_trust_region_key_list,
            selection.shrink_trust_region_key_list,
            selection.rejected_key_list,
            selection.exhausted_key_list)
    };
    const auto quarantine_transition{
        iteration_state.quarantine_state.UpdateAfterIteration(
            context,
            selection.accepted_cluster_diagnostic_list,
            selection.rejected_cluster_diagnostic_list,
            raw_iteration_result.block_activity,
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
        iteration_internal::UsesPolish(assembled_polish_provenance)
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
    result.diagnostics.trust_region_update = std::move(trust_region_update);
    iteration_state.rollback_atom_mask =
        raw_iteration_result.block_activity.BuildCombinedFixedAtomMask();
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
        GetMaximumTransformedChange(
            operator_proposal_change_summary.maximum_list)
    };

    if (selection.accepted_key_list.empty())
    {
        result.all_rejected_resolution = attempt_number >= kMaximumIterations ?
            AllRejectedResolution::MaximumIterations :
            AllRejectedResolution::BacktrackingExhausted;
        return result;
    }

    std::vector<std::size_t> group_id_by_atom_index;
    group_id_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_id_by_atom_index.emplace_back(atom_context.group_id);
    }
    const auto active_population{
        BuildActiveCoordinatePopulation(
            active_index_list,
            cluster_key_list,
            group_id_by_atom_index,
            raw_iteration_result.block_activity,
            quarantine_activity)
    };
    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            assembled_state,
            previous_state,
            iteration_state.active_index_list)
    };
    const auto accepted_active_dof_change_summary{
        SummarizeActiveDofChanges(
            assembled_state,
            previous_state,
            active_population)
    };
    auto certificate{
        SummarizeFixedPointOperator(
            raw_iteration_result.fixed_point_operator,
            previous_state,
            iteration_state.active_index_list,
            cluster_key_list,
            group_id_by_atom_index)
    };
    certificate.accepted_active_population = active_population;
    certificate.accepted_active_movement =
        accepted_active_dof_change_summary;
    certificate.solver_qualification = EvaluateSolverQualificationAudit(
        iteration_state.active_index_list,
        cluster_key_list,
        group_id_by_atom_index,
        raw_iteration_result.block_activity,
        quarantine_activity,
        raw_iteration_result.local_refit_status_by_atom,
        raw_iteration_result.health_by_key);
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
        std::ranges::any_of(
            result.diagnostics.rejected_cluster_diagnostic_list,
            [&](const auto & diagnostic)
            {
                return std::ranges::find(
                    result.diagnostics.trust_region_update.changed_key_list,
                    diagnostic.key) !=
                    result.diagnostics.trust_region_update.changed_key_list.end();
            })
    };
    const auto has_pending_quarantine_lifecycle{
        HasPendingQuarantineLifecycle(
            iteration_state.quarantine_state.state_by_target)
    };
    if (objective_domain_changed || improved_best_audit ||
        changed_rejected_trust_radius || has_pending_quarantine_lifecycle ||
        has_quarantine_transition)
    {
        iteration_state.audit_patience_count = 0;
    }
    else
    {
        iteration_state.audit_patience_count++;
    }

    result.progress.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.progress.accepted_maximum_transformed_change =
        GetMaximumTransformedChange(
            transformed_change_summary.maximum_list);
    result.transformed_change_percentile =
        accepted_active_dof_change_summary.percentile_list;
    result.audit_patience_exhausted = iteration_state.audit_patience_count >= kAuditPatience;
    certificate.blockers = ConvergenceOrthogonalBlockers{
        objective_domain_changed,
        has_quarantine_transition,
        has_suspicious_offset_fallback,
        !selection.rejected_key_list.empty()
    };
    result.converged = certificate.ProductionConverged();

    if (!options.quiet_mode && Logger::GetLogLevel() >= LogLevel::Debug)
    {
        const auto accepted_raw_change_summary{
            SummarizeTransformedChanges(
                assembled_state,
                operator_proposal_state,
                iteration_state.active_index_list)
        };
        LogConvergenceSafeguardAudit(
            options.quiet_mode,
            attempt_number,
            result.progress,
            certificate,
            accepted_active_dof_change_summary,
            active_population,
            result.diagnostics,
            raw_iteration_result.block_activity,
            raw_iteration_result.local_refit_status_by_atom,
            GetMaximumTransformedChange(
                accepted_raw_change_summary.maximum_list) == 0.0,
            assembled_uses_polish,
            has_quarantine_transition,
            objective_domain_changed);
    }

    iteration_state.previous_state = std::move(assembled_state);
    iteration_state.previous_polish_provenance = std::move(assembled_polish_provenance);
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
    const FinalDependencyPolishResult & polish_result,
    std::string_view safety_policy,
    std::string_view safety_status,
    bool applied,
    const ConvergenceCertificate * base_certificate = nullptr,
    const ConvergenceCertificate * candidate_certificate = nullptr)
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
        << ", residual-safety-policy=" << safety_policy
        << ", residual-safety=" << safety_status
        << ", applied=" << (applied ? "yes" : "no");
    const auto append_certificate = [&](
        std::string_view prefix,
        const ConvergenceCertificate & certificate)
    {
        message << ", " << prefix << "-solver-qualified="
            << (certificate.solver_qualification.solver_qualified ? "yes" : "no")
            << ", " << prefix << "-operator-complete="
            << (certificate.OperatorComplete() ? "yes" : "no")
            << ", " << prefix << "-residual-p99=";
        AppendAuditValues(
            message,
            certificate.operator_nominal_residual.percentile_list);
        message << ", " << prefix << "-residual-max=";
        AppendAuditValues(message, certificate.operator_nominal_residual.maximum_list);
    };
    if (base_certificate != nullptr)
    {
        append_certificate("base", *base_certificate);
    }
    if (candidate_certificate != nullptr)
    {
        append_certificate("candidate", *candidate_certificate);
    }
    message
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

enum class FinalPolishCertificationPolicy
{
    RequireResidualNonRegression,
    RequireStrictFixedPoint
};

enum class FinalPolishResidualSafetyStatus
{
    NotEvaluated,
    AbsolutePassed,
    RelativePassed,
    Failed,
    Error
};

struct FinalPolishResidualSafetyResult
{
    FinalPolishResidualSafetyStatus status{
        FinalPolishResidualSafetyStatus::NotEvaluated
    };
    std::optional<ConvergenceCertificate> base{};
    std::optional<ConvergenceCertificate> candidate{};

    bool Passed() const
    {
        return status == FinalPolishResidualSafetyStatus::AbsolutePassed ||
            status == FinalPolishResidualSafetyStatus::RelativePassed;
    }
};

static std::string_view GetFinalPolishCertificationPolicyText(
    FinalPolishCertificationPolicy policy)
{
    switch (policy)
    {
    case FinalPolishCertificationPolicy::RequireResidualNonRegression:
        return "non-regression";
    case FinalPolishCertificationPolicy::RequireStrictFixedPoint:
        return "strict-fixed-point";
    }
    throw std::invalid_argument("Unknown final polish certification policy.");
}

static std::string_view GetFinalPolishResidualSafetyStatusText(
    FinalPolishResidualSafetyStatus status)
{
    switch (status)
    {
    case FinalPolishResidualSafetyStatus::NotEvaluated:
        return "not-evaluated";
    case FinalPolishResidualSafetyStatus::AbsolutePassed:
        return "absolute-passed";
    case FinalPolishResidualSafetyStatus::RelativePassed:
        return "relative-passed";
    case FinalPolishResidualSafetyStatus::Failed:
        return "failed";
    case FinalPolishResidualSafetyStatus::Error:
        return "error";
    }
    throw std::invalid_argument("Unknown final polish residual safety status.");
}

static std::optional<ConvergenceCertificate> EvaluateFinalPolishCertificate(
    const SecondStageContext & context,
    const FitOptions & options,
    IterationState & iteration_state,
    const SuspiciousBlockActivity & final_block_activity,
    const FitState & candidate_state)
{
    try
    {
        const auto cluster_key_list{
            BuildGraphClusterKeyList(iteration_state.graph_partition)
        };
        auto ridge_atom_mask{ iteration_state.rollback_atom_mask };
        const auto quarantine_atom_mask{
            final_block_activity.BuildCombinedFixedAtomMask()
        };
        for (std::size_t atom_index = 0;
            atom_index < ridge_atom_mask.size(); atom_index++)
        {
            if (quarantine_atom_mask.at(atom_index) != 0)
            {
                ridge_atom_mask.at(atom_index) = 1;
            }
        }
        const auto joint_offset_ridge_multiplier_list{
            iteration_internal::BuildSuspiciousJointOffsetRidgeMultiplierList(
                ridge_atom_mask)
        };
        auto certificate_options{ options };
        certificate_options.quiet_mode = true;
        auto proposal_result{
            RunProposalIteration(
                context,
                cluster_key_list,
                candidate_state,
                certificate_options,
                joint_offset_ridge_multiplier_list,
                final_block_activity,
                iteration_state.solver_workspace_by_key)
        };
        std::vector<std::size_t> group_id_by_atom_index;
        group_id_by_atom_index.reserve(context.size());
        for (const auto & atom_context : context)
        {
            group_id_by_atom_index.emplace_back(atom_context.group_id);
        }
        auto certificate{ SummarizeFixedPointOperator(
            proposal_result.fixed_point_operator,
            candidate_state,
            iteration_state.active_index_list,
            cluster_key_list,
            group_id_by_atom_index) };
        certificate.accepted_active_population = BuildActiveCoordinatePopulation(
            iteration_state.active_index_list,
            cluster_key_list,
            group_id_by_atom_index,
            proposal_result.block_activity,
            final_block_activity);
        certificate.solver_qualification = EvaluateSolverQualificationAudit(
            iteration_state.active_index_list,
            cluster_key_list,
            group_id_by_atom_index,
            proposal_result.block_activity,
            final_block_activity,
            proposal_result.local_refit_status_by_atom,
            proposal_result.health_by_key);
        return certificate;
    }
    catch (const std::exception &)
    {
        return std::nullopt;
    }
    catch (...)
    {
        return std::nullopt;
    }
}

static bool HasComparableFinalPolishOperatorEvidence(
    const ConvergenceCertificate & certificate)
{
    const auto & percentile_list{
        certificate.operator_nominal_residual.percentile_list
    };
    return certificate.solver_qualification.solver_qualified &&
        certificate.OperatorComplete() &&
        certificate.InvariantsClear() &&
        std::ranges::all_of(
            percentile_list,
            [](double value) { return std::isfinite(value); });
}

static bool IsFinalPolishResidualNonWorsening(
    const ConvergenceCertificate & base_certificate,
    const ConvergenceCertificate & candidate_certificate)
{
    if (!HasComparableFinalPolishOperatorEvidence(base_certificate) ||
        !HasComparableFinalPolishOperatorEvidence(candidate_certificate))
    {
        return false;
    }
    const auto & base_percentile_list{
        base_certificate.operator_nominal_residual.percentile_list
    };
    const auto & candidate_percentile_list{
        candidate_certificate.operator_nominal_residual.percentile_list
    };
    for (std::size_t index = 0; index < kTransformedChangeSize; index++)
    {
        if (candidate_percentile_list.at(index) >
            std::max(base_percentile_list.at(index), kTransformedChangeTolerance))
        {
            return false;
        }
    }
    return true;
}

static FinalPolishResidualSafetyResult EvaluateFinalPolishResidualSafety(
    const SecondStageContext & context,
    const FitOptions & options,
    IterationState & iteration_state,
    const SuspiciousBlockActivity & final_block_activity,
    const FitState & base_state,
    const FitState & candidate_state,
    FinalPolishCertificationPolicy policy)
{
    FinalPolishResidualSafetyResult result;
    result.candidate = EvaluateFinalPolishCertificate(
        context,
        options,
        iteration_state,
        final_block_activity,
        candidate_state);
    if (!result.candidate.has_value())
    {
        result.status = FinalPolishResidualSafetyStatus::Error;
        return result;
    }
    if (result.candidate->StrictOperatorPassed())
    {
        result.status = FinalPolishResidualSafetyStatus::AbsolutePassed;
        return result;
    }
    if (policy == FinalPolishCertificationPolicy::RequireStrictFixedPoint)
    {
        result.status = FinalPolishResidualSafetyStatus::Failed;
        return result;
    }
    if (!HasComparableFinalPolishOperatorEvidence(*result.candidate))
    {
        result.status = FinalPolishResidualSafetyStatus::Failed;
        return result;
    }
    result.base = EvaluateFinalPolishCertificate(
        context,
        options,
        iteration_state,
        final_block_activity,
        base_state);
    if (!result.base.has_value())
    {
        result.status = FinalPolishResidualSafetyStatus::Error;
        return result;
    }
    result.status = IsFinalPolishResidualNonWorsening(
            *result.base,
            *result.candidate) ?
        FinalPolishResidualSafetyStatus::RelativePassed :
        FinalPolishResidualSafetyStatus::Failed;
    return result;
}

static void FinalizeSecondStageState(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & graph_topology,
    IterationState & iteration_state,
    bool use_best_audit_state,
    FinalPolishCertificationPolicy certification_policy,
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
    FinalPolishResidualSafetyResult residual_safety;
    if (polish_result.accepted && polish_result.objective.has_value())
    {
        residual_safety = EvaluateFinalPolishResidualSafety(
            context,
            options,
            iteration_state,
            final_block_activity,
            base_state,
            polish_result.state,
            certification_policy);
    }
    const auto polish_applied{
        polish_result.accepted &&
        polish_result.objective.has_value() &&
        residual_safety.Passed()
    };
    LogFinalDependencyPolish(
        options.quiet_mode,
        polish_result,
        GetFinalPolishCertificationPolicyText(certification_policy),
        GetFinalPolishResidualSafetyStatusText(residual_safety.status),
        polish_applied,
        residual_safety.base.has_value() ? &*residual_safety.base : nullptr,
        residual_safety.candidate.has_value() ? &*residual_safety.candidate : nullptr);
    if (polish_applied)
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

struct TerminalAuditOutcome
{
    FitState finalized_state{};
    std::optional<ObjectiveBreakdown> objective{};
    std::array<std::size_t, 3> comparison_domain_size{};
};

static TerminalAuditOutcome BuildTerminalAuditOutcome(
    const SecondStageContext & context,
    const FitState & finalized_state,
    const ObjectiveDomain & comparison_objective_domain)
{
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, finalized_state)
    };
    return TerminalAuditOutcome{
        finalized_state,
        EvaluateAuditObjective(
            comparison_objective_domain,
            SnapshotResidualEvaluator{ context, model_snapshot }),
        std::array{
            comparison_objective_domain.active_atom_count,
            comparison_objective_domain.fit_sample_count,
            comparison_objective_domain.tail_sample_count }
    };
}

static void AppendAuditObjective(
    std::ostringstream & message,
    const std::optional<ObjectiveBreakdown> & objective)
{
    if (!objective.has_value())
    {
        message << "-/-/-/-";
        return;
    }
    message << objective->fit_range_residual_objective << "/"
        << objective->GetTailValidationPenalty() << "/"
        << objective->offset_plausibility_penalty << "/"
        << objective->GetTotalObjective();
}

static void LogAuditAtomState(
    std::string_view marker,
    const SecondStageContext & context,
    const FitState & state)
{
    for (std::size_t atom_index = 0; atom_index < state.size(); atom_index++)
    {
        const auto & model{ state.at(atom_index).mdpde.GetModel() };
        std::ostringstream message;
        message << std::scientific << std::setprecision(17)
            << marker << " schema=1"
            << ", serial=" << context.at(atom_index).atom->GetSerialID()
            << ", group=" << context.at(atom_index).group_id
            << ", amplitude=" << model.GetAmplitude()
            << ", width=" << model.GetWidth()
            << ", offset=" << model.GetOffset();
        Logger::Log(LogLevel::Debug, message.str());
    }
}

static void LogSecondStageAuditTerminal(
    bool quiet_mode,
    const SecondStageContext & context,
    std::string_view reason,
    std::size_t attempt_number,
    std::size_t accepted_iteration_count,
    const TerminalAuditOutcome & outcome)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "Second-stage audit terminal: schema=2"
        << ", reason=" << reason
        << ", try=" << attempt_number
        << ", acc=" << accepted_iteration_count
        << ", fixed-domain=" << outcome.comparison_domain_size.at(0)
        << "/" << outcome.comparison_domain_size.at(1)
        << "/" << outcome.comparison_domain_size.at(2)
        << ", objective=";
    AppendAuditObjective(message, outcome.objective);
    Logger::Log(LogLevel::Debug, message.str());
    LogAuditAtomState(
        "Second-stage audit terminal atom:", context, outcome.finalized_state);
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
    const TransformedChange & transformed_change_percentile,
    const IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << iteration_state.accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_percentile.at(kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_percentile.at(kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_percentile.at(kOffsetToPeakRatioChangeIndex);
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
            iteration_internal::UsesPolish(
                iteration_state.previous_polish_provenance)
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
    SecondStageContext context;
    FitState initial_state;
    {
        auto initialization{
            BuildSecondStageInitialization(model_object, options)
        };
        StoreSecondStageNeighborCounts(
            model_object,
            initialization.context,
            initialization.neighbor_count_list);
        if (!options.quiet_mode)
        {
            Logger::Log(LogLevel::Info,
                "Run 2nd-stage local atom fitting with iterations...");
        }

        if (initialization.failure !=
            SecondStageInitializationResult::Failure::None)
        {
            if (!options.quiet_mode)
            {
                const auto unselected_seed_failure{
                    initialization.failure ==
                        SecondStageInitializationResult::Failure::
                            UnselectedSeedUnavailable
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
            initialization.selection_record_list,
            options.quiet_mode);
        LogUnselectedSecondStageSeedSelections(
            initialization.unselected_selection_record_list,
            options.quiet_mode);
        context = std::move(initialization.context);
        initial_state = std::move(initialization.state);
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

    const auto audit_comparison_objective_domain{ iteration_state.objective_domain };
    std::size_t last_attempt_number{ 0 };
    const auto log_audit_terminal = [&](
        std::string_view reason,
        bool use_best_audit_state)
    {
        const auto & finalized_state{
            use_best_audit_state && iteration_state.best_audit_state.has_value() ?
                iteration_state.best_audit_state->state :
                iteration_state.previous_state
        };
        LogSecondStageAuditTerminal(
            options.quiet_mode,
            context,
            reason,
            last_attempt_number,
            iteration_state.accepted_iteration_count,
            BuildTerminalAuditOutcome(
                context,
                finalized_state,
                audit_comparison_objective_domain));
    };

    std::string_view final_stop_reason;
    bool maximum_iterations_reached{ false };
    for (std::size_t iter = 0; iter < kMaximumIterations; iter++)
    {
        last_attempt_number = iter + 1;
        if (iteration_state.active_index_list.empty())
        {
            FinalizeSecondStageState(
                model_object,
                context,
                options,
                graph_topology,
                iteration_state,
                false,
                FinalPolishCertificationPolicy::RequireResidualNonRegression,
                performance_counters);
            log_audit_terminal("quarantine", false);
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
        LogAcceptedCandidateSearchDiagnostics(
            options.quiet_mode,
            iteration_result.diagnostics);
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
        LogTrustModelShadowDiagnostics(
            options.quiet_mode,
            iteration_result.diagnostics,
            iteration_result.progress);
#endif
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

        if (iteration_result.converged)
        {
            FinalizeSecondStageState(
                model_object,
                context,
                options,
                graph_topology,
                iteration_state,
                false,
                FinalPolishCertificationPolicy::RequireStrictFixedPoint,
                performance_counters);
            log_audit_terminal("converged", false);
            if (iteration_state.quarantine_state.HasFailures())
            {
                LogQuarantineFallback(options.quiet_mode, iteration_state);
            }
            else
            {
                LogConverged(
                    options.quiet_mode,
                    iteration_result.transformed_change_percentile,
                    iteration_state);
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state,
                "converged",
                false);
            return true;
        }

        if (iteration_result.audit_patience_exhausted)
        {
            final_stop_reason = "audit-patience";
            break;
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
            FinalPolishCertificationPolicy::RequireResidualNonRegression,
            performance_counters);
        log_audit_terminal(final_stop_reason, audit_state != nullptr);
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
        GetMaximumTransformedChange(drift_summary.maximum_list)
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

bool ConvergenceOrthogonalBlockers::Clear() const
{
    return !objective_domain_changed && !quarantine_transition &&
        !suspicious_offset_fallback && !rejected_cluster;
}

bool ConvergenceCertificate::AcceptedPercentilePassed() const
{
    return IsTransformedPercentileConverged(accepted_active_movement);
}

bool ConvergenceCertificate::OperatorPercentilePassed() const
{
    return IsTransformedPercentileConverged(operator_nominal_residual);
}

bool ConvergenceCertificate::OperatorComplete() const
{
    return std::ranges::all_of(
        operator_unavailable_count,
        [](std::size_t count) { return count == 0; });
}

bool ConvergenceCertificate::InvariantsClear() const
{
    return !MixedSharedGroup();
}

bool ConvergenceCertificate::MixedSharedGroup() const
{
    return accepted_active_population.mixed_offset_group_count != 0 ||
        solver_qualification.mixed_offset_group_count != 0;
}

bool ConvergenceCertificate::StrictOperatorPassed() const
{
    return solver_qualification.solver_qualified && OperatorComplete() &&
        InvariantsClear() && OperatorPercentilePassed();
}

bool ConvergenceCertificate::ProductionConverged() const
{
    return StrictOperatorPassed() && AcceptedPercentilePassed() &&
        blockers.Clear();
}

FixedPointResidualInterpretation EvaluateFixedPointResidualInterpretation(
    bool operator_complete,
    bool qualification_passed,
    const TransformedChangeSummary & accepted_change,
    const TransformedChangeSummary & fixed_point_residual)
{
    if (!operator_complete)
    {
        return FixedPointResidualInterpretation::Restricted;
    }
    const auto accepted_small{
        IsTransformedPercentileConverged(accepted_change)
    };
    const auto residual_small{
        IsTransformedPercentileConverged(fixed_point_residual)
    };
    const auto tail_clean = [](const TransformedChangeSummary & summary)
    {
        return std::ranges::all_of(
            summary.maximum_list,
            [](double value)
            {
                return std::isfinite(value) &&
                    value < kLegacyMaximumTransformedChangeTolerance;
            });
    };
    if (accepted_small && residual_small && !qualification_passed)
    {
        return FixedPointResidualInterpretation::UnqualifiedSmall;
    }
    if (accepted_small && !residual_small)
    {
        return FixedPointResidualInterpretation::StepLimited;
    }
    if (!accepted_small && residual_small)
    {
        return FixedPointResidualInterpretation::PostprocessedMovement;
    }
    if (accepted_small && residual_small &&
        (!tail_clean(accepted_change) || !tail_clean(fixed_point_residual)))
    {
        return FixedPointResidualInterpretation::BulkFixedPointWithTail;
    }
    if (accepted_small && residual_small)
    {
        return FixedPointResidualInterpretation::FixedPointConverged;
    }
    return FixedPointResidualInterpretation::Progressing;
}

std::string_view GetFixedPointResidualInterpretationText(
    FixedPointResidualInterpretation interpretation)
{
    switch (interpretation)
    {
    case FixedPointResidualInterpretation::Restricted:
        return "restricted";
    case FixedPointResidualInterpretation::UnqualifiedSmall:
        return "unqualified-small";
    case FixedPointResidualInterpretation::StepLimited:
        return "step-limited";
    case FixedPointResidualInterpretation::PostprocessedMovement:
        return "postprocessed-movement";
    case FixedPointResidualInterpretation::BulkFixedPointWithTail:
        return "bulk-fixed-point-with-tail";
    case FixedPointResidualInterpretation::FixedPointConverged:
        return "fixed-point-converged";
    case FixedPointResidualInterpretation::Progressing:
        return "progressing";
    }
    throw std::invalid_argument("Unknown fixed-point residual interpretation.");
}

SuspiciousUpdateMask BuildSuspiciousFailureAtomMask(
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom)
{
    const auto atom_count{ assessment_by_atom.size() };
    if (block_activity.shape_fixed_atom_mask.size() != atom_count ||
        block_activity.offset_fixed_atom_mask.size() != atom_count ||
        block_activity.hard_failure_atom_mask.size() != atom_count)
    {
        throw std::invalid_argument(
            "Suspicious failure activity and assessment sizes are inconsistent.");
    }
    SuspiciousUpdateMask result(atom_count, 0);
    for (std::size_t atom_index = 0; atom_index < atom_count; atom_index++)
    {
        const auto has_fixed_endpoint{
            block_activity.shape_fixed_atom_mask.at(atom_index) != 0 ||
            block_activity.offset_fixed_atom_mask.at(atom_index) != 0
        };
        result.at(atom_index) =
            block_activity.hard_failure_atom_mask.at(atom_index) != 0 ||
            (has_fixed_endpoint &&
                assessment_by_atom[atom_index].reason !=
                    SuspiciousGaussianReason::None) ? 1 : 0;
    }
    return result;
}

ActiveCoordinatePopulation BuildActiveCoordinatePopulation(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity)
{
    ValidateActiveCoordinateInputs(
        group_id_by_atom_index.size(),
        group_id_by_atom_index,
        block_activity,
        quarantine_activity);

    ActiveCoordinatePopulation result;
    for (const auto atom_index : atom_index_list)
    {
        if (block_activity.HasActiveShape(atom_index))
        {
            result.active_atom_index_list_by_parameter.at(kLogPeakHeightChangeIndex)
                .emplace_back(atom_index);
            result.active_atom_index_list_by_parameter.at(kLogWidthChangeIndex)
                .emplace_back(atom_index);
        }
        if (block_activity.HasActiveOffset(atom_index))
        {
            result.active_atom_index_list_by_parameter.at(kOffsetToPeakRatioChangeIndex)
                .emplace_back(atom_index);
        }
    }

    const auto offset_group_list{
        BuildOffsetGroupEntries(cluster_key_list, group_id_by_atom_index)
    };
    result.total_offset_group_count = offset_group_list.size();
    for (const auto & group : offset_group_list)
    {
        const auto active_count{ static_cast<std::size_t>(std::ranges::count_if(
            group.atom_index_list,
            [&](const auto atom_index)
            {
                return block_activity.HasActiveOffset(atom_index);
            })) };
        const auto quarantine_count{ static_cast<std::size_t>(std::ranges::count_if(
            group.atom_index_list,
            [&](const auto atom_index)
            {
                return !quarantine_activity.HasActiveOffset(atom_index);
            })) };
        if (active_count == 0)
        {
            if (quarantine_count != 0) result.quarantined_offset_group_count++;
            else result.fixed_offset_group_count++;
            continue;
        }

        const auto is_mixed{ active_count != group.atom_index_list.size() };
        if (is_mixed) result.mixed_offset_group_count++;
        result.active_offset_group_atom_index_list.emplace_back(group.atom_index_list);
        result.mixed_offset_group_mask.emplace_back(is_mixed ? 1 : 0);
    }
    return result;
}

TransformedChangeSummary SummarizeActiveDofChanges(
    const std::vector<TransformedChange> & change_list,
    const ActiveCoordinatePopulation & population)
{
    auto result{ SummarizeTransformedChangesByParameter(
        change_list,
        population.active_atom_index_list_by_parameter) };

    if (population.active_offset_group_atom_index_list.size() !=
        population.mixed_offset_group_mask.size())
    {
        throw std::invalid_argument(
            "Active-coordinate shared-offset population is inconsistent.");
    }

    std::vector<double> offset_group_change_list;
    offset_group_change_list.reserve(
        population.active_offset_group_atom_index_list.size());
    for (std::size_t group_position = 0;
        group_position < population.active_offset_group_atom_index_list.size();
        group_position++)
    {
        double maximum_change{ 0.0 };
        bool is_finite{ population.mixed_offset_group_mask.at(group_position) == 0 };
        for (const auto atom_index :
            population.active_offset_group_atom_index_list.at(group_position))
        {
            if (atom_index >= change_list.size())
            {
                throw std::invalid_argument(
                    "Active-coordinate shared-offset change input is inconsistent.");
            }
            const auto value{
                change_list.at(atom_index).at(
                    kOffsetToPeakRatioChangeIndex)
            };
            if (!std::isfinite(value))
            {
                is_finite = false;
                break;
            }
            maximum_change = std::max(maximum_change, std::abs(value));
        }
        offset_group_change_list.emplace_back(
            is_finite ? maximum_change : std::numeric_limits<double>::infinity());
    }

    result.population_size_list.at(kOffsetToPeakRatioChangeIndex) =
        offset_group_change_list.size();
    result.percentile_list.at(
        kOffsetToPeakRatioChangeIndex) = array_helper::ComputePercentile(
            offset_group_change_list,
            kConvergencePercentile);
    result.maximum_list.at(kOffsetToPeakRatioChangeIndex) =
        offset_group_change_list.empty() ? 0.0 :
            *std::ranges::max_element(offset_group_change_list);
    return result;
}

TransformedChangeSummary SummarizeActiveDofChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const ActiveCoordinatePopulation & population)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument(
            "Active-coordinate transformed state sizes are inconsistent.");
    }
    std::vector<TransformedChange> change_list;
    change_list.reserve(current_state.size());
    for (std::size_t atom_index = 0; atom_index < current_state.size(); atom_index++)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, atom_index),
            GetFitModel(previous_state, atom_index)));
    }
    return SummarizeActiveDofChanges(change_list, population);
}

SolverQualificationAudit EvaluateSolverQualificationAudit(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key)
{
    const auto atom_count{ group_id_by_atom_index.size() };
    ValidateActiveCoordinateInputs(
        atom_count,
        group_id_by_atom_index,
        block_activity,
        quarantine_activity);
    if (local_refit_status_by_atom.size() != atom_count)
    {
        throw std::invalid_argument(
            "Convergence audit qualification inputs are inconsistent.");
    }

    SolverQualificationAudit result;
    result.production_qualified = std::ranges::all_of(
        health_by_key | std::views::values,
        &ClusterHealth::production_convergence_qualified);
    for (const auto & health : health_by_key | std::views::values)
    {
        result.joint_offset_status_count.at(
            static_cast<std::size_t>(health.joint_offset_status))++;
    }
    for (const auto atom_index : atom_index_list)
    {
        if (!block_activity.HasActiveShape(atom_index))
        {
            if (!quarantine_activity.HasActiveShape(atom_index))
            {
                result.quarantined_shape_count++;
            }
            else
            {
                result.fixed_shape_count++;
            }
            continue;
        }

        result.active_shape_count++;
        const auto & status{ local_refit_status_by_atom[atom_index] };
        if (status.has_value() && IsLocalRefitStatusSolverQualified(*status))
        {
            result.qualified_shape_count++;
        }
        else if (status.has_value())
        {
            result.soft_unqualified_shape_count++;
            result.solver_qualified = false;
        }
        else
        {
            result.hard_failure_shape_count++;
            result.solver_qualified = false;
        }
    }

    for (const auto & group :
        BuildOffsetGroupEntries(cluster_key_list, group_id_by_atom_index))
    {
        const auto active_count{ static_cast<std::size_t>(std::ranges::count_if(
            group.atom_index_list,
            [&](const auto atom_index)
            {
                return block_activity.HasActiveOffset(atom_index);
            })) };
        const auto quarantine_count{ static_cast<std::size_t>(std::ranges::count_if(
            group.atom_index_list,
            [&](const auto atom_index)
            {
                return !quarantine_activity.HasActiveOffset(atom_index);
            })) };
        if (active_count == 0)
        {
            if (quarantine_count != 0) result.quarantined_offset_group_count++;
            else result.fixed_offset_group_count++;
            continue;
        }
        if (active_count != group.atom_index_list.size())
        {
            result.mixed_offset_group_count++;
            result.solver_qualified = false;
            continue;
        }

        result.active_offset_group_count++;
        const auto health_iter{ health_by_key.find(group.cluster_key) };
        if (health_iter == health_by_key.end())
        {
            result.hard_failure_offset_group_count++;
            result.solver_qualified = false;
            continue;
        }
        if (health_iter->second.joint_offset_status ==
            JointOffsetSolveStatus::Converged)
        {
            result.qualified_offset_group_count++;
        }
        else if (IsJointOffsetSolveHardFailure(
                health_iter->second.joint_offset_status))
        {
            result.hard_failure_offset_group_count++;
            result.solver_qualified = false;
        }
        else
        {
            result.soft_unqualified_offset_group_count++;
            result.solver_qualified = false;
        }
    }

    result.restricted_active_set = result.fixed_shape_count != 0 ||
        result.quarantined_shape_count != 0 ||
        result.fixed_offset_group_count != 0 ||
        result.quarantined_offset_group_count != 0 ||
        result.mixed_offset_group_count != 0;
    result.all_fixed = result.active_shape_count == 0 &&
        result.active_offset_group_count == 0 &&
        result.mixed_offset_group_count == 0;
    return result;
}

bool HasPendingQuarantineLifecycle(
    const QuarantineFailureStateMap & state_by_target)
{
    return std::ranges::any_of(
        state_by_target | std::views::values,
        [](const auto & state)
        {
            return !state.quarantined || !state.probation_exhausted;
        });
}

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
