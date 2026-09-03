#include "core/detail/IterationProcess.hpp"

#include "core/detail/Diagnosis.hpp"
#include "core/detail/PreparedLocalGaussianFit.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "data/detail/AtomClassifier.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core::detail {

namespace {

constexpr std::size_t kAuditPatience{ 3 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kConvergencePercentile{ 0.99 };

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

struct SecondStageInitializationResult
{
    SecondStageContext context{};
    FitState state{};
    std::vector<int> neighbor_count_list{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
    std::vector<UnselectedSecondStageSeedSelectionRecord>
        unselected_selection_record_list{};
    SecondStageInitializationFailure failure{
        SecondStageInitializationFailure::None
    };
};

struct QuarantineState
{
    QuarantineFailureStateMap state_by_target{};
    std::vector<QuarantineTarget> probation_target_list{};
    std::size_t entered_target_count{ 0 };
    std::size_t released_target_count{ 0 };
    std::size_t failed_probation_count{ 0 };
    bool force_probation{ false };
    QuarantineState() = default;

    explicit QuarantineState(std::size_t atom_count)
        : m_atom_count(atom_count)
    {
    }

    SuspiciousBlockActivity BeginIteration(std::size_t accepted_iteration_count);
    SuspiciousBlockActivity BuildFinalActivity() const;
    bool UpdateAfterIteration(
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
    std::size_t AtomCount() const;
    std::size_t TargetCount() const
    {
        return static_cast<std::size_t>(std::ranges::count_if(
            state_by_target,
            [](const auto & entry)
            {
                return entry.second.lifecycle != QuarantineLifecycle::Tracking;
            }));
    }
private:
    std::size_t m_atom_count{ 0 };
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

struct FixedPointOperatorEvidence
{
    FitState state{};
    std::vector<char> shape_available_atom_mask{};
    std::vector<char> offset_available_atom_mask{};
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

struct OffsetGroupEntry
{
    ClusterKey cluster_key{};
    ClusterKey atom_index_list{};
};

static void ValidateBlockActivitySize(
    std::size_t atom_count,
    const SuspiciousBlockActivity & block_activity)
{
    if (block_activity.shape_fixed_atom_mask.size() != atom_count ||
        block_activity.offset_fixed_atom_mask.size() != atom_count ||
        block_activity.hard_failure_atom_mask.size() != atom_count)
    {
        throw std::invalid_argument("Active-coordinate convergence inputs are inconsistent.");
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
            atom_index_list_by_group[group_id_by_atom_index.at(atom_index)].emplace_back(atom_index);
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
    const auto operator_nominal_population{ BuildActiveCoordinatePopulation(
        atom_index_list,
        cluster_key_list,
        group_id_by_atom_index,
        nominal_activity) };

    std::vector<TransformedChange> change_list;
    change_list.reserve(previous_state.size());
    for (std::size_t atom_index = 0;
        atom_index < previous_state.size(); atom_index++)
    {
        auto change{ CalculateTransformedChange(
            GetFitModel(evidence.state, atom_index),
            GetFitModel(previous_state, atom_index)) };
        if (evidence.shape_available_atom_mask.at(atom_index) == 0)
        {
            change.at(GaussianModel3D::LogPeakHeightCoordinateIndex()) = std::numeric_limits<double>::infinity();
            change.at(GaussianModel3D::LogWidthCoordinateIndex()) = std::numeric_limits<double>::infinity();
        }
        if (evidence.offset_available_atom_mask.at(atom_index) == 0)
        {
            change.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) = std::numeric_limits<double>::infinity();
        }
        change_list.emplace_back(std::move(change));
    }
    result.operator_nominal_residual = SummarizeActiveDofChanges(change_list, operator_nominal_population);
    result.operator_complete = std::ranges::all_of(
        atom_index_list,
        [&](const auto atom_index)
        {
            return evidence.shape_available_atom_mask.at(atom_index) != 0 &&
                evidence.offset_available_atom_mask.at(atom_index) != 0;
        });
    return result;
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

    std::unordered_map<GroupKey, std::size_t> selected_group_id_by_key;
    selected_group_id_by_key.reserve(context.size());
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        atom_index_map.emplace(context.at(i).atom, i);
    }
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto group_key{ data_internal::GetGroupKey(atom) };
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        atom_context.raw_sampling_entries = local_view.GetRawSamplingEntries(false);
        state.at(atom_index) = local_view.GetGaussianResult(FittingStage::Second);
        atom_context.alpha_r = local_view.GetAlphaR(FittingStage::Second);
        atom_context.refit_design = PreparedLocalGaussianDesign{
            atom_context.raw_sampling_entries,
            options.distance_min,
            options.distance_max
        };
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

    std::unordered_map<const AtomObject *, std::size_t> unselected_atom_index_map;
    std::vector<int> unselected_atom_serial_id_list;
    build_result.neighbor_count_list.reserve(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{
            atom->FindNeighborAtoms(kNeighborAtomSearchRange)
        };
        std::unordered_set<const AtomObject *> neighbor_atom_set;

        atom_context.neighbor_atom_sample_offset_list.reserve(atom_context.raw_sampling_entries.size() + 1);
        atom_context.neighbor_atom_sample_offset_list.emplace_back(0);
        atom_context.neighbor_atom_sample_list.reserve(
            atom_context.raw_sampling_entries.size() * neighbor_atom_list.size());
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{
                atom_context.raw_sampling_entries.at(sample_index)
            };
            for (auto * neighbor_atom : neighbor_atom_list)
            {
                if (options.exclude_hydrogen && neighbor_atom->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto distance{
                    array_helper::ComputeNorm(sample.point.position, neighbor_atom->GetPositionRef())
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
        build_result.neighbor_count_list.emplace_back(static_cast<int>(neighbor_atom_set.size()));
    }

    std::vector<GaussianModel3D> global_models;
    global_models.reserve(context.size());

    for (const auto & result : state)
    {
        global_models.emplace_back(result.mdpde.GetModel());
    }
    const auto global_median{ BuildGaussianParameterMedian(global_models) };

    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto & result{ state.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        const auto selection{ SelectSecondStageSeed(result.mdpde, global_median) };
        if (!selection.has_value())
        {
            build_result.failure = SecondStageInitializationFailure::SelectedSeedUnavailable;
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

    if (!context.unselected_atom_list.empty() && !global_median.has_value())
    {
        build_result.failure = SecondStageInitializationFailure::UnselectedSeedUnavailable;
        return build_result;
    }
    for (std::size_t i = 0; i < context.unselected_atom_list.size(); i++)
    {
        auto & unselected_atom_contributor{
            context.unselected_atom_list.at(i)
        };
        unselected_atom_contributor.initial_seed = *global_median;
        build_result.unselected_selection_record_list.emplace_back(
            UnselectedSecondStageSeedSelectionRecord{
                unselected_atom_serial_id_list.at(i),
                *global_median
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
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        analysis.SetAtomLocalNeighborCountForPeeling(
            *context.at(atom_index).atom,
            neighbor_count_list.at(atom_index));
    }
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
                assembled_polish_provenance.at(atom_index) = previous_polish_provenance.at(atom_index);
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
            assembled_polish_provenance.at(atom_index) = previous_polish_provenance.at(atom_index);
        }
    }
}

SuspiciousBlockActivity QuarantineState::BeginIteration(std::size_t accepted_iteration_count)
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
        if (state.lifecycle != QuarantineLifecycle::Quarantined) continue;
        const auto probation_due{
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
        state_by_target.at(target).lifecycle = QuarantineLifecycle::Probation;
        probation_target_list.emplace_back(target);
        selected_probation_atom_index_set.insert(
            target.atom_index_list.begin(),
            target.atom_index_list.end());
    }
    for (auto & [target, state] : state_by_target)
    {
        if (state.lifecycle == QuarantineLifecycle::Tracking ||
            state.lifecycle == QuarantineLifecycle::Probation)
        {
            continue;
        }
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
        if (state.lifecycle == QuarantineLifecycle::Tracking) continue;
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
            if (std::max(
                change.at(GaussianModel3D::LogPeakHeightCoordinateIndex()),
                change.at(GaussianModel3D::LogWidthCoordinateIndex())) >=
                kTransformedChangeTolerance)
            {
                return true;
            }
        }
        else if (target.kind == QuarantineTargetKind::OffsetGroup)
        {
            if (change.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) >= kTransformedChangeTolerance)
            {
                return true;
            }
        }
        else if (IsTransformedChangeMaterial(change, kTransformedChangeTolerance))
        {
            return true;
        }
    }
    return false;
}

static bool DoesFailureAffectTarget(
    const QuarantineTarget & failure_target,
    const QuarantineTarget & target)
{
    const auto overlaps{
        std::ranges::any_of(
            target.atom_index_list,
            [&](const auto atom_index)
            {
                return std::ranges::binary_search(failure_target.atom_index_list, atom_index);
            })
    };
    if (!overlaps) return false;
    if (target.kind == QuarantineTargetKind::HardFailureCluster) return true;
    if (target.kind == QuarantineTargetKind::OffsetGroup)
    {
        return failure_target.kind != QuarantineTargetKind::ShapeAtom;
    }
    return failure_target.kind != QuarantineTargetKind::OffsetGroup;
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

bool QuarantineState::UpdateAfterIteration(
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
    QuarantineFailureReasonMap failure_reason_by_target;
    for (const auto & [key, health] : health_by_key)
    {
        if (!IsJointOffsetSolveHardFailure(health.joint_offset_status)) continue;
        failure_reason_by_target.try_emplace(
            QuarantineTarget{ QuarantineTargetKind::HardFailureCluster, key },
            health.joint_offset_status);
    }
    const auto append_terminal_observations = [&](const auto & diagnostic_list)
    {
        for (const auto & diagnostic : diagnostic_list)
        {
            for (const auto & terminal : diagnostic.attempt.terminal_diagnostic_list)
            {
                const auto reason{ terminal.reason };
                if (reason == StabilizationTerminalReason::None) continue;
                const StabilizationTerminalFailure failure{ reason, terminal.guard_reason };
                if (reason == StabilizationTerminalReason::GuardInfeasible &&
                    terminal.guard_atom_index.has_value() &&
                    terminal.guard_mode.has_value())
                {
                    const auto atom_index{ *terminal.guard_atom_index };
                    if (*terminal.guard_mode == SuspiciousUpdateMode::OffsetOnly)
                    {
                        std::vector<std::size_t> group_atom_index_list;
                        const auto group_id{ context.at(atom_index).group_id };
                        for (std::size_t member_index = 0; member_index < context.size(); member_index++)
                        {
                            if (context.at(member_index).group_id == group_id)
                            {
                                group_atom_index_list.emplace_back(member_index);
                            }
                        }
                        failure_reason_by_target.try_emplace(
                            QuarantineTarget{
                                QuarantineTargetKind::OffsetGroup,
                                std::move(group_atom_index_list)
                            },
                            failure);
                    }
                    else
                    {
                        failure_reason_by_target.try_emplace(
                            QuarantineTarget{
                                QuarantineTargetKind::ShapeAtom,
                                { atom_index }
                            },
                            failure);
                    }
                    continue;
                }
                failure_reason_by_target.try_emplace(
                    QuarantineTarget{
                        QuarantineTargetKind::HardFailureCluster,
                        diagnostic.key
                    },
                    failure);
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
                failure_reason_by_target,
                [&](const auto & failure)
                {
                    return DoesFailureAffectTarget(failure.first, target);
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
            failure_reason_by_target,
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
    entered_target_count += transition.entered_target_list.size();
    released_target_count += transition.released_target_list.size();
    failed_probation_count += transition.failed_probation_target_list.size();
    return !transition.entered_target_list.empty() ||
        !transition.released_target_list.empty() ||
        !transition.failed_probation_target_list.empty();
}

std::size_t QuarantineState::AtomCount() const
{
    std::set<std::size_t> atom_index_set;
    for (const auto & [target, state] : state_by_target)
    {
        if (state.lifecycle == QuarantineLifecycle::Tracking) continue;
        atom_index_set.insert(target.atom_index_list.begin(), target.atom_index_list.end());
    }
    return atom_index_set.size();
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
            .reason = SuspiciousGaussianReason::InvalidModel,
            .normalized_margin = std::numeric_limits<double>::infinity()
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
    std::vector<std::optional<GaussianModel3DWithUncertainty>> result(
        context.size());
    int refit_thread_size{ options.thread_size };
#ifdef USE_OPENMP
    const bool parallel_refits{
        !IsDebugLogLevelEnabled() &&
        options.thread_size > 1 && context.size() > 1
    };
    if (parallel_refits) refit_thread_size = 1;
#pragma omp parallel for schedule(dynamic) if(parallel_refits) num_threads(options.thread_size)
#endif
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
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
    const auto is_debug_logging_enabled{ IsDebugLogLevelEnabled() };
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
    fixed_point_operator.shape_available_atom_mask.assign(context.size(), 0);
    fixed_point_operator.offset_available_atom_mask.assign(context.size(), 0);
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
                continue;
            }
            const auto proposed_offset{
                offset_result.offset(static_cast<Eigen::Index>(position))
            };
            const auto operator_model{
                previous_state.at(atom_index).mdpde.GetModel().WithOffset(proposed_offset)
            };
            if (!IsValidSecondStageGaussianModel(operator_model))
            {
                continue;
            }
            auto & operator_result{ fixed_point_operator.state.at(atom_index) };
            SetLocalResultOffset(operator_result, proposed_offset);
            fixed_point_operator.offset_available_atom_mask.at(atom_index) = 1;
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
                    current_model_snapshot.selected.at(atom_index) = accepted_model_list.at(member_position);
                    assessment_by_atom.at(atom_index) = accepted_assessment_list.at(member_position);
                }
                else
                {
                    block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                    current_model_snapshot.selected.at(atom_index) = previous_state.at(atom_index).mdpde.GetModel();
                }
            }
            if (!accepted)
            {
                health.all_local_refits_solver_qualified = false;
            }
        }
    }

    bool operator_offsets_complete{ true };
    bool operator_offsets_match_proposal{ true };
    for (const auto & key : cluster_key_list)
    {
        for (const auto atom_index : key)
        {
            if (fixed_point_operator.offset_available_atom_mask.at(atom_index) == 0)
            {
                operator_offsets_complete = false;
                operator_offsets_match_proposal = false;
                continue;
            }
            operator_offsets_match_proposal =
                operator_offsets_match_proposal &&
                fixed_point_operator.state.at(atom_index).mdpde.GetModel().GetOffset() ==
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
    std::vector<std::optional<LocalAtomRefitResult>> refit_result_list(context.size());
    std::vector<std::exception_ptr> refit_exception_list(context.size());
#ifdef USE_OPENMP
    const bool parallel_refits{
        !is_debug_logging_enabled &&
        options.thread_size > 1 &&
        context.size() > 1
    };
#else
    const bool parallel_refits{ false };
#endif
    FitOptions refit_options{ options };
    if (parallel_refits)
    {
        refit_options.thread_size = 1;
    }
    const auto run_refit = [&](std::size_t atom_index)
    {
        try
        {
            refit_result_list.at(atom_index) =
                FitAtomWithJointOffsetFallback(
                    context.at(atom_index),
                    previous_state.at(atom_index),
                    GetFitModel(refit_model_bundle.selected, atom_index),
                    refit_response_cache.at(atom_index),
                    refit_options);
        }
        catch (...)
        {
            refit_exception_list.at(atom_index) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (parallel_refits)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
        for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
        {
            run_refit(atom_index);
        }
    }
    else
#endif
    {
        for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
        {
            run_refit(atom_index);
        }
    }
    for (const auto & exception : refit_exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    for (const auto & key : cluster_key_list)
    {
        auto & health{ health_by_key.at(key) };
        for (const auto atom_index : key)
        {
            auto refit_result{ std::move(refit_result_list.at(atom_index)) };
            if (!refit_result.has_value())
            {
                health.all_local_refits_solver_qualified = false;
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
                    .reason = SuspiciousGaussianReason::InvalidModel,
                    .normalized_margin = std::numeric_limits<double>::infinity()
                };
                continue;
            }
            if (operator_offsets_complete && operator_offsets_match_proposal)
            {
                if (refit_result->unrestricted_mdpde.has_value())
                {
                    fixed_point_operator.state.at(atom_index).mdpde = *refit_result->unrestricted_mdpde;
                    fixed_point_operator.shape_available_atom_mask.at(atom_index) = 1;
                }
            }
            local_refit_status_by_atom.at(atom_index) = refit_result->attempted_refit_status;
            const auto shape_solver_qualified{
                refit_result->attempted_refit_status.has_value() &&
                IsLocalRefitStatusSolverQualified(*refit_result->attempted_refit_status)
            };
            if (!shape_solver_qualified)
            {
                health.all_local_refits_solver_qualified = false;
            }
            if (!refit_result->unrestricted_mdpde.has_value())
            {
                block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            assessment_by_atom.at(atom_index) = refit_result->assessment;
            iteration_state.at(atom_index) = std::move(refit_result->result);
        }
    }
    if (operator_offsets_complete && !operator_offsets_match_proposal)
    {
        const auto unrestricted_shape_list{
            RunUnrestrictedShapeRefits(
                context,
                fixed_point_operator.state,
                group_id_by_atom_index,
                options)
        };
        for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
        {
            if (unrestricted_shape_list.at(atom_index).has_value())
            {
                fixed_point_operator.state.at(atom_index).mdpde = *unrestricted_shape_list.at(atom_index);
                fixed_point_operator.shape_available_atom_mask.at(atom_index) = 1;
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

static bool AreGraphPartitionsEqual(const CouplingGraphPartition & lhs, const CouplingGraphPartition & rhs)
{
    return lhs.boundary_sample_dependency_list.size() == rhs.boundary_sample_dependency_list.size() &&
        lhs.boundary_sample_dependency_list == rhs.boundary_sample_dependency_list &&
        lhs.sample_id_list_by_key == rhs.sample_id_list_by_key;
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

    FinishProgressLine(options.quiet_mode);
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
    const auto elapsed_milliseconds{
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - rebuild_start).count()
    };
    performance_counters.RecordTopologyRebuild(elapsed_milliseconds, partition_changed);

    LogAdaptiveTopologyRebuild(
        options.quiet_mode,
        iteration_state.accepted_iteration_count,
        decision,
        graph_topology,
        rebuilt_topology,
        iteration_state.graph_partition,
        rebuilt_partition,
        partition_changed);

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
        probation_atom_index_set.insert(target.atom_index_list.begin(), target.atom_index_list.end());
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
        iteration_internal::BuildSuspiciousJointOffsetRidgeMultiplierList(ridge_atom_mask)
    };
    const auto iteration_phase_start{ std::chrono::steady_clock::now() };
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

    LogUnrestrictedOperatorAssessments(
        options.quiet_mode,
        raw_iteration_result.assessment_by_atom,
        raw_iteration_result.block_activity);

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
    IterationResult result;
    result.trust_region_update =
        iteration_state.trust_region_state.ApplyRadiusUpdates(
            selection.grow_trust_region_key_list,
            selection.shrink_trust_region_key_list,
            selection.rejected_key_list,
            selection.exhausted_key_list);
    const auto has_quarantine_transition{
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
    if (has_quarantine_transition) selection.final_audit_objective.reset();
    const auto assembled_uses_polish{
        iteration_internal::UsesPolish(assembled_polish_provenance)
    };
    result.accepted_cluster_diagnostic_list = std::move(selection.accepted_cluster_diagnostic_list);
    result.rejected_cluster_diagnostic_list = std::move(selection.rejected_cluster_diagnostic_list);
    result.boundary_reconciliation_diagnostic_list = std::move(selection.boundary_reconciliation_diagnostic_list);
    iteration_state.rollback_atom_mask = raw_iteration_result.block_activity.BuildCombinedFixedAtomMask();
    result.attempt_number = attempt_number;
    result.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.quarantine_atom_count = iteration_state.quarantine_state.AtomCount();
    result.active_atom_count = context.size() - result.quarantine_atom_count;
    result.polish_progress = selection.polish_progress;
    result.suspicious_atom_count = iteration_suspicious_atom_count;
    result.operator_maximum_transformed_change = std::ranges::max(operator_proposal_change_summary.maximum_list);

    if (selection.accepted_key_list.empty())
    {
        result.stop_reason = attempt_number >= kMaximumIterations ?
            SecondStageStopReason::AllRejectedAtMaximumIterations :
            SecondStageStopReason::AllRejectedBacktrackingExhausted;
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
            raw_iteration_result.block_activity)
    };
    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            assembled_state,
            previous_state,
            iteration_state.active_index_list)
    };
    auto certificate{
        SummarizeFixedPointOperator(
            raw_iteration_result.fixed_point_operator,
            previous_state,
            iteration_state.active_index_list,
            cluster_key_list,
            group_id_by_atom_index)
    };
    certificate.accepted_active_movement = SummarizeActiveDofChanges(
        assembled_state,
        previous_state,
        active_population);
    certificate.solver_qualified = AreActiveCoordinatesSolverQualified(
        iteration_state.active_index_list,
        cluster_key_list,
        group_id_by_atom_index,
        raw_iteration_result.block_activity,
        raw_iteration_result.local_refit_status_by_atom,
        raw_iteration_result.health_by_key);
    iteration_state.accepted_iteration_count++;
    iteration_state.accepted_iterations_since_topology_rebuild++;
    result.objective_domain_changed = TryRebuildAdaptiveTopology(
        context,
        options,
        assembled_state,
        assembled_uses_polish,
        graph_topology,
        iteration_state,
        performance_counters);
    if (result.objective_domain_changed)
    {
        iteration_state.quarantine_state.force_probation = true;
    }
    bool improved_best_audit{ false };
    if (!result.objective_domain_changed)
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
            result.rejected_cluster_diagnostic_list,
            [&](const auto & diagnostic)
            {
                return std::ranges::find(
                    result.trust_region_update.changed_key_list,
                    diagnostic.key) !=
                    result.trust_region_update.changed_key_list.end();
            })
    };
    const auto has_pending_quarantine_lifecycle{
        std::ranges::any_of(
            iteration_state.quarantine_state.state_by_target | std::views::values,
            [](const auto & state)
            {
                return state.lifecycle != QuarantineLifecycle::Exhausted;
            })
    };
    if (result.objective_domain_changed || improved_best_audit ||
        changed_rejected_trust_radius || has_pending_quarantine_lifecycle ||
        has_quarantine_transition)
    {
        iteration_state.audit_patience_count = 0;
    }
    else
    {
        iteration_state.audit_patience_count++;
    }

    result.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.accepted_maximum_transformed_change = std::ranges::max(transformed_change_summary.maximum_list);
    result.transformed_change_percentile = certificate.accepted_active_movement.percentile_list;
    certificate.objective_domain_changed = result.objective_domain_changed;
    certificate.quarantine_transition = has_quarantine_transition;
    certificate.suspicious_offset_fallback = has_suspicious_offset_fallback;
    certificate.rejected_cluster = !selection.rejected_key_list.empty();
    if (certificate.ProductionConverged())
    {
        result.stop_reason = SecondStageStopReason::Converged;
    }
    else if (iteration_state.audit_patience_count >= kAuditPatience)
    {
        result.stop_reason = SecondStageStopReason::AuditPatience;
    }

    LogConvergenceSafeguardAudit(options.quiet_mode, result, certificate);

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

struct FinalPolishResidualSafetyResult
{
    FinalPolishResidualSafetyStatus status{ FinalPolishResidualSafetyStatus::NotEvaluated };
    std::optional<ConvergenceCertificate> base{};
    std::optional<ConvergenceCertificate> candidate{};
};

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
        for (std::size_t atom_index = 0; atom_index < ridge_atom_mask.size(); atom_index++)
        {
            if (quarantine_atom_mask.at(atom_index) != 0)
            {
                ridge_atom_mask.at(atom_index) = 1;
            }
        }
        const auto joint_offset_ridge_multiplier_list{
            iteration_internal::BuildSuspiciousJointOffsetRidgeMultiplierList(ridge_atom_mask)
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
        certificate.solver_qualified = AreActiveCoordinatesSolverQualified(
            iteration_state.active_index_list,
            cluster_key_list,
            group_id_by_atom_index,
            proposal_result.block_activity,
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

static bool HasComparableFinalPolishOperatorEvidence(const ConvergenceCertificate & certificate)
{
    const auto & percentile_list{
        certificate.operator_nominal_residual.percentile_list
    };
    return certificate.solver_qualified && certificate.operator_complete &&
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
    for (std::size_t index = 0; index < base_percentile_list.size(); index++)
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
        (residual_safety.status == FinalPolishResidualSafetyStatus::AbsolutePassed ||
            residual_safety.status == FinalPolishResidualSafetyStatus::RelativePassed)
    };
    LogFinalDependencyPolish(
        options.quiet_mode,
        polish_result,
        certification_policy,
        residual_safety.status,
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
        LogSecondStageStart(options.quiet_mode);

        if (initialization.failure !=
            SecondStageInitializationFailure::None)
        {
            LogSecondStageInitializationFailure(
                options.quiet_mode,
                initialization.failure);
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
    IterationResult terminal_result;
    if (context.size() == 0)
    {
        terminal_result.attempt_number = 1;
        terminal_result.accepted_iteration_count =
            iteration_state.accepted_iteration_count;
        terminal_result.stop_reason = SecondStageStopReason::Quarantine;
    }
    else
    {
        for (std::size_t iter = 0; iter < kMaximumIterations; iter++)
        {
            terminal_result = RunIteration(
                context,
                graph_topology,
                options,
                iter + 1,
                iteration_state,
                performance_counters);
            if (terminal_result.objective_domain_changed)
            {
                LogObjectiveDomain(
                    iteration_state.objective_domain,
                    options.quiet_mode,
                    true);
            }
            LogAcceptedCandidateSearchDiagnostics(
                options.quiet_mode,
                terminal_result);
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
            LogTrustModelShadowDiagnostics(
                options.quiet_mode,
                terminal_result);
#endif
            LogRejectedClusterDiagnostics(
                options.quiet_mode,
                terminal_result.rejected_cluster_diagnostic_list);
            LogIterationProgress(
                options.quiet_mode,
                progress_column_widths,
                terminal_result);

            if (terminal_result.stop_reason ==
                    SecondStageStopReason::AllRejectedBacktrackingExhausted ||
                terminal_result.stop_reason ==
                    SecondStageStopReason::AllRejectedAtMaximumIterations)
            {
                LogAllRejectedResolution(
                    options.quiet_mode,
                    terminal_result);
                break;
            }
            if (terminal_result.stop_reason == SecondStageStopReason::Converged ||
                terminal_result.stop_reason == SecondStageStopReason::AuditPatience)
            {
                break;
            }
            if (iter + 1 == kMaximumIterations)
            {
                terminal_result.stop_reason =
                    SecondStageStopReason::MaximumIterations;
                break;
            }
        }
    }

    const auto converged{
        terminal_result.stop_reason == SecondStageStopReason::Converged
    };
    const auto use_best_audit_state{
        !converged &&
        terminal_result.stop_reason != SecondStageStopReason::Quarantine &&
        iteration_state.best_audit_state.has_value()
    };
    FinalizeSecondStageState(
        model_object,
        context,
        options,
        graph_topology,
        iteration_state,
        use_best_audit_state,
        converged ?
            FinalPolishCertificationPolicy::RequireStrictFixedPoint :
            FinalPolishCertificationPolicy::RequireResidualNonRegression,
        performance_counters);
    const auto & finalized_state{
        use_best_audit_state ?
            iteration_state.best_audit_state->state :
            iteration_state.previous_state
    };
    LogSecondStageAuditTerminal(
        options.quiet_mode,
        context,
        terminal_result.stop_reason,
        terminal_result.attempt_number,
        iteration_state.accepted_iteration_count,
        finalized_state,
        audit_comparison_objective_domain);

    if (terminal_result.stop_reason == SecondStageStopReason::Quarantine)
    {
        if (iteration_state.quarantine_state.TargetCount() != 0)
        {
            LogQuarantineFallback(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                iteration_state.quarantine_state.entered_target_count,
                iteration_state.quarantine_state.released_target_count,
                iteration_state.quarantine_state.failed_probation_count,
                iteration_state.quarantine_state.TargetCount(),
                finalized_state);
        }
        LogNoSelectedAtoms(options.quiet_mode);
    }
    else if (converged)
    {
        if (iteration_state.quarantine_state.TargetCount() != 0)
        {
            LogQuarantineFallback(
                options.quiet_mode,
                iteration_state.accepted_iteration_count,
                iteration_state.quarantine_state.entered_target_count,
                iteration_state.quarantine_state.released_target_count,
                iteration_state.quarantine_state.failed_probation_count,
                iteration_state.quarantine_state.TargetCount(),
                finalized_state);
        }
        else
        {
            LogConverged(
                options.quiet_mode,
                terminal_result,
                finalized_state);
        }
    }
    else if (terminal_result.stop_reason ==
        SecondStageStopReason::MaximumIterations)
    {
        LogMaximumIterations(
            options.quiet_mode,
            iteration_state.quarantine_state.entered_target_count,
            iteration_state.quarantine_state.released_target_count,
            iteration_state.quarantine_state.failed_probation_count,
            iteration_state.quarantine_state.TargetCount(),
            iteration_state.best_audit_state,
            iteration_state.previous_state);
    }
    LogSecondStageSummary(
        options.quiet_mode,
        iteration_state.accepted_iteration_count,
        iteration_state.best_audit_state,
        iteration_state.previous_polish_provenance,
        terminal_result.stop_reason,
        use_best_audit_state);
    return true;
}

} // namespace

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(
    const GaussianModel3DWithUncertainty & local_mdpde,
    const std::optional<GaussianModel3D> & global_median)
{
    if (IsValidSecondStageGaussianModel(local_mdpde.GetModel()))
    {
        return SecondStageSeedSelection{
            SecondStageSeedSource::LocalMdpde,
            local_mdpde
        };
    }
    if (global_median.has_value() &&
        IsValidSecondStageGaussianModel(*global_median))
    {
        return SecondStageSeedSelection{
            SecondStageSeedSource::GlobalMedian,
            GaussianModel3DWithUncertainty{
                *global_median,
                GaussianModel3DUncertainty{}
            }
        };
    }
    return std::nullopt;
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
    const auto maximum_transformed_drift{ std::ranges::max(drift_summary.maximum_list) };
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

bool ConvergenceCertificate::StrictOperatorPassed() const
{
    return solver_qualified && operator_complete &&
        IsTransformedPercentileConverged(operator_nominal_residual);
}

bool ConvergenceCertificate::ProductionConverged() const
{
    return StrictOperatorPassed() &&
        IsTransformedPercentileConverged(accepted_active_movement) &&
        !objective_domain_changed && !quarantine_transition &&
        !suspicious_offset_fallback && !rejected_cluster;
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
    const SuspiciousBlockActivity & block_activity)
{
    ValidateBlockActivitySize(group_id_by_atom_index.size(), block_activity);

    ActiveCoordinatePopulation result;
    for (const auto atom_index : atom_index_list)
    {
        if (block_activity.HasActiveShape(atom_index))
        {
            result.active_atom_index_list_by_parameter.at(
                GaussianModel3D::LogPeakHeightCoordinateIndex())
                .emplace_back(atom_index);
            result.active_atom_index_list_by_parameter.at(
                GaussianModel3D::LogWidthCoordinateIndex())
                .emplace_back(atom_index);
        }
        if (block_activity.HasActiveOffset(atom_index))
        {
            result.active_atom_index_list_by_parameter.at(
                GaussianModel3D::OffsetToPeakRatioCoordinateIndex())
                .emplace_back(atom_index);
        }
    }

    const auto offset_group_list{
        BuildOffsetGroupEntries(cluster_key_list, group_id_by_atom_index)
    };
    for (const auto & group : offset_group_list)
    {
        const auto active_count{ static_cast<std::size_t>(std::ranges::count_if(
            group.atom_index_list,
            [&](const auto atom_index)
            {
                return block_activity.HasActiveOffset(atom_index);
            })) };
        if (active_count == 0)
        {
            continue;
        }

        const auto is_mixed{ active_count != group.atom_index_list.size() };
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
                    GaussianModel3D::OffsetToPeakRatioCoordinateIndex())
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

    result.population_size_list.at(
        GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        offset_group_change_list.size();
    result.percentile_list.at(
        GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        array_helper::ComputePercentile(
            offset_group_change_list,
            kConvergencePercentile);
    result.maximum_list.at(
        GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
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

bool AreActiveCoordinatesSolverQualified(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key)
{
    const auto atom_count{ group_id_by_atom_index.size() };
    ValidateBlockActivitySize(atom_count, block_activity);
    if (local_refit_status_by_atom.size() != atom_count)
    {
        throw std::invalid_argument(
            "Convergence audit qualification inputs are inconsistent.");
    }

    for (const auto atom_index : atom_index_list)
    {
        if (!block_activity.HasActiveShape(atom_index)) continue;
        const auto & status{ local_refit_status_by_atom[atom_index] };
        if (!status.has_value() || !IsLocalRefitStatusSolverQualified(*status))
        {
            return false;
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
        if (active_count == 0) continue;
        if (active_count != group.atom_index_list.size())
        {
            return false;
        }

        const auto health_iter{ health_by_key.find(group.cluster_key) };
        if (health_iter == health_by_key.end() ||
            health_iter->second.joint_offset_status !=
                JointOffsetSolveStatus::Converged)
        {
            return false;
        }
    }
    return true;
}

QuarantineStateTransition UpdateQuarantineFailureState(
    const QuarantineFailureReasonMap & failure_reason_by_target,
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
        if (state.lifecycle != QuarantineLifecycle::Probation)
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
        state.probation_count++;
        state.lifecycle =
            state.probation_count >= kQuarantineMaximumProbationCount ?
                QuarantineLifecycle::Exhausted :
                QuarantineLifecycle::Quarantined;
        state.next_probation_iteration =
            accepted_iteration_count + kQuarantineProbationCooldown;
        transition.failed_probation_target_list.emplace_back(target);
        ++iter;
    }

    for (auto iter = state_by_target.begin(); iter != state_by_target.end();)
    {
        if (iter->second.lifecycle == QuarantineLifecycle::Tracking &&
            !failure_reason_by_target.contains(iter->first))
        {
            iter = state_by_target.erase(iter);
            continue;
        }
        ++iter;
    }

    for (const auto & [target, reason] : failure_reason_by_target)
    {
        auto [iter, inserted]{
            state_by_target.try_emplace(
                target,
                QuarantineFailureState{
                    reason,
                    0,
                    0,
                    0,
                    QuarantineLifecycle::Tracking
                })
        };
        auto & state{ iter->second };
        if (state.lifecycle != QuarantineLifecycle::Tracking) continue;
        if (!inserted && state.reason != reason)
        {
            state.reason = reason;
            state.stable_iteration_count = 0;
        }
        state.stable_iteration_count++;
        if (state.stable_iteration_count >= kPersistentQuarantineFailureIterationLimit)
        {
            state.lifecycle = QuarantineLifecycle::Quarantined;
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
    return detail::RunSecondStageIterations(model_object, options);
}

} // namespace rhbm_gem::core
