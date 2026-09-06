#include "core/detail/Quarantine.hpp"

#include "core/detail/CandidateSelection.hpp"

#include <algorithm>
#include <set>

namespace rhbm_gem::core::detail {

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
            if (target.kind == QuarantineTargetKind::OffsetAtom)
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

static void ApplyQuarantineTargetActivity(
    const QuarantineTarget & target,
    SuspiciousBlockActivity & activity)
{
    for (const auto atom_index : target.atom_index_list)
    {
        if (target.kind != QuarantineTargetKind::OffsetAtom)
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
                    return probation_target.kind == QuarantineTargetKind::HardFailureCluster &&
                        target.kind != QuarantineTargetKind::HardFailureCluster &&
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
        ApplyQuarantineTargetActivity(target, activity);
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
        ApplyQuarantineTargetActivity(target, activity);
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
        else if (target.kind == QuarantineTargetKind::OffsetAtom)
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
    if (target.kind == QuarantineTargetKind::OffsetAtom)
    {
        return failure_target.kind != QuarantineTargetKind::ShapeAtom;
    }
    return failure_target.kind != QuarantineTargetKind::OffsetAtom;
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
        if (target.kind != QuarantineTargetKind::OffsetAtom &&
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
                        failure_reason_by_target.try_emplace(
                            QuarantineTarget{
                                QuarantineTargetKind::OffsetAtom,
                                { atom_index }
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
    ApplyQuarantineFallbackTargets(
        transition.entered_target_list,
        previous_state,
        previous_polish_provenance,
        assembled_state,
        assembled_polish_provenance);
    ApplyQuarantineFallbackTargets(
        transition.failed_probation_target_list,
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
