#include "core/detail/Quarantine.hpp"

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/IterationProposal.hpp"

#include <algorithm>
#include <set>

namespace rhbm_gem::core::detail {

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

SuspiciousBlockActivity QuarantineState::BeginIteration(std::size_t domain_revision)
{
    SuspiciousBlockActivity activity{
        SuspiciousUpdateMask(m_atom_count, 0),
        SuspiciousUpdateMask(m_atom_count, 0),
        SuspiciousUpdateMask(m_atom_count, 0)
    };
    retry_target_list.clear();
    for (auto & [target, state] : state_by_target)
    {
        if (state.lifecycle != QuarantineLifecycle::Frozen) continue;
        if (domain_revision > state.last_domain_revision)
        {
            state.last_domain_revision = domain_revision;
            retry_target_list.emplace_back(target);
        }
        else
        {
            ApplyQuarantineTargetActivity(target, activity);
        }
    }
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
        if (state.lifecycle == QuarantineLifecycle::Active) continue;
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
    const ClusterHealthMap & health_by_key,
    const FixedPointOperatorEvidence & operator_evidence,
    const FitState & previous_state)
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
        if (atom_index >= operator_evidence.state.size() ||
            atom_index >= operator_evidence.shape_available_atom_mask.size() ||
            atom_index >= operator_evidence.offset_available_atom_mask.size()) return false;
        const auto change{ CalculateTransformedChange(operator_evidence.state.at(atom_index),
            previous_state.at(atom_index).mdpde.GetModel()) };
        if (target.kind != QuarantineTargetKind::OffsetAtom &&
            (!block_activity.HasActiveShape(atom_index) ||
                !operator_evidence.shape_available_atom_mask.at(atom_index) ||
                !(change.at(GaussianModel3D::LogPeakHeightCoordinateIndex()) < kTransformedChangeTolerance) ||
                !(change.at(GaussianModel3D::LogWidthCoordinateIndex()) < kTransformedChangeTolerance))) return false;
        if (target.kind != QuarantineTargetKind::ShapeAtom &&
            (!block_activity.HasActiveOffset(atom_index) ||
                !operator_evidence.offset_available_atom_mask.at(atom_index) ||
                !(change.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) < kTransformedChangeTolerance))) return false;
    }
    return true;
}

bool QuarantineState::UpdateAfterIteration(
    std::span<const ClusterCandidateDiagnostic> accepted_diagnostic_list,
    std::span<const ClusterCandidateDiagnostic> rejected_diagnostic_list,
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    const ClusterHealthMap & health_by_key,
    const FixedPointOperatorEvidence & operator_evidence,
    const FitState & assembled_state,
    const FitState & previous_state,
    std::size_t domain_revision)
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

    std::vector<QuarantineTarget> successful_retry_target_list;
    for (const auto & target : retry_target_list)
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
        const auto fully_active{ std::ranges::all_of(target.atom_index_list, [&](const auto atom_index)
        {
            return (target.kind == QuarantineTargetKind::OffsetAtom || block_activity.HasActiveShape(atom_index)) &&
                (target.kind == QuarantineTargetKind::ShapeAtom || block_activity.HasActiveOffset(atom_index));
        }) };
        if (!has_affecting_observation && fully_active &&
            (accepted_material_proposal ||
                IsGuardSafeNonMaterialSolverQualifiedEndpoint(
                    target,
                    block_activity,
                    assessment_by_atom,
                    health_by_key, operator_evidence, previous_state)))
        {
            successful_retry_target_list.emplace_back(target);
        }
    }
    auto transition{
        UpdateQuarantineFailureState(
            failure_reason_by_target,
            retry_target_list,
            successful_retry_target_list,
            domain_revision,
            state_by_target)
    };
    entered_target_count += transition.entered_target_list.size();
    released_target_count += transition.released_target_list.size();
    failed_retry_count += transition.failed_retry_target_list.size();
    return !transition.entered_target_list.empty() ||
        !transition.released_target_list.empty() ||
        !transition.failed_retry_target_list.empty();
}

std::size_t QuarantineState::AtomCount() const
{
    std::set<std::size_t> atom_index_set;
    for (const auto & [target, state] : state_by_target)
    {
        if (state.lifecycle == QuarantineLifecycle::Active) continue;
        atom_index_set.insert(target.atom_index_list.begin(), target.atom_index_list.end());
    }
    return atom_index_set.size();
}

QuarantineStateTransition UpdateQuarantineFailureState(
    const QuarantineFailureReasonMap & failure_reason_by_target,
    const std::vector<QuarantineTarget> & retry_target_list,
    const std::vector<QuarantineTarget> & successful_retry_target_list,
    std::size_t domain_revision,
    QuarantineFailureStateMap & state_by_target)
{
    QuarantineStateTransition transition;
    for (const auto & target : retry_target_list)
    {
        auto & state{ state_by_target.at(target) };
        if (std::ranges::find(successful_retry_target_list, target) != successful_retry_target_list.end())
        {
            transition.released_target_list.emplace_back(target);
            state_by_target.erase(target);
        }
        else
        {
            state.last_domain_revision = domain_revision;
            transition.failed_retry_target_list.emplace_back(target);
        }
    }
    for (auto iter = state_by_target.begin(); iter != state_by_target.end();)
    {
        if (iter->second.lifecycle == QuarantineLifecycle::Active &&
            !failure_reason_by_target.contains(iter->first))
        {
            iter = state_by_target.erase(iter);
            continue;
        }
        ++iter;
    }
    for (const auto & [target, reason] : failure_reason_by_target)
    {
        auto [iter, inserted]{ state_by_target.try_emplace(target,
            QuarantineFailureState{ .reason = reason }) };
        auto & state{ iter->second };
        if (state.lifecycle != QuarantineLifecycle::Active) continue;
        if (!inserted && state.reason != reason)
        {
            state.reason = reason;
            state.stable_iteration_count = 0;
        }
        state.stable_iteration_count++;
        if (state.stable_iteration_count >= kPersistentQuarantineFailureIterationLimit)
        {
            state.lifecycle = QuarantineLifecycle::Frozen;
            state.last_domain_revision = domain_revision;
            transition.entered_target_list.emplace_back(target);
        }
    }
    return transition;
}

} // namespace rhbm_gem::core::detail
