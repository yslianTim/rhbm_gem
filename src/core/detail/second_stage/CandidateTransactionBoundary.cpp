#include "core/detail/second_stage/ComponentAssembly.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"

#include "core/detail/second_stage/observation/PerformanceCounters.hpp"

#include <algorithm>
#include <chrono>
#include <iterator>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core::detail {
static ClusterKey FlattenClusterKeyList(const std::vector<ClusterKey> & key_list)
{
    ClusterKey atom_index_list;
    for (const auto & key : key_list)
    {
        atom_index_list.insert(
            atom_index_list.end(),
            key.begin(),
            key.end());
    }
    std::ranges::sort(atom_index_list);
    atom_index_list.erase(
        std::ranges::unique(atom_index_list).begin(),
        atom_index_list.end());
    return atom_index_list;
}

void CandidateTransactionBuilder::RejectSelectionKeys(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & key_list,
    bool exhausted)
{
    auto & selection{ m_selection };
    for (const auto & key : key_list)
    {
        auto & candidate{ m_candidate_by_key.at(key) };
        if (!candidate.selected) continue;
        for (const auto atom_index : key)
        {
            selection.assembled_state.at(atom_index) = inputs.previous_state.at(atom_index);
            selection.assembled_polish_provenance.at(atom_index) = inputs.previous_polish_provenance.at(atom_index);
        }
        candidate.selected = false;
        candidate.rejection_order = m_next_rejection_order++;
        if (exhausted) candidate.exhausted = true;
    }
}

static FitStatePatch BuildSelectionPatch(
    const CandidateSelection & selection,
    const std::vector<ClusterKey> & key_list)
{
    return FitStatePatch::FromState(
        selection.assembled_state,
        FlattenClusterKeyList(key_list));
}

static bool OverlayFitStatePatch(FitStatePatch & base_patch, const FitStatePatch & overlay_patch)
{
    for (std::size_t position = 0; position < overlay_patch.atom_index_list.size(); position++)
    {
        const auto atom_index{ overlay_patch.atom_index_list.at(position) };
        const auto iter{
            std::ranges::lower_bound(base_patch.atom_index_list, atom_index)
        };
        if (iter == base_patch.atom_index_list.end() || *iter != atom_index)
        {
            return false;
        }
        base_patch.mdpde_list.at(static_cast<std::size_t>(std::distance(
            base_patch.atom_index_list.begin(),
            iter))) = overlay_patch.mdpde_list.at(position);
    }
    return true;
}

std::optional<CandidateTransactionBuilder::ComponentCandidate>
CandidateTransactionBuilder::TryBoundaryJointCorrection(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown & previous_audit_objective,
    const ObjectiveBreakdown & improvement_reference_objective,
    const FitStatePatch & endpoint_patch,
    BoundaryComponentDecision & decision,
    BoundaryObservationScope & observations, BoundaryAcceptancePolicy policy)
{
    auto & observation{ observations.Trials() };
    auto & selection{ m_selection };
    if (component.halo_atom_index_list.empty()) return std::nullopt;

    std::vector<std::size_t> shape_active_atom_index_list;
    for (const auto atom_index : component.halo_atom_index_list)
    {
        if (selection.block_activity.HasActiveShape(atom_index))
        {
            shape_active_atom_index_list.emplace_back(atom_index);
        }
    }
    std::vector<std::size_t> offset_active_atom_index_list;
    for (const auto atom_index : component.halo_atom_index_list)
    {
        if (selection.block_activity.HasActiveOffset(atom_index))
        {
            offset_active_atom_index_list.emplace_back(atom_index);
        }
    }
    if (shape_active_atom_index_list.empty() &&
        offset_active_atom_index_list.empty())
    {
        return std::nullopt;
    }
    const FitStateView endpoint_state_view{
        inputs.previous_state,
        endpoint_patch
    };
    std::vector<BoundaryJointTrustRegion> trust_region_list;
    trust_region_list.reserve(component.key_list.size());
    for (const auto & key : component.key_list)
    {
        trust_region_list.emplace_back(BoundaryJointTrustRegion{
            key,
            inputs.trust_region_state.GetRadius(key)
        });
    }
    const BoundaryJointCorrectionWorkspaceKey workspace_key{
        shape_active_atom_index_list,
        offset_active_atom_index_list,
        component.affected_sample_ref_list
    };
    auto & solver{
        inputs.boundary_joint_correction_workspace_by_key.try_emplace(workspace_key).first->second
    };
    const auto start_time{ std::chrono::steady_clock::now() };
    const auto record_performance = [&]
    {
        inputs.performance_counters.RecordBoundaryJointCorrection(std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start_time).count());
    };
    auto correction_result{
        BuildBoundaryJointCorrection(
            inputs.context,
            endpoint_state_view,
            shape_active_atom_index_list,
            offset_active_atom_index_list,
            component.affected_sample_ref_list,
            inputs.ridge_multiplier_list,
            trust_region_list,
            solver)
    };
    if (correction_result.status != BoundaryJointCorrectionStatus::CandidateReady ||
        !correction_result.patch.has_value())
    {
        observations.BeginTrial(BoundaryObservationStage::Correction, correction_result.damping);
        observation.RejectCorrectionUnavailable();
        record_performance();
        return std::nullopt;
    }

    auto corrected_component_patch{ endpoint_patch };
    if (!OverlayFitStatePatch(
            corrected_component_patch,
            *correction_result.patch))
    {
        record_performance();
        return std::nullopt;
    }
    const CandidateEvaluationOverlay corrected_overlay{
        inputs.context,
        inputs.residual_baseline,
        inputs.previous_state,
        corrected_component_patch
    };
    observation.BeginBoundary(policy, BoundaryObservationStage::Correction, correction_result.damping);
    const auto correction_accepted{ EvaluateBoundaryCorrection(corrected_overlay,
        BoundaryCandidateReference{
            .policy = policy,
            .samples_by_key = inputs.partition.sample_id_list_by_key,
            .domain = inputs.objective_domain,
            .previous_objective_by_key = inputs.previous_objective_by_key,
            .best_audit = inputs.best_audit_state ? &inputs.best_audit_state->objective : nullptr,
            .counters = inputs.performance_counters,
            .component = component}, endpoint_state_view, previous_audit_objective,
        improvement_reference_objective, &observation) };
    if (!correction_accepted)
    {
        record_performance();
        return std::nullopt;
    }

    ComponentCandidate candidate{ .patch = std::move(corrected_component_patch) };
    const FitStateView candidate_state{ inputs.previous_state, candidate.patch };
    for (const auto atom_index : component.halo_atom_index_list)
    {
        const auto change{
            CalculateTransformedChange(
                candidate_state.GetModel(atom_index),
                endpoint_state_view.GetModel(atom_index))
        };
        if (IsTransformedChangeMaterial(change, kTransformedChangeTolerance))
        {
            candidate.provenance_updates.emplace_back(atom_index, 1);
        }
    }
    decision.accepted_source = BoundaryComponentAcceptedSource::JointCorrection;
    record_performance();
    return candidate;
}

std::optional<CandidateTransactionBuilder::ComponentCandidate>
CandidateTransactionBuilder::TryBacktrackBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    const FitStatePatch & endpoint_patch,
    BoundaryComponentDecision & decision,
    BoundaryObservationScope & observations, BoundaryAcceptancePolicy policy)
{
    auto & observation{ observations.Trials() };
    auto & selection{ m_selection };
    BacktrackingWorkspace backtracking_workspace{
        inputs.previous_state,
        endpoint_patch,
        kTransformedChangeTolerance
    };
    BacktrackingStep step;
    std::optional<ObjectiveBreakdown> accepted_evaluation;
    for (step = backtracking_workspace.BuildNextCandidate();
        step.status == BacktrackingStepStatus::CandidateReady;
        step = backtracking_workspace.BuildNextCandidate())
    {
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            inputs.previous_state,
            backtracking_workspace.GetCandidatePatch()
        };
        observations.BeginTrial(BoundaryObservationStage::Backtracking, step.factor, step.trial_number);
        accepted_evaluation = EvaluateBoundaryCandidate(candidate_overlay,
            BoundaryCandidateReference{
                .policy = policy,
                .samples_by_key = inputs.partition.sample_id_list_by_key,
                .domain = inputs.objective_domain,
                .previous_objective_by_key = inputs.previous_objective_by_key,
                .best_audit = inputs.best_audit_state ? &inputs.best_audit_state->objective : nullptr,
                .counters = inputs.performance_counters,
                .component = component}, previous_audit_objective, &observation);
        if (accepted_evaluation.has_value())
        {
            break;
        }
    }
    if (!accepted_evaluation.has_value())
    {
        decision.exhausted = step.status == BacktrackingStepStatus::Exhausted;
        return std::nullopt;
    }

    ComponentCandidate candidate{ .patch = backtracking_workspace.GetCandidatePatch() };
    const auto reconciled_provenance{
        backtracking_workspace.BuildCandidatePolishProvenance(
            inputs.previous_polish_provenance,
            selection.assembled_polish_provenance)
    };
    for (const auto atom_index : FlattenClusterKeyList(component.key_list))
    {
        candidate.provenance_updates.emplace_back(atom_index, reconciled_provenance.at(atom_index));
    }
    decision.accepted_factor = step.factor;
    decision.accepted_source = BoundaryComponentAcceptedSource::Backtracking;
    return candidate;
}

void CandidateTransactionBuilder::ApplyComponentCandidate(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const FitStatePatch & endpoint_patch,
    ComponentCandidate candidate,
    BoundaryComponentAcceptedSource accepted_source)
{
    auto & selection{ m_selection };
    const FitStatePatch * patch{ &candidate.patch };
    ApplyComponentPatches(selection.assembled_state, { &patch, 1 });
    for (const auto & [atom_index, provenance] : candidate.provenance_updates)
        selection.assembled_polish_provenance.at(atom_index) = provenance;

    const FitStateView endpoint_state{ inputs.previous_state, endpoint_patch };
    for (const auto & key : component.key_list)
    {
        auto & pending{ m_candidate_by_key.at(key) };
        if (pending.selected) continue;
        pending.selected = true;
        pending.exhausted = false;
        for (const auto atom_index : key)
        {
            const auto & model{ selection.assembled_state.at(atom_index).mdpde.GetModel() };
            if (!IsTransformedChangeMaterial(CalculateTransformedChange(
                    model, inputs.previous_state.at(atom_index).mdpde.GetModel()),
                    kTransformedChangeTolerance))
            {
                selection.assembled_polish_provenance.at(atom_index) = inputs.previous_polish_provenance.at(atom_index);
                continue;
            }
            const auto correction_changed_endpoint{
                accepted_source == BoundaryComponentAcceptedSource::JointCorrection &&
                IsTransformedChangeMaterial(CalculateTransformedChange(
                    model, endpoint_state.GetModel(atom_index)), kTransformedChangeTolerance)
            };
            selection.assembled_polish_provenance.at(atom_index) = correction_changed_endpoint ? 1 : 0;
        }
    }
}

bool CandidateTransactionBuilder::ReconcileBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    BoundaryAcceptancePolicy policy)
{
    auto & selection{ m_selection };
    const bool cooperative{ policy == BoundaryAcceptancePolicy::CooperativeRescue };
    std::vector<ClusterKey> cooperative_key_list;
    if (cooperative)
    {
        for (const auto & key : component.key_list)
        {
            const auto & pending{ m_candidate_by_key.at(key) };
            if (!pending.selected && pending.cooperative_patch)
                cooperative_key_list.emplace_back(key);
        }
        if (cooperative_key_list.empty()) return false;
    }

    BoundaryComponentDecision decision;
    BoundaryObservationScope observations(inputs, component, policy);
    auto & observation{ observations.Trials() };

    auto endpoint_patch{ BuildSelectionPatch(selection, component.key_list) };
    for (const auto & key : cooperative_key_list)
    {
        if (!OverlayFitStatePatch(endpoint_patch, *m_candidate_by_key.at(key).cooperative_patch))
            throw std::logic_error("Boundary rescue candidate patch does not match its component.");
    }
    const CandidateEvaluationOverlay endpoint_overlay{
        inputs.context, inputs.residual_baseline, inputs.previous_state, endpoint_patch
    };
    observations.BeginTrial(BoundaryObservationStage::Endpoint, 1.0);
    const auto endpoint_evaluation{ EvaluateBoundaryCandidate(endpoint_overlay,
        BoundaryCandidateReference{
            .policy = policy,
            .samples_by_key = inputs.partition.sample_id_list_by_key,
            .domain = inputs.objective_domain,
            .previous_objective_by_key = inputs.previous_objective_by_key,
            .best_audit = inputs.best_audit_state ? &inputs.best_audit_state->objective : nullptr,
            .counters = inputs.performance_counters,
            .component = component}, previous_audit_objective, &observation) };

    std::optional<ComponentCandidate> accepted;
    if (previous_audit_objective != nullptr)
    {
        accepted = TryBoundaryJointCorrection(inputs, component, *previous_audit_objective,
            endpoint_evaluation ? *endpoint_evaluation : *previous_audit_objective,
            endpoint_patch, decision, observations, policy);
    }
    if (!accepted && endpoint_evaluation)
    {
        accepted = ComponentCandidate{ .patch = endpoint_patch };
        decision.accepted_factor = 1.0;
        decision.accepted_source = BoundaryComponentAcceptedSource::Endpoint;
    }
    if (!accepted)
        accepted = TryBacktrackBoundaryComponent(inputs, component, previous_audit_objective,
            endpoint_patch, decision, observations, policy);
    if (accepted)
    {
        ApplyComponentCandidate(inputs, component, endpoint_patch, std::move(*accepted),
            decision.accepted_source);
    }
    else if (!cooperative)
        RejectSelectionKeys(inputs, component.key_list, decision.exhausted);
    observations.Finish(decision);
    return accepted.has_value();
}

static std::vector<BoundaryReconciliationComponent>
BuildExpandedBoundaryReconciliationComponents(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & key_list)
{
    auto component_list{
        BuildBoundaryReconciliationComponents(
            inputs.context,
            inputs.partition,
            key_list)
    };
    for (auto & component : component_list)
    {
        component = ExpandBoundaryReconciliationHalo(
            inputs.context,
            std::move(component),
            inputs.options.second_stage_boundary_halo_depth);
    }
    return component_list;
}

bool CandidateTransactionBuilder::ReconcileCooperativeComponents(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective)
{
    std::vector<ClusterKey> eligible_key_list;
    for (const auto & [key, candidate] : m_candidate_by_key)
    {
        if (candidate.selected || candidate.cooperative_patch)
            eligible_key_list.emplace_back(key);
    }
    bool accepted_any{ false };
    for (const auto & component : BuildExpandedBoundaryReconciliationComponents(inputs, eligible_key_list))
    {
        accepted_any = ReconcileBoundaryComponent(inputs, component, &previous_audit_objective,
            BoundaryAcceptancePolicy::CooperativeRescue) || accepted_any;
    }
    return accepted_any;
}

static std::optional<ObjectiveBreakdown> EvaluateFinalSelectionAudit(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const CandidateSelection & selection,
    const std::vector<ClusterKey> & selected_key_list, bool rescue_audit)
{
    if (selected_key_list.empty()) return std::nullopt;
    const auto candidate_patch{
        BuildSelectionPatch(selection, selected_key_list)
    };
    const CandidateEvaluationOverlay candidate_overlay{
        inputs.context,
        inputs.residual_baseline,
        inputs.previous_state,
        candidate_patch
    };
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(inputs.partition, selected_key_list)
    };
    const auto * best_audit_objective{
        inputs.best_audit_state.has_value() ? &inputs.best_audit_state->objective : nullptr
    };
    return EvaluateGlobalCandidate(candidate_overlay,
        GlobalCandidateReference{affected_sample_ref_list, inputs.objective_domain,
            best_audit_objective, &previous_audit_objective, inputs.performance_counters}, inputs.observation, rescue_audit);
}

static std::vector<std::pair<double, std::vector<ClusterKey>>> BuildRejectionCandidates(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const CandidateSelection & selection,
    const std::vector<ClusterKey> & selected_key_list)
{
    auto component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            selected_key_list)
    };
    for (const auto & key : selected_key_list)
    {
        if (std::ranges::any_of(
                component_list,
                [&](const auto & component)
                {
                    return std::ranges::find(component.key_list, key) !=
                        component.key_list.end();
                }))
        {
            continue;
        }
        component_list.emplace_back(BoundaryReconciliationComponent{
            .key_list = { key },
            .affected_sample_ref_list =
                inputs.partition.sample_id_list_by_key.at(key)
        });
    }

    std::vector<std::pair<double, std::vector<ClusterKey>>>
        rejection_candidate_list;
    for (const auto & component : component_list)
    {
        const auto candidate_patch{
            BuildSelectionPatch(selection, component.key_list)
        };
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            inputs.previous_state,
            candidate_patch
        };
        const auto candidate_audit_objective{
            EvaluateObjectiveDelta(
                candidate_overlay,
                component.affected_sample_ref_list,
                inputs.objective_domain,
                previous_audit_objective,
                inputs.performance_counters)
        };
        if (candidate_audit_objective.has_value() &&
            IsBetterAuditObjective(
                candidate_audit_objective->GetTotalObjective(),
                previous_audit_objective.GetTotalObjective(),
                kObjectiveStrictTolerance))
        {
            continue;
        }
        rejection_candidate_list.emplace_back(
            candidate_audit_objective.has_value() ?
                candidate_audit_objective->GetTotalObjective() :
                std::numeric_limits<double>::infinity(),
            component.key_list);
    }
    std::ranges::sort(
        rejection_candidate_list,
        [](const auto & lhs, const auto & rhs)
        {
            if (lhs.first != rhs.first)
            {
                return lhs.first > rhs.first;
            }
            return lhs.second < rhs.second;
        });

    return rejection_candidate_list;
}

void CandidateTransactionBuilder::AuditAndSalvageFinalSelection(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective, bool rescue_audit)
{
    auto & selection{ m_selection };
    const auto initial_count{ SelectedKeys().size() };
    const auto evaluate = [&]
    {
        return EvaluateFinalSelectionAudit(
            inputs, previous_audit_objective, selection, SelectedKeys(), rescue_audit);
    };
    selection.final_audit_objective = evaluate();
    if (!selection.final_audit_objective.has_value() && !SelectedKeys().empty())
    {
        const auto rejection_candidate_list{ BuildRejectionCandidates(
            inputs, previous_audit_objective, selection, SelectedKeys()) };
        for (const auto & candidate : rejection_candidate_list)
        {
            RejectSelectionKeys(inputs, candidate.second, false);
            selection.final_audit_objective = evaluate();
            if (selection.final_audit_objective.has_value() || SelectedKeys().empty()) break;
        }
    }
    if (selection.final_audit_objective.has_value() || SelectedKeys().empty())
    {
        const auto remaining{ SelectedKeys().size() };
        if (inputs.observation) inputs.observation->ObserveSelectionAudit(rescue_audit, true,
            remaining ? "passed" : "empty_after_salvage", remaining ? "" : "no-selection-remains", initial_count - remaining);
        return;
    }
    if (inputs.observation) inputs.observation->ObserveSelectionAudit(rescue_audit, true, "rejected", "salvage-failed", initial_count);

    const auto remaining_key_list{ SelectedKeys() };
    RejectSelectionKeys(
        inputs,
        remaining_key_list,
        true);
    selection.final_audit_objective.reset();
}

void CandidateTransactionBuilder::ReconcileSelectedBoundaries(
    const CandidateSelectionInputs & inputs)
{
    const auto boundary_component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            SelectedKeys())
    };
    const auto previous_audit_objective{
        EvaluateAuditObjective(inputs.objective_domain, inputs.residual_baseline)
    };
    if (inputs.observation) inputs.observation->ObserveScoreReferences(previous_audit_objective,
        inputs.best_audit_state ? &inputs.best_audit_state->objective : nullptr);
    if (!boundary_component_list.empty())
    {
        const auto boundary_reconciliation_start{ std::chrono::steady_clock::now() };
        for (const auto & component : boundary_component_list)
        {
            ReconcileBoundaryComponent(
                inputs,
                component,
                previous_audit_objective.has_value() ?
                    &*previous_audit_objective : nullptr, BoundaryAcceptancePolicy::Ordinary);
        }
        if (!previous_audit_objective.has_value())
        {
            if (inputs.observation) inputs.observation->ObserveSelectionAudit(false, false, "unavailable", "previous-objective-unavailable");
            const auto remaining_key_list{ SelectedKeys() };
            RejectSelectionKeys(
                inputs,
                remaining_key_list,
                true);
        }
        else
        {
            AuditAndSalvageFinalSelection(
                inputs,
                *previous_audit_objective);
        }
        inputs.performance_counters.RecordBoundaryReconciliation(std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - boundary_reconciliation_start).count());
    }
    if (previous_audit_objective.has_value() &&
        ReconcileCooperativeComponents(inputs, *previous_audit_objective))
    {
        AuditAndSalvageFinalSelection(
            inputs,
            *previous_audit_objective, true);
    }
    if (!previous_audit_objective)
        if (inputs.observation) inputs.observation->ObserveSelectionAudit(true, false, "unavailable", "previous-objective-unavailable");
    MaterializeSelection();
}
} // namespace rhbm_gem::core::detail
