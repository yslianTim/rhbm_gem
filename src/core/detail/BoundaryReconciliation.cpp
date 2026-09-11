#include "core/detail/ClusterHistoryObserver.hpp"
#include "core/detail/PhaseAudit.hpp"
#include "core/detail/CandidateTransaction.hpp"
#include "core/detail/CandidateEvaluation.hpp"

#include "core/detail/Diagnosis.hpp"

#include <algorithm>
#include <chrono>
#include <iterator>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core::detail {

static bool ContainsClusterKey(const std::vector<ClusterKey> & key_list, const ClusterKey & key)
{
    return std::ranges::find(key_list, key) != key_list.end();
}

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
        if (inputs.context.cluster_history) inputs.context.cluster_history->Reject(key);
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
    BoundaryComponentReconciliationDiagnostic & diagnostic,
    CandidateScope scope)
{
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
    diagnostic.shape_active_atom_count = shape_active_atom_index_list.size();
    diagnostic.offset_active_atom_count = offset_active_atom_index_list.size();
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
    const auto record_performance = [&](bool accepted)
    {
        inputs.performance_counters.RecordBoundaryJointCorrection(
            accepted,
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start_time).count());
    };
    diagnostic.joint_reference_component_objective = improvement_reference_objective.GetTotalObjective();
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
    diagnostic.joint_correction_status = correction_result.status;
    diagnostic.joint_parameter_count = correction_result.parameter_count;
    if (correction_result.status == BoundaryJointCorrectionStatus::CandidateReady)
    {
        diagnostic.joint_damping = correction_result.damping;
        diagnostic.maximum_normalized_trust_step = correction_result.maximum_normalized_trust_step;
    }
    if (correction_result.status != BoundaryJointCorrectionStatus::CandidateReady ||
        !correction_result.patch.has_value())
    {
        record_performance(false);
        return std::nullopt;
    }

    auto corrected_component_patch{ endpoint_patch };
    if (!OverlayFitStatePatch(
            corrected_component_patch,
            *correction_result.patch))
    {
        record_performance(false);
        return std::nullopt;
    }
    const CandidateEvaluationOverlay corrected_overlay{
        inputs.context,
        inputs.residual_baseline,
        inputs.previous_state,
        corrected_component_patch
    };
    const auto correction_evaluation{ EvaluateCandidate(corrected_overlay,
        (scope == CandidateScope::CooperativeRescue) ? CandidateScope::CooperativeRescue : CandidateScope::Boundary,
        BoundaryCorrectionReference{inputs, component, endpoint_state_view, previous_audit_objective,
            improvement_reference_objective, diagnostic.objective_diagnostic_list, correction_result.damping}) };
    diagnostic.suspicious_candidate_atom_count = correction_evaluation.suspicious_atom_count;
    if (diagnostic.suspicious_candidate_atom_count != 0)
    {
        ObservePhaseCorrection(inputs.context,
            (scope == CandidateScope::CooperativeRescue) ? "rescue-correction" : "boundary-correction", corrected_component_patch.atom_index_list,
            corrected_overlay.GetState(), endpoint_state_view, correction_result.damping, "rejected", "suspicious",
            inputs, component.key_list, improvement_reference_objective);
        record_performance(false);
        return std::nullopt;
    }
    if (correction_evaluation.raw_objective)
        diagnostic.joint_candidate_component_objective = correction_evaluation.raw_objective->GetTotalObjective();
    const auto & candidate_evaluation{ correction_evaluation.members };
    auto * record{ correction_evaluation.record };
    const auto is_strict_improvement{ correction_evaluation.accepted };
    ObservePhaseCorrection(inputs.context,
        (scope == CandidateScope::CooperativeRescue) ? "rescue-correction" : "boundary-correction", corrected_component_patch.atom_index_list,
        corrected_overlay.GetState(), endpoint_state_view, correction_result.damping,
        is_strict_improvement ? "accepted" : "rejected",
        candidate_evaluation && !is_strict_improvement ? "strict-improvement" : (record ? record->outcome : ""),
        inputs, component.key_list, improvement_reference_objective);
    if (!is_strict_improvement)
    {
        if (record && candidate_evaluation)
            record->outcome = "members-passed-strict-improvement-failed";
        record_performance(false);
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
    candidate.history_observation = record ? record->history_observation : 0;
    diagnostic.accepted_source = BoundaryComponentAcceptedSource::JointCorrection;
    diagnostic.candidate_component_objective = candidate_evaluation->audit_objective.GetTotalObjective();
    diagnostic.locally_deteriorated_member_count = candidate_evaluation->locally_deteriorated_member_count;
    diagnostic.maximum_local_deterioration = candidate_evaluation->maximum_local_deterioration;
    record_performance(true);
    return candidate;
}

std::optional<CandidateTransactionBuilder::ComponentCandidate>
CandidateTransactionBuilder::TryBacktrackBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    const FitStatePatch & endpoint_patch,
    BoundaryComponentReconciliationDiagnostic & diagnostic,
    CandidateScope scope)
{
    auto & selection{ m_selection };
    BacktrackingWorkspace backtracking_workspace{
        inputs.previous_state,
        endpoint_patch,
        kTransformedChangeTolerance
    };
    BacktrackingStep step;
    std::optional<BoundaryCandidateEvaluation> accepted_evaluation;
    std::size_t accepted_history_observation{ 0 };
    for (step = backtracking_workspace.BuildNextCandidate();
        step.status == BacktrackingStepStatus::CandidateReady;
        step = backtracking_workspace.BuildNextCandidate())
    {
        diagnostic.trial_count = step.trial_number;
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            inputs.previous_state,
            backtracking_workspace.GetCandidatePatch()
        };
        auto * record{ BeginJointCandidateDiagnostic(inputs.options.quiet_mode,
            diagnostic.objective_diagnostic_list,
            (scope == CandidateScope::CooperativeRescue) ? "rescue-backtracking" : "backtracking", step.factor) };
        accepted_evaluation = EvaluateCandidate(candidate_overlay, (scope == CandidateScope::CooperativeRescue) ? CandidateScope::CooperativeRescue : CandidateScope::Boundary,
            BoundaryCandidateReference{inputs, component, previous_audit_objective, record});
        ObservePhaseCandidate(inputs.context,
            (scope == CandidateScope::CooperativeRescue) ? "rescue-backtracking" : "boundary-backtracking", endpoint_patch.atom_index_list,
            candidate_overlay.GetState(), nullptr, step.factor, accepted_evaluation ? "accepted" : "rejected",
            record ? record->outcome : "", false, false);
        if (accepted_evaluation.has_value())
        {
            accepted_history_observation = record ? record->history_observation : 0;
            break;
        }
    }
    if (!accepted_evaluation.has_value())
    {
        diagnostic.exhausted = step.status == BacktrackingStepStatus::Exhausted;
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
    candidate.history_observation = accepted_history_observation;
    diagnostic.accepted_factor = step.factor;
    diagnostic.accepted_source = BoundaryComponentAcceptedSource::Backtracking;
    diagnostic.candidate_component_objective = accepted_evaluation->audit_objective.GetTotalObjective();
    diagnostic.locally_deteriorated_member_count = accepted_evaluation->locally_deteriorated_member_count;
    diagnostic.maximum_local_deterioration = accepted_evaluation->maximum_local_deterioration;
    return candidate;
}

void CandidateTransactionBuilder::ApplyComponentCandidate(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const FitStatePatch & endpoint_patch,
    ComponentCandidate candidate,
    BoundaryComponentAcceptedSource accepted_source,
    CandidateScope scope)
{
    auto & selection{ m_selection };
    candidate.patch.ApplyTo(selection.assembled_state);
    for (const auto & [atom_index, provenance] : candidate.provenance_updates)
        selection.assembled_polish_provenance.at(atom_index) = provenance;
    if (inputs.context.cluster_history)
        inputs.context.cluster_history->AcceptBoundary(candidate.history_observation);

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
        if (pending.diagnostic)
            pending.diagnostic->boundary_rescued = scope == CandidateScope::CooperativeRescue;
    }
}

bool CandidateTransactionBuilder::ReconcileBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    CandidateScope scope)
{
    auto & selection{ m_selection };
    const bool cooperative{ scope == CandidateScope::CooperativeRescue };
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

    BoundaryComponentReconciliationDiagnostic diagnostic;
    diagnostic.key_list = component.key_list;
    diagnostic.atom_count = FlattenClusterKeyList(component.key_list).size();
    diagnostic.boundary_sample_count = component.boundary_sample_count;
    diagnostic.interface_atom_count = component.interface_atom_index_list.size();
    diagnostic.shape_active_atom_count = component.halo_atom_index_list.size();
    diagnostic.is_rescue_attempt = cooperative;
    if (cooperative)
    {
        diagnostic.accepted_cluster_count = component.key_list.size() - cooperative_key_list.size();
        diagnostic.rescue_candidate_cluster_count = cooperative_key_list.size();
    }
    if (previous_audit_objective != nullptr)
        diagnostic.previous_component_objective = previous_audit_objective->GetTotalObjective();

    auto endpoint_patch{ BuildSelectionPatch(selection, component.key_list) };
    for (const auto & key : cooperative_key_list)
    {
        if (!OverlayFitStatePatch(endpoint_patch, *m_candidate_by_key.at(key).cooperative_patch))
            throw std::logic_error("Boundary rescue candidate patch does not match its component.");
    }
    const CandidateEvaluationOverlay endpoint_overlay{
        inputs.context, inputs.residual_baseline, inputs.previous_state, endpoint_patch
    };
    auto * endpoint_record{ BeginJointCandidateDiagnostic(inputs.options.quiet_mode,
        diagnostic.objective_diagnostic_list, cooperative ? "rescue-endpoint" : "endpoint", 1.0) };
    const auto endpoint_evaluation{ EvaluateCandidate(endpoint_overlay, scope,
        BoundaryCandidateReference{inputs, component, previous_audit_objective, endpoint_record}) };
    const auto endpoint_history_observation{ endpoint_record ? endpoint_record->history_observation : 0 };
    ObservePhaseCandidate(inputs.context,
        cooperative ? "rescue-endpoint" : "boundary-endpoint", endpoint_patch.atom_index_list,
        endpoint_overlay.GetState(), nullptr, 1.0, endpoint_evaluation ? "accepted" : "rejected",
        endpoint_record ? endpoint_record->outcome : "");

    std::optional<ComponentCandidate> accepted;
    if (endpoint_evaluation)
        diagnostic.endpoint_component_objective = endpoint_evaluation->audit_objective.GetTotalObjective();
    if (previous_audit_objective != nullptr)
    {
        accepted = TryBoundaryJointCorrection(inputs, component, *previous_audit_objective,
            endpoint_evaluation ? endpoint_evaluation->audit_objective : *previous_audit_objective,
            endpoint_patch, diagnostic, scope);
    }
    if (!accepted && endpoint_evaluation)
    {
        accepted = ComponentCandidate{ .patch = endpoint_patch,
            .history_observation = endpoint_history_observation };
        diagnostic.accepted_factor = 1.0;
        diagnostic.accepted_source = BoundaryComponentAcceptedSource::Endpoint;
        diagnostic.candidate_component_objective = endpoint_evaluation->audit_objective.GetTotalObjective();
        diagnostic.locally_deteriorated_member_count = endpoint_evaluation->locally_deteriorated_member_count;
        diagnostic.maximum_local_deterioration = endpoint_evaluation->maximum_local_deterioration;
    }
    if (!accepted)
        accepted = TryBacktrackBoundaryComponent(inputs, component, previous_audit_objective,
            endpoint_patch, diagnostic, scope);
    if (accepted)
    {
        ApplyComponentCandidate(inputs, component, endpoint_patch, std::move(*accepted),
            diagnostic.accepted_source, scope);
        if (cooperative)
        {
            if (diagnostic.previous_component_objective && diagnostic.candidate_component_objective)
                diagnostic.component_improvement = *diagnostic.previous_component_objective -
                    *diagnostic.candidate_component_objective;
            diagnostic.rescued_cluster_count = cooperative_key_list.size();
        }
    }
    else if (!cooperative)
        RejectSelectionKeys(inputs, component.key_list, diagnostic.exhausted);
    if (cooperative)
        inputs.performance_counters.RecordBoundaryRescue(accepted.has_value(),
            diagnostic.accepted_source != BoundaryComponentAcceptedSource::Endpoint);
    selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
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
            CandidateScope::CooperativeRescue) || accepted_any;
    }
    return accepted_any;
}

static std::optional<ObjectiveBreakdown> EvaluateFinalSelectionAudit(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const CandidateSelection & selection,
    const std::vector<ClusterKey> & selected_key_list)
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
    return EvaluateCandidate(candidate_overlay,
        GlobalCandidateReference{affected_sample_ref_list, inputs.objective_domain,
            best_audit_objective, &previous_audit_objective, inputs.performance_counters});
}

void CandidateTransactionBuilder::MarkBoundaryDiagnosticRejected(
    const std::vector<ClusterKey> & key_list,
    bool exhausted)
{
    auto & selection{ m_selection };
    auto iter{ std::ranges::find(
        selection.boundary_reconciliation_diagnostic_list | std::views::reverse,
        key_list,
        &BoundaryComponentReconciliationDiagnostic::key_list) };
    if (iter == selection.boundary_reconciliation_diagnostic_list.rend()) return;
    iter->accepted_factor.reset();
    iter->accepted_source = BoundaryComponentAcceptedSource::None;
    iter->exhausted = exhausted;
}

void CandidateTransactionBuilder::AuditAndSalvageFinalSelection(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective)
{
    auto & selection{ m_selection };
    selection.final_audit_objective = EvaluateFinalSelectionAudit(
        inputs,
        previous_audit_objective,
        selection, SelectedKeys());
    if (selection.final_audit_objective.has_value() ||
        SelectedKeys().empty())
    {
        return;
    }

    auto component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            SelectedKeys())
    };
    for (const auto & key : SelectedKeys())
    {
        if (std::ranges::any_of(
                component_list,
                [&](const auto & component)
                {
                    return ContainsClusterKey(component.key_list, key);
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

    for (const auto & rejection_candidate : rejection_candidate_list)
    {
        const auto & key_list{ rejection_candidate.second };
        MarkBoundaryDiagnosticRejected(key_list, false);
        RejectSelectionKeys(
            inputs,
            key_list,
            false);
        selection.final_audit_objective = EvaluateFinalSelectionAudit(
            inputs,
            previous_audit_objective,
            selection, SelectedKeys());
        if (selection.final_audit_objective.has_value()) return;
    }

    const auto remaining_key_list{ SelectedKeys() };
    for (const auto & component : BuildBoundaryReconciliationComponents(
        inputs.context,
        inputs.partition,
        remaining_key_list))
    {
        MarkBoundaryDiagnosticRejected(component.key_list, true);
    }
    RejectSelectionKeys(
        inputs,
        remaining_key_list,
        true);
    selection.final_audit_objective.reset();
}

void CandidateTransactionBuilder::ReconcileSelectedBoundaries(
    const CandidateSelectionInputs & inputs)
{
    auto & selection{ m_selection };
    const auto boundary_component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            SelectedKeys())
    };
    const auto previous_audit_objective{
        EvaluateAuditObjective(inputs.objective_domain, inputs.residual_baseline)
    };
    if (!boundary_component_list.empty())
    {
        const auto boundary_reconciliation_start{ std::chrono::steady_clock::now() };
        for (const auto & component : boundary_component_list)
        {
            ReconcileBoundaryComponent(
                inputs,
                component,
                previous_audit_objective.has_value() ?
                    &*previous_audit_objective : nullptr, CandidateScope::Boundary);
        }
        if (!previous_audit_objective.has_value())
        {
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
        const auto backtracked_component_count{
            std::ranges::count_if(
                selection.boundary_reconciliation_diagnostic_list,
                [](const auto & diagnostic)
                {
                    return diagnostic.accepted_source !=
                            BoundaryComponentAcceptedSource::None &&
                        diagnostic.accepted_factor.has_value() &&
                        *diagnostic.accepted_factor < 1.0;
                })
        };
        const auto rejected_component_count{
            std::ranges::count_if(
                selection.boundary_reconciliation_diagnostic_list,
                [](const auto & diagnostic)
                {
                    return diagnostic.accepted_source ==
                        BoundaryComponentAcceptedSource::None;
                })
        };
        inputs.performance_counters.RecordBoundaryReconciliation(
            selection.boundary_reconciliation_diagnostic_list.size(),
            static_cast<std::size_t>(backtracked_component_count),
            static_cast<std::size_t>(rejected_component_count),
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - boundary_reconciliation_start).count());
    }
    if (previous_audit_objective.has_value() &&
        ReconcileCooperativeComponents(inputs, *previous_audit_objective))
    {
        AuditAndSalvageFinalSelection(
            inputs,
            *previous_audit_objective);
    }
    if (previous_audit_objective.has_value() &&
        selection.final_audit_objective.has_value())
    {
        const auto global_improvement{
            previous_audit_objective->GetTotalObjective() -
            selection.final_audit_objective->GetTotalObjective()
        };
        for (auto & diagnostic : selection.boundary_reconciliation_diagnostic_list)
        {
            if (diagnostic.is_rescue_attempt &&
                diagnostic.accepted_source !=
                    BoundaryComponentAcceptedSource::None)
            {
                diagnostic.global_improvement = global_improvement;
            }
        }
    }
    MaterializeSelection();
}

} // namespace rhbm_gem::core::detail
