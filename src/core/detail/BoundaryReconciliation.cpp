#include "core/detail/BoundaryReconciliation.hpp"

#include "core/detail/Diagnosis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
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

static void RejectSelectionKeys(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & key_list,
    bool exhausted,
    CandidateSelection & selection)
{
    for (const auto & key : key_list)
    {
        if (!ContainsClusterKey(selection.accepted_key_list, key)) continue;
        for (const auto atom_index : key)
        {
            selection.assembled_state.at(atom_index) = inputs.previous_state.at(atom_index);
            selection.assembled_polish_provenance.at(atom_index) = inputs.previous_polish_provenance.at(atom_index);
        }
        std::erase(selection.accepted_key_list, key);
        std::erase(selection.grow_trust_region_key_list, key);
        selection.rejected_key_list.emplace_back(key);
        if (exhausted) selection.exhausted_key_list.emplace_back(key);
        selection.cluster_objective_state.at(key) = inputs.cluster_objective_state.at(key);

        const auto diagnostic_iter{
            std::ranges::find(
                selection.accepted_cluster_diagnostic_list,
                key,
                &ClusterCandidateDiagnostic::key)
        };
        if (diagnostic_iter != selection.accepted_cluster_diagnostic_list.end())
        {
            selection.rejected_cluster_diagnostic_list.emplace_back(std::move(*diagnostic_iter));
            selection.accepted_cluster_diagnostic_list.erase(diagnostic_iter);
        }
    }
    std::ranges::sort(selection.accepted_key_list);
    std::ranges::sort(selection.rejected_key_list);
    selection.rejected_key_list.erase(
        std::ranges::unique(selection.rejected_key_list).begin(),
        selection.rejected_key_list.end());
    std::ranges::sort(selection.exhausted_key_list);
    selection.exhausted_key_list.erase(
        std::ranges::unique(selection.exhausted_key_list).begin(),
        selection.exhausted_key_list.end());
}

struct BoundaryCandidateEvaluation
{
    ClusterObjectiveStateMap objective_state_by_key{};
    ObjectiveBreakdown audit_objective{};
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
};

static FitStatePatch BuildSelectionPatch(
    const CandidateSelection & selection,
    const std::vector<ClusterKey> & key_list)
{
    return FitStatePatch::FromState(
        selection.assembled_state,
        FlattenClusterKeyList(key_list));
}

static std::optional<BoundaryCandidateEvaluation>
EvaluateBoundaryComponentCandidate(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const CandidateEvaluationOverlay & candidate_overlay,
    const ObjectiveBreakdown * previous_audit_objective,
    bool cooperative = false)
{
    BoundaryCandidateEvaluation evaluation;
    for (const auto & key : component.key_list)
    {
        auto objective_state{ inputs.cluster_objective_state.at(key) };
        ObjectiveAttemptDiagnostic diagnostic;
        const auto & previous_objective{ inputs.previous_objective_by_key.at(key) };
        if (!cooperative)
        {
            if (!TryCommitClusterCandidate(
                    candidate_overlay,
                    key,
                    inputs.partition.sample_id_list_by_key.at(key),
                    previous_objective.has_value() ? &*previous_objective : nullptr,
                    false,
                    inputs.objective_domain,
                    objective_state,
                    diagnostic,
                    inputs.performance_counters))
            {
                return std::nullopt;
            }
        }
        else
        {
            diagnostic.candidate_objective = EvaluateObjectiveContribution(
                candidate_overlay,
                key,
                inputs.partition.sample_id_list_by_key.at(key),
                inputs.objective_domain);
            if (!diagnostic.candidate_objective.has_value() || !previous_objective.has_value())
            {
                return std::nullopt;
            }
            const auto candidate_value{
                diagnostic.candidate_objective->GetTotalObjective()
            };
            const auto previous_value{ previous_objective->GetTotalObjective() };
            if (!std::isfinite(candidate_value) ||
                IsObjectiveDeteriorated(
                    candidate_value,
                    previous_value,
                    kObjectiveProgressTolerance))
            {
                return std::nullopt;
            }
            if (candidate_value > previous_value)
            {
                evaluation.locally_deteriorated_member_count++;
                evaluation.maximum_local_deterioration = std::max(
                    evaluation.maximum_local_deterioration,
                    candidate_value - previous_value);
            }
            const auto improves_member{
                IsBetterAuditObjective(
                    candidate_value,
                    previous_value,
                    kObjectiveStrictTolerance)
            };
            if (improves_member &&
                (!objective_state.best_objective.has_value() ||
                    IsBetterAuditObjective(
                        candidate_value,
                        objective_state.best_objective->GetTotalObjective(),
                        kObjectiveStrictTolerance)))
            {
                objective_state.best_objective = diagnostic.candidate_objective;
                objective_state.best_maximum_transformed_change =
                    std::ranges::max(
                        SummarizeTransformedChanges(
                            candidate_overlay.GetState(),
                            candidate_overlay.GetBaseline().model_snapshot.node,
                            key).maximum_list);
            }
        }
        evaluation.objective_state_by_key.emplace(key, std::move(objective_state));
    }
    const auto * best_audit_objective{
        cooperative && inputs.best_audit_state.has_value() ?
            &inputs.best_audit_state->objective : nullptr
    };
    const auto audit_objective{ EvaluateCombinedObjective(
        candidate_overlay,
        component.affected_sample_ref_list,
        inputs.objective_domain,
        best_audit_objective,
        previous_audit_objective,
        inputs.performance_counters) };
    if (!audit_objective.has_value()) return std::nullopt;
    if (cooperative &&
        !IsBetterAuditObjective(
            audit_objective->GetTotalObjective(),
            previous_audit_objective->GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        return std::nullopt;
    }
    evaluation.audit_objective = *audit_objective;
    return evaluation;
}

static void CommitBoundaryObjectiveState(
    const BoundaryCandidateEvaluation & evaluation,
    ClusterObjectiveStateMap & working_objective_state)
{
    for (const auto & [key, objective_state] : evaluation.objective_state_by_key)
    {
        working_objective_state.at(key) = objective_state;
    }
}

static void RemoveTrustGrowthForKeys(
    const std::vector<ClusterKey> & key_list,
    CandidateSelection & selection)
{
    for (const auto & key : key_list)
    {
        std::erase(selection.grow_trust_region_key_list, key);
    }
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

static bool TryBoundaryJointCorrection(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown & previous_audit_objective,
    const ObjectiveBreakdown & improvement_reference_objective,
    const FitStatePatch & endpoint_patch,
    CandidateSelection & selection,
    BoundaryComponentReconciliationDiagnostic & diagnostic)
{
    if (component.halo_atom_index_list.empty()) return false;

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
        return false;
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
        return false;
    }

    auto corrected_component_patch{ endpoint_patch };
    if (!OverlayFitStatePatch(
            corrected_component_patch,
            *correction_result.patch))
    {
        record_performance(false);
        return false;
    }
    const FitStateView corrected_state_view{
        inputs.previous_state,
        corrected_component_patch
    };
    diagnostic.suspicious_candidate_atom_count =
        CountSuspiciousPolishAtoms(
            inputs.context,
            inputs.options,
            component.halo_atom_index_list,
            endpoint_state_view,
            corrected_state_view);
    if (diagnostic.suspicious_candidate_atom_count != 0)
    {
        record_performance(false);
        return false;
    }
    const CandidateEvaluationOverlay corrected_overlay{
        inputs.context,
        inputs.residual_baseline,
        corrected_state_view
    };
    const auto raw_candidate_objective{
        EvaluateObjectiveDelta(
            corrected_overlay,
            component.affected_sample_ref_list,
            inputs.objective_domain,
            previous_audit_objective,
            inputs.performance_counters)
    };
    if (raw_candidate_objective.has_value())
    {
        diagnostic.joint_candidate_component_objective = raw_candidate_objective->GetTotalObjective();
    }
    const auto candidate_evaluation{
        EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            corrected_overlay,
            &previous_audit_objective,
            diagnostic.is_rescue_attempt)
    };
    const auto is_strict_improvement{
        candidate_evaluation.has_value() &&
        IsBetterAuditObjective(
            candidate_evaluation->audit_objective.GetTotalObjective(),
            improvement_reference_objective.GetTotalObjective(),
            kObjectiveStrictTolerance)
    };
    if (!is_strict_improvement)
    {
        record_performance(false);
        return false;
    }

    corrected_component_patch.ApplyTo(selection.assembled_state);
    for (const auto atom_index : component.halo_atom_index_list)
    {
        const auto change{
            CalculateTransformedChange(
                selection.assembled_state.at(atom_index).mdpde.GetModel(),
                endpoint_state_view.GetModel(atom_index))
        };
        if (IsTransformedChangeMaterial(change, kTransformedChangeTolerance))
        {
            selection.assembled_polish_provenance.at(atom_index) = 1;
        }
    }
    CommitBoundaryObjectiveState(*candidate_evaluation, selection.cluster_objective_state);
    RemoveTrustGrowthForKeys(component.key_list, selection);
    diagnostic.accepted_source = BoundaryComponentAcceptedSource::JointCorrection;
    diagnostic.candidate_component_objective = candidate_evaluation->audit_objective.GetTotalObjective();
    diagnostic.locally_deteriorated_member_count = candidate_evaluation->locally_deteriorated_member_count;
    diagnostic.maximum_local_deterioration = candidate_evaluation->maximum_local_deterioration;
    record_performance(true);
    return true;
}

static bool TryBacktrackBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    const FitStatePatch & endpoint_patch,
    CandidateSelection & selection,
    BoundaryComponentReconciliationDiagnostic & diagnostic)
{
    BacktrackingWorkspace backtracking_workspace{
        inputs.previous_state,
        endpoint_patch,
        kTransformedChangeTolerance
    };
    BacktrackingStep step;
    std::optional<BoundaryCandidateEvaluation> accepted_evaluation;
    for (step = backtracking_workspace.BuildNextCandidate();
        step.status == BacktrackingStepStatus::CandidateReady;
        step = backtracking_workspace.BuildNextCandidate())
    {
        diagnostic.trial_count = step.trial_number;
        const FitStateView candidate_state_view{
            inputs.previous_state,
            backtracking_workspace.GetCandidatePatch()
        };
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            candidate_state_view
        };
        accepted_evaluation = EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            candidate_overlay,
            previous_audit_objective,
            diagnostic.is_rescue_attempt);
        if (accepted_evaluation.has_value()) break;
    }
    if (!accepted_evaluation.has_value())
    {
        diagnostic.exhausted = step.status == BacktrackingStepStatus::Exhausted;
        return false;
    }

    backtracking_workspace.GetCandidatePatch().ApplyTo(selection.assembled_state);
    const auto reconciled_provenance{
        backtracking_workspace.BuildCandidatePolishProvenance(
            inputs.previous_polish_provenance,
            selection.assembled_polish_provenance)
    };
    for (const auto atom_index : FlattenClusterKeyList(component.key_list))
    {
        selection.assembled_polish_provenance.at(atom_index) = reconciled_provenance.at(atom_index);
    }
    CommitBoundaryObjectiveState(*accepted_evaluation, selection.cluster_objective_state);
    diagnostic.accepted_factor = step.factor;
    diagnostic.accepted_source = BoundaryComponentAcceptedSource::Backtracking;
    diagnostic.candidate_component_objective = accepted_evaluation->audit_objective.GetTotalObjective();
    diagnostic.locally_deteriorated_member_count = accepted_evaluation->locally_deteriorated_member_count;
    diagnostic.maximum_local_deterioration = accepted_evaluation->maximum_local_deterioration;
    return true;
}

static void ReconcileBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    CandidateSelection & selection)
{
    BoundaryComponentReconciliationDiagnostic diagnostic;
    diagnostic.key_list = component.key_list;
    diagnostic.atom_count = FlattenClusterKeyList(component.key_list).size();
    diagnostic.boundary_sample_count = component.boundary_sample_count;
    diagnostic.interface_atom_count = component.interface_atom_index_list.size();
    diagnostic.shape_active_atom_count = component.halo_atom_index_list.size();
    if (previous_audit_objective != nullptr)
    {
        diagnostic.previous_component_objective = previous_audit_objective->GetTotalObjective();
    }
    const auto endpoint_patch{ BuildSelectionPatch(selection, component.key_list) };
    const FitStateView endpoint_state_view{
        inputs.previous_state,
        endpoint_patch
    };
    const CandidateEvaluationOverlay endpoint_overlay{
        inputs.context,
        inputs.residual_baseline,
        endpoint_state_view
    };
    const auto endpoint_evaluation{
        EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            endpoint_overlay,
            previous_audit_objective)
    };
    if (endpoint_evaluation.has_value())
    {
        diagnostic.endpoint_component_objective = endpoint_evaluation->audit_objective.GetTotalObjective();
        if (previous_audit_objective != nullptr &&
            TryBoundaryJointCorrection(
                inputs,
                component,
                *previous_audit_objective,
                endpoint_evaluation->audit_objective,
                endpoint_patch,
                selection,
                diagnostic))
        {
            selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
            return;
        }
        CommitBoundaryObjectiveState(*endpoint_evaluation, selection.cluster_objective_state);
        diagnostic.accepted_factor = 1.0;
        diagnostic.accepted_source = BoundaryComponentAcceptedSource::Endpoint;
        diagnostic.candidate_component_objective = endpoint_evaluation->audit_objective.GetTotalObjective();
        selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
        return;
    }

    if (previous_audit_objective != nullptr &&
        TryBoundaryJointCorrection(
            inputs,
            component,
            *previous_audit_objective,
            *previous_audit_objective,
            endpoint_patch,
            selection,
            diagnostic))
    {
        selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
        return;
    }

    if (!TryBacktrackBoundaryComponent(
            inputs,
            component,
            previous_audit_objective,
            endpoint_patch,
            selection,
            diagnostic))
    {
        RejectSelectionKeys(
            inputs,
            component.key_list,
            diagnostic.exhausted,
            selection);
        selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
        return;
    }

    RemoveTrustGrowthForKeys(component.key_list, selection);
    selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
}

static void PromoteBoundaryRescueKeys(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & rescue_key_list,
    const FitStatePatch & endpoint_patch,
    BoundaryComponentAcceptedSource accepted_source,
    CandidateSelection & selection)
{
    const FitStateView endpoint_state{ inputs.previous_state, endpoint_patch };
    for (const auto & key : rescue_key_list)
    {
        std::erase(selection.rejected_key_list, key);
        std::erase(selection.exhausted_key_list, key);
        std::erase(selection.grow_trust_region_key_list, key);
        if (!ContainsClusterKey(selection.accepted_key_list, key))
        {
            selection.accepted_key_list.emplace_back(key);
        }
        for (const auto atom_index : key)
        {
            const auto changed_from_previous{
                IsTransformedChangeMaterial(
                    CalculateTransformedChange(
                        selection.assembled_state.at(atom_index).mdpde.GetModel(),
                        inputs.previous_state.at(atom_index).mdpde.GetModel()),
                    kTransformedChangeTolerance)
            };
            if (!changed_from_previous)
            {
                selection.assembled_polish_provenance.at(atom_index) = inputs.previous_polish_provenance.at(atom_index);
                continue;
            }
            const auto correction_changed_endpoint{
                accepted_source == BoundaryComponentAcceptedSource::JointCorrection &&
                IsTransformedChangeMaterial(
                    CalculateTransformedChange(
                        selection.assembled_state.at(atom_index).mdpde.GetModel(),
                        endpoint_state.GetModel(atom_index)),
                    kTransformedChangeTolerance)
            };
            selection.assembled_polish_provenance.at(atom_index) = correction_changed_endpoint ? 1 : 0;
        }

        const auto diagnostic_iter{
            std::ranges::find(
                selection.rejected_cluster_diagnostic_list,
                key,
                &ClusterCandidateDiagnostic::key)
        };
        if (diagnostic_iter != selection.rejected_cluster_diagnostic_list.end())
        {
            diagnostic_iter->boundary_rescued = true;
            selection.accepted_cluster_diagnostic_list.emplace_back(std::move(*diagnostic_iter));
            selection.rejected_cluster_diagnostic_list.erase(diagnostic_iter);
        }
    }
    std::ranges::sort(selection.accepted_key_list);
    std::ranges::sort(selection.rejected_key_list);
    std::ranges::sort(
        selection.accepted_cluster_diagnostic_list,
        {},
        &ClusterCandidateDiagnostic::key);
}

static bool TryRescueBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown & previous_audit_objective,
    const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key,
    CandidateSelection & selection)
{
    std::vector<ClusterKey> rescue_key_list;
    for (const auto & key : component.key_list)
    {
        if (ContainsClusterKey(selection.rejected_key_list, key) &&
            rescue_patch_by_key.contains(key))
        {
            rescue_key_list.emplace_back(key);
        }
    }
    if (rescue_key_list.empty()) return false;

    BoundaryComponentReconciliationDiagnostic diagnostic;
    diagnostic.key_list = component.key_list;
    diagnostic.atom_count = FlattenClusterKeyList(component.key_list).size();
    diagnostic.boundary_sample_count = component.boundary_sample_count;
    diagnostic.interface_atom_count = component.interface_atom_index_list.size();
    diagnostic.shape_active_atom_count = component.halo_atom_index_list.size();
    diagnostic.accepted_cluster_count = component.key_list.size() - rescue_key_list.size();
    diagnostic.rescue_candidate_cluster_count = rescue_key_list.size();
    diagnostic.is_rescue_attempt = true;
    diagnostic.previous_component_objective = previous_audit_objective.GetTotalObjective();

    auto endpoint_patch{ BuildSelectionPatch(selection, component.key_list) };
    for (const auto & key : rescue_key_list)
    {
        if (!OverlayFitStatePatch(endpoint_patch, rescue_patch_by_key.at(key)))
        {
            throw std::logic_error(
                "Boundary rescue candidate patch does not match its component.");
        }
    }
    const FitStateView endpoint_state_view{ inputs.previous_state, endpoint_patch };
    const CandidateEvaluationOverlay endpoint_overlay{
        inputs.context,
        inputs.residual_baseline,
        endpoint_state_view
    };
    const auto endpoint_evaluation{
        EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            endpoint_overlay,
            &previous_audit_objective,
            true)
    };
    if (endpoint_evaluation.has_value())
    {
        diagnostic.endpoint_component_objective =
            endpoint_evaluation->audit_objective.GetTotalObjective();
        if (!TryBoundaryJointCorrection(
                inputs,
                component,
                previous_audit_objective,
                endpoint_evaluation->audit_objective,
                endpoint_patch,
                selection,
                diagnostic))
        {
            endpoint_patch.ApplyTo(selection.assembled_state);
            CommitBoundaryObjectiveState(*endpoint_evaluation, selection.cluster_objective_state);
            diagnostic.accepted_factor = 1.0;
            diagnostic.accepted_source = BoundaryComponentAcceptedSource::Endpoint;
            diagnostic.candidate_component_objective =
                endpoint_evaluation->audit_objective.GetTotalObjective();
            diagnostic.locally_deteriorated_member_count =
                endpoint_evaluation->locally_deteriorated_member_count;
            diagnostic.maximum_local_deterioration =
                endpoint_evaluation->maximum_local_deterioration;
        }
    }
    else if (!TryBoundaryJointCorrection(
            inputs,
            component,
            previous_audit_objective,
            previous_audit_objective,
            endpoint_patch,
            selection,
            diagnostic))
    {
        TryBacktrackBoundaryComponent(
            inputs,
            component,
            &previous_audit_objective,
            endpoint_patch,
            selection,
            diagnostic);
    }

    if (diagnostic.accepted_source != BoundaryComponentAcceptedSource::None)
    {
        if (diagnostic.previous_component_objective.has_value() &&
            diagnostic.candidate_component_objective.has_value())
        {
            diagnostic.component_improvement =
                *diagnostic.previous_component_objective -
                *diagnostic.candidate_component_objective;
        }
        diagnostic.rescued_cluster_count = rescue_key_list.size();
        PromoteBoundaryRescueKeys(
            inputs,
            rescue_key_list,
            endpoint_patch,
            diagnostic.accepted_source,
            selection);
    }
    inputs.performance_counters.RecordBoundaryRescue(
        diagnostic.accepted_source != BoundaryComponentAcceptedSource::None,
        diagnostic.accepted_source != BoundaryComponentAcceptedSource::Endpoint);
    selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
    return selection.boundary_reconciliation_diagnostic_list.back()
            .accepted_source != BoundaryComponentAcceptedSource::None;
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

static bool RescueRejectedBoundaryClusters(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key,
    CandidateSelection & selection)
{
    std::vector<ClusterKey> eligible_key_list;
    for (const auto & key : selection.accepted_key_list)
    {
        eligible_key_list.emplace_back(key);
    }
    for (const auto & key : selection.rejected_key_list)
    {
        if (rescue_patch_by_key.contains(key))
        {
            eligible_key_list.emplace_back(key);
        }
    }
    std::ranges::sort(eligible_key_list);
    eligible_key_list.erase(
        std::ranges::unique(eligible_key_list).begin(),
        eligible_key_list.end());

    bool rescued_any{ false };
    for (const auto & component : BuildExpandedBoundaryReconciliationComponents(
        inputs,
        eligible_key_list))
    {
        rescued_any = TryRescueBoundaryComponent(
            inputs,
            component,
            previous_audit_objective,
            rescue_patch_by_key,
            selection) || rescued_any;
    }
    return rescued_any;
}

static std::optional<ObjectiveBreakdown> EvaluateFinalSelectionAudit(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const CandidateSelection & selection)
{
    if (selection.accepted_key_list.empty()) return std::nullopt;
    const auto candidate_patch{
        BuildSelectionPatch(selection, selection.accepted_key_list)
    };
    const FitStateView candidate_state_view{
        inputs.previous_state,
        candidate_patch
    };
    const CandidateEvaluationOverlay candidate_overlay{
        inputs.context,
        inputs.residual_baseline,
        candidate_state_view
    };
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(inputs.partition, selection.accepted_key_list)
    };
    const auto * best_audit_objective{
        inputs.best_audit_state.has_value() ? &inputs.best_audit_state->objective : nullptr
    };
    return EvaluateCombinedObjective(
        candidate_overlay,
        affected_sample_ref_list,
        inputs.objective_domain,
        best_audit_objective,
        &previous_audit_objective,
        inputs.performance_counters);
}

static void MarkBoundaryDiagnosticRejected(
    const std::vector<ClusterKey> & key_list,
    bool exhausted,
    CandidateSelection & selection)
{
    auto iter{ std::ranges::find(
        selection.boundary_reconciliation_diagnostic_list | std::views::reverse,
        key_list,
        &BoundaryComponentReconciliationDiagnostic::key_list) };
    if (iter == selection.boundary_reconciliation_diagnostic_list.rend()) return;
    iter->accepted_factor.reset();
    iter->accepted_source = BoundaryComponentAcceptedSource::None;
    iter->exhausted = exhausted;
}

static void AuditAndSalvageFinalSelection(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    CandidateSelection & selection)
{
    selection.final_audit_objective = EvaluateFinalSelectionAudit(
        inputs,
        previous_audit_objective,
        selection);
    if (selection.final_audit_objective.has_value() ||
        selection.accepted_key_list.empty())
    {
        return;
    }

    auto component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            selection.accepted_key_list)
    };
    for (const auto & key : selection.accepted_key_list)
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
        const FitStateView candidate_state_view{
            inputs.previous_state,
            candidate_patch
        };
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            candidate_state_view
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
        MarkBoundaryDiagnosticRejected(key_list, false, selection);
        RejectSelectionKeys(
            inputs,
            key_list,
            false,
            selection);
        selection.final_audit_objective = EvaluateFinalSelectionAudit(
            inputs,
            previous_audit_objective,
            selection);
        if (selection.final_audit_objective.has_value()) return;
    }

    const auto remaining_key_list{ selection.accepted_key_list };
    for (const auto & component : BuildBoundaryReconciliationComponents(
        inputs.context,
        inputs.partition,
        remaining_key_list))
    {
        MarkBoundaryDiagnosticRejected(component.key_list, true, selection);
    }
    RejectSelectionKeys(
        inputs,
        remaining_key_list,
        true,
        selection);
    selection.final_audit_objective.reset();
}

void ReauditFallbackSelection(const CandidateSelectionInputs & inputs, CandidateSelection & selection)
{
    selection.cluster_objective_state = inputs.cluster_objective_state;
    const auto accepted_keys{ selection.accepted_key_list };
    for (const auto & key : accepted_keys)
    {
        auto patch{ FitStatePatch::FromState(selection.assembled_state, key) };
        const FitStateView view{ inputs.previous_state, patch };
        std::vector<GaussianModel3D> previous_models;
        std::vector<GaussianModel3D> candidate_models;
        for (const auto node : key)
        {
            previous_models.emplace_back(inputs.previous_state.at(node).mdpde.GetModel());
            candidate_models.emplace_back(view.GetModel(node));
        }
        const auto norm{ CalculateModelTrustRegionStepNorm(previous_models, candidate_models) };
        bool safe{ norm.has_value() && IsTrustRegionStepWithinRadius(*norm, inputs.trust_region_state.GetRadius(key)) &&
            !EvaluateClusterCandidateGuard(inputs.context, inputs.options, inputs.residual_baseline.model_snapshot,
                key, view, selection.block_activity).has_value() };
        if (safe)
        {
            ObjectiveAttemptDiagnostic diagnostic;
            const auto & previous{ inputs.previous_objective_by_key.at(key) };
            safe = TryCommitClusterCandidate(CandidateEvaluationOverlay{ inputs.context, inputs.residual_baseline, view },
                key, inputs.partition.sample_id_list_by_key.at(key), previous.has_value() ? &*previous : nullptr,
                false, inputs.objective_domain, selection.cluster_objective_state.at(key), diagnostic, inputs.performance_counters);
        }
        if (safe) patch.ApplyTo(selection.assembled_state);
        else RejectSelectionKeys(inputs, { key }, false, selection);
    }
    const auto previous_audit{ EvaluateAuditObjective(inputs.objective_domain, inputs.residual_baseline) };
    if (previous_audit.has_value())
        AuditAndSalvageFinalSelection(inputs, *previous_audit, selection);
    else
    {
        const auto remaining_keys{ selection.accepted_key_list };
        RejectSelectionKeys(inputs, remaining_keys, false, selection);
    }
}

void ReconcileSelectedBoundaries(
    const CandidateSelectionInputs & inputs,
    const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key,
    CandidateSelection & selection)
{
    const auto boundary_component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            selection.accepted_key_list)
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
                    &*previous_audit_objective : nullptr,
                    selection);
        }
        if (!previous_audit_objective.has_value())
        {
            const auto remaining_key_list{ selection.accepted_key_list };
            RejectSelectionKeys(
                inputs,
                remaining_key_list,
                true,
                selection);
        }
        else
        {
            AuditAndSalvageFinalSelection(
                inputs,
                *previous_audit_objective,
                selection);
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
        RescueRejectedBoundaryClusters(
            inputs,
            *previous_audit_objective,
            rescue_patch_by_key,
            selection))
    {
        AuditAndSalvageFinalSelection(
            inputs,
            *previous_audit_objective,
            selection);
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
}

} // namespace rhbm_gem::core::detail
