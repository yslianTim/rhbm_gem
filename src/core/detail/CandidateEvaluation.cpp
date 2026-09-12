#include "core/detail/SecondStageObservation.hpp"
#include "core/detail/CandidateEvaluation.hpp"
#include "core/detail/Diagnosis.hpp"
#include "core/detail/GaussianModelOperations.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <algorithm>
#include <cmath>
#include <ranges>

namespace rhbm_gem::core::detail {

static bool EvaluateLocalObjective(
    const CandidateEvaluationOverlay & candidate_overlay,
    const LocalCandidateReference & reference,
    ObjectiveAttemptDiagnostic & diagnostic)
{
    const auto & key{ reference.key };
    const auto & objective_sample_ref_list{ reference.samples };
    const auto * previous_objective{ reference.objective_reference };
    const auto & domain{ reference.domain };
    auto & performance_counters{ reference.counters };

    const auto unique_sample_count{
        domain.unique_sample_count
    };
    performance_counters.RecordObjectiveSampleEvaluation(
        CountObjectiveSamples(objective_sample_ref_list, domain),
        unique_sample_count);
    const auto domain_iter{ domain.cluster_by_key.find(key) };
    diagnostic.scale.reset();
    if (domain_iter != domain.cluster_by_key.end())
    {
        diagnostic.fit_sample_count = domain_iter->second.fit_sample_ref_list.size();
        diagnostic.tail_sample_count = domain_iter->second.tail_sample_ref_list.size();
        diagnostic.scale = domain_iter->second.scale;
    }
    diagnostic.candidate_objective =
        EvaluateObjectiveContribution(
            candidate_overlay,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.previous_objective.reset();
    if (previous_objective != nullptr)
    {
        diagnostic.previous_objective = *previous_objective;
    }
    if (!diagnostic.candidate_objective.has_value() || previous_objective == nullptr)
    {
        return false;
    }
    const auto candidate_objective_value{ diagnostic.candidate_objective->GetTotalObjective() };
    const auto previous_objective_value{ previous_objective->GetTotalObjective() };
    if (reference.policy == LocalObjectivePolicy::StrictReferenceImprovement)
    {
        const bool accepted{ std::isfinite(previous_objective_value) &&
            IsBetterAuditObjective(candidate_objective_value, previous_objective_value,
                kObjectiveStrictTolerance) };
        // Strict improvement implies non-regression; retain the old rejection classification.
        diagnostic.rejected_by_previous = !accepted && IsObjectiveDeteriorated(
            candidate_objective_value, previous_objective_value, kObjectiveProgressTolerance);
        return accepted;
    }
    diagnostic.rejected_by_previous = IsObjectiveDeteriorated(
        candidate_objective_value,
        previous_objective_value,
        kObjectiveProgressTolerance);
    return !diagnostic.rejected_by_previous;
}

CandidatePreflightEvaluation EvaluateCandidate(const CandidateEvaluationOverlay & candidate,
    const CandidatePreflightReference & reference)
{
    if (!IsTrustRegionStepWithinRadius(reference.step_norm, reference.radius))
        return {CandidateFailureStage::Trust, {}};
    const auto failure{ EvaluateClusterCandidateGuard(candidate.GetContext(), candidate.GetBaseline().model_snapshot,
        reference.key, candidate.GetState(), reference.activity) };
    return {failure ? CandidateFailureStage::Guard : CandidateFailureStage::None, failure};
}

LocalCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay & candidate_overlay,
    const LocalCandidateReference & reference)
{
    LocalCandidateEvaluation result;
    result.diagnostic = reference.diagnostic;
    result.accepted = EvaluateLocalObjective(candidate_overlay, reference, result.diagnostic);
    return result;
}

static std::optional<BoundaryCandidateEvaluation> EvaluateBoundaryCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const BoundaryCandidateReference & reference,
    const std::optional<ObjectiveBreakdown> * precomputed_objective)
{
    const auto & inputs{ reference.inputs };
    const auto & component{ reference.component };
    const auto * previous_audit_objective{ reference.previous_audit };
    const bool cooperative{ reference.policy == BoundaryAcceptancePolicy::CooperativeRescue };
    auto * record{ reference.record };

    ObserveBoundaryHistory(inputs.context, record);
    const auto observe_member = [&](const ClusterKey & key, bool accepted,
                                    const ObjectiveAttemptDiagnostic & diagnostic)
    {
        ObserveBoundaryMemberHistory(candidate_overlay, key,
            inputs.partition.sample_id_list_by_key.at(key), inputs.objective_domain,
            accepted, diagnostic, record);
    };
    BoundaryCandidateEvaluation evaluation;
    for (const auto & key : component.key_list)
    {
        ObjectiveAttemptDiagnostic diagnostic;
        if (record)
        {
            diagnostic.trial_count = record->candidate_number;
            diagnostic.accepted_factor = record->factor;
        }
        const auto & previous_objective{ inputs.previous_objective_by_key.at(key) };
        if (!cooperative)
        {
            const auto member{ EvaluateCandidate(candidate_overlay,
                LocalCandidateReference{LocalObjectivePolicy::PreviousNonRegression, key, inputs.partition.sample_id_list_by_key.at(key),
                    previous_objective ? &*previous_objective : nullptr, inputs.objective_domain,
                    diagnostic, inputs.performance_counters}) };
            diagnostic = member.diagnostic;
            if (!member.accepted)
            {
                RecordJointMemberRejection(record, key, diagnostic.previous_objective,
                    std::nullopt, diagnostic.candidate_objective, false);
                observe_member(key, false, diagnostic);
                return std::nullopt;
            }
        }
        else
        {
            diagnostic.previous_objective = previous_objective;
            diagnostic.candidate_objective = EvaluateObjectiveContribution(
                candidate_overlay,
                key,
                inputs.partition.sample_id_list_by_key.at(key),
                inputs.objective_domain);
            if (!diagnostic.candidate_objective.has_value() || !previous_objective.has_value())
            {
                RecordJointMemberRejection(record, key, previous_objective,
                    std::nullopt, diagnostic.candidate_objective, false);
                observe_member(key, false, diagnostic);
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
                RecordJointMemberRejection(record, key, previous_objective,
                    std::nullopt, diagnostic.candidate_objective, false);
                observe_member(key, false, diagnostic);
                return std::nullopt;
            }
            if (candidate_value > previous_value)
            {
                evaluation.locally_deteriorated_member_count++;
                evaluation.maximum_local_deterioration = std::max(
                    evaluation.maximum_local_deterioration,
                    candidate_value - previous_value);
            }
        }
        observe_member(key, true, diagnostic);
    }
    const auto * best_audit_objective{
        cooperative && inputs.best_audit_state.has_value() ?
            &inputs.best_audit_state->objective : nullptr
    };
    // A supplied empty optional is unavailable evidence, not a request to recompute.
    const auto audit_objective{ precomputed_objective ? *precomputed_objective :
        (previous_audit_objective ? EvaluateObjectiveDelta(candidate_overlay,
            component.affected_sample_ref_list, inputs.objective_domain,
            *previous_audit_objective, inputs.performance_counters) : std::nullopt) };
    if (!audit_objective.has_value() || previous_audit_objective == nullptr)
    {
        if (record) record->outcome = "members-passed-global-objective-rejected-or-unavailable";
        return std::nullopt;
    }
    const auto candidate_value{ audit_objective->GetTotalObjective() };
    const auto previous_value{ previous_audit_objective->GetTotalObjective() };
    if (cooperative)
    {
        if (!std::isfinite(candidate_value) || !std::isfinite(previous_value) ||
            (best_audit_objective && IsObjectiveDeteriorated(candidate_value,
                best_audit_objective->GetTotalObjective(), kObjectiveProgressTolerance)))
        {
            if (record) record->outcome = "members-passed-global-objective-rejected-or-unavailable";
            return std::nullopt;
        }
        if (!IsBetterAuditObjective(candidate_value, previous_value, kObjectiveStrictTolerance))
        {
            // The previous bound is only needed to preserve the original failure reason.
            if (record) record->outcome = IsObjectiveDeteriorated(
                candidate_value, previous_value, kObjectiveProgressTolerance) ?
                "members-passed-global-objective-rejected-or-unavailable" :
                "members-passed-strict-improvement-failed";
            return std::nullopt;
        }
    }
    else if (!IsAuditObjectiveAcceptableForProgress(candidate_value, previous_value,
        best_audit_objective, kObjectiveProgressTolerance))
    {
        if (record) record->outcome = "members-passed-global-objective-rejected-or-unavailable";
        return std::nullopt;
    }
    evaluation.audit_objective = *audit_objective;
    return evaluation;
}

std::optional<BoundaryCandidateEvaluation> EvaluateCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const BoundaryCandidateReference & reference)
{
    return EvaluateBoundaryCandidate(candidate_overlay, reference, nullptr);
}

BoundaryCorrectionEvaluation EvaluateCandidate(const CandidateEvaluationOverlay & candidate,
    const BoundaryCorrectionReference & reference)
{
    BoundaryCorrectionEvaluation result;
    const auto & inputs{ reference.inputs };
    result.suspicious_atom_count = CountSuspiciousPolishAtoms(inputs.context,
        reference.component.halo_atom_index_list, reference.endpoint, candidate.GetState());
    if (result.suspicious_atom_count != 0) return result;
    result.raw_objective = EvaluateObjectiveDelta(candidate, reference.component.affected_sample_ref_list,
        inputs.objective_domain, reference.previous_audit, inputs.performance_counters);
    result.record = BeginJointCandidateDiagnostic(inputs.options.quiet_mode, reference.records,
        BoundaryDiagnosticName(reference.policy, BoundaryObservationStage::Correction), reference.damping);
    result.members = EvaluateBoundaryCandidate(candidate,
        BoundaryCandidateReference{reference.policy, inputs, reference.component, &reference.previous_audit, result.record},
        &result.raw_objective);
    result.accepted = result.members && IsBetterAuditObjective(result.members->audit_objective.GetTotalObjective(),
        reference.improvement.GetTotalObjective(), kObjectiveStrictTolerance);
    return result;
}

std::optional<ObjectiveBreakdown> EvaluateCandidate(const CandidateEvaluationOverlay & candidate,
    const GlobalCandidateReference & reference)
{
    return EvaluateCombinedObjective(candidate, reference.samples, reference.domain,
        reference.best, reference.previous, reference.counters);
}

FinalPolishCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay & candidate_overlay,
    const FinalPolishCandidateReference & reference)
{
    FinalPolishCandidateEvaluation evaluation;
    const auto & context{ candidate_overlay.GetContext() };
    const auto & base_baseline{ candidate_overlay.GetBaseline() };
    const auto & component{ reference.component };
    const auto & partition{ reference.partition };
    const auto & objective_domain{ reference.domain };
    const auto & endpoint_state_view{ reference.endpoint };
    auto & performance_counters{ reference.counters };
    const auto has_invalid_model{
        std::ranges::any_of(
            component.atom_index_list,
            [&](const auto atom_index)
            {
                return !IsValidSecondStageGaussianModel(
                    candidate_overlay.GetState().GetModel(atom_index));
            })
    };
    if (has_invalid_model) return evaluation;

    const auto suspicious_atom_count{
        CountSuspiciousPolishAtoms(
            context,
            component.atom_index_list,
            endpoint_state_view,
            candidate_overlay.GetState())
    };
    evaluation.suspicious_atom_count = suspicious_atom_count;
    if (suspicious_atom_count != 0) return evaluation;

    auto * record{ BeginJointCandidateDiagnostic(reference.quiet,
        reference.records, "final-polish", reference.damping, reference.round) };
    const auto candidate_objective{
        EvaluateObjectiveDelta(
            candidate_overlay,
            component.affected_sample_ref_list,
            objective_domain,
            reference.base_objective,
            performance_counters)
    };
    if (!candidate_objective.has_value() ||
        !IsBetterAuditObjective(
            candidate_objective->GetTotalObjective(),
            reference.endpoint_objective.GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        if (record)
        {
            record->previous = reference.endpoint_objective;
            record->candidate = candidate_objective;
            record->outcome = "members-not-evaluated-global-improvement-failed-or-unavailable";
        }
        return evaluation;
    }

    const auto member_guard_passed{
        std::ranges::all_of(
            component.key_list,
            [&](const auto & key)
            {
                const auto sample_iter{
                    partition.sample_id_list_by_key.find(key)
                };
                if (sample_iter ==
                    partition.sample_id_list_by_key.end())
                {
                    RecordJointMemberRejection(record, key, std::nullopt, std::nullopt, std::nullopt, false);
                    if (record) record->outcome = "member-samples-unavailable";
                    return false;
                }
                auto owned_sample_ref_list{ sample_iter->second };
                owned_sample_ref_list.erase(
                    std::remove_if(
                        owned_sample_ref_list.begin(),
                        owned_sample_ref_list.end(),
                        [&](const auto & sample_ref)
                        {
                            return sample_ref.atom_index >=
                                    objective_domain.owner_key_by_atom_index.size() ||
                                objective_domain.owner_key_by_atom_index.at(
                                    sample_ref.atom_index).empty();
                        }),
                    owned_sample_ref_list.end());
                const auto base_contribution{
                    EvaluateObjectiveContribution(
                        base_baseline,
                        key,
                        owned_sample_ref_list,
                        objective_domain)
                };
                const auto candidate_contribution{
                    EvaluateObjectiveContribution(
                        candidate_overlay,
                        key,
                        owned_sample_ref_list,
                        objective_domain)
                };
                const bool passed{ base_contribution.has_value() &&
                    candidate_contribution.has_value() &&
                    !IsObjectiveDeteriorated(
                        candidate_contribution->GetTotalObjective(),
                        base_contribution->GetTotalObjective(),
                        kObjectiveProgressTolerance) };
                if (!passed) RecordJointMemberRejection(record, key, base_contribution,
                    std::nullopt, candidate_contribution, false);
                return passed;
            })
    };
    if (!member_guard_passed) return evaluation;

    evaluation.objective = candidate_objective;
    return evaluation;
}

} // namespace rhbm_gem::core::detail
