#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <algorithm>
#include <cmath>
#include <ranges>

namespace rhbm_gem::core::detail {

bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance,
    ObjectiveProgressGateEvidence * evidence)
{
    ValidateObjectiveTolerance(tolerance);
    if (evidence) *evidence = {};
    if (!std::isfinite(candidate) || !std::isfinite(previous))
    {
        if (evidence) evidence->reason = "objective-nonfinite";
        return false;
    }
    if (evidence) evidence->previous_checked = true;
    if (IsObjectiveDeteriorated(candidate, previous, tolerance))
    {
        if (evidence) evidence->reason = "previous-gate";
        return false;
    }
    if (best != nullptr)
    {
        if (evidence) evidence->best_checked = true;
        if (IsObjectiveDeteriorated(candidate, best->GetTotalObjective(), tolerance))
        {
            if (evidence) evidence->reason = "best-gate";
            return false;
        }
    }
    if (evidence) evidence->reason = "";
    return true;
}

static bool EvaluateLocalObjective(
    const CandidateEvaluationOverlay & candidate_overlay,
    const LocalCandidateReference & reference,
    CandidateDecisionEvidence & evidence)
{
    const auto & key{ reference.key };
    const auto & objective_sample_ref_list{ reference.samples };
    const auto * previous_objective{ reference.objective_reference };
    const auto & domain{ reference.domain };
    auto & performance_counters{ reference.counters };

    const auto unique_sample_count{
        domain.unique_sample_count
    };
    if (performance_counters.AuditEnabled()) performance_counters.RecordObjectiveSampleEvaluation(
        CountObjectiveSamples(objective_sample_ref_list, domain),
        unique_sample_count);
    evidence.candidate_objective =
        EvaluateObjectiveContribution(
            candidate_overlay,
            key,
            objective_sample_ref_list,
            domain);
    evidence.previous_objective.reset();
    if (previous_objective != nullptr)
    {
        evidence.previous_objective = *previous_objective;
    }
    if (!evidence.candidate_objective.has_value() || previous_objective == nullptr)
    {
        return false;
    }
    const auto candidate_objective_value{ evidence.candidate_objective->GetTotalObjective() };
    const auto previous_objective_value{ previous_objective->GetTotalObjective() };
    if (reference.policy == LocalObjectivePolicy::StrictReferenceImprovement)
    {
        const bool accepted{ std::isfinite(previous_objective_value) &&
            IsBetterAuditObjective(candidate_objective_value, previous_objective_value,
                kObjectiveStrictTolerance) };
        // Strict improvement implies non-regression; retain the old rejection classification.
        evidence.rejected_by_previous = !accepted && IsObjectiveDeteriorated(
            candidate_objective_value, previous_objective_value, kObjectiveProgressTolerance);
        return accepted;
    }
    evidence.rejected_by_previous = IsObjectiveDeteriorated(
        candidate_objective_value,
        previous_objective_value,
        kObjectiveProgressTolerance);
    if (reference.member_best)
    {
        const auto patch{ OverlayMemberBest(candidate_overlay.GetState(), *reference.member_best) };
        const CandidateEvaluationOverlay historical{ candidate_overlay.GetContext(), candidate_overlay.GetBaseline(),
            candidate_overlay.GetState().GetBaseState(), patch };
        if (performance_counters.AuditEnabled()) performance_counters.RecordObjectiveSampleEvaluation(
            CountObjectiveSamples(objective_sample_ref_list, domain), unique_sample_count);
        evidence.member_best_objective = EvaluateObjectiveContribution(historical, key, objective_sample_ref_list, domain);
        evidence.rejected_by_member_best = !evidence.member_best_objective || IsObjectiveDeteriorated(
            candidate_objective_value, evidence.member_best_objective->GetTotalObjective(), kObjectiveProgressTolerance);
    }
    return !evidence.rejected_by_previous && !evidence.rejected_by_member_best;
}

LocalCandidateEvaluation EvaluateLocalCandidate(const CandidateEvaluationOverlay & candidate_overlay,
    const LocalCandidateReference & reference)
{
    LocalCandidateEvaluation result;
    result.evidence = reference.evidence;
    result.accepted = EvaluateLocalObjective(candidate_overlay, reference, result.evidence);
    return result;
}

static std::optional<ObjectiveBreakdown> EvaluateBoundaryCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const BoundaryCandidateReference & reference,
    const ObjectiveBreakdown * previous_audit_objective,
    const std::optional<ObjectiveBreakdown> * precomputed_objective,
    JointCandidateObservation * observation)
{
    const auto & component{ reference.component };
    const bool cooperative{ reference.policy == BoundaryAcceptancePolicy::CooperativeRescue };
    const auto observe_member = [&](const ClusterKey & key, bool accepted,
                                    const CandidateDecisionEvidence & evidence)
    {
        if (observation) observation->Member(key, accepted, evidence);
    };
    for (const auto & key : component.key_list)
    {
        CandidateDecisionEvidence evidence;
        const auto & previous_objective{ reference.previous_objective_by_key.at(key) };
        if (!cooperative)
        {
            const auto member{ EvaluateLocalCandidate(candidate_overlay,
                LocalCandidateReference{LocalObjectivePolicy::PreviousNonRegression, key, reference.samples_by_key.at(key),
                    previous_objective ? &*previous_objective : nullptr, reference.domain,
                    evidence, reference.counters, reference.member_best ? &reference.member_best->at(key) : nullptr}) };
            evidence = member.evidence;
            if (!member.accepted)
            {
                observe_member(key, false, evidence);
                return std::nullopt;
            }
        }
        else
        {
            evidence.previous_objective = previous_objective;
            evidence.candidate_objective = EvaluateObjectiveContribution(
                candidate_overlay,
                key,
                reference.samples_by_key.at(key),
                reference.domain);
            if (!evidence.candidate_objective.has_value() || !previous_objective.has_value())
            {
                observe_member(key, false, evidence);
                return std::nullopt;
            }
            const auto candidate_value{
                evidence.candidate_objective->GetTotalObjective()
            };
            const auto previous_value{ previous_objective->GetTotalObjective() };
            if (!std::isfinite(candidate_value) ||
                IsObjectiveDeteriorated(
                    candidate_value,
                    previous_value,
                    kObjectiveProgressTolerance))
            {
                observe_member(key, false, evidence);
                return std::nullopt;
            }
        }
        observe_member(key, true, evidence);
    }
    const auto * best_audit_objective{
        cooperative ? reference.best_audit : nullptr
    };
    // A supplied empty optional is unavailable evidence, not a request to recompute.
    const auto audit_objective{ precomputed_objective ? *precomputed_objective :
        (previous_audit_objective ? EvaluateObjectiveDelta(candidate_overlay,
            component.affected_sample_ref_list, reference.domain,
            *previous_audit_objective, reference.counters) : std::nullopt) };
    if (observation) observation->Global(previous_audit_objective, audit_objective, best_audit_objective);
    if (!audit_objective.has_value() || previous_audit_objective == nullptr)
    {
        if (observation) observation->RejectGlobalObjective();
        return std::nullopt;
    }
    const auto candidate_value{ audit_objective->GetTotalObjective() };
    const auto previous_value{ previous_audit_objective->GetTotalObjective() };
    std::optional<ObjectiveProgressGateEvidence> gate;
    if (observation && observation->IsRecording()) gate.emplace();
    if (cooperative)
    {
        if (!std::isfinite(candidate_value) || !std::isfinite(previous_value))
        {
            if (gate) gate->reason = "nonfinite-objective";
            if (gate) { observation->Gate(*gate); observation->RejectGlobalObjective(); }
            return std::nullopt;
        }
        if (gate) gate->best_checked = best_audit_objective != nullptr;
        if (best_audit_objective && IsObjectiveDeteriorated(candidate_value,
            best_audit_objective->GetTotalObjective(), kObjectiveProgressTolerance))
        {
            if (gate) gate->reason = "best-gate";
            if (gate) { observation->Gate(*gate); observation->RejectGlobalObjective(); }
            return std::nullopt;
        }
        if (gate) gate->previous_checked = true;
        if (!IsBetterAuditObjective(candidate_value, previous_value, kObjectiveStrictTolerance))
        {
            if (gate) { observation->Gate(*gate); observation->RejectStrictImprovement(); }
            return std::nullopt;
        }
    }
    else if (!IsAuditObjectiveAcceptableForProgress(candidate_value, previous_value,
        best_audit_objective, kObjectiveProgressTolerance, gate ? &*gate : nullptr))
    {
        if (gate) { observation->Gate(*gate); observation->RejectGlobalObjective(); }
        return std::nullopt;
    }
    if (gate) observation->Gate(*gate);
    return audit_objective;
}

std::optional<ObjectiveBreakdown> EvaluateBoundaryCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const BoundaryCandidateReference & reference, const ObjectiveBreakdown * previous_audit, JointCandidateObservation * observation)
{
    return EvaluateBoundaryCandidate(candidate_overlay, reference, previous_audit, nullptr, observation);
}

bool EvaluateBoundaryCorrection(const CandidateEvaluationOverlay & candidate,
    const BoundaryCandidateReference & reference, const FitStateView & endpoint,
    const ObjectiveBreakdown & previous_audit, const ObjectiveBreakdown & improvement, JointCandidateObservation * observation)
{
    const auto suspicious_atom_count{ CountSuspiciousPolishAtoms(candidate.GetContext(),
        reference.component.halo_atom_index_list, endpoint, candidate.GetState()) };
    if (suspicious_atom_count != 0)
    {
        if (observation) observation->RejectSuspicious();
        return false;
    }
    const auto raw_objective{ EvaluateObjectiveDelta(candidate, reference.component.affected_sample_ref_list,
        reference.domain, previous_audit, reference.counters) };
    const auto members{ EvaluateBoundaryCandidate(candidate,
        reference, &previous_audit,
        &raw_objective, observation) };
    const bool accepted{ members && IsBetterAuditObjective(members->GetTotalObjective(),
        improvement.GetTotalObjective(), kObjectiveStrictTolerance) };
    if (members && !accepted && observation)
        observation->RejectBoundaryImprovement(improvement, raw_objective);
    return accepted;
}

std::optional<ObjectiveBreakdown> EvaluateGlobalCandidate(const CandidateEvaluationOverlay & candidate,
    const GlobalCandidateReference & reference, SecondStageObservationSession * observation, bool rescue_audit)
{
    if (reference.previous == nullptr)
    {
        if (observation) observation->ObserveGlobalGate(rescue_audit, nullptr, {}, reference.best, false);
        return std::nullopt;
    }
    const auto candidate_objective{
        EvaluateObjectiveDelta(
            candidate,
            reference.samples,
            reference.domain,
            *reference.previous,
            reference.counters)
    };
    std::optional<ObjectiveProgressGateEvidence> gate;
    if (observation && observation->Enabled()) gate.emplace();
    const bool accepted{ candidate_objective.has_value() &&
        IsAuditObjectiveAcceptableForProgress(candidate_objective->GetTotalObjective(),
            reference.previous->GetTotalObjective(), reference.best, kObjectiveProgressTolerance, gate ? &*gate : nullptr) };
    if (gate) observation->ObserveGlobalGate(rescue_audit, reference.previous,
        candidate_objective, reference.best, accepted, *gate);
    if (!accepted) return std::nullopt;
    return candidate_objective;
}

FinalPolishCandidateEvaluation EvaluateFinalPolishCandidate(const CandidateEvaluationOverlay & candidate_overlay,
    const FinalPolishCandidateReference & reference, JointCandidateObservation * observation)
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
    if (has_invalid_model)
    {
        if (observation) observation->RejectInvalidModel();
        return evaluation;
    }

    const auto suspicious_atom_count{
        CountSuspiciousPolishAtoms(
            context,
            component.atom_index_list,
            endpoint_state_view,
            candidate_overlay.GetState())
    };
    evaluation.suspicious_atom_count = suspicious_atom_count;
    if (suspicious_atom_count != 0)
    {
        if (observation) observation->RejectSuspicious();
        return evaluation;
    }
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
        if (observation) observation->RejectPolishImprovement(reference.endpoint_objective, candidate_objective);
        return evaluation;
    }

    if (observation && observation->IsRecording()) { observation->Global(&reference.endpoint_objective, candidate_objective, nullptr); observation->Gate({true, false, ""}); }
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
                    if (observation) observation->RejectMemberSamplesUnavailable(key);
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
                if (!passed && observation) observation->RejectPolishMember(key, base_contribution, candidate_contribution);
                return passed;
            })
    };
    if (!member_guard_passed) return evaluation;

    evaluation.objective = candidate_objective;
    return evaluation;
}

} // namespace rhbm_gem::core::detail
