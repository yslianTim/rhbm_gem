#include "core/detail/CandidateTransaction.hpp"
#include "core/detail/IterationProcess.hpp"
#include "core/detail/Diagnosis.hpp"
#include "core/detail/PhaseAudit.hpp"

#include <algorithm>

namespace rhbm_gem::core::detail {

CandidateTransaction CandidateTransactionBuilder::Finish(const CandidateSelectionInputs & inputs,
    const QuarantineState & quarantine, std::span<const SuspiciousGaussianAssessment> assessments,
    const ClusterHealthMap & health, std::size_t accepted_iteration) &&
{
    // Failure evidence is sampled before quarantine changes activity or rolls back a model.
    const auto failure_mask{ BuildSuspiciousFailureAtomMask(m_selection.block_activity, assessments) };
    const auto suspicious_count{ static_cast<std::size_t>(std::ranges::count_if(failure_mask,
        [](char value) { return value != 0; })) };
    auto next_quarantine{ quarantine };
    const auto transition{ next_quarantine.UpdateAfterIteration(
        m_selection.accepted_cluster_diagnostic_list, m_selection.rejected_cluster_diagnostic_list,
        m_selection.block_activity, assessments, health, m_selection.assembled_state,
        inputs.previous_state, inputs.previous_polish_provenance,
        m_selection.assembled_polish_provenance, accepted_iteration) };
    if (transition)
    {
        m_selection.final_audit_objective.reset();
        ReauditFallbackSelection(inputs);
    }
    ObservePhaseState(inputs.context, "final-selection",
            m_selection.accepted_key_list.empty() ? inputs.previous_state : m_selection.assembled_state);
    return CandidateTransaction(std::move(m_selection), std::move(next_quarantine), suspicious_count, transition);
}

CandidateCommitResult CandidateTransaction::Commit(const SecondStageContext & context,
    FitState & previous_state, FitState & accepted_state, PolishProvenance & provenance,
    ClusterObjectiveStateMap & history, QuarantineState & quarantine,
    TrustRegionStateSet & radii, IterationResult & result) &&
{
    const bool accepted{ !m_selection.accepted_key_list.empty() };
    quarantine = std::move(m_quarantine);
    history = std::move(m_selection.cluster_objective_state);
    LogBestObjectivePublication(context, history);
    result.trust_region_update = radii.ApplyRadiusUpdates(
        m_selection.grow_trust_region_key_list, m_selection.shrink_trust_region_key_list,
        m_selection.rejected_key_list, m_selection.exhausted_key_list);
    result.accepted_cluster_diagnostic_list = std::move(m_selection.accepted_cluster_diagnostic_list);
    result.rejected_cluster_diagnostic_list = std::move(m_selection.rejected_cluster_diagnostic_list);
    result.boundary_reconciliation_diagnostic_list = std::move(m_selection.boundary_reconciliation_diagnostic_list);
    if (accepted)
    {
        accepted_state = std::move(m_selection.assembled_state);
        provenance = std::move(m_selection.assembled_polish_provenance);
    }
    else accepted_state = std::move(previous_state);
    return {std::move(m_selection.block_activity), m_selection.final_audit_objective,
        m_selection.polish_progress, m_suspicious_atom_count, accepted,
        !m_selection.rejected_key_list.empty(), m_quarantine_transition};
}
} // namespace rhbm_gem::core::detail
