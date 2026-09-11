#include "core/detail/ClusterHistoryObserver.hpp"
#include "core/detail/CandidateTransaction.hpp"
#include "core/detail/IterationProcess.hpp"
#include "core/detail/Diagnosis.hpp"
#include "core/detail/PhaseAudit.hpp"

#include <algorithm>

namespace rhbm_gem::core::detail {

CandidateTransactionBuilder::CandidateTransactionBuilder(CandidateSelection initial)
    : m_selection(std::move(initial))
{
    for (const auto & key : m_selection.accepted_key_list)
        m_candidate_by_key[key].selected = true;
    for (const auto & key : m_selection.rejected_key_list)
        m_candidate_by_key[key].rejection_order = m_next_rejection_order++;
    for (const auto & key : m_selection.exhausted_key_list)
        m_candidate_by_key[key].exhausted = true;
    for (const auto & key : m_selection.shrink_trust_region_key_list)
        m_candidate_by_key[key].shrink_trust_region = true;
    for (auto & diagnostic : m_selection.accepted_cluster_diagnostic_list)
        m_candidate_by_key[diagnostic.key].diagnostic = std::move(diagnostic);
    for (auto & diagnostic : m_selection.rejected_cluster_diagnostic_list)
        m_candidate_by_key[diagnostic.key].diagnostic = std::move(diagnostic);
    m_selection.accepted_key_list.clear();
    m_selection.rejected_key_list.clear();
    m_selection.exhausted_key_list.clear();
    m_selection.shrink_trust_region_key_list.clear();
    m_selection.accepted_cluster_diagnostic_list.clear();
    m_selection.rejected_cluster_diagnostic_list.clear();
}

std::vector<ClusterKey> CandidateTransactionBuilder::SelectedKeys() const
{
    std::vector<ClusterKey> keys;
    for (const auto & [key, candidate] : m_candidate_by_key)
        if (candidate.selected) keys.emplace_back(key);
    return keys;
}

void CandidateTransactionBuilder::MaterializeSelection()
{
    // Classify once after component selection and global salvage. Preserve the
    // historical rejection event order without moving diagnostics between lists.
    std::vector<PendingCandidate *> rejected;
    for (auto & [key, candidate] : m_candidate_by_key)
    {
        if (candidate.selected)
        {
            m_selection.accepted_key_list.emplace_back(key);
            if (candidate.diagnostic)
                m_selection.accepted_cluster_diagnostic_list.emplace_back(std::move(*candidate.diagnostic));
        }
        else
        {
            m_selection.rejected_key_list.emplace_back(key);
            rejected.emplace_back(&candidate);
        }
        if (candidate.exhausted) m_selection.exhausted_key_list.emplace_back(key);
        if (candidate.shrink_trust_region) m_selection.shrink_trust_region_key_list.emplace_back(key);
    }
    std::ranges::sort(rejected, {}, &PendingCandidate::rejection_order);
    for (auto * candidate : rejected)
        if (candidate->diagnostic)
            m_selection.rejected_cluster_diagnostic_list.emplace_back(std::move(*candidate->diagnostic));
}

CandidateTransaction CandidateTransactionBuilder::Finish(const CandidateSelectionInputs & inputs,
    const QuarantineState & quarantine, std::span<const SuspiciousGaussianAssessment> assessments,
    const ClusterHealthMap & health, const FixedPointOperatorEvidence & operator_evidence,
    std::size_t recovery_revision) &&
{
    // Quarantine publishes only next-iteration activity, never changes the audited model.
    const auto failure_mask{ BuildSuspiciousFailureAtomMask(m_selection.block_activity, assessments) };
    const auto suspicious_count{ static_cast<std::size_t>(std::ranges::count_if(failure_mask,
        [](char value) { return value != 0; })) };
    auto next_quarantine{ quarantine };
    const auto transition{ next_quarantine.UpdateAfterIteration(
        m_selection.accepted_cluster_diagnostic_list, m_selection.rejected_cluster_diagnostic_list,
        m_selection.block_activity, assessments, health, operator_evidence,
        m_selection.accepted_key_list.empty() ? inputs.previous_state : m_selection.assembled_state,
        inputs.previous_state, recovery_revision) };
    ObservePhaseState(inputs.context, "final-selection",
            m_selection.accepted_key_list.empty() ? inputs.previous_state : m_selection.assembled_state);
    return CandidateTransaction(std::move(m_selection), std::move(next_quarantine), suspicious_count, transition);
}

CandidateCommitResult CandidateTransaction::Commit(const SecondStageContext & context,
    FitState & previous_state, FitState & accepted_state, PolishProvenance & provenance,
    QuarantineState & quarantine,
    TrustRegionStateSet & radii, IterationResult & result) &&
{
    const bool accepted{ !m_selection.accepted_key_list.empty() };
    quarantine = std::move(m_quarantine);
    result.trust_region_update = radii.ApplyRadiusUpdates(
        m_selection.shrink_trust_region_key_list,
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
    if (context.cluster_history) context.cluster_history->Publish(context);
    return {std::move(m_selection.block_activity), m_selection.final_audit_objective,
        m_selection.polish_progress, m_suspicious_atom_count, accepted,
        !m_selection.rejected_key_list.empty(), m_quarantine_transition};
}
} // namespace rhbm_gem::core::detail
