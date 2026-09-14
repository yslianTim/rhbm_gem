#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include "support/SecondStageNumericalProbe.hpp"
#else
#define RHBM_TEST_COMMIT(result, quarantine, radii) ((void)0)
#endif
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"

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
    for (auto & decision : m_selection.accepted_cluster_evidence_list)
        m_candidate_by_key[decision.key].evidence = std::move(decision);
    for (auto & decision : m_selection.rejected_cluster_evidence_list)
        m_candidate_by_key[decision.key].evidence = std::move(decision);
    m_selection.accepted_key_list.clear();
    m_selection.rejected_key_list.clear();
    m_selection.exhausted_key_list.clear();
    m_selection.shrink_trust_region_key_list.clear();
    m_selection.accepted_cluster_evidence_list.clear();
    m_selection.rejected_cluster_evidence_list.clear();
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
    // rejection event order used by terminal evidence and observation output.
    std::vector<PendingCandidate *> rejected;
    for (auto & [key, candidate] : m_candidate_by_key)
    {
        if (candidate.selected)
        {
            m_selection.accepted_key_list.emplace_back(key);
            if (candidate.evidence)
                m_selection.accepted_cluster_evidence_list.emplace_back(std::move(*candidate.evidence));
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
        if (candidate->evidence)
            m_selection.rejected_cluster_evidence_list.emplace_back(std::move(*candidate->evidence));
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
        m_selection.accepted_cluster_evidence_list, m_selection.rejected_cluster_evidence_list,
        m_selection.block_activity, assessments, health, operator_evidence,
        m_selection.accepted_key_list.empty() ? inputs.previous_state : m_selection.assembled_state,
        inputs.previous_state, recovery_revision) };
    if (inputs.observation) inputs.observation->ObserveQuarantine(quarantine, next_quarantine);
    return CandidateTransaction(std::move(m_selection), std::move(next_quarantine), suspicious_count, transition);
}

CandidateCommitResult CandidateTransaction::Commit(FitState & previous_state, FitState & accepted_state, PolishProvenance & provenance,
    QuarantineState & quarantine,
    TrustRegionStateSet & radii, SecondStageObservationSession * observation) &&
{
    const bool accepted{ !m_selection.accepted_key_list.empty() };
    quarantine = std::move(m_quarantine);
    CandidateCommitResult result;
    result.trust_region_update = radii.ApplyRadiusUpdates(
        m_selection.shrink_trust_region_key_list,
        m_selection.rejected_key_list, m_selection.exhausted_key_list);
    result.accepted_key_list = m_selection.accepted_key_list;
    result.rejected_key_list = m_selection.rejected_key_list;
    if (accepted)
    {
        accepted_state = std::move(m_selection.assembled_state);
        provenance = std::move(m_selection.assembled_polish_provenance);
    }
    else accepted_state = std::move(previous_state);
    result.block_activity = std::move(m_selection.block_activity);
    result.final_audit_objective = m_selection.final_audit_objective;
    result.polish_progress = m_selection.polish_progress;
    result.suspicious_atom_count = m_suspicious_atom_count;
    result.accepted = accepted;
    result.rejected_cluster = !m_selection.rejected_key_list.empty();
    result.quarantine_transition = m_quarantine_transition;
    if (observation) observation->ObserveCommit(m_selection, result, radii);
    RHBM_TEST_COMMIT(result, quarantine, radii);
    return result;
}
} // namespace rhbm_gem::core::detail
