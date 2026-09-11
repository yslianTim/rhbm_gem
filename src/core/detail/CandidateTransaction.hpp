#pragma once

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/Quarantine.hpp"

#include <utility>

namespace rhbm_gem::core::detail {
struct IterationResult;
enum class CandidateScope;

struct CandidateCommitResult
{
    SuspiciousBlockActivity block_activity{};
    std::optional<ObjectiveBreakdown> final_audit_objective{};
    PolishProgress polish_progress{};
    std::size_t suspicious_atom_count{ 0 };
    bool accepted{ false };
    bool rejected_cluster{ false };
    bool quarantine_transition{ false };
};

class CandidateTransaction
{
    CandidateSelection m_selection;
    QuarantineState m_quarantine;
    std::size_t m_suspicious_atom_count;
    bool m_quarantine_transition;
    friend class CandidateTransactionBuilder;
    CandidateTransaction(CandidateSelection selection, QuarantineState quarantine,
        std::size_t suspicious_count, bool transition)
        : m_selection(std::move(selection)), m_quarantine(std::move(quarantine)),
          m_suspicious_atom_count(suspicious_count), m_quarantine_transition(transition) {}
public:
    CandidateTransaction(const CandidateTransaction &) = delete;
    CandidateTransaction(CandidateTransaction &&) = default;
    CandidateCommitResult Commit(const SecondStageContext &, FitState & previous_state,
        FitState & accepted_state, PolishProvenance &,
        QuarantineState &, TrustRegionStateSet &, IterationResult &) &&;
};

class CandidateTransactionBuilder
{
    struct PendingCandidate
    {
        std::optional<ClusterCandidateDiagnostic> diagnostic{};
        std::optional<FitStatePatch> cooperative_patch{};
        bool selected{ false };
        bool exhausted{ false };
        bool shrink_trust_region{ false };
        std::size_t rejection_order{ 0 };
    };
    struct ComponentCandidate
    {
        FitStatePatch patch{};
        std::vector<std::pair<std::size_t, char>> provenance_updates{};
        std::size_t history_observation{ 0 };
    };
    CandidateSelection m_selection{};
    std::map<ClusterKey, PendingCandidate> m_candidate_by_key{};
    std::size_t m_next_rejection_order{ 0 };
    std::vector<ClusterKey> SelectedKeys() const;
    void MaterializeSelection();
    void ApplyComponentCandidate(
        const CandidateSelectionInputs &, const BoundaryReconciliationComponent &,
        const FitStatePatch & endpoint_patch, ComponentCandidate,
        BoundaryComponentAcceptedSource, CandidateScope);
    void RejectSelectionKeys(
        const CandidateSelectionInputs & inputs,
        const std::vector<ClusterKey> & key_list,
        bool exhausted);
    std::optional<ComponentCandidate> TryBoundaryJointCorrection(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown & previous_audit_objective,
        const ObjectiveBreakdown & improvement_reference_objective,
        const FitStatePatch & endpoint_patch,
        BoundaryComponentReconciliationDiagnostic & diagnostic,
        CandidateScope scope);
    std::optional<ComponentCandidate> TryBacktrackBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown * previous_audit_objective,
        const FitStatePatch & endpoint_patch,
        BoundaryComponentReconciliationDiagnostic & diagnostic,
        CandidateScope scope);
    bool ReconcileBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown * previous_audit_objective,
        CandidateScope scope);
    bool ReconcileCooperativeComponents(
        const CandidateSelectionInputs & inputs,
        const ObjectiveBreakdown & previous_audit_objective);
    void MarkBoundaryDiagnosticRejected(
        const std::vector<ClusterKey> & key_list,
        bool exhausted);
    void AuditAndSalvageFinalSelection(
        const CandidateSelectionInputs & inputs,
        const ObjectiveBreakdown & previous_audit_objective);
public:
    CandidateTransactionBuilder() = default;
    explicit CandidateTransactionBuilder(CandidateSelection initial);
    CandidateTransactionBuilder(const CandidateTransactionBuilder &) = delete;
    void Select(const CandidateSelectionInputs &);
    void ReconcileSelectedBoundaries(const CandidateSelectionInputs & inputs);
    const CandidateSelection & View() const { return m_selection; }
    CandidateTransaction Finish(const CandidateSelectionInputs &, const QuarantineState &,
        std::span<const SuspiciousGaussianAssessment>, const ClusterHealthMap &, const FixedPointOperatorEvidence &, std::size_t recovery_revision) &&;
};
} // namespace rhbm_gem::core::detail
