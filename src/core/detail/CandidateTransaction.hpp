#pragma once

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/Quarantine.hpp"

#include <utility>

namespace rhbm_gem::core::detail {
struct IterationResult;
struct CandidateEvaluation;

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
    const CandidateSelection & View() const { return m_selection; }
    CandidateCommitResult Commit(const SecondStageContext &, FitState & previous_state,
        FitState & accepted_state, PolishProvenance &, ClusterObjectiveStateMap &,
        QuarantineState &, TrustRegionStateSet &, IterationResult &) &&;
};

class CandidateTransactionBuilder
{
    CandidateSelection m_selection{};
    void RejectSelectionKeys(
        const CandidateSelectionInputs & inputs,
        const std::vector<ClusterKey> & key_list,
        bool exhausted);
    bool TryBoundaryJointCorrection(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown & previous_audit_objective,
        const ObjectiveBreakdown & improvement_reference_objective,
        const FitStatePatch & endpoint_patch,
        BoundaryComponentReconciliationDiagnostic & diagnostic);
    bool TryBacktrackBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown * previous_audit_objective,
        const FitStatePatch & endpoint_patch,
        BoundaryComponentReconciliationDiagnostic & diagnostic);
    void ReconcileBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown * previous_audit_objective);
    void PromoteBoundaryRescueKeys(
        const CandidateSelectionInputs & inputs,
        const std::vector<ClusterKey> & rescue_key_list,
        const FitStatePatch & endpoint_patch,
        BoundaryComponentAcceptedSource accepted_source);
    bool TryRescueBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown & previous_audit_objective,
        const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key);
    bool RescueRejectedBoundaryClusters(
        const CandidateSelectionInputs & inputs,
        const ObjectiveBreakdown & previous_audit_objective,
        const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key);
    void MarkBoundaryDiagnosticRejected(
        const std::vector<ClusterKey> & key_list,
        bool exhausted);
    void AuditAndSalvageFinalSelection(
        const CandidateSelectionInputs & inputs,
        const ObjectiveBreakdown & previous_audit_objective);
    void ReauditFallbackSelection(const CandidateSelectionInputs & inputs);
public:
    CandidateTransactionBuilder() = default;
    explicit CandidateTransactionBuilder(CandidateSelection initial) : m_selection(std::move(initial)) {}
    CandidateTransactionBuilder(const CandidateTransactionBuilder &) = delete;
    void Select(const CandidateSelectionInputs &);
    void ReconcileSelectedBoundaries(
        const CandidateSelectionInputs & inputs,
        const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key);
    const CandidateSelection & View() const { return m_selection; }
    CandidateTransaction Finish(const CandidateSelectionInputs &, const QuarantineState &,
        std::span<const SuspiciousGaussianAssessment>, const ClusterHealthMap &, std::size_t accepted_iteration) &&;
};
} // namespace rhbm_gem::core::detail
