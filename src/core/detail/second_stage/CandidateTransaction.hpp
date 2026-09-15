#pragma once

#include "core/detail/second_stage/CandidateState.hpp"
#include "core/detail/second_stage/Quarantine.hpp"

#include <map>
#include <utility>

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

class TrustRegionStateSet
{
    std::map<ClusterKey, unsigned int> m_shrink_level_by_key{};

public:
    void Reconcile(const std::vector<ClusterKey> & key_list);
    double GetRadius(const ClusterKey & key) const;
    void ResetToMinimum(const std::vector<ClusterKey> & key_list);
    TrustRegionRadiusUpdate ApplyRadiusUpdates(
        const std::vector<ClusterKey> & accepted_shrink_key_list,
        const std::vector<ClusterKey> & rejected_key_list,
        const std::vector<ClusterKey> & exhausted_key_list);

};

bool ShouldShrinkAcceptedTrustRegionRadius(
    std::optional<double> first_objective_evaluated_factor,
    std::optional<double> accepted_factor);

enum class BacktrackingStepStatus
{
    CandidateReady,
    InvalidCandidate,
    Exhausted
};

struct BacktrackingStep
{
    BacktrackingStepStatus status{ BacktrackingStepStatus::Exhausted };
    double factor{ 0.0 };
    std::size_t trial_number{ 0 };
};

class BacktrackingWorkspace
{
    std::size_t m_previous_state_size{ 0 };
    double m_minimum_transformed_change{ 0.0 };
    double m_next_factor{ 0.5 };
    std::size_t m_trial_number{ 1 };
    std::vector<GaussianModel3D> m_previous_model_list{};
    std::vector<GaussianModel3D> m_endpoint_model_list{};
    FitStatePatch m_candidate_patch{};

public:
    BacktrackingWorkspace(
        const FitState & previous_state,
        const FitStatePatch & endpoint_patch,
        double minimum_transformed_change);

    BacktrackingWorkspace(const BacktrackingWorkspace &) = delete;
    BacktrackingWorkspace & operator=(const BacktrackingWorkspace &) = delete;
    BacktrackingStep BuildNextCandidate();
    const FitStatePatch & GetCandidatePatch() const { return m_candidate_patch; }

    PolishProvenance BuildCandidatePolishProvenance(
        const PolishProvenance & previous_provenance,
        const PolishProvenance & endpoint_provenance) const;

private:
    bool BuildCandidate(double factor);
    double GetMaximumTransformedChange() const;

};

class SecondStageObservationSession;

struct CandidateSelectionInputs
{
    // Algorithm inputs stay unchanged; updated activity is returned in CandidateSelection.
    // Solver workspaces, counters and observation state may mutate.
    const SecondStageContext & context;
    const FitOptions & options;
    const ResidualBaseline & residual_baseline;
    const CouplingGraphPartition & partition;
    const ClusterHealthMap & health_by_key;
    const FitState & previous_state;
    const PolishProvenance & previous_polish_provenance;
    const FitState & proposal_state;
    const SuspiciousBlockActivity & block_activity;
    const std::vector<double> & ridge_multiplier_list;
    const ObjectiveDomain & objective_domain;
    const ObjectiveByKey & previous_objective_by_key;
    const BestAuditState & best_audit_state;
    const TrustRegionStateSet & trust_region_state;
    ClusterSolverWorkspaceMap & solver_workspace_by_key;
    BoundaryJointCorrectionWorkspaceMap & boundary_joint_correction_workspace_by_key;
    PerformanceCounters & performance_counters;
    SecondStageObservationSession * observation{ nullptr };
    const MemberBestState * member_best{ nullptr };
};

class BoundaryObservationScope;
enum class BoundaryAcceptancePolicy;

struct CandidateCommitResult
{
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    TrustRegionRadiusUpdate trust_region_update{};
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
    CandidateCommitResult Commit(FitState & previous_state,
        FitState & accepted_state, PolishProvenance &,
        QuarantineState &, TrustRegionStateSet &, SecondStageObservationSession * = nullptr) &&;
};

class CandidateTransactionBuilder
{
    struct PendingCandidate
    {
        std::optional<ClusterCandidateDecision> evidence{};
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
    };
    CandidateSelection m_selection{};
    std::map<ClusterKey, PendingCandidate> m_candidate_by_key{};
    std::size_t m_next_rejection_order{ 0 };
    std::vector<ClusterKey> SelectedKeys() const;
    void MaterializeSelection();
    void ApplyComponentCandidate(
        const CandidateSelectionInputs &, const BoundaryReconciliationComponent &,
        const FitStatePatch & endpoint_patch, ComponentCandidate,
        BoundaryComponentAcceptedSource);
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
        BoundaryComponentDecision & decision,
        BoundaryObservationScope & observation, BoundaryAcceptancePolicy policy);
    std::optional<ComponentCandidate> TryBacktrackBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown * previous_audit_objective,
        const FitStatePatch & endpoint_patch,
        BoundaryComponentDecision & decision,
        BoundaryObservationScope & observation, BoundaryAcceptancePolicy policy);
    bool ReconcileBoundaryComponent(
        const CandidateSelectionInputs & inputs,
        const BoundaryReconciliationComponent & component,
        const ObjectiveBreakdown * previous_audit_objective,
        BoundaryAcceptancePolicy policy);
    bool ReconcileCooperativeComponents(
        const CandidateSelectionInputs & inputs,
        const ObjectiveBreakdown & previous_audit_objective);
    void AuditAndSalvageFinalSelection(
        const CandidateSelectionInputs & inputs,
        const ObjectiveBreakdown & previous_audit_objective, bool rescue_audit = false);
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
