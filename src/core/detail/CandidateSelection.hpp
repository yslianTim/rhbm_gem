#pragma once

#include "core/detail/SecondStageFitting.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointFitting.hpp"
#include "core/detail/ObjectiveEvaluation.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <vector>

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

struct TrustRegionRadiusUpdate
{
    std::vector<ClusterKey> changed_key_list{};
    std::vector<ClusterKey> saturated_key_list{};
};

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
    const ObjectiveAttemptDiagnostic & diagnostic);

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

struct ClusterCandidateDiagnostic
{
    ClusterKey key{};
    ObjectiveAttemptDiagnostic attempt{};

    bool boundary_rescued{ false };
};

struct PolishProgress
{
    std::size_t eligible_count{ 0 };
    std::size_t accepted_count{ 0 };
    std::size_t rejected_count{ 0 };
    std::size_t skipped_count{ 0 };
};

enum class BoundaryComponentAcceptedSource
{
    None,
    Endpoint,
    JointCorrection,
    Backtracking
};

struct BoundaryComponentReconciliationDiagnostic
{
    std::vector<ClusterKey> key_list{};
    std::size_t atom_count{ 0 };
    std::size_t boundary_sample_count{ 0 };
    std::size_t trial_count{ 1 };
    std::optional<double> accepted_factor{};
    BoundaryComponentAcceptedSource accepted_source{ BoundaryComponentAcceptedSource::None };
    std::optional<BoundaryJointCorrectionStatus> joint_correction_status{};
    std::size_t interface_atom_count{ 0 };
    std::size_t shape_active_atom_count{ 0 };
    std::size_t offset_active_atom_count{ 0 };
    std::size_t suspicious_candidate_atom_count{ 0 };
    std::size_t joint_parameter_count{ 0 };
    std::optional<double> joint_damping{};
    std::optional<double> maximum_normalized_trust_step{};
    std::optional<double> previous_component_objective{};
    std::optional<double> endpoint_component_objective{};
    std::optional<double> joint_reference_component_objective{};
    std::optional<double> joint_candidate_component_objective{};
    std::optional<double> candidate_component_objective{};
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rescue_candidate_cluster_count{ 0 };
    std::size_t rescued_cluster_count{ 0 };
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
    std::optional<double> component_improvement{};
    std::optional<double> global_improvement{};
    std::vector<JointCandidateObjectiveDiagnostic> objective_diagnostic_list{};
    bool is_rescue_attempt{ false };
    bool exhausted{ false };
};

struct CandidateSelection
{
    SuspiciousBlockActivity block_activity{};
    FitState assembled_state{};
    PolishProvenance assembled_polish_provenance{};
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    std::vector<ClusterKey> shrink_trust_region_key_list{};
    std::vector<ClusterKey> exhausted_key_list{};
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<BoundaryComponentReconciliationDiagnostic> boundary_reconciliation_diagnostic_list{};
    std::optional<ObjectiveBreakdown> final_audit_objective{};
    PolishProgress polish_progress{};
};

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
};

} // namespace rhbm_gem::core::detail
