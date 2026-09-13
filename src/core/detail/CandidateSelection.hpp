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

struct ClusterCandidateDecision
{
    ClusterKey key{};
    CandidateDecisionEvidence evidence{};
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

struct BoundaryComponentDecision
{
    std::vector<ClusterKey> key_list{};
    std::optional<double> accepted_factor{};
    BoundaryComponentAcceptedSource accepted_source{ BoundaryComponentAcceptedSource::None };
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
    std::vector<ClusterCandidateDecision> accepted_cluster_evidence_list{};
    std::vector<ClusterCandidateDecision> rejected_cluster_evidence_list{};
    std::vector<BoundaryComponentDecision> boundary_decision_list{};
    std::optional<ObjectiveBreakdown> final_audit_objective{};
    PolishProgress polish_progress{};
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
};

} // namespace rhbm_gem::core::detail
