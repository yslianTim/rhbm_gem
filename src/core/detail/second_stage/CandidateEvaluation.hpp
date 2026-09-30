#pragma once

#include "core/detail/second_stage/CandidateEvidence.hpp"

namespace rhbm_gem::core::detail {
bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance);

enum class LocalObjectivePolicy
{
    PreviousNonRegression, StrictReferenceImprovement
};

enum class BoundaryAcceptancePolicy
{
    Ordinary, CooperativeRescue
};

struct LocalCandidateEvaluation
{
    bool accepted{ false };
    CandidateDecisionEvidence evidence{};
};

struct LocalCandidateReference
{
    LocalObjectivePolicy policy;
    const ClusterKey & key;
    const std::vector<SampleRef> & samples;
    const ObjectiveBreakdown * objective_reference;
    const ObjectiveDomain & domain;
    CandidateDecisionEvidence evidence;
    const FitStatePatch * member_best{ nullptr };
};

struct BoundaryCandidateReference
{
    BoundaryAcceptancePolicy policy;
    const std::map<ClusterKey, std::vector<SampleRef>> & samples_by_key;
    const ObjectiveDomain & domain;
    const ObjectiveByKey & previous_objective_by_key;
    const ObjectiveBreakdown * best_audit;
    const BoundaryReconciliationComponent & component;
    const MemberBestState * member_best{ nullptr };
};

struct GlobalCandidateReference
{
    const std::vector<SampleRef> & samples;
    const ObjectiveDomain & domain;
    const ObjectiveBreakdown * best;
    const ObjectiveBreakdown * previous;
};

struct FinalPolishCandidateReference
{
    const DependencyPolishComponent & component;
    const CouplingGraphPartition & partition;
    const ObjectiveDomain & domain;
    const FitStateView & endpoint;
    const ObjectiveBreakdown & base_objective;
    const ObjectiveBreakdown & endpoint_objective;
};

struct FinalPolishCandidateEvaluation
{
    std::optional<ObjectiveBreakdown> objective{};
    std::size_t suspicious_atom_count{ 0 };
};

LocalCandidateEvaluation EvaluateLocalCandidate(const CandidateEvaluationOverlay &, const LocalCandidateReference &);
std::optional<ObjectiveBreakdown> EvaluateBoundaryCandidate(const CandidateEvaluationOverlay &, const BoundaryCandidateReference &,
    const ObjectiveBreakdown * previous_audit);
bool EvaluateBoundaryCorrection(const CandidateEvaluationOverlay &, const BoundaryCandidateReference &,
    const FitStateView & endpoint, const ObjectiveBreakdown & previous_audit, const ObjectiveBreakdown & improvement);
std::optional<ObjectiveBreakdown> EvaluateGlobalCandidate(const CandidateEvaluationOverlay &, const GlobalCandidateReference &);
FinalPolishCandidateEvaluation EvaluateFinalPolishCandidate(const CandidateEvaluationOverlay &, const FinalPolishCandidateReference &);

} // namespace rhbm_gem::core::detail
