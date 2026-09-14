#pragma once

#include "core/detail/second_stage/CandidateEvidence.hpp"

namespace rhbm_gem::core::detail {
class SecondStageObservationSession;

class JointCandidateObservation;

bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance,
    ObjectiveProgressGateEvidence * evidence = nullptr);

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
    PerformanceCounters & counters;
};

struct BoundaryCandidateReference
{
    BoundaryAcceptancePolicy policy;
    const std::map<ClusterKey, std::vector<SampleRef>> & samples_by_key;
    const ObjectiveDomain & domain;
    const ObjectiveByKey & previous_objective_by_key;
    const ObjectiveBreakdown * best_audit;
    PerformanceCounters & counters;
    const BoundaryReconciliationComponent & component;
};

struct GlobalCandidateReference
{
    const std::vector<SampleRef> & samples;
    const ObjectiveDomain & domain;
    const ObjectiveBreakdown * best;
    const ObjectiveBreakdown * previous;
    PerformanceCounters & counters;
};

struct FinalPolishCandidateReference
{
    const DependencyPolishComponent & component;
    const CouplingGraphPartition & partition;
    const ObjectiveDomain & domain;
    const FitStateView & endpoint;
    const ObjectiveBreakdown & base_objective;
    const ObjectiveBreakdown & endpoint_objective;
    PerformanceCounters & counters;
};

struct FinalPolishCandidateEvaluation
{
    std::optional<ObjectiveBreakdown> objective{};
    std::size_t suspicious_atom_count{ 0 };
};

LocalCandidateEvaluation EvaluateLocalCandidate(const CandidateEvaluationOverlay &, const LocalCandidateReference &);
std::optional<ObjectiveBreakdown> EvaluateBoundaryCandidate(const CandidateEvaluationOverlay &, const BoundaryCandidateReference &,
    const ObjectiveBreakdown * previous_audit, JointCandidateObservation * observation = nullptr);
bool EvaluateBoundaryCorrection(const CandidateEvaluationOverlay &, const BoundaryCandidateReference &,
    const FitStateView & endpoint, const ObjectiveBreakdown & previous_audit, const ObjectiveBreakdown & improvement,
    JointCandidateObservation * observation = nullptr);
std::optional<ObjectiveBreakdown> EvaluateGlobalCandidate(const CandidateEvaluationOverlay &, const GlobalCandidateReference &,
    SecondStageObservationSession * observation = nullptr, bool rescue_audit = false);
FinalPolishCandidateEvaluation EvaluateFinalPolishCandidate(const CandidateEvaluationOverlay &, const FinalPolishCandidateReference &, JointCandidateObservation * observation = nullptr);

} // namespace rhbm_gem::core::detail
