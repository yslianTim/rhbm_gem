#pragma once

#include "core/detail/second_stage/CandidateEvidence.hpp"

namespace rhbm_gem::core::detail {
class SecondStageObservationSession;

class JointCandidateObservation;

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
    const ObjectiveBreakdown * previous_audit;
};

struct BoundaryCorrectionReference
{
    BoundaryAcceptancePolicy policy;
    const std::map<ClusterKey, std::vector<SampleRef>> & samples_by_key;
    const ObjectiveDomain & domain;
    const ObjectiveByKey & previous_objective_by_key;
    const ObjectiveBreakdown * best_audit;
    PerformanceCounters & counters;
    const BoundaryReconciliationComponent & component;
    const FitStateView & endpoint;
    const ObjectiveBreakdown & previous_audit;
    const ObjectiveBreakdown & improvement;
    double damping;
};

bool EvaluateBoundaryCorrection(const CandidateEvaluationOverlay &, const BoundaryCorrectionReference &, JointCandidateObservation * observation = nullptr);

struct GlobalCandidateReference
{
    const std::vector<SampleRef> & samples;
    const ObjectiveDomain & domain;
    const ObjectiveBreakdown * best;
    const ObjectiveBreakdown * previous;
    PerformanceCounters & counters;
    SecondStageObservationSession * observation{ nullptr };
    bool rescue_audit{ false };
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
    double damping;
    std::size_t round;
};

struct FinalPolishCandidateEvaluation
{
    std::optional<ObjectiveBreakdown> objective{};
    std::size_t suspicious_atom_count{ 0 };
};

LocalCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const LocalCandidateReference &);
std::optional<ObjectiveBreakdown> EvaluateBoundaryCandidate(const CandidateEvaluationOverlay &, const BoundaryCandidateReference &, JointCandidateObservation * observation = nullptr);
std::optional<ObjectiveBreakdown> EvaluateCandidate(const CandidateEvaluationOverlay &, const GlobalCandidateReference &);
FinalPolishCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const FinalPolishCandidateReference &, JointCandidateObservation * observation = nullptr);

} // namespace rhbm_gem::core::detail
