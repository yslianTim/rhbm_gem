#pragma once

#include "core/detail/second_stage/CandidateState.hpp"

namespace rhbm_gem::core::detail {

class JointCandidateObservation;
struct CandidateSelectionInputs;

enum class LocalObjectivePolicy
{
    PreviousNonRegression, StrictReferenceImprovement
};

enum class BoundaryAcceptancePolicy
{
    Ordinary, CooperativeRescue
};

enum class CandidateFailureStage
{
    None, Trust, Guard
};

struct LocalCandidateEvaluation
{
    bool accepted{ false };
    CandidateDecisionEvidence evidence{};
};

struct BoundaryCandidateEvaluation
{
    ObjectiveBreakdown audit_objective{};
};

struct CandidatePreflightReference
{
    const ClusterKey & key;
    const SuspiciousBlockActivity & activity;
    double step_norm;
    double radius;
};

struct CandidatePreflightEvaluation
{
    CandidateFailureStage failure_stage{ CandidateFailureStage::None };
    std::optional<StabilizationTerminalEvidence> guard_failure{};
};

CandidatePreflightEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const CandidatePreflightReference &);

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
    const CandidateSelectionInputs & inputs;
    const BoundaryReconciliationComponent & component;
    const ObjectiveBreakdown * previous_audit;
};

struct BoundaryCorrectionReference
{
    BoundaryAcceptancePolicy policy;
    const CandidateSelectionInputs & inputs;
    const BoundaryReconciliationComponent & component;
    const FitStateView & endpoint;
    const ObjectiveBreakdown & previous_audit;
    const ObjectiveBreakdown & improvement;
    double damping;
};

struct BoundaryCorrectionEvaluation
{
    std::size_t suspicious_atom_count{ 0 };
    std::optional<ObjectiveBreakdown> raw_objective{};
    std::optional<BoundaryCandidateEvaluation> members{};
    bool accepted{ false };
};

BoundaryCorrectionEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const BoundaryCorrectionReference &, JointCandidateObservation * observation = nullptr);

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
    double damping;
    std::size_t round;
};

struct FinalPolishCandidateEvaluation
{
    std::optional<ObjectiveBreakdown> objective{};
    std::size_t suspicious_atom_count{ 0 };
};

LocalCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const LocalCandidateReference &);
std::optional<BoundaryCandidateEvaluation> EvaluateCandidate(const CandidateEvaluationOverlay &, const BoundaryCandidateReference &, JointCandidateObservation * observation = nullptr);
std::optional<ObjectiveBreakdown> EvaluateCandidate(const CandidateEvaluationOverlay &, const GlobalCandidateReference &);
FinalPolishCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const FinalPolishCandidateReference &, JointCandidateObservation * observation = nullptr);

} // namespace rhbm_gem::core::detail
