#pragma once

#include "core/detail/CandidateSelection.hpp"

namespace rhbm_gem::core::detail {

enum class CandidateScope
{
    LocalSearch, LocalPolish, Boundary,
    CooperativeRescue
};

enum class CandidateFailureStage
{
    None, Trust, Guard
};

struct LocalCandidateEvaluation
{
    bool accepted{ false };
    ObjectiveAttemptDiagnostic diagnostic{};
    std::optional<ClusterObjectiveState> objective_state{};
};

struct BoundaryCandidateEvaluation
{
    ClusterObjectiveStateMap objective_state_by_key{};
    ObjectiveBreakdown audit_objective{};
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
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
    std::optional<StabilizationTerminalDiagnostic> guard_failure{};
};

CandidatePreflightEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const CandidatePreflightReference &);

struct LocalCandidateReference
{
    const ClusterKey & key;
    const std::vector<SampleRef> & samples;
    const ObjectiveBreakdown * previous;
    const ObjectiveDomain & domain;
    const ClusterObjectiveState & history;
    ObjectiveAttemptDiagnostic diagnostic;
    PerformanceCounters & counters;
    std::string_view source{ "local-candidate" };
};

struct BoundaryCandidateReference
{
    const CandidateSelectionInputs & inputs;
    const BoundaryReconciliationComponent & component;
    const ObjectiveBreakdown * previous_audit;
    JointCandidateObjectiveDiagnostic * record;
};

struct BoundaryCorrectionReference
{
    const CandidateSelectionInputs & inputs;
    const BoundaryReconciliationComponent & component;
    const FitStateView & endpoint;
    const ObjectiveBreakdown & previous_audit;
    const ObjectiveBreakdown & improvement;
    std::vector<JointCandidateObjectiveDiagnostic> & records;
    double damping;
};

struct BoundaryCorrectionEvaluation
{
    std::size_t suspicious_atom_count{ 0 };
    std::optional<ObjectiveBreakdown> raw_objective{};
    std::optional<BoundaryCandidateEvaluation> members{};
    JointCandidateObjectiveDiagnostic * record{ nullptr };
    bool accepted{ false };
};

BoundaryCorrectionEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, CandidateScope, const BoundaryCorrectionReference &);

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
    bool quiet;
    std::vector<JointCandidateObjectiveDiagnostic> & records;
    double damping;
    std::size_t round;
};

struct FinalPolishCandidateEvaluation
{
    std::optional<ObjectiveBreakdown> objective{};
    std::size_t suspicious_atom_count{ 0 };
};

LocalCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, CandidateScope, const LocalCandidateReference &);
std::optional<BoundaryCandidateEvaluation> EvaluateCandidate(const CandidateEvaluationOverlay &, CandidateScope, const BoundaryCandidateReference &);
std::optional<ObjectiveBreakdown> EvaluateCandidate(const CandidateEvaluationOverlay &, const GlobalCandidateReference &);
FinalPolishCandidateEvaluation EvaluateCandidate(const CandidateEvaluationOverlay &, const FinalPolishCandidateReference &);

} // namespace rhbm_gem::core::detail
