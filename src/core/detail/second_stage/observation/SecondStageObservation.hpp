#pragma once

#include "core/detail/second_stage/observation/SecondStageDiagnostics.hpp"
#include "utils/hrl/EstimationAudit.hpp"
#include <array>

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

bool IsDebugLogLevelEnabled();

void RecordJointMemberRejection(
    JointCandidateObjectiveDiagnostic * record,
    const ClusterKey & key,
    const std::optional<ObjectiveBreakdown> & previous,
    const std::optional<ObjectiveBreakdown> & best,
    const std::optional<ObjectiveBreakdown> & candidate,
    bool best_checked);

JointCandidateObjectiveDiagnostic * BeginJointCandidateDiagnostic(
    bool quiet_mode,
    std::vector<JointCandidateObjectiveDiagnostic> & records,
    std::string_view source,
    std::optional<double> factor = std::nullopt,
    std::size_t round = 0);

struct BestObjectiveTraceEnvironment;
class ClusterHistoryObserver;
class PhaseAudit;
class TrustModelAudit;

class SecondStageObservationSession
{
public:
    std::shared_ptr<BestObjectiveTraceEnvironment> best_trace{};
    std::shared_ptr<ClusterHistoryObserver> cluster_history{};
    std::shared_ptr<PhaseAudit> phase_audit{};
    std::shared_ptr<TrustModelAudit> trust_model_audit{};
    IterationObservation iteration{};
    FinalDependencyPolishDiagnostic final_polish{};

    SecondStageObservationSession() = default;
    SecondStageObservationSession(const SecondStageObservationSession &) = delete;
    SecondStageObservationSession & operator=(const SecondStageObservationSession &) = delete;
};

struct CandidateSelectionInputs;
struct CandidateSelection;
struct BoundaryComponentDecision;
struct IterationResult;
struct IterationProposalResult;
struct TrustModelTrialRecord;
enum class BoundaryAcceptancePolicy;

// Observation entry points preserve the existing lifecycle and never choose candidates.
void BeginClusterHistoryObserver(SecondStageObservationSession &, bool quiet) noexcept;
void ObserveHistoryPartition(SecondStageObservationSession * observation, const SecondStageContext &, const CouplingGraphPartition &,
    const ObjectiveDomain &, const FitState &) noexcept;
void ObserveHistoryBackground(SecondStageObservationSession * observation, const SecondStageContext &, const std::shared_ptr<const FrozenBackground> &,
    const CouplingGraphPartition &, const ObjectiveDomain &, const FitState &) noexcept;
void ObserveHistoryAttempt(SecondStageObservationSession * observation, SecondStageContext &, const ObjectiveByKey &, const FitState &,
    const CouplingGraphPartition &, const ObjectiveDomain &, std::size_t attempt, std::size_t accepted) noexcept;
void ObserveHistorySearch(SecondStageObservationSession * observation, const ClusterKey &) noexcept;
void ObserveLocalHistory(SecondStageObservationSession * observation, const CandidateEvaluationOverlay &, const ClusterKey &,
    const std::vector<SampleRef> &, const ObjectiveDomain &, std::string_view source,
    bool accepted, ObjectiveAttemptDiagnostic &) noexcept;
void ObserveBoundaryHistory(SecondStageObservationSession * observation, JointCandidateObjectiveDiagnostic *) noexcept;
void ObserveBoundaryMemberHistory(SecondStageObservationSession * observation, const CandidateEvaluationOverlay &, const ClusterKey &,
    const std::vector<SampleRef> &, const ObjectiveDomain &, bool accepted,
    const CandidateDecisionEvidence &, JointCandidateObjectiveDiagnostic *) noexcept;
void ObserveBoundaryHistoryAccepted(SecondStageObservationSession * observation, std::size_t token) noexcept;
void ObserveHistoryRejected(SecondStageObservationSession * observation, const ClusterKey &) noexcept;
void ObserveHistoryPublication(SecondStageObservationSession * observation) noexcept;

std::shared_ptr<PhaseAudit> BeginPhaseAudit(const SecondStageContext &, bool quiet,
    const ObjectiveDomain &, const FitState &, const std::vector<ClusterKey> &,
    std::size_t attempt, std::size_t domain_id) noexcept;

void BeginPhaseObservation(SecondStageObservationSession * observation, SecondStageContext &, bool quiet, const ObjectiveDomain &, const FitState &,
    const std::vector<ClusterKey> &, std::size_t attempt, std::size_t domain_id) noexcept;

class ProductionObservationScope
{
    estimation_audit::Scope m_scope;
public:
    ProductionObservationScope(const SecondStageObservationSession *, std::size_t attempt);
    ProductionObservationScope(const ProductionObservationScope &) = delete;
    ProductionObservationScope & operator=(const ProductionObservationScope &) = delete;
};

enum class BoundaryObservationStage { Endpoint, Correction, Backtracking };
std::string_view BoundaryDiagnosticName(BoundaryAcceptancePolicy, BoundaryObservationStage) noexcept;
std::string_view BoundaryPhaseName(BoundaryAcceptancePolicy, BoundaryObservationStage) noexcept;

// Record indices and history publication tokens never enter numerical results.
class JointCandidateObservation
{
    SecondStageObservationSession * m_session;
    bool m_quiet;
    std::vector<JointCandidateObjectiveDiagnostic> & m_records;
    std::optional<std::size_t> m_current{};
    std::array<std::optional<std::size_t>, 3> m_record_by_stage{};
public:
    JointCandidateObservation(SecondStageObservationSession * session, bool quiet,
        std::vector<JointCandidateObjectiveDiagnostic> & records)
        : m_session(session), m_quiet(quiet), m_records(records) {}
    void Begin(BoundaryObservationStage, std::string_view, double factor, std::size_t round = 0);
    JointCandidateObjectiveDiagnostic * Record();
    void Accept(BoundaryComponentAcceptedSource, BoundaryComponentReconciliationDiagnostic &) noexcept;
};

class BoundaryObservationScope
{
    BoundaryComponentReconciliationDiagnostic m_unobserved{};
    BoundaryComponentReconciliationDiagnostic & m_diagnostic;
    JointCandidateObservation m_trials;
public:
    BoundaryObservationScope(const CandidateSelectionInputs &, const BoundaryReconciliationComponent &);
    BoundaryComponentReconciliationDiagnostic & Diagnostic() { return m_diagnostic; }
    JointCandidateObservation & Trials() { return m_trials; }
    void Finish(const BoundaryComponentDecision &);
};

void BeginCandidateObservation(const CandidateSelectionInputs &, const std::vector<ClusterKey> &);
void ObserveCandidateDecision(const CandidateSelectionInputs &, const ClusterKey &, const CandidateDecisionEvidence &);
void ObserveCandidateSelection(SecondStageObservationSession *, const CandidateSelection &);
void ObserveBoundaryRescue(SecondStageObservationSession *, const ClusterKey &);
void ObserveBoundaryGlobalImprovement(SecondStageObservationSession *, double improvement);
void ObserveBoundaryRejected(SecondStageObservationSession *, const std::vector<ClusterKey> &, bool exhausted);

class LocalSearchObservation
{
    const CandidateSelectionInputs & m_inputs;
    const ClusterKey & m_key;
    const std::vector<SampleRef> & m_samples;
    ObjectiveAttemptDiagnostic * m_diagnostic;
public:
    LocalSearchObservation(const CandidateSelectionInputs &, const ClusterKey &, const std::vector<SampleRef> &);
    void BeginSearch(double radius);
    void Generated();
    void Step(double norm);
    void TrustSkipped();
    void Trial(const CandidateEvaluationOverlay &, const CandidateDecisionEvidence &, bool accepted, bool polish = false);
    void Nonmaterial();
};

class TrustModelTrialObserver
{
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    const CandidateSelectionInputs & inputs;
    const ClusterKey & key;
    const std::vector<SampleRef> & samples;
    TrustModelTrialRecord * record;
    std::optional<std::size_t> final_trial{};
    std::size_t search_pass{ 0 };
    std::size_t trial_number{ 0 };
public:
    TrustModelTrialObserver(const CandidateSelectionInputs &, const ClusterKey &, const std::vector<SampleRef> &);
    void BeginSearch() { ++search_pass; trial_number = 0; }
    void Generated();
    void Invalid();
    void Nonmaterial();
    void TrustSkipped();
    void GuardRejected();
    void Trial(const FitStatePatch &, const CandidateDecisionEvidence &, bool polish, double factor, bool accepted);
    void Finish(bool shrink_trust_region, std::optional<double>, const CandidateDecisionEvidence &);
#else
public:
    TrustModelTrialObserver(const CandidateSelectionInputs &, const ClusterKey &, const std::vector<SampleRef> &) {}
    void BeginSearch() {}
    void Generated() {}
    void Invalid() {}
    void Nonmaterial() {}
    void TrustSkipped() {}
    void GuardRejected() {}
    void Trial(const FitStatePatch &, const CandidateDecisionEvidence &, bool, double, bool) {}
    void Finish(bool, std::optional<double>, const CandidateDecisionEvidence &) {}
#endif
};

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
void BeginTrustModelAudit(SecondStageObservationSession *, const std::vector<ClusterKey> &);
void FinalizeTrustModelAudit(SecondStageObservationSession *, const CandidateSelection &);
void LogTrustModelAudit(const SecondStageObservationSession *, bool, const IterationResult &);
#else
inline void BeginTrustModelAudit(SecondStageObservationSession *, const std::vector<ClusterKey> &) {}
inline void FinalizeTrustModelAudit(SecondStageObservationSession *, const CandidateSelection &) {}
inline void LogTrustModelAudit(const SecondStageObservationSession *, bool, const IterationResult &) {}
#endif

#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE

void ObservePhaseMissing(SecondStageObservationSession * observation, std::string_view, const ClusterKey &, std::string_view) noexcept;
void ObservePhaseState(SecondStageObservationSession * observation, std::string_view, const FitState &, bool probe = false) noexcept;
void ObservePhaseCandidate(SecondStageObservationSession * observation, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView * parent = nullptr, double factor = 1.0, std::string_view disposition = "observed",
    std::string_view reason = "", bool probe = false, bool recertify = true) noexcept;
void ObservePhaseCorrection(SecondStageObservationSession * observation, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, std::string_view, std::string_view, const CandidateSelectionInputs &,
    const std::vector<ClusterKey> &, const ObjectiveBreakdown &) noexcept;
void ObservePhaseLocalPolish(SecondStageObservationSession * observation, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, bool, const CandidateDecisionEvidence &) noexcept;
void ObservePhaseSearchAssembly(SecondStageObservationSession * observation, const FitState &) noexcept;
void ObservePhaseProposal(SecondStageObservationSession * observation, const IterationProposalResult &) noexcept;
void ObservePhaseFinish(SecondStageObservationSession * observation, const FitOptions &, const std::vector<double> &,
    const SuspiciousBlockActivity &, const IterationProposalResult &, const FitState &) noexcept;
void ObservePhaseIntermediate(SecondStageObservationSession * observation, const FittedGaussianSnapshot &) noexcept;
#else

inline void ObservePhaseMissing(SecondStageObservationSession *, std::string_view, const ClusterKey &, std::string_view) noexcept {}
inline void ObservePhaseState(SecondStageObservationSession *, std::string_view, const FitState &, bool = false) noexcept {}
inline void ObservePhaseCandidate(SecondStageObservationSession *, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView * = nullptr, double = 1.0, std::string_view = "observed",
    std::string_view = "", bool = false, bool = true) noexcept {}
inline void ObservePhaseCorrection(SecondStageObservationSession *, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, std::string_view, std::string_view, const CandidateSelectionInputs &,
    const std::vector<ClusterKey> &, const ObjectiveBreakdown &) noexcept {}
inline void ObservePhaseLocalPolish(SecondStageObservationSession *, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, bool, const CandidateDecisionEvidence &) noexcept {}
inline void ObservePhaseSearchAssembly(SecondStageObservationSession *, const FitState &) noexcept {}
inline void ObservePhaseProposal(SecondStageObservationSession *, const IterationProposalResult &) noexcept {}
inline void ObservePhaseFinish(SecondStageObservationSession *, const FitOptions &, const std::vector<double> &,
    const SuspiciousBlockActivity &, const IterationProposalResult &, const FitState &) noexcept {}
inline void ObservePhaseIntermediate(SecondStageObservationSession *, const FittedGaussianSnapshot &) noexcept {}
#endif

} // namespace rhbm_gem::core::detail
