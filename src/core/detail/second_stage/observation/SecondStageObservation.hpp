#pragma once

#include "core/detail/second_stage/observation/SecondStageDiagnostics.hpp"
#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <span>

namespace rhbm_gem::core { struct FitOptions; }
namespace rhbm_gem::core::detail {

class TrustRegionStateSet;
struct CandidateSelectionInputs;
struct CandidateSelection;
struct CandidateCommitResult;
struct BoundaryComponentDecision;
struct IterationResult;
struct IterationProposalResult;
struct QuarantineState;
struct FinalDependencyPolishResult;
enum class SecondStageStopReason;
enum class BoundaryAcceptancePolicy;

// Scheduling uses the raw log level; never replace it with audit enablement.
bool IsDebugLogLevelEnabled();
bool IsSecondStageAuditEnabled(bool quiet) noexcept;

struct SecondStageAuditData
{
    AuditBatch batch{};
    bool polish_attempted{ false }, polish_accepted{ false }, polish_applied{ false };
    FinalPolishResidualSafetyStatus polish_status{ FinalPolishResidualSafetyStatus::NotEvaluated };
    std::optional<ConvergenceAssessment> polish_certificate{};
    std::optional<ObjectiveBreakdown> final_objective{};
    std::array<SelectionAuditDiagnostic, 2> selection{};
    std::optional<ConvergenceAssessment> convergence{};
    std::optional<ObjectiveBreakdown> previous{}, candidate{}, best{};
    std::string_view score_source{ "not_evaluated" };
    std::size_t attempt{ 0 }, objective_revision{ 0 }, recovery_revision{ 0 };
    std::size_t background_revision{ 0 }, partition_revision{ 0 };
    std::size_t entered{ 0 }, retried{ 0 }, released{ 0 }, failed_retry{ 0 };
};

class SecondStageObservationSession
{
    std::unique_ptr<SecondStageAuditData> m_audit;
    std::atomic<bool> m_enabled{ false };
    std::mutex m_mutex;
    std::chrono::steady_clock::time_point m_start{ std::chrono::steady_clock::now() };
public:
    using Writer = void (*)(std::string_view);
    explicit SecondStageObservationSession(bool quiet = true, Writer writer = nullptr) noexcept;
    SecondStageObservationSession(const SecondStageObservationSession &) = delete;
    SecondStageObservationSession & operator=(const SecondStageObservationSession &) = delete;
    IterationDiagnostics iteration{};
    FinalDependencyPolishDiagnostic final_polish{};
    Writer writer;
    bool Enabled() const noexcept { return m_enabled.load(std::memory_order_relaxed); }
    const SecondStageAuditData * Audit() const noexcept { return Enabled() ? m_audit.get() : nullptr; }
    void Disable() noexcept { m_enabled.store(false, std::memory_order_relaxed); }
    void Merge(const AuditBatch &) noexcept;
    void Record(AuditEvent) noexcept;
    void Write(std::string_view) noexcept;
    double ElapsedMilliseconds() const noexcept;
    void BeginAttempt(std::size_t attempt, std::size_t objective_revision, std::size_t recovery_revision,
        bool background_changed, bool partition_changed) noexcept;
    void ObserveProposal(const IterationProposalResult &) noexcept;
    void ObserveCommit(const CandidateSelection &, const CandidateCommitResult &, const TrustRegionStateSet &) noexcept;
    void ObserveQuarantine(const QuarantineState &, const QuarantineState &) noexcept;
    void ObserveRetries(const QuarantineState &) noexcept;
    void ObserveSelectionAudit(bool rescue, bool executed,
        std::string_view result, std::string_view reason, std::size_t removed = 0) noexcept;
    void ObserveGlobalGate(bool rescue, const ObjectiveBreakdown *,
        const std::optional<ObjectiveBreakdown> &, const ObjectiveBreakdown *, bool accepted,
        const ObjectiveProgressGateEvidence & = {}) noexcept;
    void ObserveScoreReferences(const std::optional<ObjectiveBreakdown> & previous, const ObjectiveBreakdown * best) noexcept;
    void ObserveCandidateScoreSource(bool from_selection) noexcept;
    void ObserveScores(const std::optional<ObjectiveBreakdown> & previous,
        const std::optional<ObjectiveBreakdown> & candidate, const ObjectiveBreakdown * best) noexcept;
    void ObserveConvergence(const ConvergenceAssessment &) noexcept;
    void BeginFinalization(const ObjectiveBreakdown * selected_best) noexcept;
    void ObserveFinalPolishAttempt() noexcept;
    void ObserveFinalCertification(const FinalDependencyPolishResult &, FinalPolishResidualSafetyStatus,
        const std::optional<ConvergenceAssessment> &, bool applied) noexcept;
    void ObserveFinalPolishCorrectionFailure(std::span<const std::size_t> key, std::size_t round) noexcept;
    void ObserveFinalPolishComponentFailure(std::span<const std::size_t> key) noexcept;
};

// Each local worker owns five detail slots; only the bounded merge takes a lock.
class LocalSearchObservation
{
    SecondStageObservationSession * m_session;
    const ClusterKey & m_key;
    std::optional<AuditBatch> m_batch;
    double m_radius{ 0.0 };
    std::size_t m_trial{ 0 };
    void Failure(AuditCategory, std::string_view reason) noexcept;
public:
    LocalSearchObservation(const CandidateSelectionInputs &, const ClusterKey &) noexcept;
    ~LocalSearchObservation();
    void BeginSearch(double radius) noexcept { if (m_batch) m_radius = radius; }
    void Generated() noexcept { if (m_batch) ++m_trial; }
    void TrustSkipped() noexcept;
    void Trial(const CandidateDecisionEvidence &, bool accepted, bool polish = false) noexcept;
    void Nonmaterial() noexcept;
    void InvalidCandidate() noexcept;
    void GuardRejected() noexcept;
};

enum class BoundaryObservationStage { Endpoint, Correction, Backtracking };
class JointCandidateObservation
{
    SecondStageObservationSession * m_session;
    std::optional<AuditEvent> m_event;
    bool m_pending{ false };
    std::size_t m_first_atom{ 0 }, m_atom_count{ 0 };
    AuditEvent * Record() noexcept;
    void Begin(BoundaryObservationStage, std::string_view, double factor, std::size_t round) noexcept;
public:
    explicit JointCandidateObservation(SecondStageObservationSession * session, std::span<const std::size_t> key = {}) noexcept;
    ~JointCandidateObservation();
    bool IsRecording() const noexcept;
    void BeginFinalPolish(double factor, std::size_t round) noexcept;
    void Flush() noexcept;
    void BeginBoundary(BoundaryAcceptancePolicy, BoundaryObservationStage, double factor, std::size_t trial = 0) noexcept;
    void Member(const ClusterKey &, bool accepted, const CandidateDecisionEvidence &) noexcept;
    void Gate(const ObjectiveProgressGateEvidence &) noexcept;
    void RejectGlobalObjective() noexcept;
    void RejectStrictImprovement() noexcept;
    void Global(const ObjectiveBreakdown *, const std::optional<ObjectiveBreakdown> &, const ObjectiveBreakdown *) noexcept;
    void RejectInvalidModel() noexcept;
    void RejectSuspicious() noexcept;
    void RejectCorrectionUnavailable() noexcept;
    void RejectBoundaryImprovement(const ObjectiveBreakdown &, const std::optional<ObjectiveBreakdown> &) noexcept;
    void RejectPolishImprovement(const ObjectiveBreakdown &, const std::optional<ObjectiveBreakdown> &) noexcept;
    void RejectMemberSamplesUnavailable(const ClusterKey &) noexcept;
    void RejectPolishMember(const ClusterKey &, const std::optional<ObjectiveBreakdown> & previous,
        const std::optional<ObjectiveBreakdown> & candidate) noexcept;
};

class BoundaryObservationScope
{
    SecondStageObservationSession * m_session;
    BoundaryAcceptancePolicy m_policy;
    const BoundaryReconciliationComponent & m_component;
    JointCandidateObservation m_trials;
public:
    BoundaryObservationScope(const CandidateSelectionInputs &, const BoundaryReconciliationComponent &,
        BoundaryAcceptancePolicy) noexcept;
    JointCandidateObservation & Trials() noexcept { return m_trials; }
    void BeginTrial(BoundaryObservationStage, double factor, std::size_t trial_number = 1) noexcept;
    void Finish(const BoundaryComponentDecision &) noexcept;
};

} // namespace rhbm_gem::core::detail
