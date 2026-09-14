#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include "support/SecondStageNumericalProbe.hpp"
#else
#define RHBM_TEST_AUDIT_FAULT(point) ((void)0)
#endif
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/IterationResult.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <cmath>
#include <tuple>

namespace rhbm_gem::core::detail {
namespace {
bool Before(const AuditEvent & a, const AuditEvent & b) noexcept
{
    return std::tie(a.stage, a.first_atom, a.atom_count, a.trial, a.category, a.reason) <
        std::tie(b.stage, b.first_atom, b.atom_count, b.trial, b.category, b.reason);
}
void KeepDetail(AuditBatch & batch, const AuditEvent & event) noexcept
{
    auto end = batch.details.begin() + static_cast<std::ptrdiff_t>(batch.detail_count);
    auto position = std::lower_bound(batch.details.begin(), end, event, Before);
    if (position == batch.details.end()) return;
    if (batch.detail_count < kAuditDetailLimit) ++batch.detail_count;
    end = batch.details.begin() + static_cast<std::ptrdiff_t>(batch.detail_count);
    std::move_backward(position, end - 1, end);
    *position = event;
}
void SetKey(AuditEvent & event, std::span<const std::size_t> key) noexcept
{
    event.first_atom = key.empty() ? 0 : key.front();
    event.atom_count = key.size();
}
void DefaultWriter(std::string_view text) { Logger::Log(LogLevel::Debug, std::string(text)); }
}

void AuditBatch::Add(const AuditEvent & event) noexcept
{
    auto & count = stages[static_cast<std::size_t>(event.stage)];
    ++count.total;
    if (event.outcome == "accepted") ++count.accepted;
    else if (event.outcome == "rejected") ++count.rejected;
    else ++count.skipped;
    if (event.category == AuditCategory::None) return;
    ++categories[static_cast<std::size_t>(event.category)];
    ++abnormal_count;
    KeepDetail(*this, event);
}
void AuditBatch::Merge(const AuditBatch & other) noexcept
{
    for (std::size_t i = 0; i < stages.size(); ++i)
    {
        stages[i].total += other.stages[i].total;
        stages[i].accepted += other.stages[i].accepted;
        stages[i].rejected += other.stages[i].rejected;
        stages[i].skipped += other.stages[i].skipped;
    }
    for (std::size_t i = 0; i < categories.size(); ++i) categories[i] += other.categories[i];
    abnormal_count += other.abnormal_count;
    for (std::size_t i = 0; i < other.detail_count; ++i) KeepDetail(*this, other.details[i]);
}

bool IsDebugLogLevelEnabled() { return Logger::GetLogLevel() >= LogLevel::Debug; }
bool IsSecondStageAuditEnabled(bool quiet) noexcept
{
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT
    return !quiet && IsDebugLogLevelEnabled();
#else
    (void)quiet;
    return false;
#endif
}
SecondStageObservationSession::SecondStageObservationSession(bool quiet, Writer output) noexcept
    : writer(output ? output : DefaultWriter)
{
    if (!IsSecondStageAuditEnabled(quiet)) return;
    try { RHBM_TEST_AUDIT_FAULT(Allocation); m_audit = std::make_unique<SecondStageAuditData>(); m_enabled = true; }
    catch (...) { Disable(); }
}
void SecondStageObservationSession::Merge(const AuditBatch & batch) noexcept
{
    if (!Enabled()) return;
    try { RHBM_TEST_AUDIT_FAULT(Collection); std::lock_guard lock(m_mutex); if (Enabled()) m_audit->batch.Merge(batch); }
    catch (...) { Disable(); }
}
void SecondStageObservationSession::Record(AuditEvent event) noexcept
{
    if (!Enabled()) return;
    AuditBatch batch;
    batch.Add(event);
    Merge(batch);
}
void SecondStageObservationSession::Write(std::string_view text) noexcept
{
    if (!Enabled()) return;
    try { RHBM_TEST_AUDIT_FAULT(Writer); writer(text); } catch (...) { Disable(); }
}
double SecondStageObservationSession::ElapsedMilliseconds() const noexcept
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - m_start).count();
}
void SecondStageObservationSession::BeginAttempt(std::size_t attempt, std::size_t objective_revision,
    std::size_t recovery_revision, bool background_changed, bool partition_changed) noexcept
{
    auto * data = Audit();
    if (!data) return;
    const auto background = data->background_revision + static_cast<std::size_t>(background_changed || attempt == 1);
    const auto partition = data->partition_revision + static_cast<std::size_t>(partition_changed || attempt == 1);
    *data = SecondStageAuditData{};
    data->selection[1].reason = "no-accepted-rescue";
    data->attempt = attempt;
    data->objective_revision = objective_revision;
    data->recovery_revision = recovery_revision;
    data->background_revision = background;
    data->partition_revision = partition;
    if (partition_changed) Record({.stage=AuditStage::Partition, .category=AuditCategory::Partition,
        .outcome="changed", .reason="partition-applied"});
}

void ObserveProposal(SecondStageObservationSession * session, const IterationProposalResult & proposal) noexcept
{
    if (!session || !session->Enabled()) return;
    for (const auto & [key, health] : proposal.health_by_key)
    {
        AuditEvent event{ .stage=AuditStage::Proposal };
        SetKey(event, key);
        if (!health.IsSolverQualified())
        {
            event.category = AuditCategory::Solver;
            event.outcome = "rejected";
            event.reason = IsJointOffsetSolveHardFailure(health.joint_offset_status) ? "solver-hard-failure" : "solver-unqualified";
        }
        session->Record(event);
    }
}
void ObserveCommit(SecondStageObservationSession * session, const CandidateSelection & selection,
    const CandidateCommitResult & result, const TrustRegionStateSet & radii) noexcept
{
    if (!session || !session->Enabled()) return;
    if (!result.accepted) { session->Audit()->candidate=session->Audit()->previous; session->Audit()->score_source="restored_previous"; }
    for (const auto & key : result.accepted_key_list)
    {
        AuditEvent event{ .stage=AuditStage::Commit }; SetKey(event,key); session->Record(event);
    }
    for (const auto & key : result.rejected_key_list)
    {
        AuditEvent event{ .stage=AuditStage::Commit, .category=AuditCategory::Rejected, .outcome="rejected", .reason="final-selection-rejected" };
        SetKey(event,key); session->Record(event);
    }
    for (const auto & key : result.trust_region_update.changed_key_list)
    {
        AuditEvent event{ .stage=AuditStage::Commit, .category=AuditCategory::Shrink, .outcome="changed", .reason="radius-shrunk" };
        event.radius = radii.GetRadius(key); SetKey(event,key); session->Record(event);
    }
    for (const auto & key : selection.exhausted_key_list)
    {
        AuditEvent event{ .stage=AuditStage::Commit, .category=AuditCategory::Exhausted, .outcome="rejected", .reason="search-exhausted" };
        SetKey(event,key); session->Record(event);
    }
}
void ObserveQuarantine(SecondStageObservationSession * session, const QuarantineState & before,
    const QuarantineState & after) noexcept
{
    if (!session || !session->Enabled()) return;
    auto * data = session->Audit();
    data->entered += after.entered_target_count - before.entered_target_count;
    data->released += after.released_target_count - before.released_target_count;
    data->failed_retry += after.failed_retry_count - before.failed_retry_count;
    for (const auto & [target, state] : after.state_by_target)
    {
        const auto previous = before.state_by_target.find(target);
        if (state.lifecycle == QuarantineLifecycle::Frozen &&
            (previous == before.state_by_target.end() || previous->second.lifecycle != QuarantineLifecycle::Frozen))
        {
            AuditEvent event{ .stage=AuditStage::Quarantine, .category=AuditCategory::Enter, .outcome="changed", .reason="quarantine-enter" };
            SetKey(event,target.atom_index_list); session->Record(event);
        }
    }
    for (const auto & [target, state] : before.state_by_target)
        if (state.lifecycle == QuarantineLifecycle::Frozen && !after.state_by_target.contains(target))
        {
            AuditEvent event{ .stage=AuditStage::Quarantine, .category=AuditCategory::Release, .outcome="changed", .reason="quarantine-release" };
            SetKey(event,target.atom_index_list); session->Record(event);
        }
}
void ObserveSelectionAudit(SecondStageObservationSession * session, bool rescue, bool executed,
    std::string_view result, std::string_view reason, std::size_t removed) noexcept
{
    if (!session || !session->Enabled()) return;
    auto & entry = session->Audit()->selection[rescue ? 1 : 0];
    entry.executed = executed; entry.result = result; entry.reason = reason; entry.removed_clusters = removed;
    AuditEvent event{ .stage=rescue ? AuditStage::SelectionRescue : AuditStage::SelectionOrdinary,
        .category=removed ? AuditCategory::Salvage : (result == "unavailable" ? AuditCategory::Unavailable :
            (result == "rejected" || result == "empty_after_salvage" ? AuditCategory::Rejected : AuditCategory::None)),
        .outcome=result == "passed" ? "accepted" : (executed ? "rejected" : "skipped"), .reason=reason, .scope="global" };
    session->Record(event);
}
void ObserveGlobalGate(SecondStageObservationSession * session, bool rescue, const ObjectiveBreakdown * previous,
    const std::optional<ObjectiveBreakdown> & candidate, const ObjectiveBreakdown * best, bool accepted,
    const ObjectiveProgressGateEvidence & gate) noexcept
{
    if (!session || !session->Enabled()) return;
    auto & entry = session->Audit()->selection[rescue ? 1 : 0];
    ++entry.evaluations;
    entry.previous = previous ? std::optional{*previous} : std::nullopt;
    entry.candidate = candidate; entry.best = best ? std::optional{*best} : std::nullopt;
    AuditEvent event{ .stage=rescue ? AuditStage::SelectionRescue : AuditStage::SelectionOrdinary,
        .category=accepted ? AuditCategory::None : (candidate && previous ? AuditCategory::Rejected : AuditCategory::Unavailable),
        .trial=entry.evaluations, .outcome=accepted ? "accepted" : "rejected",
        .reason=accepted ? "" : gate.reason,
        .scope="global", .previous=entry.previous, .candidate=candidate, .best=entry.best,
        .previous_checked=gate.previous_checked, .best_checked=gate.best_checked };
    session->Record(event);
}

LocalSearchObservation::LocalSearchObservation(const CandidateSelectionInputs & inputs, const ClusterKey & key) noexcept : m_session(inputs.observation), m_key(key)
{
    if (m_session && m_session->Enabled()) m_batch.emplace();
}
LocalSearchObservation::~LocalSearchObservation() { if (m_batch) m_session->Merge(*m_batch); }
void LocalSearchObservation::Failure(AuditCategory category, std::string_view reason) noexcept
{
    if (!m_batch) return;
    AuditEvent event{ .category=category, .trial=m_trial, .outcome="rejected", .reason=reason, .radius=m_radius };
    SetKey(event,m_key); m_batch->Add(event);
}
void LocalSearchObservation::TrustSkipped() noexcept { Failure(AuditCategory::Trust,"outside-radius"); }
void LocalSearchObservation::Nonmaterial() noexcept
{
    if (!m_batch) return;
    AuditEvent event{ .trial=m_trial, .outcome="skipped", .reason="nonmaterial" }; SetKey(event,m_key); m_batch->Add(event);
}
void LocalSearchObservation::Trial(const CandidateDecisionEvidence & evidence,
    bool accepted, bool polish) noexcept
{
    if (!m_batch) return;
    AuditEvent event{ .stage=polish ? AuditStage::LocalPolish : AuditStage::LocalSearch,
        .category=accepted ? AuditCategory::None : (evidence.candidate_objective && evidence.previous_objective ? AuditCategory::Rejected : AuditCategory::Unavailable),
        .trial=m_trial, .outcome=accepted ? "accepted" : "rejected",
        .reason=accepted ? "" : (evidence.rejected_by_previous ? "previous-gate" :
            (evidence.candidate_objective && evidence.previous_objective ? "strict-improvement" : "objective-unavailable")),
        .reference=polish ? "local_search_candidate" : "iteration_previous",
        .previous=evidence.previous_objective, .candidate=evidence.candidate_objective,
        .factor=evidence.accepted_factor, .radius=m_radius,
        .previous_checked=evidence.previous_objective.has_value() && evidence.candidate_objective.has_value() };
    SetKey(event,m_key); m_batch->Add(event);
}

JointCandidateObservation::JointCandidateObservation(SecondStageObservationSession * session,
    std::span<const std::size_t> key) noexcept : m_session(session) { m_first_atom = key.empty() ? 0 : key.front(); m_atom_count = key.size(); }
JointCandidateObservation::~JointCandidateObservation() { Flush(); }
AuditEvent * JointCandidateObservation::Record() noexcept { return m_session && m_session->Enabled() && m_pending ? &*m_event : nullptr; }
void JointCandidateObservation::Flush() noexcept
{
    if (auto * record = Record()) m_session->Record(*record);
    m_pending = false;
}
void JointCandidateObservation::Begin(BoundaryObservationStage stage, std::string_view source, double factor, std::size_t round) noexcept
{
    Flush();
    if (!m_session || !m_session->Enabled()) return;
    m_event.emplace(); m_event->first_atom = m_first_atom; m_event->atom_count = m_atom_count;
    m_event->stage = source == "final-polish" ? AuditStage::FinalPolish :
        static_cast<AuditStage>(static_cast<int>(AuditStage::BoundaryEndpoint) + static_cast<int>(stage));
    m_event->reference = source == "final-polish" ? "final_polish_endpoint" : "iteration_previous";
    m_event->scope = "global"; m_event->factor = factor; m_event->trial = round;
    m_pending = true;
}
void JointCandidateObservation::BeginBoundary(BoundaryAcceptancePolicy policy, BoundaryObservationStage stage, double factor) noexcept
{
    Begin(stage,"boundary",factor);
    if (m_pending && policy == BoundaryAcceptancePolicy::CooperativeRescue)
        m_event->stage = static_cast<AuditStage>(static_cast<int>(AuditStage::RescueEndpoint) + static_cast<int>(stage));
}
void RecordJointMemberRejection(AuditEvent * record, const ClusterKey & key,
    const std::optional<ObjectiveBreakdown> & previous, const std::optional<ObjectiveBreakdown> & best,
    const std::optional<ObjectiveBreakdown> & candidate, bool best_checked) noexcept
{
    if (!record) return;
    SetKey(*record,key); record->scope="cluster"; record->previous=previous; record->candidate=candidate; record->best=best;
    record->previous_checked=previous.has_value() && candidate.has_value(); record->best_checked=best_checked;
    record->outcome="rejected";
    record->category=previous && candidate ? AuditCategory::Rejected : AuditCategory::Unavailable;
    record->reason=previous && candidate ? "member-gate" : "member-objective-unavailable";
}
void JointCandidateObservation::Member(const ClusterKey & key, bool accepted,
    const CandidateDecisionEvidence & evidence) noexcept
{
    if (!accepted) RecordJointMemberRejection(Record(),key,evidence.previous_objective,{},evidence.candidate_objective,false);
}
void JointCandidateObservation::Global(const ObjectiveBreakdown * previous,
    const std::optional<ObjectiveBreakdown> & candidate, const ObjectiveBreakdown * best) noexcept
{
    if (auto * record = Record())
    {
        record->previous=previous ? std::optional{*previous} : std::nullopt;
        record->candidate=candidate; record->best=best ? std::optional{*best} : std::nullopt;
        record->previous_checked=false; record->best_checked=false;
    }
}
void JointCandidateObservation::Gate(const ObjectiveProgressGateEvidence & gate) noexcept
{
    if (auto * record=Record()) { record->previous_checked=gate.previous_checked; record->best_checked=gate.best_checked; record->reason=gate.reason; }
}
void JointCandidateObservation::RejectGlobalObjective() noexcept
{
    if (auto * record=Record())
    {
        record->outcome="rejected"; record->category=record->candidate && record->previous ? AuditCategory::Rejected : AuditCategory::Unavailable;
        if (record->reason.empty()) record->reason="global-objective-unavailable";
    }
}
void JointCandidateObservation::RejectStrictImprovement() noexcept
{
    if (auto * record=Record()) { record->outcome="rejected"; record->category=AuditCategory::Rejected; record->reason="strict-improvement"; }
}
BoundaryObservationScope::BoundaryObservationScope(const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component, BoundaryAcceptancePolicy policy) noexcept
    : m_session(inputs.observation), m_policy(policy), m_component(component), m_trials(inputs.observation, component.key_list.empty() ? std::span<const std::size_t>{} : std::span<const std::size_t>{component.key_list.front()}) {}
void BoundaryObservationScope::BeginTrial(BoundaryObservationStage stage, double factor, std::size_t trial) noexcept
{
    m_trials.BeginBoundary(m_policy,stage,factor);
    if (auto * record=m_trials.Record()) record->trial=trial;
}
void BoundaryObservationScope::Finish(const BoundaryComponentDecision & decision) noexcept
{
    m_trials.Flush();
    if (!m_session || !m_session->Enabled()) return;
    AuditEvent event{ .stage=m_policy == BoundaryAcceptancePolicy::CooperativeRescue ? AuditStage::RescueEndpoint : AuditStage::BoundaryEndpoint,
        .category=decision.accepted_source == BoundaryComponentAcceptedSource::None ? (decision.exhausted ? AuditCategory::Exhausted : AuditCategory::Rejected) :
            (m_policy == BoundaryAcceptancePolicy::CooperativeRescue ? AuditCategory::Rescue : AuditCategory::None),
        .outcome=decision.accepted_source == BoundaryComponentAcceptedSource::None ? "rejected" : "accepted",
        .reason=decision.accepted_source == BoundaryComponentAcceptedSource::None ? (decision.exhausted ? "component-exhausted" : "component-rejected") : "component-selected",
        .scope="component", .factor=decision.accepted_factor };
    SetKey(event, m_component.key_list.empty() ? std::span<const std::size_t>{} : std::span<const std::size_t>{m_component.key_list.front()}); m_session->Record(event);
}
} // namespace rhbm_gem::core::detail
