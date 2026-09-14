#include "core/detail/second_stage/CandidateTransaction.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/observation/ClusterHistoryObserver.hpp"
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE
#include "core/detail/second_stage/observation/PhaseAudit.hpp"
#endif
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <cmath>
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <algorithm>
#include <ranges>

namespace rhbm_gem::core::detail {

std::shared_ptr<PhaseAudit> BeginPhaseAudit(const SecondStageContext & context, bool quiet,
    const ObjectiveDomain & domain, const FitState & baseline, const std::vector<ClusterKey> & keys,
    std::size_t attempt, std::size_t domain_id) noexcept
{
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE
    if (!quiet && Logger::GetLogLevel() >= LogLevel::Debug)
        try { return std::make_shared<PhaseAudit>(context, domain, baseline, keys, attempt, domain_id); }
        catch (...) { Logger::Log(LogLevel::Debug, "Second-stage phase audit error: capture initialization failed"); }
#else
    (void)context; (void)quiet; (void)domain; (void)baseline; (void)keys; (void)attempt; (void)domain_id;
#endif
    return {};
}

bool IsDebugLogLevelEnabled()
{
    return Logger::GetLogLevel() >= LogLevel::Debug;
}

void RecordJointMemberRejection(
    JointCandidateObjectiveDiagnostic * record,
    const ClusterKey & key,
    const std::optional<ObjectiveBreakdown> & previous,
    const std::optional<ObjectiveBreakdown> & best,
    const std::optional<ObjectiveBreakdown> & candidate,
    bool best_checked)
{
    if (record == nullptr) return;
    record->member_key = key;
    record->previous = previous;
    record->best = best;
    record->candidate = candidate;
    record->best_checked = best_checked;
    if (best_checked && !best)
        record->outcome = "best-reference-unavailable";
    else if (!previous || !candidate)
        record->outcome = "member-objective-unavailable";
    else if (!std::isfinite(previous->GetTotalObjective()) ||
        !std::isfinite(candidate->GetTotalObjective()) ||
        (best_checked && best && !std::isfinite(best->GetTotalObjective())))
        record->outcome = "member-objective-nonfinite";
    else
    {
        const bool previous_failed{ IsObjectiveDeteriorated(candidate->GetTotalObjective(),
            previous->GetTotalObjective(), kObjectiveProgressTolerance) };
        const bool best_failed{ best_checked && best && IsObjectiveDeteriorated(
            candidate->GetTotalObjective(), best->GetTotalObjective(), kObjectiveProgressTolerance) };
        record->outcome = previous_failed ? (best_failed ? "previous+best" : "previous") :
            (best_failed ? "best" : "member-check-failed");
    }
}

JointCandidateObjectiveDiagnostic * BeginJointCandidateDiagnostic(
    bool quiet_mode,
    std::vector<JointCandidateObjectiveDiagnostic> & records,
    std::string_view source,
    std::optional<double> factor,
    std::size_t round)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return nullptr;
    return &records.emplace_back(JointCandidateObjectiveDiagnostic{
        .source = source, .round = round, .candidate_number = records.size() + 1, .factor = factor });
}

void ObserveHistoryPartition(SecondStageObservationSession * observation, const SecondStageContext & context, const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain, const FitState & state) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->ResetPartition(context, partition, domain, state);
}

void ObserveHistoryBackground(SecondStageObservationSession * observation, const SecondStageContext & context,
    const std::shared_ptr<const FrozenBackground> & background, const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain, const FitState & state) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->ResetBackground(context, background, partition, domain, state);
}

void ObserveHistoryAttempt(SecondStageObservationSession * observation, SecondStageContext & context, const ObjectiveByKey & previous,
    const FitState & state, const CouplingGraphPartition & partition, const ObjectiveDomain & domain,
    std::size_t attempt, std::size_t accepted) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->BeginAttempt(context, previous, state, partition, domain, attempt, accepted);
}

void ObserveHistorySearch(SecondStageObservationSession * observation, const ClusterKey & key) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->BeginSearch(key);
}

void ObserveLocalHistory(SecondStageObservationSession * observation, const CandidateEvaluationOverlay & candidate, const ClusterKey & key,
    const std::vector<SampleRef> & samples, const ObjectiveDomain & domain, std::string_view source,
    bool accepted, ObjectiveAttemptDiagnostic & diagnostic) noexcept
{
    if (observation && observation->cluster_history)
        diagnostic.history = observation->cluster_history->Local(
            candidate, key, samples, domain, source, accepted, diagnostic);
}

void ObserveBoundaryHistoryAccepted(SecondStageObservationSession * observation, std::size_t token) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->AcceptBoundary(token);
}

void ObserveHistoryRejected(SecondStageObservationSession * observation, const ClusterKey & key) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->Reject(key);
}

void ObserveHistoryPublication(SecondStageObservationSession * observation) noexcept
{
    if (observation && observation->cluster_history) observation->cluster_history->Publish();
}

void BeginPhaseObservation(SecondStageObservationSession * observation, SecondStageContext & context, bool quiet, const ObjectiveDomain & domain,
    const FitState & state, const std::vector<ClusterKey> & keys, std::size_t attempt, std::size_t domain_id) noexcept
{
    if (observation) observation->phase_audit = BeginPhaseAudit(context, quiet, domain, state, keys, attempt, domain_id);
}

ProductionObservationScope::ProductionObservationScope(const SecondStageObservationSession * observation, std::size_t attempt)
    : m_scope(observation && observation->phase_audit ? attempt : 0, "production") {}

static std::string_view BoundaryDiagnosticName(BoundaryAcceptancePolicy policy, BoundaryObservationStage stage) noexcept
{
    const bool cooperative{ policy == BoundaryAcceptancePolicy::CooperativeRescue };
    switch (stage)
    {
    case BoundaryObservationStage::Endpoint: return cooperative ? "rescue-endpoint" : "endpoint";
    case BoundaryObservationStage::Correction: return cooperative ? "rescue-joint-correction" : "joint-correction";
    case BoundaryObservationStage::Backtracking: return cooperative ? "rescue-backtracking" : "backtracking";
    }
    return {};
}

static std::string_view BoundaryPhaseName(BoundaryAcceptancePolicy policy, BoundaryObservationStage stage) noexcept
{
    const bool cooperative{ policy == BoundaryAcceptancePolicy::CooperativeRescue };
    switch (stage)
    {
    case BoundaryObservationStage::Endpoint: return cooperative ? "rescue-endpoint" : "boundary-endpoint";
    case BoundaryObservationStage::Correction: return cooperative ? "rescue-correction" : "boundary-correction";
    case BoundaryObservationStage::Backtracking: return cooperative ? "rescue-backtracking" : "boundary-backtracking";
    }
    return {};
}

void JointCandidateObservation::Begin(BoundaryObservationStage stage, std::string_view source,
    double factor, std::size_t round)
{
    m_current.reset();
    if (!m_session) return;
    if (BeginJointCandidateDiagnostic(m_quiet, m_records, source, factor, round))
        m_current = m_records.size() - 1;
    m_record_by_stage.at(static_cast<std::size_t>(stage)) = m_current;
}

JointCandidateObjectiveDiagnostic * JointCandidateObservation::Record()
{
    return m_current ? &m_records.at(*m_current) : nullptr;
}

void JointCandidateObservation::BeginBoundary(BoundaryAcceptancePolicy policy, BoundaryObservationStage stage, double factor)
{
    Begin(stage, BoundaryDiagnosticName(policy, stage), factor);
}

void JointCandidateObservation::BeginMembers()
{
    auto * record{ Record() };
    if (m_session && m_session->cluster_history) m_session->cluster_history->BeginBoundary(record);
}

void JointCandidateObservation::Member(const CandidateEvaluationOverlay & candidate, const ClusterKey & key,
    const std::vector<SampleRef> & samples, const ObjectiveDomain & domain,
    BoundaryAcceptancePolicy policy, bool accepted, const CandidateDecisionEvidence & evidence)
{
    auto * record{ Record() };
    if (!accepted)
        RecordJointMemberRejection(record, key, evidence.previous_objective, std::nullopt, evidence.candidate_objective, false);
    else if (record && policy == BoundaryAcceptancePolicy::CooperativeRescue)
    {
        const auto candidate_value{ evidence.candidate_objective->GetTotalObjective() };
        const auto previous_value{ evidence.previous_objective->GetTotalObjective() };
        if (candidate_value > previous_value)
        {
            record->locally_deteriorated_member_count++;
            record->maximum_local_deterioration = std::max(record->maximum_local_deterioration, candidate_value - previous_value);
        }
    }
    // Keep history payload copying inside the original noexcept observation boundary.
    [&]() noexcept
    {
        if (m_session && m_session->cluster_history)
        {
            ObjectiveAttemptDiagnostic diagnostic;
            static_cast<CandidateDecisionEvidence &>(diagnostic) = evidence;
            if (record)
            {
                diagnostic.trial_count = record->candidate_number;
                diagnostic.accepted_factor = record->factor;
            }
            m_session->cluster_history->BoundaryMember(candidate, key, samples, domain, accepted, diagnostic, record);
        }
    }();
}

void JointCandidateObservation::RejectGlobalObjective()
{
    if (auto * record{ Record() }) record->outcome = "members-passed-global-objective-rejected-or-unavailable";
}

void JointCandidateObservation::RejectStrictImprovement(double candidate, double previous)
{
    if (auto * record{ Record() })
        record->outcome = IsObjectiveDeteriorated(candidate, previous, kObjectiveProgressTolerance) ?
            "members-passed-global-objective-rejected-or-unavailable" : "members-passed-strict-improvement-failed";
}

void JointCandidateObservation::Accept(BoundaryComponentAcceptedSource source,
    BoundaryComponentReconciliationDiagnostic & diagnostic) noexcept
{
    if (source == BoundaryComponentAcceptedSource::None) return;
    const auto stage{ source == BoundaryComponentAcceptedSource::Endpoint ? BoundaryObservationStage::Endpoint :
        source == BoundaryComponentAcceptedSource::JointCorrection ? BoundaryObservationStage::Correction :
        BoundaryObservationStage::Backtracking };
    const auto index{ m_record_by_stage.at(static_cast<std::size_t>(stage)) };
    if (!index) return;
    const auto & record{ m_records.at(*index) };
    diagnostic.locally_deteriorated_member_count = record.locally_deteriorated_member_count;
    diagnostic.maximum_local_deterioration = record.maximum_local_deterioration;
    ObserveBoundaryHistoryAccepted(m_session, record.history_observation);
}

BoundaryObservationScope::BoundaryObservationScope(const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component, BoundaryAcceptancePolicy policy,
    const ObjectiveBreakdown * previous, std::size_t rescue_candidate_count)
    : m_inputs(inputs), m_component(component), m_policy(policy),
      m_diagnostic(inputs.observation ?
        inputs.observation->iteration.boundary_reconciliation_diagnostic_list.emplace_back() : m_unobserved),
      m_trials(inputs.observation, inputs.options.quiet_mode, m_diagnostic.objective_diagnostic_list)
{
    m_diagnostic.key_list = component.key_list;
    ClusterKey atoms;
    for (const auto & key : component.key_list) atoms.insert(atoms.end(), key.begin(), key.end());
    std::ranges::sort(atoms);
    atoms.erase(std::ranges::unique(atoms).begin(), atoms.end());
    m_diagnostic.atom_count = atoms.size();
    m_diagnostic.boundary_sample_count = component.boundary_sample_count;
    m_diagnostic.interface_atom_count = component.interface_atom_index_list.size();
    m_diagnostic.shape_active_atom_count = component.halo_atom_index_list.size();
    m_diagnostic.is_rescue_attempt = policy == BoundaryAcceptancePolicy::CooperativeRescue;
    if (m_diagnostic.is_rescue_attempt)
    {
        m_diagnostic.accepted_cluster_count = component.key_list.size() - rescue_candidate_count;
        m_diagnostic.rescue_candidate_cluster_count = rescue_candidate_count;
    }
    if (previous) m_diagnostic.previous_component_objective = previous->GetTotalObjective();
}

void BoundaryObservationScope::BeginTrial(BoundaryObservationStage stage, double factor, std::size_t trial_number)
{
    if (stage == BoundaryObservationStage::Backtracking) m_diagnostic.trial_count = trial_number;
    m_trials.BeginBoundary(m_policy, stage, factor);
}

void BoundaryObservationScope::CandidateEvaluated(BoundaryObservationStage stage,
    const FitStatePatch & patch, const FitStateView & state, double factor,
    const std::optional<BoundaryCandidateEvaluation> & evaluation)
{
    const auto * record{ m_trials.Record() };
    ObservePhaseCandidate(m_inputs.observation, BoundaryPhaseName(m_policy, stage), patch.atom_index_list,
        state, nullptr, factor, evaluation ? "accepted" : "rejected", record ? record->outcome : "",
        false, stage != BoundaryObservationStage::Backtracking);
    if (stage == BoundaryObservationStage::Endpoint && evaluation)
        m_diagnostic.endpoint_component_objective = evaluation->audit_objective.GetTotalObjective();
}

void BoundaryObservationScope::CorrectionActivity(std::size_t shape_count, std::size_t offset_count)
{
    m_diagnostic.shape_active_atom_count = shape_count;
    m_diagnostic.offset_active_atom_count = offset_count;
}

void BoundaryObservationScope::BeginCorrectionSolve(const ObjectiveBreakdown & reference)
{
    m_diagnostic.joint_reference_component_objective = reference.GetTotalObjective();
}

void BoundaryObservationScope::CorrectionSolved(const BoundaryJointCorrectionResult & result)
{
    m_diagnostic.joint_correction_status = result.status;
    m_diagnostic.joint_parameter_count = result.parameter_count;
    if (result.status == BoundaryJointCorrectionStatus::CandidateReady)
    {
        m_diagnostic.joint_damping = result.damping;
        m_diagnostic.maximum_normalized_trust_step = result.maximum_normalized_trust_step;
    }
}

void BoundaryObservationScope::CorrectionEvaluated(const FitStatePatch & patch, const FitStateView & state,
    const FitStateView & endpoint, double damping, const ObjectiveBreakdown & reference,
    const BoundaryCorrectionEvaluation & evaluation)
{
    m_diagnostic.suspicious_candidate_atom_count = evaluation.suspicious_atom_count;
    if (evaluation.suspicious_atom_count != 0)
    {
        ObservePhaseCorrection(m_inputs.observation, BoundaryPhaseName(m_policy, BoundaryObservationStage::Correction),
            patch.atom_index_list, state, endpoint, damping, "rejected", "suspicious", m_inputs, m_component.key_list, reference);
        return;
    }
    if (evaluation.raw_objective)
        m_diagnostic.joint_candidate_component_objective = evaluation.raw_objective->GetTotalObjective();
    auto * record{ m_trials.Record() };
    ObservePhaseCorrection(m_inputs.observation, BoundaryPhaseName(m_policy, BoundaryObservationStage::Correction),
        patch.atom_index_list, state, endpoint, damping, evaluation.accepted ? "accepted" : "rejected",
        evaluation.members && !evaluation.accepted ? "strict-improvement" : (record ? record->outcome : ""),
        m_inputs, m_component.key_list, reference);
    if (!evaluation.accepted && record && evaluation.members)
        record->outcome = "members-passed-strict-improvement-failed";
}

void BoundaryObservationScope::AcceptedCandidate(const BoundaryCandidateEvaluation & evaluation)
{
    m_diagnostic.candidate_component_objective = evaluation.audit_objective.GetTotalObjective();
}

void BoundaryObservationScope::Accept(const BoundaryComponentDecision & decision)
{
    m_trials.Accept(decision.accepted_source, m_diagnostic);
    if (m_diagnostic.is_rescue_attempt)
        m_diagnostic.rescued_cluster_count = m_diagnostic.rescue_candidate_cluster_count;
}

void BoundaryObservationScope::Finish(const BoundaryComponentDecision & decision)
{
    m_diagnostic.accepted_source = decision.accepted_source;
    m_diagnostic.accepted_factor = decision.accepted_factor;
    m_diagnostic.exhausted = decision.exhausted;
    if (m_diagnostic.is_rescue_attempt && decision.accepted_source != BoundaryComponentAcceptedSource::None &&
        m_diagnostic.previous_component_objective && m_diagnostic.candidate_component_objective)
        m_diagnostic.component_improvement = *m_diagnostic.previous_component_objective - *m_diagnostic.candidate_component_objective;
}

void ObserveBoundaryGlobalImprovement(SecondStageObservationSession * session, double improvement)
{
    if (!session) return;
    for (auto & diagnostic : session->iteration.boundary_reconciliation_diagnostic_list)
        if (diagnostic.is_rescue_attempt && diagnostic.accepted_source != BoundaryComponentAcceptedSource::None)
            diagnostic.global_improvement = improvement;
}

void BeginCandidateObservation(const CandidateSelectionInputs & inputs, const std::vector<ClusterKey> & keys)
{
    if (!inputs.observation) return;
    // Allocate entries before workers start; each key then has one writer.
    auto & records{ inputs.observation->iteration.candidate_by_key };
    records.clear();
    if (inputs.options.quiet_mode || !IsDebugLogLevelEnabled()) return;
    for (const auto & key : keys) records.emplace(key, ClusterCandidateDiagnostic{ .key = key });
}

void ObserveCandidateDecision(const CandidateSelectionInputs & inputs, const ClusterKey & key,
    const CandidateDecisionEvidence & evidence)
{
    if (!inputs.observation) return;
    auto & records{ inputs.observation->iteration.candidate_by_key };
    const auto iter{ records.find(key) };
    if (iter != records.end()) static_cast<CandidateDecisionEvidence &>(iter->second.attempt) = evidence;
}

void ObserveCandidateSelection(SecondStageObservationSession * session, const CandidateSelection & selection)
{
    if (!session) return;
    auto & output{ session->iteration };
    const auto append = [&](const auto & decisions, auto & records)
    {
        records.clear();
        for (const auto & decision : decisions)
        {
            const auto iter{ output.candidate_by_key.find(decision.key) };
            if (iter != output.candidate_by_key.end()) records.emplace_back(iter->second);
        }
    };
    append(selection.accepted_cluster_evidence_list, output.accepted_cluster_diagnostic_list);
    append(selection.rejected_cluster_evidence_list, output.rejected_cluster_diagnostic_list);
}

void ObserveBoundaryRescue(SecondStageObservationSession * session, const ClusterKey & key)
{
    if (!session) return;
    const auto iter{ session->iteration.candidate_by_key.find(key) };
    if (iter != session->iteration.candidate_by_key.end()) iter->second.boundary_rescued = true;
}

void ObserveBoundaryRejected(SecondStageObservationSession * session,
    const std::vector<ClusterKey> & keys, bool exhausted)
{
    if (!session) return;
    auto & records{ session->iteration.boundary_reconciliation_diagnostic_list };
    const auto iter{ std::ranges::find(records | std::views::reverse, keys,
        &BoundaryComponentReconciliationDiagnostic::key_list) };
    if (iter == records.rend()) return;
    iter->accepted_source = BoundaryComponentAcceptedSource::None;
    iter->accepted_factor.reset();
    iter->exhausted = exhausted;
}

LocalSearchObservation::LocalSearchObservation(const CandidateSelectionInputs & inputs,
    const ClusterKey & key, const std::vector<SampleRef> & samples)
    : m_inputs(inputs), m_key(key), m_samples(samples), m_diagnostic(nullptr)
{
    if (!inputs.observation) return;
    const auto iter{ inputs.observation->iteration.candidate_by_key.find(key) };
    if (iter != inputs.observation->iteration.candidate_by_key.end()) m_diagnostic = &iter->second.attempt;
}

void LocalSearchObservation::BeginSearch(double radius)
{
    ObserveHistorySearch(m_inputs.observation, m_key);
    if (!m_diagnostic) return;
    *m_diagnostic = ObjectiveAttemptDiagnostic{};
    m_diagnostic->trust_region_radius = radius;
}
void LocalSearchObservation::Generated() { if (m_diagnostic) ++m_diagnostic->trial_count; }
void LocalSearchObservation::Step(double norm)
{
    if (!m_diagnostic) return;
    m_diagnostic->pre_objective_attempted_step_norm = norm;
    m_diagnostic->trust_region_step_norm = norm;
}
void LocalSearchObservation::TrustSkipped() { if (m_diagnostic) ++m_diagnostic->trust_skipped_trial_count; }
void LocalSearchObservation::Nonmaterial() { if (m_diagnostic) m_diagnostic->trust_region_step_norm = 0.0; }
void LocalSearchObservation::Trial(const CandidateEvaluationOverlay & candidate,
    const CandidateDecisionEvidence & evidence, bool accepted, bool polish)
{
    if (!m_diagnostic) return;
    ObjectiveAttemptDiagnostic polish_diagnostic;
    auto & diagnostic{ polish ? polish_diagnostic : *m_diagnostic };
    static_cast<CandidateDecisionEvidence &>(diagnostic) = evidence;
    diagnostic.pre_objective_attempted_step_norm.reset();
    diagnostic.scale.reset();
    const auto domain_iter{ m_inputs.objective_domain.cluster_by_key.find(m_key) };
    if (domain_iter != m_inputs.objective_domain.cluster_by_key.end())
    {
        diagnostic.fit_sample_count = domain_iter->second.fit_sample_ref_list.size();
        diagnostic.tail_sample_count = domain_iter->second.tail_sample_ref_list.size();
        diagnostic.scale = domain_iter->second.scale;
    }
    ObserveLocalHistory(m_inputs.observation, candidate, m_key, m_samples, m_inputs.objective_domain,
        polish ? "local-polish" : "local-candidate", accepted, diagnostic);
}

} // namespace rhbm_gem::core::detail
