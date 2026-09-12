#include "core/detail/SecondStageObservation.hpp"
#include "core/detail/CandidateEvaluation.hpp"
#include "core/detail/ClusterHistoryObserver.hpp"
#include "core/detail/PhaseAudit.hpp"

namespace rhbm_gem::core::detail {

void ObserveHistoryPartition(const SecondStageContext & context, const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain, const FitState & state) noexcept
{
    if (context.cluster_history) context.cluster_history->ResetPartition(context, partition, domain, state);
}

void ObserveHistoryBackground(const SecondStageContext & context,
    const std::shared_ptr<const FrozenBackground> & background, const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain, const FitState & state) noexcept
{
    if (context.cluster_history) context.cluster_history->ResetBackground(context, background, partition, domain, state);
}

void ObserveHistoryAttempt(SecondStageContext & context, const ObjectiveByKey & previous,
    const FitState & state, const CouplingGraphPartition & partition, const ObjectiveDomain & domain,
    std::size_t attempt, std::size_t accepted) noexcept
{
    if (context.cluster_history) context.cluster_history->BeginAttempt(context, previous, state, partition, domain, attempt, accepted);
}

void ObserveHistorySearch(const SecondStageContext & context, const ClusterKey & key) noexcept
{
    if (context.cluster_history) context.cluster_history->BeginSearch(key);
}

void ObserveLocalHistory(const CandidateEvaluationOverlay & candidate, const ClusterKey & key,
    const std::vector<SampleRef> & samples, const ObjectiveDomain & domain, std::string_view source,
    bool accepted, ObjectiveAttemptDiagnostic & diagnostic) noexcept
{
    if (candidate.GetContext().cluster_history)
        diagnostic.history = candidate.GetContext().cluster_history->Local(
            candidate, key, samples, domain, source, accepted, diagnostic);
}

void ObserveBoundaryHistory(const SecondStageContext & context, JointCandidateObjectiveDiagnostic * record) noexcept
{
    if (context.cluster_history) context.cluster_history->BeginBoundary(record);
}

void ObserveBoundaryMemberHistory(const CandidateEvaluationOverlay & candidate, const ClusterKey & key,
    const std::vector<SampleRef> & samples, const ObjectiveDomain & domain, bool accepted,
    const ObjectiveAttemptDiagnostic & diagnostic, JointCandidateObjectiveDiagnostic * record) noexcept
{
    if (candidate.GetContext().cluster_history)
        candidate.GetContext().cluster_history->BoundaryMember(candidate, key, samples, domain, accepted, diagnostic, record);
}

void ObserveBoundaryHistoryAccepted(const SecondStageContext & context, std::size_t observation) noexcept
{
    if (context.cluster_history) context.cluster_history->AcceptBoundary(observation);
}

void ObserveHistoryRejected(const SecondStageContext & context, const ClusterKey & key) noexcept
{
    if (context.cluster_history) context.cluster_history->Reject(key);
}

void ObserveHistoryPublication(const SecondStageContext & context) noexcept
{
    if (context.cluster_history) context.cluster_history->Publish(context);
}

void BeginPhaseObservation(SecondStageContext & context, bool quiet, const ObjectiveDomain & domain,
    const FitState & state, const std::vector<ClusterKey> & keys, std::size_t attempt, std::size_t domain_id) noexcept
{
    context.phase_audit = BeginPhaseAudit(context, quiet, domain, state, keys, attempt, domain_id);
}

ProductionObservationScope::ProductionObservationScope(const SecondStageContext & context, std::size_t attempt)
    : m_scope(context.phase_audit ? attempt : 0, "production") {}

std::string_view BoundaryDiagnosticName(BoundaryAcceptancePolicy policy, BoundaryObservationStage stage) noexcept
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

std::string_view BoundaryPhaseName(BoundaryAcceptancePolicy policy, BoundaryObservationStage stage) noexcept
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

} // namespace rhbm_gem::core::detail
