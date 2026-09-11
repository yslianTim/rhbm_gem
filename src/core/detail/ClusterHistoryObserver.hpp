#pragma once

#include "core/detail/ObjectiveEvaluation.hpp"
#include <atomic>
#include <mutex>

namespace rhbm_gem::core::detail {

struct ClusterHistoryDiagnostic
{
    std::optional<ObjectiveBreakdown> best_objective{};
    std::optional<ObjectiveBreakdown> stored_best_objective{};
    bool best_reference_unavailable{ false };
};

struct ClusterHistoryCounters
{
    std::atomic<std::size_t> objective_samples{ 0 };
    void RecordObjectiveSampleEvaluation(std::size_t samples, std::size_t) { objective_samples += samples; }
};

struct BestObjectiveSource
{
    std::string id{}, predecessor_id{};
    ClusterKey key{};
    std::size_t attempt{ 0 }, accepted_iteration{ 0 }, sequence{ 0 }, candidate_number{ 0 };
    std::string_view source{}, reason{};
    std::optional<double> factor{};
    std::optional<ObjectiveBreakdown> before{}, objective{};
    double before_step{ 0.0 }, step{ 0.0 };
    SecondStageModelSnapshot snapshot{};
    std::shared_ptr<const ObjectiveDomain> domain{};
    std::vector<SampleRef> sample_refs{};
};

struct BestObjectiveTraceEnvironment
{
    std::size_t attempt{ 0 }, accepted_iteration{ 0 };
    std::shared_ptr<const ObjectiveDomain> domain{};
    std::mutex mutex{};
    std::map<ClusterKey, std::size_t> sequence_by_key{};
    std::vector<std::shared_ptr<const BestObjectiveSource>> events{};
};

struct ClusterObjectiveState
{
    // Observation only: reevaluated with candidate neighbors and the current domain.
    std::optional<ObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
    FitStatePatch best_parameters{};
    std::shared_ptr<const BestObjectiveSource> best_source{};
    std::string_view best_reset_reason{ "initialize" };
    std::shared_ptr<const BestObjectiveSource> reset_source{};
};

using ClusterObjectiveStateMap = std::map<ClusterKey, ClusterObjectiveState>;
void ReconcileClusterObjectiveState(
    const ObjectiveByKey & previous_objective_by_key,
    const FitState & accepted_state,
    ClusterObjectiveStateMap & state_by_key);

FitStatePatch CaptureClusterParameters(const FitStateView & state, const ClusterKey & key);

std::optional<ObjectiveBreakdown> EvaluateBestObjectiveReference(
    const CandidateEvaluationOverlay & candidate,
    const ClusterKey & key,
    const std::vector<SampleRef> & samples,
    const ObjectiveDomain & domain,
    const ClusterObjectiveState & state,
    ClusterHistoryCounters & counters);

void BeginBestObjectiveTrace(
    SecondStageContext & context, bool quiet_mode, const ObjectiveDomain & domain,
    std::size_t attempt, std::size_t accepted_iteration);
void CaptureBestObjectiveSource(
    const SecondStageContext & context, const ClusterKey & key,
    SecondStageModelSnapshot snapshot, const std::vector<SampleRef> & sample_refs,
    ClusterObjectiveState & state, const std::optional<ObjectiveBreakdown> & before,
    double before_step, std::string_view source, std::string_view reason,
    std::size_t candidate_number = 0, std::optional<double> factor = std::nullopt);
void LogBestObjectivePublication(const SecondStageContext & context, const ClusterObjectiveStateMap & states);
void DiagnoseBestObjectiveComparison(
    JointCandidateObjectiveDiagnostic * record, const CandidateEvaluationOverlay & candidate,
    const ClusterKey & key, const std::vector<SampleRef> & samples,
    const ObjectiveDomain & domain, const ClusterObjectiveState & state);


class ClusterHistoryObserver
{
    std::atomic<bool> m_disabled{ false };
    ClusterObjectiveStateMap m_baseline, m_staged;
    std::map<std::size_t, ClusterObjectiveStateMap> m_boundary;
    std::size_t m_next_observation{ 0 };
    ClusterHistoryCounters m_counters;
    void Disable() noexcept;
public:
    void BeginAttempt(SecondStageContext &, const ObjectiveByKey &, const FitState &,
        const CouplingGraphPartition &, const ObjectiveDomain &, std::size_t, std::size_t) noexcept;
    void ResetPartition(const SecondStageContext &, const CouplingGraphPartition &,
        const ObjectiveDomain &, const FitState &) noexcept;
    void ResetBackground(const SecondStageContext &, const std::shared_ptr<const FrozenBackground> &,
        const CouplingGraphPartition &, const ObjectiveDomain &, const FitState &) noexcept;
    void BeginSearch(const ClusterKey &) noexcept;
    std::shared_ptr<const ClusterHistoryDiagnostic> Local(const CandidateEvaluationOverlay &,
        const ClusterKey &, const std::vector<SampleRef> &, const ObjectiveDomain &,
        std::string_view, bool accepted, const ObjectiveAttemptDiagnostic &) noexcept;
    void BeginBoundary(JointCandidateObjectiveDiagnostic *) noexcept;
    void BoundaryMember(const CandidateEvaluationOverlay &, const ClusterKey &,
        const std::vector<SampleRef> &, const ObjectiveDomain &, bool accepted,
        const ObjectiveAttemptDiagnostic &, JointCandidateObjectiveDiagnostic *) noexcept;
    void AcceptBoundary(std::size_t observation) noexcept;
    void Reject(const ClusterKey &) noexcept;
    void Publish(const SecondStageContext &) noexcept;
    std::optional<ClusterObjectiveState> Snapshot(const ClusterKey &) noexcept;
    std::optional<ClusterObjectiveState> BaselineSnapshot(const ClusterKey &) noexcept;
};

void BeginClusterHistoryObserver(SecondStageContext &, bool quiet) noexcept;

} // namespace rhbm_gem::core::detail
