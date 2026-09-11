#pragma once

#include "core/detail/CouplingGraph.hpp"
#include "core/detail/SuspiciousUpdate.hpp"

#include <map>
#include <mutex>
#include <string>
#include <string_view>

namespace rhbm_gem::core::detail {

class PerformanceCounters;

inline constexpr double kObjectiveRobustLossCutoffMultiplier{ 1.345 };
inline constexpr double kFitRangeWeight{ 1.0 };
inline constexpr double kOffsetPlausibilityPenaltyWeight{ 1.0e-2 };

constexpr double kTailValidationWeight{ 0.25 };

struct ObjectiveBreakdown
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };

    constexpr double GetTailValidationPenalty() const noexcept
    {
        return kTailValidationWeight * tail_validation_loss;
    }

    constexpr double GetTotalObjective() const noexcept
    {
        return fit_range_residual_objective + GetTailValidationPenalty() + offset_plausibility_penalty;
    }
};

struct ObjectiveTolerance
{
    double absolute_tolerance{ 0.0 };
    double relative_tolerance{ 0.0 };
};

inline constexpr ObjectiveTolerance kObjectiveProgressTolerance{ 1.0e-8, 1.0e-3 };
inline constexpr ObjectiveTolerance kObjectiveStrictTolerance{ 1.0e-10, 1.0e-8 };

struct AuditedState
{
    ObjectiveBreakdown objective{};
    FitState state{};
    bool uses_polish{ false };
    std::size_t source_iteration{ 0 };
};

using BestAuditState = std::optional<AuditedState>;

class CandidateEvaluationOverlay
{
    const SecondStageContext & m_context;
    const ResidualBaseline & m_baseline;
    FitStateView m_candidate_state;

public:
    CandidateEvaluationOverlay(
        const SecondStageContext & context,
        const ResidualBaseline & baseline,
        const FitState & base_state,
        const FitStatePatch & candidate_patch);

    std::optional<ResidualSample> operator()(const SampleRef & sample_ref) const;
    const FitStateView & GetState() const { return m_candidate_state; }
    const SecondStageContext & GetContext() const { return m_context; }
    const ResidualBaseline & GetBaseline() const { return m_baseline; }
};

enum class PreObjectiveFailureReason
{
    None,
    InvalidModel,
    NoCandidateWithinTrustRegion
};

struct ObjectiveScale
{
    double fit{ 0.0 };
    double tail{ 0.0 };
};

struct JointCandidateObjectiveDiagnostic
{
    std::string_view source{};
    std::size_t round{ 0 };
    std::size_t candidate_number{ 0 };
    std::optional<double> factor{};
    ClusterKey member_key{};
    std::optional<ObjectiveBreakdown> previous{};
    std::optional<ObjectiveBreakdown> best{};
    std::optional<ObjectiveBreakdown> candidate{};
    std::optional<ObjectiveBreakdown> stored_best{};
    bool best_checked{ false };
    std::string_view outcome{ "accepted" };
    std::string best_source_id{};
    std::vector<std::string> best_comparison_lines{};
};

struct ObjectiveAttemptDiagnostic
{
    std::optional<double> accepted_factor{};
    PreObjectiveFailureReason pre_objective_failure_reason{ PreObjectiveFailureReason::None };
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<ObjectiveScale> scale{};
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
    std::optional<ObjectiveBreakdown> candidate_objective{};
    std::optional<ObjectiveBreakdown> previous_objective{};
    std::optional<ObjectiveBreakdown> best_objective{};
    std::optional<ObjectiveBreakdown> stored_best_objective{};
    bool best_reference_unavailable{ false };
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    bool rejected_by_previous{ false };
    std::size_t trial_count{ 0 };
    std::size_t invalid_trial_count{ 0 };
    std::size_t trust_skipped_trial_count{ 0 };
    std::size_t guard_rejected_trial_count{ 0 };
    std::size_t objective_rejected_trial_count{ 0 };
    std::vector<StabilizationTerminalDiagnostic> terminal_diagnostic_list{};
};

double CalculateClusterAtomWeight(std::size_t cluster_atom_count, std::size_t active_atom_count);

std::optional<ObjectiveBreakdown> BuildObjectiveBreakdown(
    double fit_range_residual_objective,
    double tail_validation_loss,
    double offset_plausibility_penalty);

bool IsBetterAuditObjective(double candidate, double best, ObjectiveTolerance tolerance);

bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance);

struct ObjectiveClusterDomain
{
    // Unique union of the independently selected fit and tail samples.
    std::vector<SampleRef> sample_ref_list{};
    std::vector<SampleRef> fit_sample_ref_list{};
    std::vector<SampleRef> tail_sample_ref_list{};
    std::optional<ObjectiveScale> scale{};
    std::size_t selected_atom_count{ 0 };
};

struct ObjectiveDomain
{
    std::map<ClusterKey, ObjectiveClusterDomain> cluster_by_key{};
    std::vector<ClusterKey> owner_key_by_atom_index{};
    std::vector<std::vector<char>> fit_sample_mask_by_atom{};
    std::vector<std::vector<char>> tail_sample_mask_by_atom{};
    std::size_t unique_sample_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
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
    // Historical value for provenance; acceptance reevaluates best_parameters
    // with the candidate neighbors and current objective domain.
    std::optional<ObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
    FitStatePatch best_parameters{};
    std::shared_ptr<const BestObjectiveSource> best_source{};
    std::string_view best_reset_reason{ "initialize" };
    std::shared_ptr<const BestObjectiveSource> reset_source{};
};

using ClusterObjectiveStateMap = std::map<ClusterKey, ClusterObjectiveState>;
using ObjectiveByKey = std::map<ClusterKey, std::optional<ObjectiveBreakdown>>;

std::size_t CountObjectiveSamples(const std::vector<SampleRef> &, const ObjectiveDomain &);

ObjectiveDomain BuildObjectiveDomain(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<ClusterKey> & cluster_key_list);


std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator);

std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot);

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator);

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot);

bool TryUpdateBestAuditState(
    const FitState & candidate_state,
    bool candidate_uses_polish,
    std::size_t source_iteration,
    const ObjectiveBreakdown & candidate_objective,
    BestAuditState & audit_state);

void ReevaluateBestAuditState(
    const SecondStageContext & context,
    const ObjectiveDomain & domain,
    BestAuditState & audit_state);

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
    PerformanceCounters & performance_counters);

std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const ResidualBaseline & evaluator,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain);

std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const CandidateEvaluationOverlay & evaluator,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain);

double CalculateObjectiveTolerance(
    double reference,
    ObjectiveTolerance tolerance);

bool IsObjectiveDeteriorated(
    double candidate,
    double reference,
    ObjectiveTolerance tolerance);

std::optional<ObjectiveBreakdown> EvaluateObjectiveDelta(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown & baseline,
    PerformanceCounters & performance_counters);

std::optional<ObjectiveBreakdown> EvaluateCombinedObjective(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown * best_objective,
    const ObjectiveBreakdown * previous_objective,
    PerformanceCounters & performance_counters);

} // namespace rhbm_gem::core::detail
