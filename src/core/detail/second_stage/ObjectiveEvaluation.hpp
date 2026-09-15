#pragma once

#include "core/detail/second_stage/CouplingGraph.hpp"

#include <map>

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
using MemberBestState = std::map<ClusterKey, FitStatePatch>;

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

FitStatePatch OverlayMemberBest(const FitStateView &, const FitStatePatch &);

struct ObjectiveScale
{
    double fit{ 0.0 };
    double tail{ 0.0 };
};

double CalculateClusterAtomWeight(std::size_t cluster_atom_count, std::size_t active_atom_count);

std::optional<ObjectiveBreakdown> BuildObjectiveBreakdown(
    double fit_range_residual_objective,
    double tail_validation_loss,
    double offset_plausibility_penalty);

bool IsBetterAuditObjective(double candidate, double best, ObjectiveTolerance tolerance);

void ValidateObjectiveTolerance(ObjectiveTolerance tolerance);

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

void UpdateMemberBestState(const SecondStageContext &, const ObjectiveDomain &,
    const FitState &, const std::vector<ClusterKey> &, MemberBestState &);

} // namespace rhbm_gem::core::detail
