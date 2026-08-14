#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "core/detail/CouplingGraph.hpp"
#include "core/detail/CandidateEvaluationOverlay.hpp"
#include "core/detail/PerformanceCounters.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/RobustScale.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/TransformedChange.hpp"

namespace rhbm_gem::core::detail {

struct ObjectiveBreakdown
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double tail_validation_penalty{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    double total_objective{ 0.0 };
};

struct ObjectiveTolerance
{
    double absolute_tolerance{ 0.0 };
    double relative_tolerance{ 0.0 };
};

constexpr double kObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kFitRangeWeight{ 1.0 };
constexpr double kTailValidationWeight{ 0.25 };
constexpr double kOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kOffsetPeakRatioMax{ 1.0 };
constexpr double kObjectiveRobustLossCutoffMultiplier{ 1.345 };
constexpr ObjectiveTolerance
kObjectiveStrictTolerance{ 1.0e-10, 1.0e-8 };
constexpr ObjectiveTolerance
kObjectiveProgressTolerance{ 1.0e-8, 1.0e-3 };

struct AuditedState
{
    ObjectiveBreakdown objective{};
    FitState state{};
    PolishProvenance polish_provenance{};
    std::optional<std::size_t> accepted_iteration{};
};

struct BestAuditState
{
    std::optional<AuditedState> best{};
};

enum class PreObjectiveFailureReason
{
    None,
    InvalidModel,
    PreviousSharedOffsetProjectionOutsideTrustRegion,
    NoCandidateWithinTrustRegion
};

struct ObjectiveAttemptDiagnostic
{
    double effective_damping{ 1.0 };
    bool is_invalid_model{ false };
    PreObjectiveFailureReason pre_objective_failure_reason{
        PreObjectiveFailureReason::None
    };
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<double> fit_scale{};
    std::optional<double> tail_scale{};
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
    std::optional<ObjectiveBreakdown> candidate_objective{};
    std::optional<ObjectiveBreakdown> previous_objective{};
    std::optional<ObjectiveBreakdown> best_objective{};
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    std::size_t backtracking_trial_count{ 0 };
    std::optional<double> accepted_backtracking_factor{};
    bool backtracking_exhausted{ false };
};

inline double CalculateClusterAtomWeight(
    std::size_t cluster_atom_count,
    std::size_t active_atom_count)
{
    if (cluster_atom_count == 0 || active_atom_count == 0 ||
        cluster_atom_count > active_atom_count)
    {
        throw std::invalid_argument(
            "Local fitting cluster atom counts are invalid.");
    }
    return static_cast<double>(cluster_atom_count) /
        static_cast<double>(active_atom_count);
}

inline void ValidateObjectiveTolerance(
    const ObjectiveTolerance & tolerance)
{
    if (!std::isfinite(tolerance.absolute_tolerance) ||
        tolerance.absolute_tolerance < 0.0 ||
        !std::isfinite(tolerance.relative_tolerance) ||
        tolerance.relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerances must be finite and "
            "non-negative.");
    }
}

inline double CalculateObjectiveTolerance(
    double reference,
    const ObjectiveTolerance & tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(reference))
    {
        throw std::invalid_argument(
            "Local fitting audit objective reference must be finite.");
    }
    return tolerance.absolute_tolerance +
        tolerance.relative_tolerance * std::abs(reference);
}

inline std::optional<ObjectiveBreakdown>
BuildObjectiveBreakdown(
    double fit_range_residual_objective,
    double tail_validation_loss,
    double offset_plausibility_penalty)
{
    if (!std::isfinite(fit_range_residual_objective) ||
        !std::isfinite(tail_validation_loss) ||
        !std::isfinite(offset_plausibility_penalty))
    {
        return std::nullopt;
    }

    ObjectiveBreakdown breakdown;
    breakdown.fit_range_residual_objective = fit_range_residual_objective;
    breakdown.tail_validation_loss = tail_validation_loss;
    breakdown.tail_validation_penalty =
        kTailValidationWeight * tail_validation_loss;
    breakdown.offset_plausibility_penalty = offset_plausibility_penalty;
    breakdown.total_objective =
        breakdown.fit_range_residual_objective +
        breakdown.tail_validation_penalty +
        breakdown.offset_plausibility_penalty;
    if (!std::isfinite(breakdown.total_objective)) return std::nullopt;
    return breakdown;
}

inline bool IsBetterAuditObjective(
    double candidate,
    double best,
    const ObjectiveTolerance & tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate)) return false;
    if (!std::isfinite(best)) return true;
    return candidate < best -
        CalculateObjectiveTolerance(best, tolerance);
}

inline bool IsAuditObjectiveAcceptableForProgress(
    const std::optional<double> & candidate,
    const std::optional<double> & previous,
    const std::optional<double> & best,
    const ObjectiveTolerance & tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!candidate.has_value() || !previous.has_value() ||
        !std::isfinite(*candidate) || !std::isfinite(*previous))
    {
        return false;
    }
    const auto is_deteriorated = [&](double reference)
    {
        if (!std::isfinite(reference)) return true;
        return *candidate > reference +
            CalculateObjectiveTolerance(reference, tolerance);
    };
    return !is_deteriorated(*previous) &&
        (!best.has_value() || !is_deteriorated(*best));
}

struct ObjectiveScale
{
    double fit{ 0.0 };
    std::optional<double> tail{};
};

struct ObjectiveClusterDomain
{
    std::vector<ObjectiveSampleRef> fit_sample_ref_list{};
    std::vector<ObjectiveSampleRef> tail_sample_ref_list{};
    std::optional<ObjectiveScale> scale{};
};

struct ObjectiveDomain
{
    std::map<ClusterKey, ObjectiveClusterDomain> cluster_by_key{};
    std::vector<ClusterKey> owner_key_by_atom_index{};
    std::vector<std::vector<char>> fit_sample_mask_by_atom{};
    std::size_t active_atom_count{ 0 };
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
};

struct ClusterObjectiveState
{
    std::optional<ObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
};

using ClusterObjectiveStateMap = std::map<ClusterKey, ClusterObjectiveState>;
using ObjectiveByKey =
    std::map<ClusterKey, std::optional<ObjectiveBreakdown>>;

struct CombinedObjectiveCheck
{
    bool accepted{ false };
    std::optional<ObjectiveBreakdown> candidate_objective{};
};

struct CombinedCandidateObjectiveCheck
{
    bool guard_required{ false };
    bool accepted{ true };
    std::optional<ObjectiveBreakdown> previous_objective{};
    std::optional<ObjectiveBreakdown> candidate_objective{};
};

inline FitStatePatch BuildStatePatch(const FitState & state, ClusterKey atom_index_list)
{
    std::sort(atom_index_list.begin(), atom_index_list.end());
    atom_index_list.erase(
        std::unique(atom_index_list.begin(), atom_index_list.end()),
        atom_index_list.end());
    FitStatePatch patch;
    patch.atom_index_list = std::move(atom_index_list);
    patch.mdpde_list.reserve(patch.atom_index_list.size());
    for (const auto atom_index : patch.atom_index_list)
    {
        patch.mdpde_list.emplace_back(state.at(atom_index).mdpde);
    }
    return patch;
}


inline std::optional<double> BuildFixedObjectiveScale(
    const std::vector<double> & residual_list,
    const std::vector<double> & adjusted_response_list)
{
    if (residual_list.empty() ||
        residual_list.size() != adjusted_response_list.size())
    {
        return std::nullopt;
    }
    const auto scale{
        std::max({
            CalculateMedianAbsoluteDeviationScale(residual_list),
            kObjectiveResidualScaleFloorRatio *
                CalculateMedianAbsoluteDeviationScale(adjusted_response_list),
            kRobustScaleMin
        })
    };
    return numeric_validation::IsFinitePositive(scale) ?
        std::optional<double>{ scale } : std::nullopt;
}


inline ObjectiveDomain BuildObjectiveDomain(
    const SecondStageContext & context,
    const FitState & initial_state,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitOptions & options)
{
    ObjectiveDomain domain;
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, initial_state)
    };
    domain.owner_key_by_atom_index.resize(context.size());
    domain.fit_sample_mask_by_atom.resize(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        domain.fit_sample_mask_by_atom.at(atom_index).resize(
            context.at(atom_index).raw_sampling_entries.size(), 0);
    }
    for (const auto & key : cluster_key_list)
    {
        auto & cluster_domain{ domain.cluster_by_key[key] };
        domain.active_atom_count += key.size();
        std::vector<double> fit_residual_list;
        std::vector<double> fit_response_list;
        std::vector<double> tail_residual_list;
        std::vector<double> tail_response_list;
        for (const auto atom_index : key)
        {
            domain.owner_key_by_atom_index.at(atom_index) = key;
            const auto & raw_sampling_entries{
                context.at(atom_index).raw_sampling_entries
            };
            for (std::size_t sample_index = 0; sample_index < raw_sampling_entries.size(); sample_index++)
            {
                const ObjectiveSampleRef sample_ref{
                    atom_index,
                    sample_index
                };
                const auto residual_sample{
                    EvaluateResidualSample(context, initial_state, sample_ref, model_snapshot)
                };
                const auto distance{
                    static_cast<double>(raw_sampling_entries.at(sample_index).point.distance)
                };
                const auto is_fit_range{
                    distance >= options.distance_min && distance <= options.distance_max
                };
                domain.fit_sample_mask_by_atom.at(atom_index).at(sample_index) =
                    is_fit_range ? 1 : 0;
                auto & sample_ref_list{ is_fit_range ?
                    cluster_domain.fit_sample_ref_list : cluster_domain.tail_sample_ref_list
                };
                sample_ref_list.emplace_back(sample_ref);
                if (!residual_sample.has_value()) continue;
                auto & residual_list{
                    is_fit_range ? fit_residual_list : tail_residual_list
                };
                auto & response_list{
                    is_fit_range ? fit_response_list : tail_response_list
                };
                residual_list.emplace_back(residual_sample->residual);
                response_list.emplace_back(residual_sample->adjusted_response);
            }
        }
        domain.fit_sample_count += cluster_domain.fit_sample_ref_list.size();
        domain.tail_sample_count += cluster_domain.tail_sample_ref_list.size();
        const auto fit_scale{
            BuildFixedObjectiveScale(fit_residual_list, fit_response_list)
        };
        if (!fit_scale.has_value() ||
            fit_residual_list.size() != cluster_domain.fit_sample_ref_list.size())
        {
            continue;
        }
        ObjectiveScale scale;
        scale.fit = *fit_scale;
        if (!cluster_domain.tail_sample_ref_list.empty())
        {
            scale.tail = BuildFixedObjectiveScale(tail_residual_list, tail_response_list);
            if (!scale.tail.has_value() ||
                tail_residual_list.size() != cluster_domain.tail_sample_ref_list.size())
            {
                continue;
            }
        }
        cluster_domain.scale = scale;
    }
    return domain;
}


template <typename State, typename ResidualEvaluator>
inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContributionImpl(
    const State & state,
    const ClusterKey & changed_key,
    const std::vector<ObjectiveSampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    const ResidualEvaluator & residual_evaluator,
    bool include_offset_penalty = true)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    ObjectiveBreakdown breakdown;
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & owner_key{
            domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
        const auto owner_iter{ domain.cluster_by_key.find(owner_key) };
        if (owner_iter == domain.cluster_by_key.end() || !owner_iter->second.scale.has_value())
        {
            return std::nullopt;
        }
        const auto residual_sample{ residual_evaluator(sample_ref) };
        if (!residual_sample.has_value()) return std::nullopt;
        const auto is_fit_range{
            domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        const auto sample_count{ is_fit_range ?
            owner_iter->second.fit_sample_ref_list.size() : owner_iter->second.tail_sample_ref_list.size()
        };
        if (sample_count == 0) return std::nullopt;
        const auto scale{ is_fit_range ?
            std::optional<double>{ owner_iter->second.scale->fit } : owner_iter->second.scale->tail
        };
        if (!scale.has_value()) return std::nullopt;
        const auto loss{
            algorithm::CalculateCauchyLoss(
                residual_sample->residual / *scale,
                kObjectiveRobustLossCutoffMultiplier)
        };
        const auto coefficient{
            CalculateClusterAtomWeight(owner_key.size(), domain.active_atom_count) /
            static_cast<double>(sample_count)
        };
        if (is_fit_range)
        {
            breakdown.fit_range_residual_objective += kFitRangeWeight * coefficient * loss;
        }
        else
        {
            breakdown.tail_validation_loss += coefficient * loss;
        }
    }
    if (include_offset_penalty)
    {
        for (const auto atom_index : changed_key)
        {
            const auto owner_iter{
                domain.cluster_by_key.find(domain.owner_key_by_atom_index.at(atom_index))
            };
            if (owner_iter == domain.cluster_by_key.end() || !owner_iter->second.scale.has_value())
            {
                return std::nullopt;
            }
            const auto & model{ GetFitModel(state, atom_index) };
            if (!IsValidSecondStageGaussianModel(model)) return std::nullopt;
            const auto peak_signal{ model.SignalAtDistance(0.0) };
            const auto offset_peak{ model.GetOffset() * model.OffsetBasisAtDistance(0.0) };
            if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak))
            {
                return std::nullopt;
            }
            const auto offset_ratio{
                std::abs(offset_peak) /
                std::max({
                    std::abs(peak_signal),
                    owner_iter->second.scale->fit, kRobustScaleMin
                })
            };
            const auto offset_excess{
                std::max(0.0, offset_ratio - kOffsetPeakRatioMax)
            };
            breakdown.offset_plausibility_penalty +=
                kOffsetPlausibilityPenaltyWeight * offset_excess * offset_excess /
                static_cast<double>(domain.active_atom_count);
        }
    }
    return BuildObjectiveBreakdown(
        breakdown.fit_range_residual_objective,
        breakdown.tail_validation_loss,
        breakdown.offset_plausibility_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const ClusterKey & changed_key,
    const std::vector<ObjectiveSampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    bool include_offset_penalty = true)
{
    return EvaluateObjectiveContributionImpl(
        model_snapshot.selected,
        changed_key,
        sample_ref_list,
        domain,
        [&](const ObjectiveSampleRef & sample_ref)
        {
            return EvaluateResidualSample(
                context,
                model_snapshot.selected,
                sample_ref,
                model_snapshot);
        },
        include_offset_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const ResidualBaseline & baseline,
    const ClusterKey & changed_key,
    const std::vector<ObjectiveSampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    bool include_offset_penalty = true)
{
    return EvaluateObjectiveContributionImpl(
        baseline.model_snapshot.selected,
        changed_key,
        sample_ref_list,
        domain,
        [&](const ObjectiveSampleRef & sample_ref)
        {
            return baseline.sample_list.at(sample_ref.atom_index).at(
                sample_ref.sample_index);
        },
        include_offset_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const CandidateEvaluationOverlay & overlay,
    const ClusterKey & changed_key,
    const std::vector<ObjectiveSampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    bool include_offset_penalty = true)
{
    return EvaluateObjectiveContributionImpl(
        overlay.GetCandidateState(),
        changed_key,
        sample_ref_list,
        domain,
        [&](const ObjectiveSampleRef & sample_ref)
        {
            return overlay.Evaluate(sample_ref);
        },
        include_offset_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const SecondStageContext & context,
    const ObjectiveDomain & domain,
    const SecondStageModelSnapshot & model_snapshot)
{
    ObjectiveBreakdown total;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        const auto fit_contribution{
            EvaluateObjectiveContribution(
                context,
                model_snapshot,
                key,
                cluster_domain.fit_sample_ref_list,
                domain,
                false)
        };
        if (!fit_contribution.has_value()) return std::nullopt;
        const auto tail_contribution{
            EvaluateObjectiveContribution(
                context,
                model_snapshot,
                key,
                cluster_domain.tail_sample_ref_list,
                domain,
                false)
        };
        if (!tail_contribution.has_value()) return std::nullopt;
        const std::vector<ObjectiveSampleRef> empty_sample_ref_list;
        const auto offset_contribution{
            EvaluateObjectiveContribution(
                context,
                model_snapshot,
                key,
                empty_sample_ref_list,
                domain)
        };
        if (!offset_contribution.has_value()) return std::nullopt;
        total.fit_range_residual_objective += fit_contribution->fit_range_residual_objective;
        total.fit_range_residual_objective += tail_contribution->fit_range_residual_objective;
        total.tail_validation_loss += fit_contribution->tail_validation_loss;
        total.tail_validation_loss += tail_contribution->tail_validation_loss;
        total.tail_validation_penalty += fit_contribution->tail_validation_penalty;
        total.tail_validation_penalty += tail_contribution->tail_validation_penalty;
        total.offset_plausibility_penalty += offset_contribution->offset_plausibility_penalty;
    }
    total.total_objective =
        total.fit_range_residual_objective +
        total.tail_validation_penalty +
        total.offset_plausibility_penalty;
    return std::isfinite(total.total_objective) ?
        std::optional<ObjectiveBreakdown>{ total } : std::nullopt;
}

inline std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const ResidualBaseline & baseline)
{
    ObjectiveBreakdown total;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        std::vector<ObjectiveSampleRef> sample_ref_list{
            cluster_domain.fit_sample_ref_list
        };
        sample_ref_list.insert(
            sample_ref_list.end(),
            cluster_domain.tail_sample_ref_list.begin(),
            cluster_domain.tail_sample_ref_list.end());
        const auto contribution{
            EvaluateObjectiveContribution(
                baseline,
                key,
                sample_ref_list,
                domain)
        };
        if (!contribution.has_value()) return std::nullopt;
        total.fit_range_residual_objective += contribution->fit_range_residual_objective;
        total.tail_validation_loss += contribution->tail_validation_loss;
        total.offset_plausibility_penalty += contribution->offset_plausibility_penalty;
    }
    return BuildObjectiveBreakdown(
        total.fit_range_residual_objective,
        total.tail_validation_loss,
        total.offset_plausibility_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveDelta(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & changed_key,
    const std::vector<ObjectiveSampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const std::optional<ObjectiveBreakdown> & baseline,
    PerformanceCounters & performance_counters)
{
    if (!baseline.has_value()) return std::nullopt;
    const auto unique_sample_count{
        domain.fit_sample_count + domain.tail_sample_count
    };
    performance_counters.RecordObjectiveSampleEvaluation(
        affected_sample_ref_list.size(),
        unique_sample_count);
    const auto candidate_changed{
        EvaluateObjectiveContribution(
            candidate_overlay,
            changed_key,
            affected_sample_ref_list,
            domain)
    };
    const auto previous_changed{
        EvaluateObjectiveContribution(
            candidate_overlay.GetBaseline(),
            changed_key,
            affected_sample_ref_list,
            domain)
    };
    if (!candidate_changed.has_value() || !previous_changed.has_value())
    {
        return std::nullopt;
    }
    return BuildObjectiveBreakdown(
        baseline->fit_range_residual_objective +
            candidate_changed->fit_range_residual_objective -
            previous_changed->fit_range_residual_objective,
        baseline->tail_validation_loss +
            candidate_changed->tail_validation_loss -
            previous_changed->tail_validation_loss,
        baseline->offset_plausibility_penalty +
            candidate_changed->offset_plausibility_penalty -
            previous_changed->offset_plausibility_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const SecondStageContext & context,
    const FitState & state,
    const ObjectiveDomain & domain)
{
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, state)
    };
    return EvaluateAuditObjective(context, domain, model_snapshot);
}

inline ObjectiveByKey BuildObjectiveByKey(
    const SecondStageContext & context,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const SecondStageModelSnapshot & model_snapshot)
{
    ObjectiveByKey objective_by_key;
    for (const auto & [key, sample_ref_list] : partition.sample_id_list_by_key)
    {
        objective_by_key.emplace(
            key,
            EvaluateObjectiveContribution(
                context,
                model_snapshot,
                key,
                sample_ref_list,
                domain));
    }
    return objective_by_key;
}

inline ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const ResidualBaseline & baseline)
{
    ObjectiveByKey objective_by_key;
    for (const auto & [key, sample_ref_list] : partition.sample_id_list_by_key)
    {
        objective_by_key.emplace(
            key,
            EvaluateObjectiveContribution(
                baseline,
                key,
                sample_ref_list,
                domain));
    }
    return objective_by_key;
}

inline CombinedObjectiveCheck EvaluateCombinedObjective(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & changed_key,
    const std::vector<ObjectiveSampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const BestAuditState & audit_state,
    const std::optional<ObjectiveBreakdown> & previous_objective,
    PerformanceCounters & performance_counters)
{
    const auto candidate_objective{
        EvaluateObjectiveDelta(
            candidate_overlay,
            changed_key,
            affected_sample_ref_list,
            domain,
            previous_objective,
            performance_counters)
    };
    return CombinedObjectiveCheck{
        IsAuditObjectiveAcceptableForProgress(
            candidate_objective.has_value() ?
                std::optional<double>{ candidate_objective->total_objective } :
                std::nullopt,
            previous_objective.has_value() ?
                std::optional<double>{ previous_objective->total_objective } :
                std::nullopt,
            audit_state.best.has_value() ?
                std::optional<double>{ audit_state.best->objective.total_objective } :
                std::nullopt,
            kObjectiveProgressTolerance),
        candidate_objective
    };
}

inline CombinedCandidateObjectiveCheck EvaluateCombinedCandidateObjective(
    const SecondStageContext & context,
    const ResidualBaseline & baseline,
    const CouplingGraphPartition & partition,
    const FitState & previous_state,
    const FitState & candidate_state,
    const std::vector<ClusterKey> & accepted_key_list,
    const ObjectiveDomain & objective_domain,
    const BestAuditState & best_audit_state,
    PerformanceCounters & performance_counters)
{
    CombinedCandidateObjectiveCheck result;
    result.guard_required = partition.boundary_sample_count > 0 && !accepted_key_list.empty();
    if (!result.guard_required) return result;

    result.previous_objective = EvaluateAuditObjective(
        objective_domain,
        baseline);

    ClusterKey changed_atom_index_list;
    for (const auto & key : accepted_key_list)
    {
        changed_atom_index_list.insert(
            changed_atom_index_list.end(),
            key.begin(),
            key.end());
    }
    const auto combined_patch{
        BuildStatePatch(
            candidate_state,
            std::move(changed_atom_index_list))
    };
    const FitStateView combined_state_view{
        previous_state,
        combined_patch
    };
    const CandidateEvaluationOverlay combined_overlay{
        context,
        baseline,
        combined_state_view
    };
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(partition, accepted_key_list)
    };
    const auto combined_check{
        EvaluateCombinedObjective(
            combined_overlay,
            combined_patch.atom_index_list,
            affected_sample_ref_list,
            objective_domain,
            best_audit_state,
            result.previous_objective,
            performance_counters)
    };
    result.accepted = combined_check.accepted;
    result.candidate_objective = combined_check.candidate_objective;
    return result;
}


inline bool TryUpdateBestAuditState(
    const SecondStageContext & context,
    const FitState & candidate_state,
    const PolishProvenance & candidate_polish_provenance,
    const std::optional<std::size_t> & accepted_iteration,
    const ObjectiveDomain & domain,
    BestAuditState & audit_state,
    const std::optional<ObjectiveBreakdown> & precomputed_objective = std::nullopt)
{
    const auto candidate_objective{ precomputed_objective.has_value() ?
        precomputed_objective :
        EvaluateAuditObjective(context, candidate_state, domain) };
    if (!candidate_objective.has_value()) return false;
    if (audit_state.best.has_value() &&
        !IsBetterAuditObjective(
            candidate_objective->total_objective,
            audit_state.best->objective.total_objective,
            kObjectiveStrictTolerance))
    {
        return false;
    }
    audit_state.best = AuditedState{
        *candidate_objective,
        candidate_state,
        candidate_polish_provenance,
        accepted_iteration
    };
    return true;
}

inline BestAuditState BuildInitialBestAuditState(
    const SecondStageContext & context,
    const FitState & initial_state,
    const PolishProvenance & initial_polish_provenance,
    const std::optional<std::size_t> & accepted_iteration,
    const ObjectiveDomain & domain)
{
    BestAuditState audit_state;
    static_cast<void>(TryUpdateBestAuditState(
        context,
        initial_state,
        initial_polish_provenance,
        accepted_iteration,
        domain,
        audit_state));
    return audit_state;
}

inline void ResetBestAuditAfterObjectiveDomainChange(
    const SecondStageContext & context,
    const FitState & validated_state,
    const PolishProvenance & validated_polish_provenance,
    std::size_t accepted_iteration,
    const ObjectiveDomain & domain,
    BestAuditState & audit_state)
{
    audit_state = BuildInitialBestAuditState(
        context,
        validated_state,
        validated_polish_provenance,
        accepted_iteration,
        domain);
}

inline ClusterObjectiveState BuildInitialClusterObjectiveState(
    const ObjectiveByKey & previous_objective_by_key,
    const ClusterKey & key)
{
    ClusterObjectiveState state;
    const auto objective_iter{ previous_objective_by_key.find(key) };
    if (objective_iter != previous_objective_by_key.end())
    {
        state.best_objective = objective_iter->second;
    }
    return state;
}

inline void ReconcileClusterObjectiveState(
    const std::vector<ClusterKey> & cluster_key_list,
    const ObjectiveByKey & previous_objective_by_key,
    ClusterObjectiveStateMap & state_by_key)
{
    ClusterObjectiveStateMap next_state_by_key;
    for (const auto & key : cluster_key_list)
    {
        auto state_iter{ state_by_key.find(key) };
        if (state_iter != state_by_key.end())
        {
            next_state_by_key.emplace(key, std::move(state_iter->second));
            continue;
        }
        next_state_by_key.emplace(
            key,
            BuildInitialClusterObjectiveState(previous_objective_by_key, key));
    }
    state_by_key = std::move(next_state_by_key);
}

inline bool TryCommitClusterCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & key,
    const std::vector<ObjectiveSampleRef> & objective_sample_ref_list,
    const std::optional<ObjectiveBreakdown> & previous_objective,
    bool requires_strict_improvement,
    const ObjectiveDomain & domain,
    ClusterObjectiveStateMap & cluster_objective_state,
    ObjectiveAttemptDiagnostic & diagnostic,
    PerformanceCounters & performance_counters)
{
    const auto unique_sample_count{
        domain.fit_sample_count + domain.tail_sample_count
    };
    performance_counters.RecordObjectiveSampleEvaluation(
        objective_sample_ref_list.size(),
        unique_sample_count);
    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            candidate_overlay.GetCandidateState(),
            candidate_overlay.GetBaselineState(),
            key)
    };
    const auto maximum_transformed_change{
        GetMaximumTransformedPercentileChange(transformed_change_summary)
    };
    const auto domain_iter{ domain.cluster_by_key.find(key) };
    if (domain_iter != domain.cluster_by_key.end())
    {
        diagnostic.fit_sample_count = domain_iter->second.fit_sample_ref_list.size();
        diagnostic.tail_sample_count = domain_iter->second.tail_sample_ref_list.size();
        if (domain_iter->second.scale.has_value())
        {
            diagnostic.fit_scale = domain_iter->second.scale->fit;
            diagnostic.tail_scale = domain_iter->second.scale->tail;
        }
    }
    auto & state{ cluster_objective_state.at(key) };
    diagnostic.candidate_objective =
        EvaluateObjectiveContribution(
            candidate_overlay,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.previous_objective = previous_objective;
    diagnostic.best_objective = state.best_objective;

    const auto is_objective_deteriorated = [](
        const std::optional<ObjectiveBreakdown> & candidate,
        const std::optional<double> & reference)
    {
        if (!reference.has_value()) return false;
        if (!candidate.has_value()) return true;
        return candidate->total_objective > *reference +
            CalculateObjectiveTolerance(*reference, kObjectiveProgressTolerance);
    };
    if (!diagnostic.candidate_objective.has_value() || !diagnostic.previous_objective.has_value())
    {
        return false;
    }
    diagnostic.rejected_by_previous = is_objective_deteriorated(
        diagnostic.candidate_objective,
        diagnostic.previous_objective->total_objective);
    diagnostic.rejected_by_best = is_objective_deteriorated(
        diagnostic.candidate_objective,
        state.best_objective.has_value() ?
            std::optional<double>{ state.best_objective->total_objective } : std::nullopt);
    if (diagnostic.rejected_by_previous || diagnostic.rejected_by_best)
    {
        return false;
    }
    if (requires_strict_improvement &&
        (!diagnostic.candidate_objective.has_value() ||
            !diagnostic.previous_objective.has_value() ||
            !IsBetterAuditObjective(
                diagnostic.candidate_objective->total_objective,
                diagnostic.previous_objective->total_objective,
                kObjectiveStrictTolerance)))
    {
        return false;
    }

    auto is_better_than_best{ !state.best_objective.has_value() };
    if (state.best_objective.has_value())
    {
        const auto candidate_objective_value{
            diagnostic.candidate_objective->total_objective
        };
        const auto best_objective_value{ state.best_objective->total_objective };
        if (IsBetterAuditObjective(
                candidate_objective_value,
                best_objective_value,
                kObjectiveStrictTolerance))
        {
            is_better_than_best = true;
        }
        else if (IsBetterAuditObjective(
                     best_objective_value,
                     candidate_objective_value,
                     kObjectiveStrictTolerance))
        {
            is_better_than_best = false;
        }
        else
        {
            is_better_than_best = maximum_transformed_change < state.best_maximum_transformed_change;
        }
    }
    if (is_better_than_best)
    {
        state.best_objective = diagnostic.candidate_objective;
        state.best_maximum_transformed_change = maximum_transformed_change;
    }
    return true;
}

} // namespace rhbm_gem::core::detail
