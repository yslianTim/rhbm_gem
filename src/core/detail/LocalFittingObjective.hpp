#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "core/detail/CouplingGraph.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingCandidateEvaluationOverlay.hpp"
#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/LocalFittingObjectiveAttemptDiagnostic.hpp"
#include "core/detail/LocalFittingPerformanceCounters.hpp"
#include "core/detail/LocalFittingResidualEvaluation.hpp"
#include "core/detail/LocalFittingRobustScale.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/LocalFittingStateView.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/LocalFittingTrustRegion.hpp"

namespace rhbm_gem::core::detail {

constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kLocalFittingFitRangeWeight{ 1.0 };
constexpr double kLocalFittingTailValidationWeight{ 0.25 };
constexpr double kLocalFittingOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingOffsetPeakRatioMax{ 1.0 };
constexpr double kLocalFittingObjectiveRobustLossCutoffMultiplier{ 1.345 };
constexpr LocalFittingObjectiveTolerance
kLocalFittingObjectiveStrictTolerance{ 1.0e-10, 1.0e-8 };
constexpr LocalFittingObjectiveTolerance
kLocalFittingObjectiveProgressTolerance{ 1.0e-8, 1.0e-3 };

struct LocalFittingObjectiveScale
{
    double fit{ 0.0 };
    std::optional<double> tail{};
};

struct LocalFittingObjectiveClusterDomain
{
    std::vector<LocalFittingObjectiveSampleRef> fit_sample_ref_list{};
    std::vector<LocalFittingObjectiveSampleRef> tail_sample_ref_list{};
    std::optional<LocalFittingObjectiveScale> scale{};
};

struct LocalFittingObjectiveDomain
{
    std::map<LocalFittingClusterKey, LocalFittingObjectiveClusterDomain> cluster_by_key{};
    std::vector<LocalFittingClusterKey> owner_key_by_atom_index{};
    std::vector<std::vector<char>> fit_sample_mask_by_atom{};
    std::size_t active_atom_count{ 0 };
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
};

struct LocalFittingClusterObjectiveState
{
    std::optional<LocalFittingObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
};

using LocalFittingClusterObjectiveStateMap = std::map<LocalFittingClusterKey, LocalFittingClusterObjectiveState>;
using LocalFittingObjectiveByKey =
    std::map<LocalFittingClusterKey, std::optional<LocalFittingObjectiveBreakdown>>;

struct LocalFittingCombinedObjectiveCheck
{
    bool accepted{ false };
    std::optional<LocalFittingObjectiveBreakdown> candidate_objective{};
};


inline std::optional<double> BuildFixedLocalFittingObjectiveScale(
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
            CalculateLocalFittingMedianAbsoluteDeviationScale(residual_list),
            kLocalFittingObjectiveResidualScaleFloorRatio *
                CalculateLocalFittingMedianAbsoluteDeviationScale(adjusted_response_list),
            kLocalFittingRobustScaleMin
        })
    };
    return numeric_validation::IsFinitePositive(scale) ?
        std::optional<double>{ scale } : std::nullopt;
}


inline LocalFittingObjectiveDomain BuildLocalFittingObjectiveDomain(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state,
    const CouplingGraphPartition & partition,
    const FitOptions & options)
{
    LocalFittingObjectiveDomain domain;
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(initial_state))
    };
    domain.owner_key_by_atom_index.resize(context.size());
    domain.fit_sample_mask_by_atom.resize(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        domain.fit_sample_mask_by_atom.at(atom_index).resize(
            context.at(atom_index).raw_sampling_entries.size(),
            0);
    }
    for (const auto & [key, affected_sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        static_cast<void>(affected_sample_ref_list);
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
            for (std::size_t sample_index = 0;
                sample_index < raw_sampling_entries.size();
                sample_index++)
            {
                const LocalFittingObjectiveSampleRef sample_ref{
                    atom_index,
                    sample_index
                };
                const auto residual_sample{
                    EvaluateLocalFittingResidualSample(
                        context,
                        initial_state,
                        sample_ref,
                        model_snapshot)
                };
                const auto distance{
                    static_cast<double>(
                        raw_sampling_entries.at(sample_index).point.distance)
                };
                const auto is_fit_range{
                    distance >= options.distance_min &&
                    distance <= options.distance_max
                };
                domain.fit_sample_mask_by_atom.at(atom_index).at(sample_index) =
                    is_fit_range ? 1 : 0;
                auto & sample_ref_list{
                    is_fit_range ?
                        cluster_domain.fit_sample_ref_list :
                        cluster_domain.tail_sample_ref_list
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
            BuildFixedLocalFittingObjectiveScale(
                fit_residual_list,
                fit_response_list)
        };
        if (!fit_scale.has_value() ||
            fit_residual_list.size() != cluster_domain.fit_sample_ref_list.size())
        {
            continue;
        }
        LocalFittingObjectiveScale scale;
        scale.fit = *fit_scale;
        if (!cluster_domain.tail_sample_ref_list.empty())
        {
            scale.tail = BuildFixedLocalFittingObjectiveScale(
                tail_residual_list,
                tail_response_list);
            if (!scale.tail.has_value() ||
                tail_residual_list.size() !=
                    cluster_domain.tail_sample_ref_list.size())
            {
                continue;
            }
        }
        cluster_domain.scale = scale;
    }
    return domain;
}


template <typename State, typename ResidualEvaluator>
inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingObjectiveContributionImpl(
    const SecondStageLocalFittingContext & context,
    const State & state,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const LocalFittingObjectiveDomain & domain,
    const ResidualEvaluator & residual_evaluator,
    bool include_offset_penalty = true)
{
    static_cast<void>(context);
    if (domain.active_atom_count == 0) return std::nullopt;
    LocalFittingObjectiveBreakdown breakdown;
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & owner_key{
            domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
        const auto owner_iter{ domain.cluster_by_key.find(owner_key) };
        if (owner_iter == domain.cluster_by_key.end() ||
            !owner_iter->second.scale.has_value())
        {
            return std::nullopt;
        }
        const auto residual_sample{
            residual_evaluator(sample_ref)
        };
        if (!residual_sample.has_value()) return std::nullopt;
        const auto is_fit_range{
            domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(
                sample_ref.sample_index) != 0
        };
        const auto sample_count{
            is_fit_range ?
                owner_iter->second.fit_sample_ref_list.size() :
                owner_iter->second.tail_sample_ref_list.size()
        };
        if (sample_count == 0) return std::nullopt;
        const auto scale{
            is_fit_range ?
                std::optional<double>{ owner_iter->second.scale->fit } :
                owner_iter->second.scale->tail
        };
        if (!scale.has_value()) return std::nullopt;
        const auto loss{
            algorithm::CalculateCauchyLoss(
                residual_sample->residual / *scale,
                kLocalFittingObjectiveRobustLossCutoffMultiplier)
        };
        const auto coefficient{
            CalculateLocalFittingClusterAtomWeight(
                owner_key.size(),
                domain.active_atom_count) /
            static_cast<double>(sample_count)
        };
        if (is_fit_range)
        {
            breakdown.fit_range_residual_objective +=
                kLocalFittingFitRangeWeight * coefficient * loss;
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
                domain.cluster_by_key.find(
                    domain.owner_key_by_atom_index.at(atom_index))
            };
            if (owner_iter == domain.cluster_by_key.end() ||
                !owner_iter->second.scale.has_value())
            {
                return std::nullopt;
            }
            const auto & model{ GetLocalFittingModel(state, atom_index) };
            if (!IsValidSecondStageGaussianModel(model)) return std::nullopt;
            const auto peak_signal{ model.SignalAtDistance(0.0) };
            const auto offset_peak{
                model.GetOffset() * model.OffsetBasisAtDistance(0.0)
            };
            if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak))
            {
                return std::nullopt;
            }
            const auto offset_ratio{
                std::abs(offset_peak) /
                std::max({
                    std::abs(peak_signal),
                    owner_iter->second.scale->fit,
            kLocalFittingRobustScaleMin
                })
            };
            const auto offset_excess{
                std::max(0.0, offset_ratio - kLocalFittingOffsetPeakRatioMax)
            };
            breakdown.offset_plausibility_penalty +=
                kLocalFittingOffsetPlausibilityPenaltyWeight *
                offset_excess * offset_excess /
                static_cast<double>(domain.active_atom_count);
        }
    }
    return BuildLocalFittingObjectiveBreakdown(
        breakdown.fit_range_residual_objective,
        breakdown.tail_validation_loss,
        breakdown.offset_plausibility_penalty,
        kLocalFittingTailValidationWeight);
}

template <typename State>
inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingObjectiveContribution(
    const SecondStageLocalFittingContext & context,
    const State & state,
    const SecondStageModelSnapshot & model_snapshot,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const LocalFittingObjectiveDomain & domain,
    bool include_offset_penalty = true)
{
    return EvaluateLocalFittingObjectiveContributionImpl(
        context,
        state,
        changed_key,
        sample_ref_list,
        domain,
        [&](const LocalFittingObjectiveSampleRef & sample_ref)
        {
            return EvaluateLocalFittingResidualSample(
                context,
                state,
                sample_ref,
                model_snapshot);
        },
        include_offset_penalty);
}

inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingObjectiveContribution(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingResidualBaseline & residual_baseline,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const LocalFittingObjectiveDomain & domain,
    bool include_offset_penalty = true)
{
    return EvaluateLocalFittingObjectiveContributionImpl(
        context,
        state,
        changed_key,
        sample_ref_list,
        domain,
        [&](const LocalFittingObjectiveSampleRef & sample_ref)
        {
            return residual_baseline.at(sample_ref.atom_index).at(
                sample_ref.sample_index);
        },
        include_offset_penalty);
}

inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingObjectiveContribution(
    const SecondStageLocalFittingContext & context,
    const LocalFittingStateView & state,
    const LocalFittingCandidateEvaluationOverlay & overlay,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const LocalFittingObjectiveDomain & domain,
    bool include_offset_penalty = true)
{
    return EvaluateLocalFittingObjectiveContributionImpl(
        context,
        state,
        changed_key,
        sample_ref_list,
        domain,
        [&](const LocalFittingObjectiveSampleRef & sample_ref)
        {
            return overlay.Evaluate(state, sample_ref);
        },
        include_offset_penalty);
}

inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingAuditObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingObjectiveDomain & domain,
    const SecondStageModelSnapshot & model_snapshot)
{
    LocalFittingObjectiveBreakdown total;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        const auto fit_contribution{
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                model_snapshot,
                key,
                cluster_domain.fit_sample_ref_list,
                domain,
                false)
        };
        if (!fit_contribution.has_value()) return std::nullopt;
        const auto tail_contribution{
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                model_snapshot,
                key,
                cluster_domain.tail_sample_ref_list,
                domain,
                false)
        };
        if (!tail_contribution.has_value()) return std::nullopt;
        const std::vector<LocalFittingObjectiveSampleRef> empty_sample_ref_list;
        const auto offset_contribution{
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                model_snapshot,
                key,
                empty_sample_ref_list,
                domain)
        };
        if (!offset_contribution.has_value()) return std::nullopt;
        total.fit_range_residual_objective +=
            fit_contribution->fit_range_residual_objective;
        total.fit_range_residual_objective +=
            tail_contribution->fit_range_residual_objective;
        total.tail_validation_loss += fit_contribution->tail_validation_loss;
        total.tail_validation_loss += tail_contribution->tail_validation_loss;
        total.tail_validation_penalty +=
            fit_contribution->tail_validation_penalty;
        total.tail_validation_penalty +=
            tail_contribution->tail_validation_penalty;
        total.offset_plausibility_penalty +=
            offset_contribution->offset_plausibility_penalty;
    }
    total.total_objective =
        total.fit_range_residual_objective +
        total.tail_validation_penalty +
        total.offset_plausibility_penalty;
    return std::isfinite(total.total_objective) ?
        std::optional<LocalFittingObjectiveBreakdown>{ total } :
        std::nullopt;
}

inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingAuditObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingObjectiveDomain & domain,
    const LocalFittingResidualBaseline & residual_baseline)
{
    LocalFittingObjectiveBreakdown total;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        std::vector<LocalFittingObjectiveSampleRef> sample_ref_list{
            cluster_domain.fit_sample_ref_list
        };
        sample_ref_list.insert(
            sample_ref_list.end(),
            cluster_domain.tail_sample_ref_list.begin(),
            cluster_domain.tail_sample_ref_list.end());
        const auto contribution{
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                residual_baseline,
                key,
                sample_ref_list,
                domain)
        };
        if (!contribution.has_value()) return std::nullopt;
        total.fit_range_residual_objective +=
            contribution->fit_range_residual_objective;
        total.tail_validation_loss += contribution->tail_validation_loss;
        total.offset_plausibility_penalty +=
            contribution->offset_plausibility_penalty;
    }
    return BuildLocalFittingObjectiveBreakdown(
        total.fit_range_residual_objective,
        total.tail_validation_loss,
        total.offset_plausibility_penalty,
        kLocalFittingTailValidationWeight);
}

inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingObjectiveDelta(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const LocalFittingStateView & candidate_state,
    const LocalFittingResidualBaseline & residual_baseline,
    const LocalFittingCandidateEvaluationOverlay & candidate_overlay,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & affected_sample_ref_list,
    const LocalFittingObjectiveDomain & domain,
    const std::optional<LocalFittingObjectiveBreakdown> & baseline,
    LocalFittingPerformanceCounters * performance_counters = nullptr)
{
    if (!baseline.has_value()) return std::nullopt;
    if (performance_counters != nullptr)
    {
        const auto unique_sample_count{
            domain.fit_sample_count + domain.tail_sample_count
        };
        performance_counters->objective_recomputed_sample_count.fetch_add(
            affected_sample_ref_list.size(),
            std::memory_order_relaxed);
        performance_counters->objective_reused_sample_count.fetch_add(
            unique_sample_count > affected_sample_ref_list.size() ?
                unique_sample_count - affected_sample_ref_list.size() : 0,
            std::memory_order_relaxed);
    }
    const auto candidate_changed{
        EvaluateLocalFittingObjectiveContribution(
            context,
            candidate_state,
            candidate_overlay,
            changed_key,
            affected_sample_ref_list,
            domain)
    };
    const auto previous_changed{
        EvaluateLocalFittingObjectiveContribution(
            context,
            previous_state,
            residual_baseline,
            changed_key,
            affected_sample_ref_list,
            domain)
    };
    if (!candidate_changed.has_value() || !previous_changed.has_value())
    {
        return std::nullopt;
    }
    return BuildLocalFittingObjectiveBreakdown(
        baseline->fit_range_residual_objective +
            candidate_changed->fit_range_residual_objective -
            previous_changed->fit_range_residual_objective,
        baseline->tail_validation_loss +
            candidate_changed->tail_validation_loss -
            previous_changed->tail_validation_loss,
        baseline->offset_plausibility_penalty +
            candidate_changed->offset_plausibility_penalty -
            previous_changed->offset_plausibility_penalty,
        kLocalFittingTailValidationWeight);
}

inline std::optional<LocalFittingObjectiveBreakdown>
EvaluateLocalFittingAuditObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingObjectiveDomain & domain)
{
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(state))
    };
    return EvaluateLocalFittingAuditObjective(
        context,
        state,
        domain,
        model_snapshot);
}

inline LocalFittingObjectiveByKey BuildLocalFittingObjectiveByKey(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const CouplingGraphPartition & partition,
    const LocalFittingObjectiveDomain & domain,
    const SecondStageModelSnapshot & model_snapshot)
{
    LocalFittingObjectiveByKey objective_by_key;
    for (const auto & [key, sample_ref_list] : partition.sample_id_list_by_key)
    {
        objective_by_key.emplace(
            key,
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                model_snapshot,
                key,
                sample_ref_list,
                domain));
    }
    return objective_by_key;
}

inline LocalFittingObjectiveByKey BuildLocalFittingObjectiveByKey(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const CouplingGraphPartition & partition,
    const LocalFittingObjectiveDomain & domain,
    const LocalFittingResidualBaseline & residual_baseline)
{
    LocalFittingObjectiveByKey objective_by_key;
    for (const auto & [key, sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        objective_by_key.emplace(
            key,
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                residual_baseline,
                key,
                sample_ref_list,
                domain));
    }
    return objective_by_key;
}

[[maybe_unused]] inline LocalFittingCombinedObjectiveCheck
EvaluateLocalFittingCombinedObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingState & previous_state,
    const LocalFittingObjectiveDomain & domain,
    const LocalFittingBestAuditState & audit_state,
    const std::optional<LocalFittingObjectiveBreakdown> & previous_objective = std::nullopt)
{
    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(context, candidate_state, domain)
    };
    const auto resolved_previous_objective{
        previous_objective.has_value() ?
            previous_objective :
            EvaluateLocalFittingAuditObjective(context, previous_state, domain)
    };
    return LocalFittingCombinedObjectiveCheck{
        IsLocalFittingAuditObjectiveAcceptableForProgress(
            candidate_objective.has_value() ?
                std::optional<double>{ candidate_objective->total_objective } :
                std::nullopt,
            resolved_previous_objective.has_value() ?
                std::optional<double>{ resolved_previous_objective->total_objective } :
                std::nullopt,
            audit_state.best.has_value() ?
                std::optional<double>{ audit_state.best->objective.total_objective } :
                std::nullopt,
            kLocalFittingObjectiveProgressTolerance),
        candidate_objective
    };
}

inline LocalFittingCombinedObjectiveCheck EvaluateLocalFittingCombinedObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const LocalFittingStateView & candidate_state,
    const LocalFittingResidualBaseline & residual_baseline,
    const LocalFittingCandidateEvaluationOverlay & candidate_overlay,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & affected_sample_ref_list,
    const LocalFittingObjectiveDomain & domain,
    const LocalFittingBestAuditState & audit_state,
    const std::optional<LocalFittingObjectiveBreakdown> & previous_objective,
    LocalFittingPerformanceCounters * performance_counters = nullptr)
{
    const auto candidate_objective{
        EvaluateLocalFittingObjectiveDelta(
            context,
            previous_state,
            candidate_state,
            residual_baseline,
            candidate_overlay,
            changed_key,
            affected_sample_ref_list,
            domain,
            previous_objective,
            performance_counters)
    };
    return LocalFittingCombinedObjectiveCheck{
        IsLocalFittingAuditObjectiveAcceptableForProgress(
            candidate_objective.has_value() ?
                std::optional<double>{ candidate_objective->total_objective } :
                std::nullopt,
            previous_objective.has_value() ?
                std::optional<double>{ previous_objective->total_objective } :
                std::nullopt,
            audit_state.best.has_value() ?
                std::optional<double>{ audit_state.best->objective.total_objective } :
                std::nullopt,
            kLocalFittingObjectiveProgressTolerance),
        candidate_objective
    };
}


inline bool TryUpdateLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingPolishProvenance & candidate_polish_provenance,
    const std::optional<std::size_t> & accepted_iteration,
    const LocalFittingObjectiveDomain & domain,
    LocalFittingBestAuditState & audit_state,
    const std::optional<LocalFittingObjectiveBreakdown> & precomputed_objective = std::nullopt)
{
    const auto candidate_objective{ precomputed_objective.has_value() ?
        precomputed_objective :
        EvaluateLocalFittingAuditObjective(context, candidate_state, domain) };
    if (!candidate_objective.has_value()) return false;
    if (audit_state.best.has_value() &&
        !IsBetterLocalFittingAuditObjective(
            candidate_objective->total_objective,
            audit_state.best->objective.total_objective,
            kLocalFittingObjectiveStrictTolerance))
    {
        return false;
    }
    audit_state.best = LocalFittingAuditedState{
        *candidate_objective,
        candidate_state,
        candidate_polish_provenance,
        accepted_iteration
    };
    return true;
}

inline LocalFittingBestAuditState BuildInitialLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state,
    const LocalFittingPolishProvenance & initial_polish_provenance,
    const std::optional<std::size_t> & accepted_iteration,
    const LocalFittingObjectiveDomain & domain)
{
    LocalFittingBestAuditState audit_state;
    static_cast<void>(TryUpdateLocalFittingBestAuditState(
        context,
        initial_state,
        initial_polish_provenance,
        accepted_iteration,
        domain,
        audit_state));
    return audit_state;
}

inline void ResetLocalFittingBestAuditAfterObjectiveDomainChange(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & validated_state,
    const LocalFittingPolishProvenance & validated_polish_provenance,
    std::size_t accepted_iteration,
    const LocalFittingObjectiveDomain & domain,
    LocalFittingBestAuditState & audit_state)
{
    audit_state = BuildInitialLocalFittingBestAuditState(
        context,
        validated_state,
        validated_polish_provenance,
        accepted_iteration,
        domain);
}


inline LocalFittingClusterObjectiveState
BuildInitialLocalFittingClusterObjectiveState(
    const LocalFittingObjectiveByKey & previous_objective_by_key,
    const LocalFittingClusterKey & key)
{
    LocalFittingClusterObjectiveState state;
    const auto objective_iter{ previous_objective_by_key.find(key) };
    if (objective_iter != previous_objective_by_key.end())
    {
        state.best_objective = objective_iter->second;
    }
    return state;
}

inline void ReconcileLocalFittingClusterObjectiveState(
    const CouplingGraphPartition & partition,
    const LocalFittingObjectiveByKey & previous_objective_by_key,
    LocalFittingClusterObjectiveStateMap & state_by_key)
{
    LocalFittingClusterObjectiveStateMap next_state_by_key;
    for (const auto & [key, objective_sample_ref_list] : partition.sample_id_list_by_key)
    {
        static_cast<void>(objective_sample_ref_list);
        auto state_iter{ state_by_key.find(key) };
        if (state_iter != state_by_key.end())
        {
            next_state_by_key.emplace(key, std::move(state_iter->second));
            continue;
        }
        next_state_by_key.emplace(
            key,
            BuildInitialLocalFittingClusterObjectiveState(
                previous_objective_by_key,
                key));
    }
    state_by_key = std::move(next_state_by_key);
}

inline bool TryCommitLocalFittingClusterCandidate(
    const SecondStageLocalFittingContext & context,
    const LocalFittingStateView & candidate_state,
    const LocalFittingCandidateEvaluationOverlay & candidate_overlay,
    const LocalFittingStateView & previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & objective_sample_ref_list,
    const std::optional<LocalFittingObjectiveBreakdown> & previous_objective,
    bool requires_strict_improvement,
    const LocalFittingObjectiveDomain & domain,
    LocalFittingClusterObjectiveStateMap & cluster_objective_state,
    LocalFittingObjectiveAttemptDiagnostic & diagnostic,
    LocalFittingPerformanceCounters * performance_counters = nullptr)
{
    if (performance_counters != nullptr)
    {
        const auto unique_sample_count{
            domain.fit_sample_count + domain.tail_sample_count
        };
        performance_counters->objective_recomputed_sample_count.fetch_add(
            objective_sample_ref_list.size(),
            std::memory_order_relaxed);
        performance_counters->objective_reused_sample_count.fetch_add(
            unique_sample_count > objective_sample_ref_list.size() ?
                unique_sample_count - objective_sample_ref_list.size() : 0,
            std::memory_order_relaxed);
    }
    const auto maximum_transformed_change{
        algorithm::GetMaximumParameterChange(
            SummarizeLocalFittingTransformedChanges(
                candidate_state,
                previous_state,
                key).percentile_stats)
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
        EvaluateLocalFittingObjectiveContribution(
            context,
            candidate_state,
            candidate_overlay,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.previous_objective = previous_objective;
    diagnostic.best_objective = state.best_objective;

    const auto is_objective_deteriorated = [](
        const std::optional<LocalFittingObjectiveBreakdown> & candidate,
        const std::optional<double> & reference)
    {
        if (!reference.has_value()) return false;
        if (!candidate.has_value()) return true;
        return candidate->total_objective > *reference +
            CalculateLocalFittingObjectiveTolerance(
                *reference,
                kLocalFittingObjectiveProgressTolerance);
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
            std::optional<double>{ state.best_objective->total_objective } :
            std::nullopt);
    if (diagnostic.rejected_by_previous || diagnostic.rejected_by_best)
    {
        return false;
    }
    if (requires_strict_improvement &&
        (!diagnostic.candidate_objective.has_value() ||
            !diagnostic.previous_objective.has_value() ||
            !IsBetterLocalFittingAuditObjective(
                diagnostic.candidate_objective->total_objective,
                diagnostic.previous_objective->total_objective,
                kLocalFittingObjectiveStrictTolerance)))
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
        if (IsBetterLocalFittingAuditObjective(
                candidate_objective_value,
                best_objective_value,
                kLocalFittingObjectiveStrictTolerance))
        {
            is_better_than_best = true;
        }
        else if (IsBetterLocalFittingAuditObjective(
                     best_objective_value,
                     candidate_objective_value,
                     kLocalFittingObjectiveStrictTolerance))
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
