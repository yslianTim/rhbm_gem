#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "core/detail/CouplingGraph.hpp"
#include "core/detail/CandidateEvaluationOverlay.hpp"
#include "core/detail/PerformanceCounters.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/SecondStageContext.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/TransformedChange.hpp"

namespace rhbm_gem::core::detail {

constexpr double kTailValidationWeight{ 0.25 };
constexpr double kObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kObjectiveResidualScaleMin{ 1.0e-12 };
constexpr double kFitRangeWeight{ 1.0 };
constexpr double kOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kOffsetPeakRatioMax{ 1.0 };
constexpr double kObjectiveRobustLossCutoffMultiplier{ 1.345 };

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

constexpr ObjectiveTolerance kObjectiveStrictTolerance{ 1.0e-10, 1.0e-8 };
constexpr ObjectiveTolerance kObjectiveProgressTolerance{ 1.0e-8, 1.0e-3 };

struct AuditedState
{
    ObjectiveBreakdown objective{};
    FitState state{};
    bool uses_polish{ false };
    std::size_t source_iteration{ 0 };
};

using BestAuditState = std::optional<AuditedState>;

struct FinalStateSelection
{
    const FitState & state;
    bool uses_polish{ false };
    const AuditedState * audit_state{ nullptr };
};

inline FinalStateSelection SelectFinalState(
    const FitState & latest_validated_state,
    bool latest_validated_uses_polish,
    const BestAuditState & audited_state)
{
    if (audited_state.has_value())
    {
        return FinalStateSelection{
            audited_state->state,
            audited_state->uses_polish,
            &*audited_state
        };
    }
    return FinalStateSelection{
        latest_validated_state,
        latest_validated_uses_polish,
        nullptr
    };
}

enum class PreObjectiveFailureReason
{
    None,
    InvalidModel,
    PreviousSharedOffsetProjectionOutsideTrustRegion,
    NoCandidateWithinTrustRegion
};

struct ObjectiveScale
{
    double fit{ 0.0 };
    double tail{ 0.0 };
};

struct ObjectiveAttemptDiagnostic
{
    double effective_damping{ 1.0 };
    bool is_invalid_model{ false };
    PreObjectiveFailureReason pre_objective_failure_reason{ PreObjectiveFailureReason::None };
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<ObjectiveScale> scale{};
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

inline double CalculateClusterAtomWeight(std::size_t cluster_atom_count, std::size_t active_atom_count)
{
    if (cluster_atom_count == 0 || active_atom_count == 0 || cluster_atom_count > active_atom_count)
    {
        throw std::invalid_argument("Local fitting cluster atom counts are invalid.");
    }
    return static_cast<double>(cluster_atom_count) / static_cast<double>(active_atom_count);
}

inline void ValidateObjectiveTolerance(ObjectiveTolerance tolerance)
{
    if (!std::isfinite(tolerance.absolute_tolerance) || tolerance.absolute_tolerance < 0.0 ||
        !std::isfinite(tolerance.relative_tolerance) || tolerance.relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerances must be finite and non-negative.");
    }
}

inline double CalculateObjectiveTolerance(double reference, ObjectiveTolerance tolerance)
{
    return tolerance.absolute_tolerance + tolerance.relative_tolerance * std::abs(reference);
}

inline bool IsObjectiveDeteriorated(double candidate, double reference, ObjectiveTolerance tolerance)
{
    if (!std::isfinite(reference)) return true;
    return candidate > reference + CalculateObjectiveTolerance(reference, tolerance);
}

inline std::optional<ObjectiveBreakdown> BuildObjectiveBreakdown(
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

    const ObjectiveBreakdown breakdown{
        fit_range_residual_objective,
        tail_validation_loss,
        offset_plausibility_penalty
    };
    if (!std::isfinite(breakdown.GetTotalObjective())) return std::nullopt;
    return breakdown;
}

inline bool IsBetterAuditObjective(double candidate, double best, ObjectiveTolerance tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate)) return false;
    if (!std::isfinite(best)) return true;
    return candidate < best - CalculateObjectiveTolerance(best, tolerance);
}

inline bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate) || !std::isfinite(previous)) return false;
    return !IsObjectiveDeteriorated(candidate, previous, tolerance) &&
        (best == nullptr ||
            !IsObjectiveDeteriorated(candidate, best->GetTotalObjective(), tolerance));
}

struct ObjectiveClusterDomain
{
    std::vector<SampleRef> fit_sample_ref_list{};
    std::vector<SampleRef> tail_sample_ref_list{};
    std::optional<ObjectiveScale> scale{};
};

struct ResidualObjectiveContribution
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
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

inline void AppendObjectiveScaleSummary(std::ostringstream & message, const std::vector<double> & scale_list)
{
    if (scale_list.empty())
    {
        message << "unavailable";
        return;
    }
    message
        << array_helper::ComputePercentile(scale_list, 0.5) << "/"
        << array_helper::ComputePercentile(scale_list, 0.99) << "/"
        << *std::max_element(scale_list.begin(), scale_list.end());
}

inline void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset = false)
{
    if (quiet_mode) return;
    std::vector<double> fit_scale_list;
    std::vector<double> tail_scale_list;
    for (const auto & entry : domain.cluster_by_key)
    {
        const auto & cluster_domain{ entry.second };
        if (!cluster_domain.scale.has_value()) continue;
        fit_scale_list.emplace_back(cluster_domain.scale->fit);
        if (!cluster_domain.tail_sample_ref_list.empty())
        {
            tail_scale_list.emplace_back(cluster_domain.scale->tail);
        }
    }
    std::ostringstream message;
    message
        << (is_terminal_reset ?
            "Reset second-stage objective domain" : "Initialize second-stage objective domain")
        << ": fit/tail/offset weights = "
        << kFitRangeWeight << "/" << kTailValidationWeight << "/" << kOffsetPlausibilityPenaltyWeight
        << ", clusters = " << domain.cluster_by_key.size()
        << ", active atoms = " << domain.active_atom_count
        << ", unique fit/tail samples = " << domain.fit_sample_count << "/" << domain.tail_sample_count
        << ", fixed fit scale median/p99/max = ";
    AppendObjectiveScaleSummary(message, fit_scale_list);
    message << ", fixed tail scale median/p99/max = ";
    AppendObjectiveScaleSummary(message, tail_scale_list);
    message << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message.str());
}

struct ClusterObjectiveState
{
    std::optional<ObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
};

using ClusterObjectiveStateMap = std::map<ClusterKey, ClusterObjectiveState>;
using ObjectiveByKey = std::map<ClusterKey, std::optional<ObjectiveBreakdown>>;

struct CombinedCandidateObjectiveCheck
{
    bool accepted{ true };
    std::optional<ObjectiveBreakdown> previous_objective{};
    std::optional<ObjectiveBreakdown> candidate_objective{};
};

inline std::optional<double> BuildFixedObjectiveScale(
    const std::vector<double> & residual_list,
    const std::vector<double> & adjusted_response_list)
{
    if (residual_list.empty() || residual_list.size() != adjusted_response_list.size())
    {
        return std::nullopt;
    }
    const auto scale{
        std::max({
            array_helper::ComputeMedianAbsoluteDeviationScale(residual_list),
            kObjectiveResidualScaleFloorRatio *
                array_helper::ComputeMedianAbsoluteDeviationScale(adjusted_response_list),
            kObjectiveResidualScaleMin
        })
    };
    return numeric_validation::IsFinitePositive(scale) ? std::optional<double>{ scale } : std::nullopt;
}

inline ObjectiveDomain BuildObjectiveDomain(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<ClusterKey> & cluster_key_list,
    double distance_min,
    double distance_max)
{
    ObjectiveDomain domain;
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
            const auto & raw_sampling_entries{ context.at(atom_index).raw_sampling_entries };
            for (std::size_t sample_index = 0; sample_index < raw_sampling_entries.size(); sample_index++)
            {
                const SampleRef sample_ref{ atom_index, sample_index };
                const auto residual_sample{
                    EvaluateResidualSample(context, model_snapshot.selected, sample_ref, model_snapshot)
                };
                const auto distance{
                    static_cast<double>(raw_sampling_entries.at(sample_index).point.distance)
                };
                const auto is_fit_range{ distance >= distance_min && distance <= distance_max };
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
            const auto tail_scale{
                BuildFixedObjectiveScale(tail_residual_list, tail_response_list)
            };
            if (!tail_scale.has_value() ||
                tail_residual_list.size() != cluster_domain.tail_sample_ref_list.size())
            {
                continue;
            }
            scale.tail = *tail_scale;
        }
        cluster_domain.scale = scale;
    }
    return domain;
}


template <typename ResidualEvaluator>
inline std::optional<ResidualObjectiveContribution> EvaluateResidualObjectiveContributionImpl(
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    const ResidualEvaluator & residual_evaluator)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    ResidualObjectiveContribution contribution;
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
        double scale;
        if (is_fit_range)
        {
            scale = owner_iter->second.scale->fit;
        }
        else
        {
            scale = owner_iter->second.scale->tail;
        }
        const auto loss{
            algorithm::CalculateCauchyLoss(
                residual_sample->residual / scale,
                kObjectiveRobustLossCutoffMultiplier)
        };
        const auto coefficient{
            CalculateClusterAtomWeight(owner_key.size(), domain.active_atom_count) /
            static_cast<double>(sample_count)
        };
        if (is_fit_range)
        {
            contribution.fit_range_residual_objective += kFitRangeWeight * coefficient * loss;
        }
        else
        {
            contribution.tail_validation_loss += coefficient * loss;
        }
    }
    if (!std::isfinite(contribution.fit_range_residual_objective) ||
        !std::isfinite(contribution.tail_validation_loss))
    {
        return std::nullopt;
    }
    return contribution;
}

template <typename State>
inline std::optional<double> EvaluateOffsetPlausibilityPenalty(
    const State & state,
    const ClusterKey & changed_key,
    const ObjectiveDomain & domain)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    double penalty{ 0.0 };
    for (const auto atom_index : changed_key)
    {
        const auto owner_iter{
            domain.cluster_by_key.find(domain.owner_key_by_atom_index.at(atom_index))
        };
        if (owner_iter == domain.cluster_by_key.end() || !owner_iter->second.scale.has_value()) return std::nullopt;
        const auto & model{ GetFitModel(state, atom_index) };
        if (!IsValidSecondStageGaussianModel(model)) return std::nullopt;
        const auto peak_signal{ model.SignalAtDistance(0.0) };
        const auto offset_peak{ model.GetOffset() * model.OffsetBasisAtDistance(0.0) };
        if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak)) return std::nullopt;
        const auto offset_ratio{
            std::abs(offset_peak) /
            std::max({ std::abs(peak_signal), owner_iter->second.scale->fit, kObjectiveResidualScaleMin })
        };
        const auto offset_excess{ std::max(0.0, offset_ratio - kOffsetPeakRatioMax) };
        penalty +=
            kOffsetPlausibilityPenaltyWeight * offset_excess * offset_excess /
            static_cast<double>(domain.active_atom_count);
    }
    return std::isfinite(penalty) ? std::optional<double>{ penalty } : std::nullopt;
}

template <typename State, typename ResidualEvaluator>
inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContributionImpl(
    const State & state,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    const ResidualEvaluator & residual_evaluator)
{
    const auto residual_contribution{
        EvaluateResidualObjectiveContributionImpl(sample_ref_list, domain, residual_evaluator)
    };
    if (!residual_contribution.has_value()) return std::nullopt;
    const auto offset_penalty{ EvaluateOffsetPlausibilityPenalty(state, changed_key, domain) };
    if (!offset_penalty.has_value()) return std::nullopt;
    return BuildObjectiveBreakdown(
        residual_contribution->fit_range_residual_objective,
        residual_contribution->tail_validation_loss,
        *offset_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    return EvaluateObjectiveContributionImpl(
        model_snapshot.selected,
        changed_key,
        sample_ref_list,
        domain,
        [&](const SampleRef & sample_ref)
        {
            return EvaluateResidualSample(
                context,
                model_snapshot.selected,
                sample_ref,
                model_snapshot);
        });
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const ResidualBaseline & baseline,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    return EvaluateObjectiveContributionImpl(
        baseline.model_snapshot.selected,
        changed_key,
        sample_ref_list,
        domain,
        [&](const SampleRef & sample_ref)
        {
            return baseline.sample_list.at(sample_ref.atom_index).at(sample_ref.sample_index);
        });
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const CandidateEvaluationOverlay & overlay,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    return EvaluateObjectiveContributionImpl(
        overlay.GetCandidateState(),
        changed_key,
        sample_ref_list,
        domain,
        [&](const SampleRef & sample_ref)
        {
            return overlay.Evaluate(sample_ref);
        });
}

template <typename ResidualEvaluator>
inline std::optional<ObjectiveBreakdown> EvaluateAuditObjectiveImpl(
    const ObjectiveDomain & domain,
    const FittedGaussianSnapshot & selected_state,
    const ResidualEvaluator & residual_evaluator)
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        const auto fit_contribution{
            EvaluateResidualObjectiveContributionImpl(
                cluster_domain.fit_sample_ref_list,
                domain,
                residual_evaluator)
        };
        if (!fit_contribution.has_value()) return std::nullopt;
        const auto tail_contribution{
            EvaluateResidualObjectiveContributionImpl(
                cluster_domain.tail_sample_ref_list,
                domain,
                residual_evaluator)
        };
        if (!tail_contribution.has_value()) return std::nullopt;
        const auto offset_contribution{
            EvaluateOffsetPlausibilityPenalty(selected_state, key, domain)
        };
        if (!offset_contribution.has_value()) return std::nullopt;
        fit_range_residual_objective += fit_contribution->fit_range_residual_objective;
        fit_range_residual_objective += tail_contribution->fit_range_residual_objective;
        tail_validation_loss += fit_contribution->tail_validation_loss;
        tail_validation_loss += tail_contribution->tail_validation_loss;
        offset_plausibility_penalty += *offset_contribution;
    }
    return BuildObjectiveBreakdown(
        fit_range_residual_objective,
        tail_validation_loss,
        offset_plausibility_penalty);
}

inline std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const SecondStageContext & context,
    const ObjectiveDomain & domain,
    const SecondStageModelSnapshot & model_snapshot)
{
    return EvaluateAuditObjectiveImpl(
        domain,
        model_snapshot.selected,
        [&](const SampleRef & sample_ref)
        {
            return EvaluateResidualSample(
                context,
                model_snapshot.selected,
                sample_ref,
                model_snapshot);
        });
}

inline std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const ResidualBaseline & baseline)
{
    return EvaluateAuditObjectiveImpl(
        domain,
        baseline.model_snapshot.selected,
        [&](const SampleRef & sample_ref)
        {
            return baseline.sample_list.at(sample_ref.atom_index).at(sample_ref.sample_index);
        });
}

inline std::optional<ObjectiveBreakdown> EvaluateObjectiveDelta(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown & baseline,
    PerformanceCounters & performance_counters)
{
    const auto & changed_key{
        candidate_overlay.GetCandidateState().GetOverrideAtomIndexList()
    };
    const auto unique_sample_count{ domain.fit_sample_count + domain.tail_sample_count };
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
    if (!candidate_changed.has_value() || !previous_changed.has_value()) return std::nullopt;
    return BuildObjectiveBreakdown(
        baseline.fit_range_residual_objective +
            candidate_changed->fit_range_residual_objective -
            previous_changed->fit_range_residual_objective,
        baseline.tail_validation_loss +
            candidate_changed->tail_validation_loss -
            previous_changed->tail_validation_loss,
        baseline.offset_plausibility_penalty +
            candidate_changed->offset_plausibility_penalty -
            previous_changed->offset_plausibility_penalty);
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
            EvaluateObjectiveContribution(context, model_snapshot, key, sample_ref_list, domain));
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
            EvaluateObjectiveContribution(baseline, key, sample_ref_list, domain));
    }
    return objective_by_key;
}

inline std::optional<ObjectiveBreakdown> EvaluateCombinedObjective(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown * best_objective,
    const ObjectiveBreakdown * previous_objective,
    PerformanceCounters & performance_counters)
{
    if (previous_objective == nullptr) return std::nullopt;
    const auto candidate_objective{
        EvaluateObjectiveDelta(
            candidate_overlay,
            affected_sample_ref_list,
            domain,
            *previous_objective,
            performance_counters)
    };
    if (!candidate_objective.has_value() ||
        !IsAuditObjectiveAcceptableForProgress(
            candidate_objective->GetTotalObjective(),
            previous_objective->GetTotalObjective(),
            best_objective,
            kObjectiveProgressTolerance))
    {
        return std::nullopt;
    }
    return candidate_objective;
}

inline CombinedCandidateObjectiveCheck EvaluateCombinedCandidateObjective(
    const SecondStageContext & context,
    const ResidualBaseline & baseline,
    const CouplingGraphPartition & partition,
    const FitState & previous_state,
    const FitState & candidate_state,
    const std::vector<ClusterKey> & accepted_key_list,
    const ObjectiveDomain & objective_domain,
    const ObjectiveBreakdown * best_objective,
    PerformanceCounters & performance_counters)
{
    CombinedCandidateObjectiveCheck result;
    if (partition.boundary_sample_count == 0 || accepted_key_list.empty()) return result;

    result.previous_objective = EvaluateAuditObjective(objective_domain, baseline);

    ClusterKey changed_atom_index_list;
    for (const auto & key : accepted_key_list)
    {
        changed_atom_index_list.insert(changed_atom_index_list.end(), key.begin(), key.end());
    }
    const auto combined_patch{
        FitStatePatch::FromState(candidate_state, std::move(changed_atom_index_list))
    };
    const FitStateView combined_state_view{ previous_state, combined_patch };
    const CandidateEvaluationOverlay combined_overlay{ context, baseline, combined_state_view };
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(partition, accepted_key_list)
    };
    const auto combined_check{
        EvaluateCombinedObjective(
            combined_overlay,
            affected_sample_ref_list,
            objective_domain,
            best_objective,
            result.previous_objective.has_value() ? &*result.previous_objective : nullptr,
            performance_counters)
    };
    result.accepted = combined_check.has_value();
    result.candidate_objective = combined_check;
    return result;
}

inline bool TryUpdateBestAuditState(
    const FitState & candidate_state,
    bool candidate_uses_polish,
    std::size_t source_iteration,
    const ObjectiveBreakdown & candidate_objective,
    BestAuditState & audit_state)
{
    if (audit_state.has_value() &&
        !IsBetterAuditObjective(
            candidate_objective.GetTotalObjective(),
            audit_state->objective.GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        return false;
    }
    audit_state = AuditedState{
        candidate_objective,
        candidate_state,
        candidate_uses_polish,
        source_iteration
    };
    return true;
}

inline void ReconcileClusterObjectiveState(
    const ObjectiveByKey & previous_objective_by_key,
    ClusterObjectiveStateMap & state_by_key)
{
    ClusterObjectiveStateMap next_state_by_key;
    for (const auto & [key, previous_objective] : previous_objective_by_key)
    {
        auto state_iter{ state_by_key.find(key) };
        if (state_iter != state_by_key.end())
        {
            next_state_by_key.emplace(key, std::move(state_iter->second));
            continue;
        }
        ClusterObjectiveState state;
        state.best_objective = previous_objective;
        next_state_by_key.emplace(key, std::move(state));
    }
    state_by_key = std::move(next_state_by_key);
}

inline bool TryCommitClusterCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveBreakdown * previous_objective,
    bool requires_strict_improvement,
    const ObjectiveDomain & domain,
    ClusterObjectiveState & objective_state,
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
            candidate_overlay.GetBaseline().model_snapshot.selected,
            key)
    };
    const auto maximum_transformed_change{
        GetMaximumTransformedChange(transformed_change_summary.percentile_stats.percentile_list)
    };
    const auto domain_iter{ domain.cluster_by_key.find(key) };
    diagnostic.scale.reset();
    if (domain_iter != domain.cluster_by_key.end())
    {
        diagnostic.fit_sample_count = domain_iter->second.fit_sample_ref_list.size();
        diagnostic.tail_sample_count = domain_iter->second.tail_sample_ref_list.size();
        diagnostic.scale = domain_iter->second.scale;
    }
    diagnostic.candidate_objective =
        EvaluateObjectiveContribution(
            candidate_overlay,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.previous_objective.reset();
    if (previous_objective != nullptr)
    {
        diagnostic.previous_objective = *previous_objective;
    }
    diagnostic.best_objective = objective_state.best_objective;

    if (!diagnostic.candidate_objective.has_value() || previous_objective == nullptr)
    {
        return false;
    }
    const auto candidate_objective_value{ diagnostic.candidate_objective->GetTotalObjective() };
    const auto previous_objective_value{ previous_objective->GetTotalObjective() };
    diagnostic.rejected_by_previous = IsObjectiveDeteriorated(
        candidate_objective_value,
        previous_objective_value,
        kObjectiveProgressTolerance);
    diagnostic.rejected_by_best = objective_state.best_objective.has_value() &&
        IsObjectiveDeteriorated(
            candidate_objective_value,
            objective_state.best_objective->GetTotalObjective(),
            kObjectiveProgressTolerance);
    if (diagnostic.rejected_by_previous || diagnostic.rejected_by_best) return false;
    if (requires_strict_improvement &&
        !IsBetterAuditObjective(
            candidate_objective_value,
            previous_objective_value,
            kObjectiveStrictTolerance))
    {
        return false;
    }

    auto is_better_than_best{ !objective_state.best_objective.has_value() };
    if (objective_state.best_objective.has_value())
    {
        const auto best_objective_value{ objective_state.best_objective->GetTotalObjective() };
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
            is_better_than_best = maximum_transformed_change < objective_state.best_maximum_transformed_change;
        }
    }
    if (is_better_than_best)
    {
        objective_state.best_objective = diagnostic.candidate_objective;
        objective_state.best_maximum_transformed_change = maximum_transformed_change;
    }
    return true;
}

} // namespace rhbm_gem::core::detail
