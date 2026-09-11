#include "core/detail/ObjectiveEvaluation.hpp"
#include "core/detail/CandidateEvaluation.hpp"

#include "core/detail/Diagnosis.hpp"
#include "core/detail/FittingRanges.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core::detail {

std::size_t CountObjectiveSamples(
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    return static_cast<std::size_t>(std::ranges::count_if(sample_ref_list, [&](const SampleRef & sample_ref)
    {
        return !domain.owner_key_by_atom_index.at(sample_ref.atom_index).empty() &&
            (domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0 ||
                domain.tail_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0);
    }));
}

namespace {

constexpr double kObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kObjectiveResidualScaleMin{ 1.0e-12 };
constexpr double kOffsetPeakRatioMax{ 1.0 };

template<typename ResidualEvaluator>
std::optional<ObjectiveBreakdown> EvaluateResidualObjectiveContribution(
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    const ResidualEvaluator & residual_evaluator)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    ObjectiveBreakdown contribution;
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & owner_key{
            domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
        const auto in_fit{
            domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        const auto in_tail{
            domain.tail_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        if (!in_fit && !in_tail) continue;
        const auto owner_iter{ domain.cluster_by_key.find(owner_key) };
        if (owner_iter == domain.cluster_by_key.end() ||
            !owner_iter->second.scale.has_value())
        {
            return std::nullopt;
        }
        const auto residual_sample{ residual_evaluator(sample_ref) };
        if (!residual_sample.has_value()) return std::nullopt;
        for (const bool is_fit_range : { true, false })
        {
            if (!(is_fit_range ? in_fit : in_tail)) continue;
            const auto sample_count{ is_fit_range ?
                owner_iter->second.fit_sample_ref_list.size() :
                owner_iter->second.tail_sample_ref_list.size()
            };
            if (sample_count == 0) return std::nullopt;
            const auto scale{ is_fit_range ? owner_iter->second.scale->fit : owner_iter->second.scale->tail };
            const auto loss{
                algorithm::CalculateCauchyLoss(
                    residual_sample->residual / scale,
                    kObjectiveRobustLossCutoffMultiplier)
            };
            const auto coefficient{
                CalculateClusterAtomWeight(
                    owner_iter->second.selected_atom_count,
                    domain.active_atom_count) / static_cast<double>(sample_count)
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
    }
    if (!std::isfinite(contribution.fit_range_residual_objective) ||
        !std::isfinite(contribution.tail_validation_loss))
    {
        return std::nullopt;
    }
    return contribution;
}

template<typename State>
std::optional<double> EvaluateOffsetPlausibilityPenalty(
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
                owner_iter->second.scale->fit,
                kObjectiveResidualScaleMin
            })
        };
        const auto offset_excess{ std::max(0.0, offset_ratio - kOffsetPeakRatioMax) };
        penalty +=
            kOffsetPlausibilityPenaltyWeight * offset_excess * offset_excess /
            static_cast<double>(domain.active_atom_count);
    }
    return std::isfinite(penalty) ? std::optional<double>{ penalty } : std::nullopt;
}

template<typename State, typename Evaluator>
std::optional<ObjectiveBreakdown> EvaluateObjectiveContributionImpl(
    const State & state,
    const Evaluator & evaluator,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    const auto residual_contribution{
        EvaluateResidualObjectiveContribution(sample_ref_list, domain, evaluator)
    };
    if (!residual_contribution.has_value()) return std::nullopt;
    const auto offset_penalty{
        EvaluateOffsetPlausibilityPenalty(state, changed_key, domain)
    };
    if (!offset_penalty.has_value()) return std::nullopt;
    return BuildObjectiveBreakdown(
        residual_contribution->fit_range_residual_objective,
        residual_contribution->tail_validation_loss,
        *offset_penalty);
}

template<typename State, typename Evaluator>
std::optional<ObjectiveBreakdown> EvaluateAuditObjectiveImpl(
    const ObjectiveDomain & domain,
    const State & state,
    const Evaluator & evaluator)
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        const auto residual_contribution{
            EvaluateResidualObjectiveContribution(
                cluster_domain.sample_ref_list,
                domain,
                evaluator)
        };
        if (!residual_contribution.has_value()) return std::nullopt;
        const auto offset_contribution{
            EvaluateOffsetPlausibilityPenalty(
                state,
                key,
                domain)
        };
        if (!offset_contribution.has_value()) return std::nullopt;
        fit_range_residual_objective += residual_contribution->fit_range_residual_objective;
        tail_validation_loss += residual_contribution->tail_validation_loss;
        offset_plausibility_penalty += *offset_contribution;
    }
    return BuildObjectiveBreakdown(
        fit_range_residual_objective,
        tail_validation_loss,
        offset_plausibility_penalty);
}

template<typename State, typename Evaluator>
ObjectiveByKey BuildObjectiveByKeyImpl(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const State & state,
    const Evaluator & evaluator)
{
    ObjectiveByKey objective_by_key;
    for (const auto & [key, sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        objective_by_key.emplace(
            key,
            EvaluateObjectiveContributionImpl(
                state,
                evaluator,
                key,
                sample_ref_list,
                domain));
    }
    return objective_by_key;
}

} // namespace

std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const ResidualBaseline & evaluator,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    return EvaluateObjectiveContributionImpl(evaluator.GetState(), evaluator, changed_key, sample_ref_list, domain);
}

std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const CandidateEvaluationOverlay & evaluator,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    return EvaluateObjectiveContributionImpl(evaluator.GetState(), evaluator, changed_key, sample_ref_list, domain);
}

CandidateEvaluationOverlay::CandidateEvaluationOverlay(
    const SecondStageContext & context,
    const ResidualBaseline & baseline,
    const FitState & base_state,
    const FitStatePatch & candidate_patch)
    : m_context{ context },
      m_baseline{ baseline },
      m_candidate_state{ base_state, candidate_patch }
{
}

std::optional<ResidualSample> CandidateEvaluationOverlay::operator()(const SampleRef & sample_ref) const
{
    const auto & baseline{
        m_baseline.sample_list.at(sample_ref.atom_index).at(sample_ref.sample_index)
    };
    if (!baseline.has_value()) return std::nullopt;
    const auto & atom_context{ m_context.atom_list.at(sample_ref.atom_index) };
    const auto & sample{
        atom_context.raw_sampling_entries.at(sample_ref.sample_index)
    };
    auto adjusted_response{ baseline->adjusted_response };
    for (const auto & neighbor_atom_sample :
        atom_context.Neighbors(sample_ref.sample_index))
    {
        if (m_candidate_state.FindOverride(neighbor_atom_sample.atom_index) == nullptr)
        {
            continue;
        }
        adjusted_response +=
            GetFitModel(
                m_baseline.model_snapshot.node,
                neighbor_atom_sample.atom_index).ResponseAtDistance(
                    neighbor_atom_sample.distance) -
            m_candidate_state.GetModel(
                neighbor_atom_sample.atom_index).ResponseAtDistance(
                    neighbor_atom_sample.distance);
    }
    const auto expected_response{
        m_candidate_state.FindOverride(sample_ref.atom_index) == nullptr ?
            baseline->adjusted_response - baseline->residual :
            m_candidate_state.GetModel(sample_ref.atom_index).ResponseAtDistance(
                sample.point.distance)
    };
    const auto residual{ adjusted_response - expected_response };
    if (!std::isfinite(adjusted_response) || !std::isfinite(residual))
    {
        return std::nullopt;
    }
    return ResidualSample{ adjusted_response, residual };
}

double CalculateClusterAtomWeight(std::size_t cluster_atom_count, std::size_t active_atom_count)
{
    if (cluster_atom_count == 0 || active_atom_count == 0 || cluster_atom_count > active_atom_count)
    {
        throw std::invalid_argument("Local fitting cluster atom counts are invalid.");
    }
    return static_cast<double>(cluster_atom_count) / static_cast<double>(active_atom_count);
}

static void ValidateObjectiveTolerance(ObjectiveTolerance tolerance)
{
    if (!std::isfinite(tolerance.absolute_tolerance) || tolerance.absolute_tolerance < 0.0 ||
        !std::isfinite(tolerance.relative_tolerance) || tolerance.relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerances must be finite and non-negative.");
    }
}

double CalculateObjectiveTolerance(
    double reference,
    ObjectiveTolerance tolerance)
{
    return tolerance.absolute_tolerance + tolerance.relative_tolerance * std::abs(reference);
}

bool IsObjectiveDeteriorated(
    double candidate,
    double reference,
    ObjectiveTolerance tolerance)
{
    if (!std::isfinite(reference)) return true;
    return candidate > reference + CalculateObjectiveTolerance(reference, tolerance);
}

std::optional<ObjectiveBreakdown> BuildObjectiveBreakdown(
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

bool IsBetterAuditObjective(double candidate, double best, ObjectiveTolerance tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate)) return false;
    if (!std::isfinite(best)) return true;
    return candidate < best - CalculateObjectiveTolerance(best, tolerance);
}

bool IsAuditObjectiveAcceptableForProgress(
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

static std::optional<double> BuildFixedObjectiveScale(
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

ObjectiveDomain BuildObjectiveDomain(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<ClusterKey> & cluster_key_list)
{
    ObjectiveDomain domain;
    domain.owner_key_by_atom_index.resize(context.atom_list.size());
    domain.fit_sample_mask_by_atom.resize(context.atom_list.size());
    domain.tail_sample_mask_by_atom.resize(context.atom_list.size());
    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        domain.fit_sample_mask_by_atom.at(atom_index).resize(
            context.atom_list.at(atom_index).raw_sampling_entries.size(), 0);
        domain.tail_sample_mask_by_atom.at(atom_index).resize(
            context.atom_list.at(atom_index).raw_sampling_entries.size(), 0);
    }
    for (const auto & key : cluster_key_list)
    {
        auto & cluster_domain{ domain.cluster_by_key[key] };
        cluster_domain.selected_atom_count = key.size();
        domain.active_atom_count += cluster_domain.selected_atom_count;
        std::vector<double> fit_residual_list;
        std::vector<double> fit_response_list;
        std::vector<double> tail_residual_list;
        std::vector<double> tail_response_list;
        for (const auto atom_index : key)
        {
            domain.owner_key_by_atom_index.at(atom_index) = key;
            const auto & raw_sampling_entries{ context.atom_list.at(atom_index).raw_sampling_entries };
            for (std::size_t sample_index = 0; sample_index < raw_sampling_entries.size(); sample_index++)
            {
                const SampleRef sample_ref{ atom_index, sample_index };
                const auto distance{ raw_sampling_entries.at(sample_index).point.distance };
                const auto in_fit{ IsSignalDistance(distance) };
                const auto in_tail{ IsTailDistance(distance) };
                domain.fit_sample_mask_by_atom.at(atom_index).at(sample_index) = in_fit;
                domain.tail_sample_mask_by_atom.at(atom_index).at(sample_index) = in_tail;
                if (!in_fit && !in_tail) continue;
                cluster_domain.sample_ref_list.emplace_back(sample_ref);
                const auto residual_sample{
                    EvaluateResidualSample(context, sample_ref, model_snapshot)
                };
                for (const bool is_fit_range : { true, false })
                {
                    if (!(is_fit_range ? in_fit : in_tail)) continue;
                    auto & sample_ref_list{ is_fit_range ?
                        cluster_domain.fit_sample_ref_list : cluster_domain.tail_sample_ref_list
                    };
                    sample_ref_list.emplace_back(sample_ref);
                    if (!residual_sample.has_value()) continue;
                    auto & residual_list{ is_fit_range ? fit_residual_list : tail_residual_list };
                    auto & response_list{ is_fit_range ? fit_response_list : tail_response_list };
                    residual_list.emplace_back(residual_sample->residual);
                    response_list.emplace_back(residual_sample->adjusted_response);
                }
            }
        }
        domain.unique_sample_count += cluster_domain.sample_ref_list.size();
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

std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator)
{
    return EvaluateAuditObjectiveImpl(domain, evaluator.GetState(), evaluator);
}

std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot)
{
    return EvaluateAuditObjectiveImpl(domain, model_snapshot.node,
        [&](const SampleRef & sample_ref)
        {
            return EvaluateResidualSample(context, sample_ref, model_snapshot);
        });
}

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator)
{
    return BuildObjectiveByKeyImpl(partition, domain, evaluator.GetState(), evaluator);
}

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot)
{
    return BuildObjectiveByKeyImpl(partition, domain, model_snapshot.node,
        [&](const SampleRef & sample_ref)
        {
            return EvaluateResidualSample(context, sample_ref, model_snapshot);
        });
}

std::optional<ObjectiveBreakdown> EvaluateObjectiveDelta(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown & baseline,
    PerformanceCounters & performance_counters)
{
    const auto & changed_key{
        candidate_overlay.GetState().GetOverrideAtomIndexList()
    };
    const auto unique_sample_count{ domain.unique_sample_count };
    performance_counters.RecordObjectiveSampleEvaluation(
        CountObjectiveSamples(affected_sample_ref_list, domain),
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

std::optional<ObjectiveBreakdown> EvaluateCombinedObjective(
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

bool TryUpdateBestAuditState(
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

void ReevaluateBestAuditState(
    const SecondStageContext & context,
    const ObjectiveDomain & domain,
    BestAuditState & audit_state)
{
    if (!audit_state.has_value()) return;
    const auto snapshot{ BuildSecondStageModelSnapshot(context, audit_state->state) };
    const auto objective{ EvaluateAuditObjective(domain, context, snapshot) };
    if (objective.has_value()) audit_state->objective = *objective;
    else audit_state.reset();
}

void ReconcileClusterObjectiveState(
    const ObjectiveByKey & previous_objective_by_key,
    const FitState & accepted_state,
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
        next_state_by_key.emplace(
            key,
            ClusterObjectiveState{ .best_objective = previous_objective,
                .best_parameters = FitStatePatch::FromState(accepted_state, key) });
    }
    state_by_key = std::move(next_state_by_key);
}

FitStatePatch CaptureClusterParameters(const FitStateView & state, const ClusterKey & key)
{
    FitStatePatch patch{ .atom_index_list = key };
    patch.mdpde_list.reserve(key.size());
    for (const auto atom : key) patch.mdpde_list.emplace_back(state.GetMdpde(atom));
    return patch;
}

std::optional<ObjectiveBreakdown> EvaluateBestObjectiveReference(
    const CandidateEvaluationOverlay & candidate,
    const ClusterKey & key,
    const std::vector<SampleRef> & samples,
    const ObjectiveDomain & domain,
    const ClusterObjectiveState & state,
    PerformanceCounters & performance_counters)
{
    if (!state.best_objective) return std::nullopt;
    if (state.best_parameters.atom_index_list != key ||
        state.best_parameters.mdpde_list.size() != key.size())
    {
        throw std::logic_error("Cluster best objective parameters are inconsistent.");
    }
    // Preserve every candidate neighbor; replace only this cluster's parameters.
    ClusterKey merged_key;
    std::ranges::set_union(candidate.GetState().GetOverrideAtomIndexList(), key,
        std::back_inserter(merged_key));
    auto patch{ CaptureClusterParameters(candidate.GetState(), merged_key) };
    for (std::size_t i = 0; i < patch.atom_index_list.size(); i++)
    {
        if (const auto * best = state.best_parameters.Find(patch.atom_index_list.at(i)))
            patch.mdpde_list.at(i) = *best;
    }
    const CandidateEvaluationOverlay reference{
        candidate.GetContext(), candidate.GetBaseline(), candidate.GetState().GetBaseState(), patch };
    performance_counters.RecordObjectiveSampleEvaluation(
        CountObjectiveSamples(samples, domain), domain.unique_sample_count);
    return EvaluateObjectiveContribution(reference, key, samples, domain);
}

bool TryCommitClusterCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveBreakdown * previous_objective,
    bool requires_strict_improvement,
    const ObjectiveDomain & domain,
    ClusterObjectiveState & objective_state,
    ObjectiveAttemptDiagnostic & diagnostic,
    PerformanceCounters & performance_counters,
    std::string_view source)
{
    const auto evaluation{ EvaluateCandidate(candidate_overlay,
        requires_strict_improvement ? CandidateScope::LocalPolish : CandidateScope::LocalSearch,
        LocalCandidateReference{key, objective_sample_ref_list, previous_objective, domain,
            objective_state, diagnostic, performance_counters, source}) };
    diagnostic = evaluation.diagnostic;
    if (evaluation.accepted) objective_state = *evaluation.objective_state;
    return evaluation.accepted;
}

} // namespace rhbm_gem::core::detail
