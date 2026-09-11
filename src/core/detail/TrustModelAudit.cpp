#include "core/detail/TrustModelAudit.hpp"
#include "core/detail/IterationProcess.hpp"
#include "core/detail/GaussianModelOperations.hpp"
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace rhbm_gem::core::detail {
namespace {
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
std::optional<double> EvaluateTrustModelResponseDirection(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    double distance)
{
    const auto previous_coordinates{ previous_model.ToTransformedCoordinates() };
    const auto candidate_coordinates{ candidate_model.ToTransformedCoordinates() };
    if (!previous_coordinates.has_value() || !candidate_coordinates.has_value())
    {
        return std::nullopt;
    }
    const auto direction{ *candidate_coordinates - *previous_coordinates };
    if (!direction.allFinite()) return std::nullopt;
    if (direction.isZero()) return 0.0;
    const auto invariants{ BuildTransformedModelInvariants(previous_model) };
    if (!invariants.has_value()) return std::nullopt;
    const auto jacobian{ EvaluateTransformedJacobian(*invariants, distance) };
    if (!jacobian.has_value()) return std::nullopt;
    const auto response_direction{ jacobian->dot(direction) };
    return std::isfinite(response_direction) ?
        std::optional<double>{ response_direction } : std::nullopt;
}
#endif
}
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
TrustRegionRadiusAction DetermineTrustModelShadowAction(
    const TrustModelShadowDiagnostic & diagnostic)
{
    if (diagnostic.objective_backtracked)
    {
        return TrustRegionRadiusAction::Shrink;
    }
    if (diagnostic.status != TrustModelPredictionStatus::Available ||
        !diagnostic.rho.has_value())
    {
        return diagnostic.current_action;
    }
    if (*diagnostic.rho < 0.25)
    {
        return TrustRegionRadiusAction::Shrink;
    }
    if (*diagnostic.rho > 0.75 &&
        std::isfinite(diagnostic.boundary_utilization) &&
        diagnostic.boundary_utilization >= 0.8)
    {
        return TrustRegionRadiusAction::Grow;
    }
    return TrustRegionRadiusAction::Keep;
}

TrustModelShadowDiagnostic EvaluateTrustModelShadow(
    const SecondStageContext & context,
    const ResidualBaseline & residual_baseline,
    const FitState & previous_state,
    const FitStatePatch & candidate_patch,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveDomain & objective_domain,
    const std::optional<ObjectiveBreakdown> & previous_objective,
    const std::optional<ObjectiveBreakdown> & candidate_objective,
    double trust_region_radius,
    TrustRegionRadiusAction current_action,
    TrustModelCandidateSource candidate_source,
    bool objective_backtracked)
{
    TrustModelShadowDiagnostic result{
        .candidate_source = candidate_source,
        .current_action = current_action,
        .objective_backtracked = objective_backtracked
    };
    if (!previous_objective.has_value() || !candidate_objective.has_value())
    {
        result.status = TrustModelPredictionStatus::ObjectiveUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    const auto previous_objective_value{ previous_objective->GetTotalObjective() };
    const auto candidate_objective_value{ candidate_objective->GetTotalObjective() };
    if (!std::isfinite(previous_objective_value) ||
        !std::isfinite(candidate_objective_value))
    {
        result.status = TrustModelPredictionStatus::Nonfinite;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.actual_reduction = previous_objective_value - candidate_objective_value;

    const FitStateView candidate_state{ previous_state, candidate_patch };
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> candidate_model_list;
    previous_model_list.reserve(key.size());
    candidate_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
        candidate_model_list.emplace_back(candidate_state.GetModel(atom_index));
    }
    const auto step_norm{
        CalculateModelTrustRegionStepNorm(previous_model_list, candidate_model_list)
    };
    if (!step_norm.has_value())
    {
        result.status = TrustModelPredictionStatus::ModelUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.step_norm = *step_norm;
    result.boundary_utilization =
        std::isfinite(trust_region_radius) && trust_region_radius > 0.0 ?
            *step_norm / trust_region_radius : 0.0;
    if (*step_norm < kTransformedChangeTolerance)
    {
        result.status = TrustModelPredictionStatus::NonmaterialStep;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    if (objective_domain.active_atom_count == 0)
    {
        result.status = TrustModelPredictionStatus::ResidualUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }

    SecondStageModelSnapshot candidate_snapshot;
    try
    {
        candidate_snapshot = BuildSecondStageModelSnapshot(context, candidate_state);
    }
    catch (const std::exception &)
    {
        result.status = TrustModelPredictionStatus::ModelUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }

    double predicted_residual_reduction{ 0.0 };
    for (const auto & sample_ref : objective_sample_ref_list)
    {
        const auto & owner_key{
            objective_domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
        const auto in_fit{
            objective_domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        const auto in_tail{
            objective_domain.tail_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        if (!in_fit && !in_tail) continue;
        const auto owner_iter{ objective_domain.cluster_by_key.find(owner_key) };
        const auto previous_residual{ residual_baseline(sample_ref) };
        if (owner_iter == objective_domain.cluster_by_key.end() ||
            !owner_iter->second.scale.has_value() ||
            !previous_residual.has_value())
        {
            result.status = TrustModelPredictionStatus::ResidualUnavailable;
            result.shadow_action = DetermineTrustModelShadowAction(result);
            return result;
        }
        const auto & atom_context{ context.atom_list.at(sample_ref.atom_index) };
        const auto target_direction{
            EvaluateTrustModelResponseDirection(
                residual_baseline.model_snapshot.node.at(sample_ref.atom_index),
                candidate_snapshot.node.at(sample_ref.atom_index),
                atom_context.raw_sampling_entries.at(sample_ref.sample_index)
                    .point.distance)
        };
        if (!target_direction.has_value())
        {
            result.status = TrustModelPredictionStatus::ModelUnavailable;
            result.shadow_action = DetermineTrustModelShadowAction(result);
            return result;
        }
        double residual_direction{ -*target_direction };
        for (const auto & neighbor : atom_context.Neighbors(sample_ref.sample_index))
        {
            const auto & previous_neighbor{
                GetFitModel(residual_baseline.model_snapshot.node, neighbor.atom_index)
            };
            const auto & candidate_neighbor{
                GetFitModel(candidate_snapshot.node, neighbor.atom_index)
            };
            const auto neighbor_direction{
                EvaluateTrustModelResponseDirection(
                    previous_neighbor,
                    candidate_neighbor,
                    neighbor.distance)
            };
            if (!neighbor_direction.has_value())
            {
                result.status = TrustModelPredictionStatus::ModelUnavailable;
                result.shadow_action = DetermineTrustModelShadowAction(result);
                return result;
            }
            residual_direction -= *neighbor_direction;
        }
        const auto linearized_residual{
            previous_residual->residual + residual_direction
        };
        for (const bool is_fit_range : { true, false })
        {
            if (!(is_fit_range ? in_fit : in_tail)) continue;
            const auto sample_count{ is_fit_range ?
                owner_iter->second.fit_sample_ref_list.size() :
                owner_iter->second.tail_sample_ref_list.size()
            };
            const auto scale{ is_fit_range ?
                owner_iter->second.scale->fit : owner_iter->second.scale->tail
            };
            if (sample_count == 0 || !std::isfinite(scale) || scale <= 0.0)
            {
                result.status = TrustModelPredictionStatus::ResidualUnavailable;
                result.shadow_action = DetermineTrustModelShadowAction(result);
                return result;
            }

            const auto weight{
                algorithm::CalculateCauchyWeight(
                    previous_residual->residual,
                    scale,
                    kObjectiveRobustLossCutoffMultiplier)
            };
            const auto coefficient{
                CalculateClusterAtomWeight(
                    owner_iter->second.selected_atom_count,
                    objective_domain.active_atom_count) /
                static_cast<double>(sample_count)
            };
            const auto range_weight{ is_fit_range ? kFitRangeWeight : kTailValidationWeight };
            const auto previous_normalized{ previous_residual->residual / scale };
            const auto linearized_normalized{ linearized_residual / scale };
            const auto contribution{
                0.5 * range_weight * coefficient * weight *
                (previous_normalized * previous_normalized -
                    linearized_normalized * linearized_normalized)
            };
            if (!std::isfinite(contribution))
            {
                result.status = TrustModelPredictionStatus::Nonfinite;
                result.shadow_action = DetermineTrustModelShadowAction(result);
                return result;
            }
            predicted_residual_reduction += contribution;
        }
    }
    result.predicted_residual_reduction = predicted_residual_reduction;
    result.predicted_penalty_reduction =
        previous_objective->offset_plausibility_penalty -
        candidate_objective->offset_plausibility_penalty;
    const auto predicted_reduction{
        *result.predicted_residual_reduction +
        *result.predicted_penalty_reduction
    };
    if (!std::isfinite(predicted_reduction) ||
        !result.actual_reduction.has_value() ||
        !std::isfinite(*result.actual_reduction))
    {
        result.status = TrustModelPredictionStatus::Nonfinite;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.predicted_reduction = predicted_reduction;
    if (predicted_reduction <= 0.0)
    {
        result.status = TrustModelPredictionStatus::NonpositivePrediction;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    if (predicted_reduction <= CalculateObjectiveTolerance(
            previous_objective_value,
            kObjectiveProgressTolerance))
    {
        result.status = TrustModelPredictionStatus::NonmaterialPrediction;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    const auto rho{ *result.actual_reduction / predicted_reduction };
    if (!std::isfinite(rho))
    {
        result.status = TrustModelPredictionStatus::Nonfinite;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.status = TrustModelPredictionStatus::Available;
    result.rho = rho;
    result.shadow_action = DetermineTrustModelShadowAction(result);
    return result;
}
#endif
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
namespace {

std::string_view GetTrustModelPredictionStatusText(
    TrustModelPredictionStatus status)
{
    switch (status)
    {
    case TrustModelPredictionStatus::Available: return "available";
    case TrustModelPredictionStatus::NonmaterialStep: return "nonmaterial-step";
    case TrustModelPredictionStatus::ObjectiveUnavailable: return "objective-unavailable";
    case TrustModelPredictionStatus::ModelUnavailable: return "model-unavailable";
    case TrustModelPredictionStatus::ResidualUnavailable: return "residual-unavailable";
    case TrustModelPredictionStatus::Nonfinite: return "nonfinite";
    case TrustModelPredictionStatus::NonpositivePrediction: return "nonpositive-prediction";
    case TrustModelPredictionStatus::NonmaterialPrediction: return "nonmaterial-prediction";
    }
    return "model-unavailable";
}

std::string_view GetTrustModelCandidateSourceText(
    TrustModelCandidateSource source)
{
    return source == TrustModelCandidateSource::Polish ? "polish" : "base";
}

std::string_view GetTrustModelTrialDispositionText(
    TrustModelTrialDisposition disposition)
{
    return disposition == TrustModelTrialDisposition::Accepted ?
        "accepted" : "objective-rejected";
}

std::string_view GetTrustRegionRadiusActionText(
    TrustRegionRadiusAction action)
{
    switch (action)
    {
    case TrustRegionRadiusAction::Keep: return "keep";
    case TrustRegionRadiusAction::Grow: return "grow";
    case TrustRegionRadiusAction::Shrink: return "shrink";
    }
    return "keep";
}

void AppendTrustModelOptionalValue(
    std::ostringstream & stream,
    const std::optional<double> & value)
{
    if (value.has_value()) stream << *value;
    else stream << "-";
}

} // namespace

void TrustModelAudit::Log(
    bool quiet_mode,
    const IterationResult & iteration_result) const
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    const auto log_records = [&](
        const auto & diagnostic_list,
        std::string_view disposition)
    {
        for (const auto & cluster_diagnostic : diagnostic_list)
        {
            const auto & observed{ records.at(cluster_diagnostic.key) };
            const auto & funnel{ observed.funnel };
            std::ostringstream funnel_message;
            funnel_message
                << "Trust-model funnel: schema=1"
                << ", try=" << iteration_result.attempt_number
                << ", acc=" << iteration_result.accepted_iteration_count
                << ", atoms=" << cluster_diagnostic.key.size()
                << ", key-first=" << cluster_diagnostic.key.front()
                << ", key-last=" << cluster_diagnostic.key.back()
                << ", disposition=" << disposition
                << ", generated=" << funnel.generated_count
                << ", invalid=" << funnel.invalid_count
                << ", trust-skipped=" << funnel.trust_skipped_count
                << ", guard-rejected=" << funnel.guard_rejected_count
                << ", nonmaterial=" << funnel.nonmaterial_count
                << ", objective-evaluated=" << funnel.objective_evaluated_count
                << ", polish-objective-evaluated="
                << funnel.polish_objective_evaluated_count;
            Logger::Log(LogLevel::Debug, funnel_message.str());

            for (const auto & diagnostic :
                observed.trials)
            {
                std::ostringstream message;
                message << std::scientific << std::setprecision(17)
                    << "Trust-model shadow: schema=3"
                    << ", try=" << iteration_result.attempt_number
                    << ", acc=" << iteration_result.accepted_iteration_count
                    << ", atoms=" << cluster_diagnostic.key.size()
                    << ", key-first=" << cluster_diagnostic.key.front()
                    << ", key-last=" << cluster_diagnostic.key.back()
                    << ", disposition=" << disposition
                    << ", boundary-touched=" << observed.boundary_touched
                    << ", boundary-rescued=" << cluster_diagnostic.boundary_rescued
                    << ", readiness-eligible=" << diagnostic.readiness_eligible
                    << ", final-local-candidate=" << diagnostic.final_local_candidate
                    << ", status=" << GetTrustModelPredictionStatusText(diagnostic.status)
                    << ", source=" << GetTrustModelCandidateSourceText(
                        diagnostic.candidate_source)
                    << ", search-pass=" << diagnostic.search_pass
                    << ", trial=" << diagnostic.trial_number
                    << ", factor=" << diagnostic.factor
                    << ", trial-disposition=" << GetTrustModelTrialDispositionText(
                        diagnostic.trial_disposition)
                    << ", rejected-by-previous=" << diagnostic.rejected_by_previous
                    << ", rejected-by-strict-polish="
                    << diagnostic.rejected_by_strict_polish
                    << ", step-norm=" << diagnostic.step_norm
                    << ", actual-reduction=";
                AppendTrustModelOptionalValue(message, diagnostic.actual_reduction);
                message << ", polish-reduction=";
                AppendTrustModelOptionalValue(message, diagnostic.polish_reduction);
                message << ", predicted-residual-reduction=";
                AppendTrustModelOptionalValue(
                    message, diagnostic.predicted_residual_reduction);
                message << ", predicted-penalty-reduction=";
                AppendTrustModelOptionalValue(
                    message, diagnostic.predicted_penalty_reduction);
                message << ", predicted-reduction=";
                AppendTrustModelOptionalValue(message, diagnostic.predicted_reduction);
                message << ", rho=";
                AppendTrustModelOptionalValue(message, diagnostic.rho);
                message
                    << ", boundary-utilization=" << diagnostic.boundary_utilization
                    << ", current-action="
                    << GetTrustRegionRadiusActionText(diagnostic.current_action)
                    << ", shadow-action=";
                if (diagnostic.shadow_action.has_value())
                {
                    message << GetTrustRegionRadiusActionText(
                        *diagnostic.shadow_action);
                }
                else
                {
                    message << "suppressed";
                }
                message
                    << ", objective-backtracked=" << diagnostic.objective_backtracked
                    << ", unselected-dependencies=0"
                    << ", elapsed-ms=" << diagnostic.elapsed_milliseconds;
                Logger::Log(LogLevel::Debug, message.str());
            }
        }
    };
    Logger::FinishProgressLine();
    log_records(iteration_result.accepted_cluster_diagnostic_list, "accepted");
    log_records(iteration_result.rejected_cluster_diagnostic_list, "rejected");
}
#endif
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
TrustModelAudit::TrustModelAudit(const std::vector<ClusterKey> & keys)
{
    for (const auto & key : keys) records.try_emplace(key);
}

TrustModelTrialObserver::TrustModelTrialObserver(const CandidateSelectionInputs & inputs,
    const ClusterKey & key, const std::vector<SampleRef> & samples)
    : inputs(inputs), key(key), samples(samples), record(inputs.context.trust_model_audit->records.at(key)) {}

void TrustModelTrialObserver::Trial(const FitStatePatch & patch,
    const ObjectiveAttemptDiagnostic & diagnostic, bool polish, double factor, bool accepted)
{
    if (polish) ++record.funnel.polish_objective_evaluated_count;
    else ++record.funnel.objective_evaluated_count;
    const auto start{ std::chrono::steady_clock::now() };
    auto shadow{ EvaluateTrustModelShadow(inputs.context, inputs.residual_baseline,
        inputs.previous_state, patch, key, samples, inputs.objective_domain,
        inputs.previous_objective_by_key.at(key), diagnostic.candidate_objective,
        inputs.trust_region_state.GetRadius(key), TrustRegionRadiusAction::Keep,
        polish ? TrustModelCandidateSource::Polish : TrustModelCandidateSource::Base, false) };
    shadow.search_pass = search_pass;
    shadow.trial_number = polish ? 1 : diagnostic.trial_count;
    shadow.factor = factor;
    shadow.trial_disposition = accepted ? TrustModelTrialDisposition::Accepted : TrustModelTrialDisposition::ObjectiveRejected;
    shadow.rejected_by_previous = diagnostic.rejected_by_previous;
    shadow.rejected_by_strict_polish = polish && !accepted && diagnostic.candidate_objective &&
        diagnostic.previous_objective && !diagnostic.rejected_by_previous;
    if (polish && diagnostic.previous_objective && diagnostic.candidate_objective)
        shadow.polish_reduction = diagnostic.previous_objective->GetTotalObjective() - diagnostic.candidate_objective->GetTotalObjective();
    shadow.shadow_action.reset();
    shadow.elapsed_milliseconds = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    record.trials.emplace_back(std::move(shadow));
    if (accepted) final_trial = record.trials.size() - 1;
}

void TrustModelTrialObserver::Finish(bool shrink_trust_region,
    std::optional<double> first_factor, const ObjectiveAttemptDiagnostic & diagnostic)
{
    if (!final_trial) return;
    auto & shadow{ record.trials.at(*final_trial) };
    shadow.final_local_candidate = true;
    shadow.readiness_eligible = true;
    shadow.current_action = shrink_trust_region ?
        TrustRegionRadiusAction::Shrink : TrustRegionRadiusAction::Keep;
    shadow.objective_backtracked = first_factor && diagnostic.accepted_factor && *diagnostic.accepted_factor < *first_factor;
    shadow.shadow_action = DetermineTrustModelShadowAction(shadow);
}

void TrustModelAudit::Finalize(const CandidateSelection & selection)
{
    const auto update = [&](const ClusterCandidateDiagnostic & diagnostic, bool accepted)
    {
        auto & observed{ records.at(diagnostic.key) };
        observed.boundary_touched = std::ranges::any_of(selection.boundary_reconciliation_diagnostic_list,
            [&](const auto & boundary) { return std::ranges::find(boundary.key_list, diagnostic.key) != boundary.key_list.end(); });
        for (auto & shadow : observed.trials)
        {
            shadow.readiness_eligible = shadow.final_local_candidate && accepted && !observed.boundary_touched && !diagnostic.boundary_rescued;
            if (!shadow.readiness_eligible) shadow.shadow_action.reset();
        }
    };
    for (const auto & diagnostic : selection.accepted_cluster_diagnostic_list) update(diagnostic, true);
    for (const auto & diagnostic : selection.rejected_cluster_diagnostic_list) update(diagnostic, false);
}

void BeginTrustModelAudit(SecondStageContext & context, const std::vector<ClusterKey> & keys)
{
    context.trust_model_audit = std::make_shared<TrustModelAudit>(keys);
}
void FinalizeTrustModelAudit(const SecondStageContext & context, const CandidateSelection & selection)
{
    context.trust_model_audit->Finalize(selection);
}
void LogTrustModelAudit(const SecondStageContext & context, bool quiet, const IterationResult & result)
{
    context.trust_model_audit->Log(quiet, result);
}
#endif
} // namespace rhbm_gem::core::detail
