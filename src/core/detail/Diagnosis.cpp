#include "core/detail/Diagnosis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr std::array<std::string_view, 6> kProgressHeaderList{
    "Try/Acc",
    "Atom A/Q",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/G"
};

using GraphEdgeSet = std::set<std::pair<std::size_t, std::size_t>>;

const char * GetSecondStageSeedSourceText(SecondStageSeedSource source)
{
    switch (source)
    {
    case SecondStageSeedSource::GroupPosterior:
        return "group-posterior";
    case SecondStageSeedSource::GroupPrior:
        return "group-prior";
    case SecondStageSeedSource::GroupMedian:
        return "group-median";
    case SecondStageSeedSource::GlobalMedian:
        return "global-median";
    }
    throw std::logic_error("Unknown second-stage seed source.");
}

std::string_view GetSecondStageStopReasonText(SecondStageStopReason reason)
{
    switch (reason)
    {
    case SecondStageStopReason::None:
        return "";
    case SecondStageStopReason::Quarantine:
        return "quarantine";
    case SecondStageStopReason::Converged:
        return "converged";
    case SecondStageStopReason::AuditPatience:
        return "audit-patience";
    case SecondStageStopReason::AllRejectedBacktrackingExhausted:
        return "all-rejected-backtracking-exhausted";
    case SecondStageStopReason::AllRejectedAtMaximumIterations:
    case SecondStageStopReason::MaximumIterations:
        return "maximum-iterations";
    }
    return "";
}

std::string_view GetAdaptiveTopologyTriggerText(
    AdaptiveTopologyRebuildTrigger trigger)
{
    switch (trigger)
    {
    case AdaptiveTopologyRebuildTrigger::Drift: return "drift";
    case AdaptiveTopologyRebuildTrigger::Interval: return "interval";
    case AdaptiveTopologyRebuildTrigger::None: return "none";
    }
    return "unknown";
}

std::string_view GetFinalPolishCertificationPolicyText(
    FinalPolishCertificationPolicy policy)
{
    switch (policy)
    {
    case FinalPolishCertificationPolicy::RequireResidualNonRegression:
        return "non-regression";
    case FinalPolishCertificationPolicy::RequireStrictFixedPoint:
        return "strict-fixed-point";
    }
    throw std::invalid_argument("Unknown final polish certification policy.");
}

std::string_view GetFinalPolishResidualSafetyStatusText(
    FinalPolishResidualSafetyStatus status)
{
    switch (status)
    {
    case FinalPolishResidualSafetyStatus::NotEvaluated:
        return "not-evaluated";
    case FinalPolishResidualSafetyStatus::AbsolutePassed:
        return "absolute-passed";
    case FinalPolishResidualSafetyStatus::RelativePassed:
        return "relative-passed";
    case FinalPolishResidualSafetyStatus::Failed:
        return "failed";
    case FinalPolishResidualSafetyStatus::Error:
        return "error";
    }
    throw std::invalid_argument("Unknown final polish residual safety status.");
}

void AppendQuarantineSummary(
    std::ostream & stream,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count)
{
    stream << "; quarantine entered/released/probation-failed/unresolved = "
        << entered_target_count << "/"
        << released_target_count << "/"
        << failed_probation_count << "/"
        << unresolved_target_count;
}

void AppendObjectiveBreakdown(
    std::ostringstream & stream,
    const std::optional<ObjectiveBreakdown> & breakdown,
    std::string_view unavailable_text = "unavailable")
{
    if (!breakdown.has_value())
    {
        stream << unavailable_text;
        return;
    }
    stream
        << breakdown->fit_range_residual_objective << "/"
        << breakdown->GetTailValidationPenalty() << "/"
        << breakdown->offset_plausibility_penalty << "/"
        << breakdown->GetTotalObjective();
}

std::string_view GetPreObjectiveFailureReasonText(
    PreObjectiveFailureReason reason)
{
    switch (reason)
    {
    case PreObjectiveFailureReason::None:
        return "none";
    case PreObjectiveFailureReason::InvalidModel:
        return "invalid-model";
    case PreObjectiveFailureReason::NoCandidateWithinTrustRegion:
        return "no-candidate-within-trust-region";
    }
    return "unknown";
}

std::string_view GetSuspiciousGaussianReasonText(
    SuspiciousGaussianReason reason)
{
    switch (reason)
    {
    case SuspiciousGaussianReason::None: return "none";
    case SuspiciousGaussianReason::InvalidModel: return "invalid-model";
    case SuspiciousGaussianReason::NonFiniteResponse: return "non-finite-response";
    case SuspiciousGaussianReason::OffsetMagnitude: return "offset-magnitude";
    case SuspiciousGaussianReason::CenterSignFlip: return "center-sign-flip";
    case SuspiciousGaussianReason::RadialRebound: return "radial-rebound";
    case SuspiciousGaussianReason::WidthGrowth: return "width-growth";
    case SuspiciousGaussianReason::AmplitudeOffsetCompensation:
        return "amplitude-offset-compensation";
    }
    return "none";
}

std::string_view GetStabilizationTerminalReasonText(
    StabilizationTerminalReason reason)
{
    switch (reason)
    {
    case StabilizationTerminalReason::None: return "none";
    case StabilizationTerminalReason::GuardInfeasible: return "guard-infeasible";
    case StabilizationTerminalReason::ObjectiveExhausted: return "objective-exhausted";
    case StabilizationTerminalReason::InvalidCandidate: return "invalid-candidate";
    }
    return "unknown";
}

void AppendAuditValues(
    std::ostringstream & message,
    const TransformedChange & value_list)
{
    for (std::size_t i = 0; i < value_list.size(); i++)
    {
        if (i != 0) message << "/";
        message << value_list.at(i);
    }
}

void AppendAuditPopulation(
    std::ostringstream & message,
    const std::array<
        std::size_t,
        GaussianModel3D::TransformedCoordinateSize()> & population_size_list)
{
    for (std::size_t i = 0; i < population_size_list.size(); i++)
    {
        if (i != 0) message << "/";
        message << population_size_list.at(i);
    }
}

GraphEdgeSet BuildGraphEdgeSet(const GraphTopology & topology)
{
    GraphEdgeSet edge_set;
    for (std::size_t atom_index = 0;
        atom_index < topology.adjacency_list.size(); atom_index++)
    {
        for (const auto neighbor_index : topology.adjacency_list.at(atom_index))
        {
            if (atom_index < neighbor_index)
            {
                edge_set.emplace(atom_index, neighbor_index);
            }
        }
    }
    return edge_set;
}

std::size_t CountGraphEdgeDifference(
    const GraphEdgeSet & source,
    const GraphEdgeSet & destination)
{
    return static_cast<std::size_t>(
        std::ranges::count_if(
            source,
            [&](const auto & edge)
            {
                return !destination.contains(edge);
            }));
}

void AppendOffsetSummary(std::ostringstream & stream, const FitState & state)
{
    std::size_t finite_count{ 0 };
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(state.size());
    for (const auto & result : state)
    {
        const auto offset{ result.mdpde.GetModel().GetOffset() };
        if (!std::isfinite(offset)) continue;
        finite_count++;
        absolute_offset_list.emplace_back(std::abs(offset));
    }
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
    if (!absolute_offset_list.empty())
    {
        median_absolute_offset = array_helper::ComputeMedian(absolute_offset_list);
        percentile_absolute_offset =
            array_helper::ComputePercentile(absolute_offset_list, 0.99);
        maximum_absolute_offset = std::ranges::max(absolute_offset_list);
    }
    stream << std::scientific << std::setprecision(2)
        << "; offsets finite = " << finite_count << " of " << state.size()
        << ", |C| median/p99/max = "
        << median_absolute_offset << "/"
        << percentile_absolute_offset << "/"
        << maximum_absolute_offset;
}

void AppendAuditSummary(
    std::ostringstream & stream,
    const AuditedState & audited_state)
{
    const auto & objective{ audited_state.objective };
    stream << "; audit best source = ";
    if (audited_state.source_iteration != 0)
    {
        stream << "accepted iteration " << audited_state.source_iteration;
    }
    else
    {
        stream << "initial";
    }
    stream
        << std::scientific << std::setprecision(2)
        << ", fixed audit objective fit/tail-weighted/offset/total = "
        << objective.fit_range_residual_objective << "/"
        << objective.GetTailValidationPenalty() << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.GetTotalObjective()
        << ", tail raw/weight = " << objective.tail_validation_loss << "/"
        << kTailValidationWeight;
}

void LogAuditAtomState(
    std::string_view marker,
    const SecondStageContext & context,
    const FitState & state)
{
    for (std::size_t atom_index = 0; atom_index < state.size(); atom_index++)
    {
        const auto & model{ state.at(atom_index).mdpde.GetModel() };
        std::ostringstream message;
        message << std::scientific << std::setprecision(17)
            << marker << " schema=1"
            << ", serial=" << context.at(atom_index).atom->GetSerialID()
            << ", group=" << context.at(atom_index).group_id
            << ", amplitude=" << model.GetAmplitude()
            << ", width=" << model.GetWidth()
            << ", offset=" << model.GetOffset();
        Logger::Log(LogLevel::Debug, message.str());
    }
}

} // namespace

bool IsDebugLogLevelEnabled()
{
    return Logger::GetLogLevel() >= LogLevel::Debug;
}

void FinishProgressLine(bool quiet_mode)
{
    if (!quiet_mode) Logger::FinishProgressLine();
}

void LogSecondStageStart(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Info,
        "Run 2nd-stage local atom fitting with iterations...");
}

void LogSecondStageInitializationFailure(
    bool quiet_mode,
    SecondStageInitializationFailure failure)
{
    if (quiet_mode || failure == SecondStageInitializationFailure::None) return;
    const auto unselected_seed_failure{
        failure == SecondStageInitializationFailure::UnselectedSeedUnavailable
    };
    Logger::Log(LogLevel::Warning,
        unselected_seed_failure ?
            "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
            "is available for every unselected neighbor atom." :
            "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
            "is available for every selected atom.");
    Logger::Log(LogLevel::Info,
        "Second-stage local fitting summary: accepted_iterations=0, "
        "best_iteration=unavailable, stop_reason=" +
        std::string(unselected_seed_failure ?
            "no-valid-unselected-neighbor-seed" : "no-valid-seed") +
        ", best_audit_objective=unavailable, final_uses_polish=unavailable, "
        "final_state_source=unavailable.");
}

void LogNoSelectedAtoms(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Skip 2nd-stage local atom fitting because no atoms are selected.");
}

void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    std::array<std::size_t, kSecondStageSeedSourceList.size()> source_count{};
    for (const auto & record : selection_record_list)
    {
        source_count.at(static_cast<std::size_t>(record.source))++;
    }

    std::ostringstream summary;
    summary << "Selected second-stage initial seeds = "
        << selection_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < kSecondStageSeedSourceList.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedSourceText(kSecondStageSeedSourceList.at(i))
            << ":" << source_count.at(i);
    }
    summary << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (std::size_t atom_index = 0;
        atom_index < selection_record_list.size(); atom_index++)
    {
        const auto & record{ selection_record_list.at(atom_index) };
        std::ostringstream detail_message;
        detail_message << "Second-stage seed selection: atom index = "
            << atom_index
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", original MDPDE A/B/C = "
            << record.original_model.GetAmplitude() << "/"
            << record.original_model.GetWidth() << "/"
            << record.original_model.GetOffset()
            << ", selected A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

void LogUnselectedSecondStageSeedSelections(
    const std::vector<UnselectedSecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    std::size_t group_median_count{ 0 };
    std::size_t global_median_count{ 0 };
    for (const auto & record : selection_record_list)
    {
        if (record.source == SecondStageSeedSource::GroupMedian)
        {
            group_median_count++;
        }
        else if (record.source == SecondStageSeedSource::GlobalMedian)
        {
            global_median_count++;
        }
    }

    std::ostringstream summary;
    summary << "Unselected second-stage neighbor seeds = "
        << selection_record_list.size()
        << ", sources = group-median:" << group_median_count
        << ", global-median:" << global_median_count << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        std::ostringstream detail_message;
        detail_message
            << "Unselected second-stage neighbor seed selection: serial ID = "
            << record.atom_serial_id
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", seed A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

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

void LogTrustModelShadowDiagnostics(
    bool quiet_mode,
    const IterationResult & iteration_result)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    const auto log_records = [&](
        const auto & diagnostic_list,
        std::string_view disposition)
    {
        for (const auto & cluster_diagnostic : diagnostic_list)
        {
            const auto & funnel{ cluster_diagnostic.trust_model_candidate_funnel };
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
                cluster_diagnostic.trust_model_shadow_trial_list)
            {
                std::ostringstream message;
                message << std::scientific << std::setprecision(17)
                    << "Trust-model shadow: schema=2"
                    << ", try=" << iteration_result.attempt_number
                    << ", acc=" << iteration_result.accepted_iteration_count
                    << ", atoms=" << cluster_diagnostic.key.size()
                    << ", key-first=" << cluster_diagnostic.key.front()
                    << ", key-last=" << cluster_diagnostic.key.back()
                    << ", disposition=" << disposition
                    << ", boundary-touched=" << cluster_diagnostic.boundary_touched
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
                    << ", rejected-by-best=" << diagnostic.rejected_by_best
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
                    << ", unselected-dependencies="
                    << diagnostic.unselected_dependency_count
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

void LogRejectedClusterDiagnostics(
    bool quiet_mode,
    const std::vector<ClusterCandidateDiagnostic> & diagnostic_list)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug ||
        diagnostic_list.empty())
    {
        return;
    }

    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : diagnostic_list)
    {
        std::ostringstream header;
        header
            << "Rejected local fitting cluster diagnostics: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/"
            << cluster_diagnostic.key.back()
            << ", breakdown order = fit/tail-weighted/offset/total";
        Logger::Log(LogLevel::Debug, header.str());

        const auto & diagnostic{ cluster_diagnostic.attempt };
        std::ostringstream message;
        message << std::scientific << std::setprecision(2)
            << "  unified accepted factor = ";
        if (diagnostic.accepted_factor.has_value())
        {
            message << *diagnostic.accepted_factor;
        }
        else
        {
            message << "-";
        }
        message
            << ", trials[total/invalid/trust-skipped/guard-rejected/objective-rejected] = "
            << diagnostic.trial_count << "/"
            << diagnostic.invalid_trial_count << "/"
            << diagnostic.trust_skipped_trial_count << "/"
            << diagnostic.guard_rejected_trial_count << "/"
            << diagnostic.objective_rejected_trial_count
            << ", terminal = "
            << (diagnostic.terminal_diagnostic_list.empty() ? "none" :
                GetStabilizationTerminalReasonText(
                    diagnostic.terminal_diagnostic_list.back().reason))
            << ", trust radius/step norm = "
            << diagnostic.trust_region_radius << "/";

        if (diagnostic.pre_objective_failure_reason !=
            PreObjectiveFailureReason::None)
        {
            if (diagnostic.pre_objective_attempted_step_norm.has_value())
            {
                message << *diagnostic.pre_objective_attempted_step_norm;
            }
            else
            {
                message << "unavailable";
            }
            message
                << ", status = "
                << GetPreObjectiveFailureReasonText(
                    diagnostic.pre_objective_failure_reason)
                << ", objective = not-evaluated";
            Logger::Log(LogLevel::Debug, message.str());
            continue;
        }

        message << diagnostic.trust_region_step_norm;
        message << ", fit/tail scales = ";
        if (diagnostic.scale.has_value())
        {
            message << diagnostic.scale->fit;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.scale.has_value() && diagnostic.tail_sample_count > 0)
        {
            message << diagnostic.scale->tail;
        }
        else
        {
            message << "empty";
        }
        message
            << ", fit/tail samples = "
            << diagnostic.fit_sample_count << "/"
            << diagnostic.tail_sample_count
            << ", tail raw/weight = ";
        if (diagnostic.candidate_objective.has_value())
        {
            message << diagnostic.candidate_objective->tail_validation_loss;
        }
        else
        {
            message << "unavailable";
        }
        message << "/" << kTailValidationWeight;
        message << ", candidate = ";
        AppendObjectiveBreakdown(message, diagnostic.candidate_objective);
        message << ", previous = ";
        AppendObjectiveBreakdown(message, diagnostic.previous_objective);
        message << ", best = ";
        AppendObjectiveBreakdown(message, diagnostic.best_objective);
        message << ", rejected-by = ";
        if (!diagnostic.candidate_objective.has_value())
        {
            message << "objective-unavailable";
        }
        else if (diagnostic.rejected_by_previous && diagnostic.rejected_by_best)
        {
            message << "previous+best";
        }
        else if (diagnostic.rejected_by_previous)
        {
            message << "previous";
        }
        else if (diagnostic.rejected_by_best)
        {
            message << "best";
        }
        else
        {
            message << "none";
        }
        Logger::Log(LogLevel::Debug, message.str());
    }
}

void LogAllRejectedResolution(
    bool quiet_mode,
    const IterationResult & iteration_result)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "All-rejected local fitting resolution: outcome = "
        << GetSecondStageStopReasonText(iteration_result.stop_reason)
        << ", radius-changed/radius-saturated = "
        << iteration_result.trust_region_update.changed_key_list.size() << "/"
        << iteration_result.trust_region_update.saturated_key_list.size() << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

void LogAcceptedCandidateSearchDiagnostics(
    bool quiet_mode,
    const IterationResult & iteration_result)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    const auto has_local_search{
        std::any_of(
            iteration_result.accepted_cluster_diagnostic_list.begin(),
            iteration_result.accepted_cluster_diagnostic_list.end(),
            [](const ClusterCandidateDiagnostic & diagnostic)
            {
                return diagnostic.attempt.trial_count > 1;
            })
    };
    const auto has_boundary_reconciliation_diagnostic{
        !iteration_result.boundary_reconciliation_diagnostic_list.empty()
    };
    if (!has_local_search && !has_boundary_reconciliation_diagnostic) return;
    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic :
        iteration_result.accepted_cluster_diagnostic_list)
    {
        const auto & diagnostic{ cluster_diagnostic.attempt };
        if (cluster_diagnostic.boundary_rescued)
        {
            std::ostringstream rescue_message;
            rescue_message
                << "Accepted local fitting cluster after boundary rescue: atoms = "
                << cluster_diagnostic.key.size()
                << ", key first/last = "
                << cluster_diagnostic.key.front() << "/"
                << cluster_diagnostic.key.back() << ".";
            Logger::Log(LogLevel::Debug, rescue_message.str());
        }
        if (diagnostic.trial_count <= 1) continue;
        std::ostringstream message;
        message
            << "Accepted local fitting candidate search: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/"
            << cluster_diagnostic.key.back()
            << ", trials/factor = " << diagnostic.trial_count << "/";
        if (diagnostic.accepted_factor.has_value())
        {
            message << *diagnostic.accepted_factor;
        }
        else
        {
            message << "-";
        }
        message << ", fixed fit/tail scales = ";
        if (diagnostic.scale.has_value())
        {
            message << diagnostic.scale->fit;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.scale.has_value() && diagnostic.tail_sample_count > 0)
        {
            message << diagnostic.scale->tail;
        }
        else
        {
            message << "empty";
        }
        message << ".";
        Logger::Log(LogLevel::Debug, message.str());
    }
    for (const auto & diagnostic :
        iteration_result.boundary_reconciliation_diagnostic_list)
    {
        const auto accepted_source_text = [&]() -> std::string_view
        {
            switch (diagnostic.accepted_source)
            {
            case BoundaryComponentAcceptedSource::None:
                return "none";
            case BoundaryComponentAcceptedSource::Endpoint:
                return "endpoint";
            case BoundaryComponentAcceptedSource::JointCorrection:
                return "joint-correction";
            case BoundaryComponentAcceptedSource::Backtracking:
                return "backtracking";
            }
            return "none";
        };
        std::ostringstream message;
        message
            << "Boundary-component reconciliation: clusters/atoms/boundary-samples = "
            << diagnostic.key_list.size() << "/"
            << diagnostic.atom_count << "/"
            << diagnostic.boundary_sample_count
            << ", accepted/rescue-candidates/rescued = "
            << diagnostic.accepted_cluster_count << "/"
            << diagnostic.rescue_candidate_cluster_count << "/"
            << diagnostic.rescued_cluster_count
            << ", locally-deteriorated/max-delta="
            << diagnostic.locally_deteriorated_member_count << "/"
            << diagnostic.maximum_local_deterioration
            << ", component/global-improvement=";
        if (diagnostic.component_improvement.has_value())
        {
            message << *diagnostic.component_improvement;
        }
        else
        {
            message << "-";
        }
        message << "/";
        if (diagnostic.global_improvement.has_value())
        {
            message << *diagnostic.global_improvement;
        }
        else
        {
            message << "-";
        }
        message
            << ", mode="
            << (diagnostic.is_rescue_attempt ? "rescue" : "accepted-only")
            << ", trials/factor/accepted/exhausted = "
            << diagnostic.trial_count << "/";
        if (diagnostic.accepted_factor.has_value())
        {
            message << *diagnostic.accepted_factor;
        }
        else
        {
            message << "-";
        }
        message
            << "/"
            << (diagnostic.accepted_source != BoundaryComponentAcceptedSource::None ?
                "yes" : "no")
            << "/" << (diagnostic.exhausted ? "yes" : "no")
            << ", accepted_source=" << accepted_source_text()
            << ", objectives previous/endpoint/final=";
        const auto append_objective = [&](const std::optional<double> & objective)
        {
            if (objective.has_value()) message << *objective;
            else message << "-";
        };
        append_objective(diagnostic.previous_component_objective);
        message << "/";
        append_objective(diagnostic.endpoint_component_objective);
        message << "/";
        append_objective(diagnostic.candidate_component_objective);
        message << ".";
        Logger::Log(LogLevel::Debug, message.str());
        if (!diagnostic.joint_correction_status.has_value()) continue;
        std::ostringstream correction_message;
        correction_message << std::scientific << std::setprecision(2)
            << "Boundary-interface joint correction: direct-interface/shape-active/offset-active/closure/parameters = "
            << diagnostic.interface_atom_count << "/"
            << diagnostic.shape_active_atom_count << "/"
            << diagnostic.offset_active_atom_count << "/"
            << diagnostic.offset_closure_atom_count << "/"
            << diagnostic.joint_parameter_count
            << ", status="
            << GetBoundaryJointCorrectionStatusText(
                *diagnostic.joint_correction_status)
            << ", damping/trust=";
        if (diagnostic.joint_damping.has_value())
        {
            correction_message << *diagnostic.joint_damping;
        }
        else
        {
            correction_message << "-";
        }
        correction_message << "/";
        if (diagnostic.maximum_normalized_trust_step.has_value())
        {
            correction_message << *diagnostic.maximum_normalized_trust_step;
        }
        else
        {
            correction_message << "-";
        }
        correction_message << ", reference/candidate=";
        if (diagnostic.joint_reference_component_objective.has_value())
        {
            correction_message << *diagnostic.joint_reference_component_objective;
        }
        else
        {
            correction_message << "-";
        }
        correction_message << "/";
        if (diagnostic.joint_candidate_component_objective.has_value())
        {
            correction_message << *diagnostic.joint_candidate_component_objective;
        }
        else
        {
            correction_message << "-";
        }
        const auto correction_accepted{
            diagnostic.accepted_source ==
                BoundaryComponentAcceptedSource::JointCorrection
        };
        std::string_view correction_outcome{ "failed" };
        if (correction_accepted)
        {
            correction_outcome =
                diagnostic.endpoint_component_objective.has_value() ?
                    "accepted-over-endpoint" : "accepted-over-previous";
        }
        else if (diagnostic.accepted_source ==
            BoundaryComponentAcceptedSource::Endpoint)
        {
            correction_outcome = "fallback-endpoint";
        }
        correction_message << ", suspicious="
            << diagnostic.suspicious_candidate_atom_count
            << ", accepted=" << (correction_accepted ? "yes" : "no")
            << ", outcome=" << correction_outcome << ".";
        Logger::Log(LogLevel::Debug, correction_message.str());
    }
}

namespace {

std::string FormatProgressMaximum(double value)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

template<typename CellList>
std::string FormatProgressRow(
    const ProgressColumnWidths & column_widths,
    const CellList & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream << std::left
            << std::setw(static_cast<int>(column_widths.at(i)))
            << cell_list.at(i);
    }
    return stream.str();
}

} // namespace

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count)
{
    const auto maximum_iteration_text{ std::to_string(kMaximumIterations) };
    const auto maximum_atom_text{ std::to_string(atom_count) };
    const auto maximum_change_text{
        FormatProgressMaximum(std::numeric_limits<double>::max())
    };
    const std::array<std::string, 6> maximum_cell_list{
        maximum_iteration_text + "/" + maximum_iteration_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text + "/" +
            maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text,
        maximum_change_text + "/" + maximum_change_text
    };
    ProgressColumnWidths column_widths;
    for (std::size_t i = 0; i < column_widths.size(); i++)
    {
        column_widths.at(i) = std::max(
            kProgressHeaderList.at(i).size(),
            maximum_cell_list.at(i).size());
    }
    return column_widths;
}

void LogProgressHeader(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths)
{
    if (quiet_mode) return;
    Logger::Log(
        LogLevel::Info,
        FormatProgressRow(column_widths, kProgressHeaderList));
}

void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result)
{
    if (quiet_mode) return;
    const std::array<std::string, 6> cell_list{
        std::to_string(iteration_result.attempt_number) + "/" +
            std::to_string(iteration_result.accepted_iteration_count),
        std::to_string(iteration_result.active_atom_count) + "/" +
            std::to_string(iteration_result.quarantine_atom_count),
        std::to_string(
            iteration_result.accepted_cluster_diagnostic_list.size()) + "/" +
            std::to_string(
                iteration_result.rejected_cluster_diagnostic_list.size()),
        std::to_string(iteration_result.polish_progress.eligible_count) + "/" +
            std::to_string(iteration_result.polish_progress.accepted_count) + "/" +
            std::to_string(iteration_result.polish_progress.rejected_count) + "/" +
            std::to_string(iteration_result.polish_progress.skipped_count),
        std::to_string(iteration_result.suspicious_atom_count),
        (iteration_result.accepted_maximum_transformed_change.has_value() ?
            FormatProgressMaximum(
                *iteration_result.accepted_maximum_transformed_change) :
            std::string{ "-" }) + "/" +
            FormatProgressMaximum(
                iteration_result.operator_maximum_transformed_change)
    };
    Logger::ProgressLine(FormatProgressRow(column_widths, cell_list));
}

void LogUnrestrictedOperatorAssessments(
    bool quiet_mode,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    const SuspiciousBlockActivity & block_activity)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    for (std::size_t atom_index = 0;
        atom_index < assessment_by_atom.size(); atom_index++)
    {
        const auto & assessment{ assessment_by_atom[atom_index] };
        const auto has_fixed_block{
            block_activity.shape_fixed_atom_mask.at(atom_index) != 0 ||
            block_activity.offset_fixed_atom_mask.at(atom_index) != 0 ||
            block_activity.hard_failure_atom_mask.at(atom_index) != 0
        };
        if (!assessment.IsSuspicious() && !has_fixed_block) continue;
        std::ostringstream message;
        message << std::scientific << std::setprecision(2)
            << "Unrestricted operator assessment: atom=" << atom_index
            << ", reason=" << GetSuspiciousGaussianReasonText(assessment.reason)
            << ", margin=" << assessment.normalized_margin
            << ", shape-fixed/offset-fixed/hard-fixed="
            << (block_activity.shape_fixed_atom_mask.at(atom_index) != 0 ?
                "yes" : "no")
            << "/"
            << (block_activity.offset_fixed_atom_mask.at(atom_index) != 0 ?
                "yes" : "no")
            << "/"
            << (block_activity.hard_failure_atom_mask.at(atom_index) != 0 ?
                "yes" : "no")
            << ".";
        Logger::Log(LogLevel::Debug, message.str());
    }
}

void LogConvergenceSafeguardAudit(
    bool quiet_mode,
    const IterationResult & iteration_result,
    const ConvergenceCertificate & certificate)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;

    const auto & accepted_production_change{
        certificate.accepted_active_movement
    };
    const auto accepted_percentile_passed{
        IsTransformedPercentileConverged(accepted_production_change)
    };
    const auto operator_percentile_passed{
        IsTransformedPercentileConverged(
            certificate.operator_nominal_residual)
    };
    const auto blockers_clear{
        !certificate.objective_domain_changed &&
        !certificate.quarantine_transition &&
        !certificate.suspicious_offset_fallback &&
        !certificate.rejected_cluster
    };
    const auto selected_atom_count{
        iteration_result.active_atom_count +
            iteration_result.quarantine_atom_count
    };
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "Convergence safeguard audit: schema=10"
        << ", try="
        << iteration_result.attempt_number
        << ", acc=" << iteration_result.accepted_iteration_count
        << ", atoms=" << selected_atom_count
        << ", quarantine=" << iteration_result.quarantine_atom_count
        << ", accepted-active-population=";
    AppendAuditPopulation(message, accepted_production_change.population_size_list);
    message << ", operator-nominal-population=";
    AppendAuditPopulation(
        message,
        certificate.operator_nominal_residual.population_size_list);
    message
        << ", certificate[solver/accepted-p99/operator-complete/operator-p99/blockers/production]="
        << certificate.solver_qualified << "/"
        << accepted_percentile_passed << "/"
        << certificate.operator_complete << "/"
        << operator_percentile_passed << "/"
        << blockers_clear << "/"
        << certificate.ProductionConverged()
        << ", accepted-active-p99=";
    AppendAuditValues(message, accepted_production_change.percentile_list);
    message << ", accepted-active-max=";
    AppendAuditValues(message, accepted_production_change.maximum_list);
    message << ", operator-nominal-residual-p99=";
    AppendAuditValues(
        message,
        certificate.operator_nominal_residual.percentile_list);
    message << ", operator-nominal-residual-max=";
    AppendAuditValues(
        message,
        certificate.operator_nominal_residual.maximum_list);
    message
        << ", blockers[objective-domain/quarantine-transition/suspicious-offset/rejected-cluster]="
        << certificate.objective_domain_changed << "/"
        << certificate.quarantine_transition << "/"
        << certificate.suspicious_offset_fallback << "/"
        << certificate.rejected_cluster << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Debug, message.str());
}

void LogAdaptiveTopologyRebuild(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const AdaptiveTopologyRebuildDecision & decision,
    const GraphTopology & previous_topology,
    const GraphTopology & rebuilt_topology,
    const CouplingGraphPartition & previous_partition,
    const CouplingGraphPartition & rebuilt_partition,
    bool partition_changed)
{
    if (quiet_mode) return;
    const auto previous_edge_set{ BuildGraphEdgeSet(previous_topology) };
    const auto rebuilt_edge_set{ BuildGraphEdgeSet(rebuilt_topology) };
    const auto removed_edge_count{
        CountGraphEdgeDifference(previous_edge_set, rebuilt_edge_set)
    };
    const auto added_edge_count{
        CountGraphEdgeDifference(rebuilt_edge_set, previous_edge_set)
    };
    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Adaptive local-fitting topology rebuild: accepted_iteration="
        << accepted_iteration_count
        << ", trigger=" << GetAdaptiveTopologyTriggerText(decision.trigger)
        << std::scientific << std::setprecision(2)
        << ", drift=" << decision.maximum_transformed_drift
        << ", clusters="
        << previous_partition.sample_id_list_by_key.size() << "/"
        << rebuilt_partition.sample_id_list_by_key.size()
        << ", boundary_samples="
        << previous_partition.boundary_sample_dependency_list.size() << "/"
        << rebuilt_partition.boundary_sample_dependency_list.size()
        << ", edges_added/removed="
        << added_edge_count << "/" << removed_edge_count
        << ", partition_changed="
        << (partition_changed ? "yes" : "no")
        << ", objective_domain_reset="
        << (partition_changed ? "yes" : "no") << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogFinalDependencyPolish(
    bool quiet_mode,
    const FinalDependencyPolishResult & polish_result,
    FinalPolishCertificationPolicy certification_policy,
    FinalPolishResidualSafetyStatus safety_status,
    bool applied,
    const ConvergenceCertificate * base_certificate,
    const ConvergenceCertificate * candidate_certificate)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    const auto & diagnostic{ polish_result.diagnostic };
    std::ostringstream message;
    message << std::scientific << std::setprecision(2)
        << "Final dependency polish: components/attempted/accepted/fallback="
        << diagnostic.component_count << "/"
        << diagnostic.attempted_component_count << "/"
        << diagnostic.accepted_component_count << "/"
        << diagnostic.component_count - diagnostic.accepted_component_count
        << ", atoms/parameters/rounds="
        << diagnostic.atom_count << "/"
        << diagnostic.parameter_count << "/"
        << diagnostic.round_count
        << ", objective before/after=";
    if (diagnostic.objective_before.has_value())
    {
        message << *diagnostic.objective_before;
    }
    else
    {
        message << "-";
    }
    message << "/";
    if (diagnostic.objective_after.has_value())
    {
        message << *diagnostic.objective_after;
    }
    else
    {
        message << "-";
    }
    message
        << ", accepted=" << (polish_result.accepted ? "yes" : "no")
        << ", residual-safety-policy="
        << GetFinalPolishCertificationPolicyText(certification_policy)
        << ", residual-safety="
        << GetFinalPolishResidualSafetyStatusText(safety_status)
        << ", applied=" << (applied ? "yes" : "no");
    const auto append_certificate = [&message](
        std::string_view prefix,
        const ConvergenceCertificate & certificate)
    {
        message << ", " << prefix << "-solver-qualified="
            << (certificate.solver_qualified ? "yes" : "no")
            << ", " << prefix << "-operator-complete="
            << (certificate.operator_complete ? "yes" : "no")
            << ", " << prefix << "-residual-p99=";
        AppendAuditValues(
            message,
            certificate.operator_nominal_residual.percentile_list);
        message << ", " << prefix << "-residual-max=";
        AppendAuditValues(
            message,
            certificate.operator_nominal_residual.maximum_list);
    };
    if (base_certificate != nullptr)
    {
        append_certificate("base", *base_certificate);
    }
    if (candidate_certificate != nullptr)
    {
        append_certificate("candidate", *candidate_certificate);
    }
    message
        << ", elapsed_ms=" << std::fixed << std::setprecision(3)
        << diagnostic.elapsed_milliseconds << ".";
    Logger::Log(LogLevel::Info, message.str());

    if (Logger::GetLogLevel() < LogLevel::Debug) return;
    for (std::size_t position = 0;
        position < diagnostic.component_list.size(); position++)
    {
        const auto & component{ diagnostic.component_list.at(position) };
        std::ostringstream component_message;
        component_message << std::scientific << std::setprecision(2)
            << "Final dependency polish component " << position + 1
            << ": clusters/atoms/parameters/rounds="
            << component.key_list.size() << "/"
            << component.atom_count << "/"
            << component.parameter_count << "/"
            << component.round_count
            << ", suspicious/symbolic="
            << component.suspicious_candidate_atom_count << "/"
            << component.symbolic_analysis_count
            << ", objective before/after=";
        if (component.objective_before.has_value())
        {
            component_message << *component.objective_before;
        }
        else
        {
            component_message << "-";
        }
        component_message << "/";
        if (component.objective_after.has_value())
        {
            component_message << *component.objective_after;
        }
        else
        {
            component_message << "-";
        }
        component_message
            << ", accepted/fallback="
            << (component.accepted ? "yes" : "no") << "/"
            << (!component.accepted ? "yes" : "no")
            << ", elapsed_ms=" << std::fixed << std::setprecision(3)
            << component.elapsed_milliseconds << ".";
        Logger::Log(LogLevel::Debug, component_message.str());
    }
}

void LogSecondStageAuditTerminal(
    bool quiet_mode,
    const SecondStageContext & context,
    SecondStageStopReason reason,
    std::size_t attempt_number,
    std::size_t accepted_iteration_count,
    const FitState & finalized_state,
    const ObjectiveDomain & comparison_objective_domain)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, finalized_state)
    };
    const auto objective{
        EvaluateAuditObjective(
            comparison_objective_domain,
            SnapshotResidualEvaluator{ context, model_snapshot })
    };
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "Second-stage audit terminal: schema=2"
        << ", reason=" << GetSecondStageStopReasonText(reason)
        << ", try=" << attempt_number
        << ", acc=" << accepted_iteration_count
        << ", fixed-domain=" << comparison_objective_domain.active_atom_count
        << "/" << comparison_objective_domain.fit_sample_count
        << "/" << comparison_objective_domain.tail_sample_count
        << ", objective=";
    AppendObjectiveBreakdown(message, objective, "-/-/-/-");
    Logger::Log(LogLevel::Debug, message.str());
    LogAuditAtomState(
        "Second-stage audit terminal atom:",
        context,
        finalized_state);
}

void LogQuarantineFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const FitState & finalized_state)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after "
        << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendQuarantineSummary(
        warning_message,
        entered_target_count,
        released_target_count,
        failed_probation_count,
        unresolved_target_count);
    AppendOffsetSummary(warning_message, finalized_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogConverged(
    bool quiet_mode,
    const IterationResult & iteration_result,
    const FitState & finalized_state)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    const auto & transformed_change_percentile{
        iteration_result.transformed_change_percentile
    };
    std::ostringstream message;
    message
        << "Converged after " << iteration_result.accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_percentile.at(
            GaussianModel3D::LogPeakHeightCoordinateIndex())
        << ", percentile log-width change = "
        << transformed_change_percentile.at(
            GaussianModel3D::LogWidthCoordinateIndex())
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_percentile.at(
            GaussianModel3D::OffsetToPeakRatioCoordinateIndex());
    AppendOffsetSummary(message, finalized_state);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogMaximumIterations(
    bool quiet_mode,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const BestAuditState & best_audit_state,
    const FitState & latest_state)
{
    if (quiet_mode) return;
    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendQuarantineSummary(
        warning_message,
        entered_target_count,
        released_target_count,
        failed_probation_count,
        unresolved_target_count);
    const auto * audit_state{
        best_audit_state.has_value() ? &*best_audit_state : nullptr
    };
    if (audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendAuditSummary(warning_message, *audit_state);
    }
    else
    {
        warning_message << "; applying latest validated state";
    }
    AppendOffsetSummary(
        warning_message,
        audit_state != nullptr ? audit_state->state : latest_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const BestAuditState & best_audit_state,
    const PolishProvenance & latest_polish_provenance,
    SecondStageStopReason stop_reason,
    bool final_uses_best_audit)
{
    if (quiet_mode) return;
    const auto final_uses_polish{
        final_uses_best_audit && best_audit_state.has_value() ?
            best_audit_state->uses_polish :
            std::ranges::any_of(
                latest_polish_provenance,
                [](char value) { return value != 0; })
    };
    Logger::FinishProgressLine();
    std::ostringstream message;
    message << " Second-Stage Local Fitting Summary : \n"
        << " - accepted_iterations = " << accepted_iteration_count << "\n"
        << " - best_iteration = ";
    if (!best_audit_state.has_value())
    {
        message << "unavailable\n";
    }
    else if (best_audit_state->source_iteration != 0)
    {
        message << best_audit_state->source_iteration << "\n";
    }
    else
    {
        message << "initial\n";
    }
    message << " - stop_reason = "
        << GetSecondStageStopReasonText(stop_reason) << "\n"
        << " - best_audit_objective = ";
    if (best_audit_state.has_value())
    {
        message << std::scientific << std::setprecision(2)
            << best_audit_state->objective.GetTotalObjective() << "\n";
    }
    else
    {
        message << "unavailable\n";
    }
    message << " - final_uses_polish = "
        << (final_uses_polish ? "yes" : "no") << "\n";
    message << " - final_state_source = "
        << (final_uses_best_audit ? "best-audit" : "latest-validated")
        << "\n";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace rhbm_gem::core::detail
