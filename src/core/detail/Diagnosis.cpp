#include "core/detail/Diagnosis.hpp"

#include "core/detail/IterationProcess.hpp"
#include "core/detail/DependencyPolish.hpp"

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
    case SecondStageSeedSource::LocalMdpde:
        return "local-mdpde";
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
            << marker << " schema=2"
            << ", serial=" << context.atom_list.at(atom_index).atom->GetSerialID()
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

void LogSecondStageInitializationFailure(bool quiet_mode)
{
    if (quiet_mode) return;
    Logger::Log(LogLevel::Warning,
        "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
        "is available for every selected atom.");
    Logger::Log(LogLevel::Info,
        "Second-stage local fitting summary: accepted_iterations=0, "
        "best_iteration=unavailable, stop_reason=no-valid-seed, "
        "best_audit_objective=unavailable, final_uses_polish=unavailable, "
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

    std::size_t local_mdpde_count{ 0 };
    std::size_t global_median_count{ 0 };
    for (const auto & record : selection_record_list)
    {
        if (record.source == SecondStageSeedSource::LocalMdpde)
        {
            local_mdpde_count++;
        }
        else
        {
            global_median_count++;
        }
    }

    std::ostringstream summary;
    summary << "Selected second-stage initial seeds = "
        << selection_record_list.size()
        << ", sources = local-mdpde:" << local_mdpde_count
        << ", global-median:" << global_median_count << ".";
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

void LogFrozenBackground(const SecondStageContext & context, bool quiet_mode)
{
    if (quiet_mode || !IsDebugLogLevelEnabled() || !context.frozen_background) return;
    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto & atom{ context.atom_list.at(atom_index) };
        const auto affected_samples{ std::ranges::count_if(atom.unselected_distance_list_by_sample,
            [](const auto & contributors) { return !contributors.empty(); }) };
        if (affected_samples == 0) continue;
        const auto & model{ context.frozen_background->model_by_atom.at(atom_index) };
        std::ostringstream message;
        message << std::scientific << std::setprecision(17)
            << "Second-stage frozen background: target=" << atom.atom->GetSerialID()
            << ", A/W/C=" << model.GetAmplitude() << "/" << model.GetWidth() << "/" << model.GetOffset()
            << ", samples=" << affected_samples << ".";
        Logger::Log(LogLevel::Debug, message.str());
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

namespace {

std::string BestTraceKey(const ClusterKey & key)
{
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < key.size(); i++)
    {
        if (i) out << ",";
        out << key.at(i);
    }
    out << "]";
    return out.str();
}

void AppendBestTraceObjective(std::ostream & out, const std::optional<ObjectiveBreakdown> & value)
{
    if (!value) { out << "unavailable"; return; }
    out << value->fit_range_residual_objective << "/" << value->GetTailValidationPenalty()
        << "/" << value->offset_plausibility_penalty << "/" << value->GetTotalObjective();
}

void AppendBestTraceModel(std::ostream & out, const GaussianModel3D & model)
{
    out << model.GetAmplitude() << "/" << model.GetWidth() << "/" << model.GetOffset();
}

bool SameBestTraceModel(const GaussianModel3D & a, const GaussianModel3D & b)
{
    return a.GetAmplitude() == b.GetAmplitude() && a.GetWidth() == b.GetWidth() && a.GetOffset() == b.GetOffset();
}

} // namespace

void BeginBestObjectiveTrace(
    SecondStageContext & context, bool quiet_mode, const ObjectiveDomain & domain,
    std::size_t attempt, std::size_t accepted_iteration)
{
    context.best_trace.reset();
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return;
    context.best_trace = std::make_shared<BestObjectiveTraceEnvironment>();
    context.best_trace->attempt = attempt;
    context.best_trace->accepted_iteration = accepted_iteration;
    context.best_trace->domain = std::make_shared<const ObjectiveDomain>(domain);
}

void CaptureBestObjectiveSource(
    const SecondStageContext & context, const ClusterKey & key,
    SecondStageModelSnapshot snapshot, const std::vector<SampleRef> & sample_refs,
    ClusterObjectiveState & state, const std::optional<ObjectiveBreakdown> & before,
    double before_step, std::string_view source, std::string_view reason,
    std::size_t candidate_number, std::optional<double> factor)
{
    if (!context.best_trace) return;
    auto & trace{ *context.best_trace };
    auto event{ std::make_shared<BestObjectiveSource>() };
    event->key = key;
    const auto predecessor{ state.reset_source ? state.reset_source : state.best_source };
    if (predecessor) event->predecessor_id = predecessor->id;
    event->attempt = trace.attempt;
    event->accepted_iteration = trace.accepted_iteration;
    event->source = source;
    event->reason = reason;
    event->candidate_number = candidate_number;
    event->factor = factor;
    event->before = before;
    event->before_step = before_step;
    event->objective = state.best_objective;
    event->step = state.best_maximum_transformed_change;
    event->snapshot = std::move(snapshot);
    event->domain = trace.domain;
    event->sample_refs = sample_refs;
    // One worker evaluates each key. Boundary evaluation follows worker completion.
    // The lock protects the shared container; IDs use per-key execution order.
    std::lock_guard lock{ trace.mutex };
    event->sequence = ++trace.sequence_by_key[key];
    event->id = std::to_string(trace.attempt) + "/" + BestTraceKey(key) + "/" + std::to_string(event->sequence);
    state.best_source = event;
    state.reset_source.reset();
    trace.events.emplace_back(std::move(event));
}

void LogBestObjectivePublication(const SecondStageContext & context, const ClusterObjectiveStateMap & states)
{
    if (!context.best_trace) return;
    auto events{ context.best_trace->events };
    std::ranges::sort(events, [](const auto & a, const auto & b)
    {
        return a->key == b->key ? a->sequence < b->sequence : a->key < b->key;
    });
    for (const auto & event : events)
    {
        const auto iter{ states.find(event->key) };
        const bool retained{ iter != states.end() && iter->second.best_source == event };
        std::ostringstream out;
        out << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "Cluster best source: schema=1, id=" << event->id << ", key=" << BestTraceKey(event->key)
            << ", try=" << event->attempt << ", acc-before=" << event->accepted_iteration
            << ", source=" << event->source << ", candidate=" << event->candidate_number
            << ", predecessor=" << (event->predecessor_id.empty() ? "unavailable" : event->predecessor_id)
            << ", reason=" << event->reason << ", retained=" << (retained ? "yes" : "no")
            << ", factor=";
        if (event->factor) out << *event->factor; else out << "unavailable";
        out << ", before="; AppendBestTraceObjective(out, event->before);
        out << ", best="; AppendBestTraceObjective(out, event->objective);
        out << ", step-before/after=" << event->before_step << "/" << event->step;
        Logger::Log(LogLevel::Debug, out.str());
    }
    for (const auto & [key, state] : states)
    {
        std::ostringstream out;
        out << "Cluster best publication: schema=1, try=" << context.best_trace->attempt
            << ", acc-before=" << context.best_trace->accepted_iteration << ", key=" << BestTraceKey(key)
            << ", best-source=" << (state.best_source ? state.best_source->id : "unavailable");
        Logger::Log(LogLevel::Debug, out.str());
    }
}

void DiagnoseBestObjectiveComparison(
    JointCandidateObjectiveDiagnostic * record, const CandidateEvaluationOverlay & candidate,
    const ClusterKey & key, const std::vector<SampleRef> & samples,
    const ObjectiveDomain & domain, const ClusterObjectiveState & state)
{
    if (!record || !candidate.GetContext().best_trace || !state.best_objective) return;
    if (!state.best_source)
    {
        record->best_comparison_lines.emplace_back("Cluster best comparison: schema=1, status=unavailable, diagnostic-only=yes");
        return;
    }
    const auto & origin{ *state.best_source };
    record->best_source_id = origin.id;
    if (origin.key != key || !origin.domain || !domain.cluster_by_key.contains(key) ||
        origin.snapshot.node.size() != candidate.GetState().size())
    {
        record->best_comparison_lines.emplace_back("Cluster best comparison: schema=1, status=not-comparable, diagnostic-only=yes");
        return;
    }
    const auto & context{ candidate.GetContext() };
    const auto & previous{ candidate.GetBaseline().model_snapshot };
    const auto proposed{ BuildSecondStageModelSnapshot(context, candidate.GetState()) };
    std::size_t evaluations{ 0 }, residual_evaluations{ 0 };
    const auto evaluate = [&](const SecondStageModelSnapshot & snapshot, const ObjectiveDomain & eval_domain,
                              const std::vector<SampleRef> & refs)
    {
        // Independent diagnostic baseline. Immutable atom samples are read from context;
        // the supplied snapshot owns the background used by EvaluateResidualSample.
        ResidualBaseline baseline{ snapshot, {} };
        baseline.sample_list.resize(context.atom_list.size());
        for (const auto & ref : refs)
        {
            auto & list{ baseline.sample_list.at(ref.atom_index) };
            if (list.empty()) list.resize(context.atom_list.at(ref.atom_index).raw_sampling_entries.size());
            list.at(ref.sample_index) = EvaluateResidualSample(context, ref, snapshot);
            residual_evaluations++;
        }
        evaluations++;
        return EvaluateObjectiveContribution(baseline, key, refs, eval_domain);
    };
    const auto historical{ evaluate(origin.snapshot, *origin.domain, origin.sample_refs) };
    const auto domain_only{ evaluate(origin.snapshot, domain, samples) };
    auto background_snapshot{ origin.snapshot };
    background_snapshot.frozen_background = proposed.frozen_background;
    const auto background_only{ evaluate(background_snapshot, domain, samples) };
    auto previous_environment{ previous };
    auto candidate_environment{ proposed };
    for (const auto atom : key)
    {
        previous_environment.node.at(atom) = origin.snapshot.node.at(atom);
        candidate_environment.node.at(atom) = origin.snapshot.node.at(atom);
    }
    const auto current_previous{ evaluate(previous_environment, domain, samples) };
    const auto current_candidate{ evaluate(candidate_environment, domain, samples) };
    std::ostringstream out;
    out << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "Cluster best comparison: schema=1, best-source=" << origin.id << ", key=" << BestTraceKey(key)
        << ", diagnostic-only=yes, evaluations=" << evaluations << ", residual-evaluations=" << residual_evaluations;
    const auto append = [&](std::string_view label, const std::optional<ObjectiveBreakdown> & objective)
    {
        out << ", " << label << "="; AppendBestTraceObjective(out, objective);
    };
    append("stored", state.best_objective);
    append("historical", historical);
    append("domain-only", domain_only);
    append("background-after-domain", background_only);
    append("previous-environment", current_previous);
    append("candidate-environment", current_candidate);
    append("previous", record->previous);
    append("candidate", record->candidate);
    const auto delta = [&](std::string_view label, const auto & a, const auto & b)
    {
        out << ", " << label << "=";
        if (!a || !b) { out << "unavailable"; return; }
        out << b->fit_range_residual_objective - a->fit_range_residual_objective << "/"
            << b->GetTailValidationPenalty() - a->GetTailValidationPenalty() << "/"
            << b->offset_plausibility_penalty - a->offset_plausibility_penalty << "/"
            << b->GetTotalObjective() - a->GetTotalObjective();
    };
    delta("historical-minus-stored", state.best_objective, historical);
    delta("domain-delta", historical, domain_only);
    delta("background-delta", domain_only, background_only);
    delta("previous-neighbor-delta", background_only, current_previous);
    delta("candidate-neighbor-delta", background_only, current_candidate);
    const auto gate = [&](std::string_view label, const std::optional<ObjectiveBreakdown> & reference)
    {
        out << ", " << label << "-gate=";
        if (!reference || !record->candidate) { out << "unavailable"; return; }
        const auto value{ reference->GetTotalObjective() };
        const auto tolerance{ CalculateObjectiveTolerance(value, kObjectiveProgressTolerance) };
        out << (IsObjectiveDeteriorated(record->candidate->GetTotalObjective(), value, kObjectiveProgressTolerance) ? "fail" : "pass")
            << ", " << label << "-reference=" << value
            << ", " << label << "-delta=" << record->candidate->GetTotalObjective() - value
            << ", " << label << "-absolute=" << kObjectiveProgressTolerance.absolute_tolerance
            << ", " << label << "-relative=" << kObjectiveProgressTolerance.relative_tolerance
            << ", " << label << "-tolerance=" << tolerance << ", " << label << "-limit=" << value + tolerance;
    };
    gate("stored", state.best_objective);
    gate("previous-environment", current_previous);
    gate("candidate-environment", current_candidate);
    record->best_comparison_lines.emplace_back(out.str());

    std::set<SampleRef> all_samples(origin.sample_refs.begin(), origin.sample_refs.end());
    all_samples.insert(samples.begin(), samples.end());
    std::set<std::size_t> contributors(key.begin(), key.end());
    bool owner_changed{ false }, mask_changed{ false }, scale_changed{ false }, normalization_changed{ false }, background_changed{ false };
    normalization_changed = origin.domain->active_atom_count != domain.active_atom_count ||
        origin.domain->unique_sample_count != domain.unique_sample_count ||
        origin.domain->fit_sample_count != domain.fit_sample_count || origin.domain->tail_sample_count != domain.tail_sample_count;
    for (const auto & ref : all_samples)
    {
        contributors.insert(ref.atom_index);
        for (const auto & neighbor : context.atom_list.at(ref.atom_index).Neighbors(ref.sample_index))
            contributors.insert(neighbor.atom_index);
        const auto & old_owner{ origin.domain->owner_key_by_atom_index.at(ref.atom_index) };
        const auto & new_owner{ domain.owner_key_by_atom_index.at(ref.atom_index) };
        owner_changed |= old_owner != new_owner;
        mask_changed |= origin.domain->fit_sample_mask_by_atom.at(ref.atom_index).at(ref.sample_index) != domain.fit_sample_mask_by_atom.at(ref.atom_index).at(ref.sample_index) ||
            origin.domain->tail_sample_mask_by_atom.at(ref.atom_index).at(ref.sample_index) != domain.tail_sample_mask_by_atom.at(ref.atom_index).at(ref.sample_index);
        const auto old_iter{ origin.domain->cluster_by_key.find(old_owner) };
        const auto new_iter{ domain.cluster_by_key.find(new_owner) };
        if (old_iter == origin.domain->cluster_by_key.end() || new_iter == domain.cluster_by_key.end())
        {
            scale_changed |= old_iter != origin.domain->cluster_by_key.end() || new_iter != domain.cluster_by_key.end();
            normalization_changed |= scale_changed;
        }
        else
        {
            const auto & a{ old_iter->second }; const auto & b{ new_iter->second };
            scale_changed |= a.scale.has_value() != b.scale.has_value() ||
                (a.scale && b.scale && (a.scale->fit != b.scale->fit || a.scale->tail != b.scale->tail));
            normalization_changed |= a.selected_atom_count != b.selected_atom_count ||
                a.fit_sample_ref_list.size() != b.fit_sample_ref_list.size() || a.tail_sample_ref_list.size() != b.tail_sample_ref_list.size();
        }
        const auto response = [&](const auto & snapshot)
        {
            return snapshot.frozen_background ? snapshot.frozen_background->response_by_atom.at(ref.atom_index).at(ref.sample_index) : 0.0;
        };
        background_changed |= response(origin.snapshot) != response(proposed);
    }
    std::ostringstream changes;
    changes << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "Cluster best environment: schema=1, best-source=" << origin.id << ", key=" << BestTraceKey(key)
        << ", samples-changed=" << (origin.sample_refs != samples) << ", owner-changed=" << owner_changed
        << ", mask-changed=" << mask_changed << ", scale-changed=" << scale_changed
        << ", normalization-changed=" << normalization_changed << ", background-response-changed=" << background_changed
        << ", models[atom:historical/previous/candidate]=";
    for (const auto atom : contributors)
    {
        const auto & a{ origin.snapshot.node.at(atom) }; const auto & b{ previous.node.at(atom) }; const auto & c{ proposed.node.at(atom) };
        if (std::ranges::find(key, atom) == key.end() && SameBestTraceModel(a, b) && SameBestTraceModel(a, c)) continue;
        changes << " {" << atom << (std::ranges::find(key, atom) != key.end() ? ":member:" : ":contributor:");
        AppendBestTraceModel(changes, a); changes << "|"; AppendBestTraceModel(changes, b); changes << "|"; AppendBestTraceModel(changes, c); changes << "}";
    }
    record->best_comparison_lines.emplace_back(changes.str());
}

JointCandidateObjectiveDiagnostic * BeginJointCandidateDiagnostic(
    bool quiet_mode,
    std::vector<JointCandidateObjectiveDiagnostic> & records,
    std::string_view source,
    std::optional<double> factor,
    std::size_t round)
{
    if (quiet_mode || Logger::GetLogLevel() < LogLevel::Debug) return nullptr;
    return &records.emplace_back(JointCandidateObjectiveDiagnostic{
        .source = source, .round = round, .candidate_number = records.size() + 1, .factor = factor });
}

void RecordJointMemberRejection(
    JointCandidateObjectiveDiagnostic * record,
    const ClusterKey & key,
    const std::optional<ObjectiveBreakdown> & previous,
    const std::optional<ObjectiveBreakdown> & best,
    const std::optional<ObjectiveBreakdown> & candidate,
    bool best_checked)
{
    if (record == nullptr) return;
    record->member_key = key;
    record->previous = previous;
    record->best = best;
    record->candidate = candidate;
    record->best_checked = best_checked;
    if (!previous || !candidate)
        record->outcome = "member-objective-unavailable";
    else if (!std::isfinite(previous->GetTotalObjective()) ||
        !std::isfinite(candidate->GetTotalObjective()) ||
        (best_checked && best && !std::isfinite(best->GetTotalObjective())))
        record->outcome = "member-objective-nonfinite";
    else
    {
        const bool previous_failed{ IsObjectiveDeteriorated(candidate->GetTotalObjective(),
            previous->GetTotalObjective(), kObjectiveProgressTolerance) };
        const bool best_failed{ best_checked && best && IsObjectiveDeteriorated(
            candidate->GetTotalObjective(), best->GetTotalObjective(), kObjectiveProgressTolerance) };
        record->outcome = previous_failed ? (best_failed ? "previous+best" : "previous") :
            (best_failed ? "best" : "member-check-failed");
    }
}

static void LogJointCandidateDiagnostics(
    const std::vector<ClusterKey> & component,
    const std::vector<JointCandidateObjectiveDiagnostic> & records)
{
    for (std::size_t index = 0; index < records.size(); index++)
    {
        const auto & record{ records.at(index) };
        if (record.outcome == "accepted") continue;
        std::ostringstream message;
        message << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "Joint candidate objective rejection: schema=1, source=" << record.source
            << ", candidate=" << index + 1 << ", round=" << record.round << ", factor=";
        if (record.factor) message << *record.factor;
        else message << "unavailable";
        const auto append_key = [&](const ClusterKey & key)
        {
            message << "[";
            for (std::size_t i = 0; i < key.size(); i++)
            {
                if (i != 0) message << ",";
                message << key.at(i);
            }
            message << "]";
        };
        message << ", component=";
        for (const auto & key : component) append_key(key);
        message << ", member=";
        if (record.member_key.empty()) message << "none";
        else append_key(record.member_key);
        message << ", atoms=" << record.member_key.size() << ", outcome=" << record.outcome;
        if (!record.best_source_id.empty()) message << ", best-source=" << record.best_source_id;
        const auto append_objective = [&](std::string_view label, const std::optional<ObjectiveBreakdown> & value)
        {
            message << ", " << label << "=";
            if (!value) { message << "unavailable"; return; }
            message << value->fit_range_residual_objective << "/" << value->GetTailValidationPenalty() << "/" << value->offset_plausibility_penalty << "/" << value->GetTotalObjective();
        };
        append_objective("previous", record.previous);
        append_objective("best", record.best);
        append_objective("candidate-objective", record.candidate);
        const auto append_gate = [&](std::string_view label, const std::optional<ObjectiveBreakdown> & reference, bool checked)
        {
            message << ", " << label << "-gate=";
            if (!checked) { message << "not-checked"; return; }
            if (!reference || !record.candidate) { message << "unavailable"; return; }
            const auto value{ reference->GetTotalObjective() };
            const auto tolerance{ CalculateObjectiveTolerance(value, kObjectiveProgressTolerance) };
            message << "candidate<=reference+tolerance"
                << ", " << label << "-reference=" << value
                << ", " << label << "-delta=" << record.candidate->GetTotalObjective() - value
                << ", " << label << "-absolute=" << kObjectiveProgressTolerance.absolute_tolerance
                << ", " << label << "-relative=" << kObjectiveProgressTolerance.relative_tolerance
                << ", " << label << "-tolerance=" << tolerance
                << ", " << label << "-limit=" << value + tolerance;
        };
        append_gate("previous", record.previous, !record.member_key.empty());
        append_gate("best", record.best, record.best_checked);
        Logger::Log(LogLevel::Debug, message.str());
        for (const auto & line : record.best_comparison_lines) Logger::Log(LogLevel::Debug, line);
    }
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
        LogJointCandidateDiagnostics(diagnostic.key_list, diagnostic.objective_diagnostic_list);
        if (!diagnostic.joint_correction_status.has_value()) continue;
        std::ostringstream correction_message;
        correction_message << std::scientific << std::setprecision(2)
            << "Boundary-interface joint correction: direct-interface/shape-active/offset-active/parameters = "
            << diagnostic.interface_atom_count << "/"
            << diagnostic.shape_active_atom_count << "/"
            << diagnostic.offset_active_atom_count << "/"
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
                iteration_result.proposal_maximum_transformed_change)
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
    double maximum_transformed_drift,
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
        << ", trigger=drift"
        << std::scientific << std::setprecision(2)
        << ", drift=" << maximum_transformed_drift
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
        << ", partition_pending=" << (partition_changed ? "yes" : "no")
        << ", objective_domain_reset=no.";
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
        LogJointCandidateDiagnostics(component.key_list, component.objective_diagnostic_list);
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
            context, model_snapshot)
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

PerformanceCounters::PerformanceCounters(
    bool quiet_mode,
    const SecondStageContext & context,
    const ClusterSolverWorkspaceMap & solver_workspace_by_key,
    const BoundaryJointCorrectionWorkspaceMap &
        boundary_joint_correction_workspace_by_key)
    : m_quiet_mode{ quiet_mode },
      m_solver_workspace_by_key{ solver_workspace_by_key },
      m_boundary_joint_correction_workspace_by_key{ boundary_joint_correction_workspace_by_key },
      m_start_time{ std::chrono::steady_clock::now() },
      m_cached_sample_count{ CountRawSamplingEntries(context) }
{
}

PerformanceCounters::~PerformanceCounters()
{
    if (m_quiet_mode) return;

    const auto symbolic_analysis_count{
        m_retired_solver_symbolic_analysis_count + CountCurrentSolverSymbolicAnalyses()
    };
    const auto total_milliseconds{
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - m_start_time).count()
    };
    std::ostringstream message_info;
    message_info
        << " Second-Stage Local Fitting Performance :\n"
        << " - boundary_reconciliation_ms = "
        << std::fixed << std::setprecision(3) << m_boundary_reconciliation_milliseconds << "\n"
        << " - boundary_joint_correction_ms = " << m_boundary_joint_correction_milliseconds << "\n"
        << " - dependency_polish_ms = " << m_dependency_polish_milliseconds << "\n"
        << " - iteration/candidate/topology/total_ms = "
        << std::fixed << std::setprecision(3)
        << m_iteration_phase_milliseconds << "/"
        << m_candidate_phase_milliseconds << "/"
        << m_topology_rebuild_milliseconds << "/"
        << total_milliseconds << "\n";

    std::ostringstream message_debug;
    message_debug
        << " - full_state_materializations = " << m_full_state_materialization_count.load() << "\n"
        << " - gaussian_cache_hit/miss = "
        << m_gaussian_cache_hit_count.load() << "/" << m_gaussian_cache_miss_count.load() << "\n"
        << " - objective_recomputed/reused_samples = "
        << m_objective_recomputed_sample_count.load() << "/"
        << m_objective_reused_sample_count.load() <<"\n"
        << " - solver_symbolic_analyses = " << symbolic_analysis_count << "\n"
        << " - topology_rebuilds/partition_changes = "
        << m_topology_rebuild_attempt_count << "/" << m_topology_partition_change_count<< "\n"
        << " - boundary_reconciliations/backtracked/rejected = "
        << m_boundary_reconciliation_attempt_count << "/"
        << m_boundary_reconciliation_backtracked_count << "/"
        << m_boundary_reconciliation_rejected_count << "\n"
        << " - boundary_joint_correction_attempts/accepted/fallback = "
        << m_boundary_joint_correction_attempt_count << "/"
        << m_boundary_joint_correction_accepted_count << "/"
        << m_boundary_joint_correction_fallback_count << "\n"
        << " - boundary_rescues/accepted/fallback/rejected = "
        << m_boundary_rescue_attempt_count << "/"
        << m_boundary_rescue_accepted_count << "/"
        << m_boundary_rescue_fallback_count << "/"
        << m_boundary_rescue_rejected_count << "\n"
        << " - boundary_rescue_exclusions_hard/invalid/no-objective = "
        << m_boundary_rescue_hard_failure_exclusion_count << "/"
        << m_boundary_rescue_invalid_proposal_exclusion_count << "/"
        << m_boundary_rescue_objective_unavailable_exclusion_count << "\n"
        << " - dependency_polish_components/attempted/accepted/fallback = "
        << m_dependency_polish_component_count << "/"
        << m_dependency_polish_attempt_count << "/"
        << m_dependency_polish_accepted_count << "/"
        << m_dependency_polish_fallback_count << "\n"
        << " - dependency_polish_atoms/parameters/rounds = "
        << m_dependency_polish_atom_count << "/"
        << m_dependency_polish_parameter_count << "/"
        << m_dependency_polish_round_count << "\n";

    Logger::Log(LogLevel::Info, message_info.str());
    Logger::Log(LogLevel::Debug, message_debug.str());
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    std::ostringstream audit_performance;
    audit_performance << std::scientific << std::setprecision(17)
        << "Trust-model performance: schema=1"
        << ", candidate-ms=" << m_candidate_phase_milliseconds
        << ", total-ms=" << total_milliseconds;
    Logger::Log(LogLevel::Debug, audit_performance.str());
#endif
}

void PerformanceCounters::RecordFullStateMaterialization()
{
    m_full_state_materialization_count.fetch_add(1, std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheMisses()
{
    m_gaussian_cache_miss_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheHits()
{
    m_gaussian_cache_hit_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
}

void PerformanceCounters::RecordObjectiveSampleEvaluation(
    std::size_t recomputed_sample_count,
    std::size_t total_sample_count)
{
    m_objective_recomputed_sample_count.fetch_add(
        recomputed_sample_count,
        std::memory_order_relaxed);
    m_objective_reused_sample_count.fetch_add(
        total_sample_count > recomputed_sample_count ?
            total_sample_count - recomputed_sample_count : 0,
        std::memory_order_relaxed);
}

void PerformanceCounters::FinishIterationPhase(std::chrono::steady_clock::time_point start_time)
{
    m_iteration_phase_milliseconds += CalculateElapsedMilliseconds(start_time);
}

void PerformanceCounters::FinishCandidatePhase(std::chrono::steady_clock::time_point start_time)
{
    m_candidate_phase_milliseconds += CalculateElapsedMilliseconds(start_time);
}

void PerformanceCounters::RecordSolverWorkspaceReset()
{
    m_retired_solver_symbolic_analysis_count += CountCurrentSolverSymbolicAnalyses();
}

void PerformanceCounters::RecordTopologyRebuild(double elapsed_milliseconds, bool partition_changed)
{
    m_topology_rebuild_attempt_count++;
    if (partition_changed) m_topology_partition_change_count++;
    m_topology_rebuild_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryReconciliation(
    std::size_t attempt_count,
    std::size_t backtracked_count,
    std::size_t rejected_count,
    double elapsed_milliseconds)
{
    m_boundary_reconciliation_attempt_count += attempt_count;
    m_boundary_reconciliation_backtracked_count += backtracked_count;
    m_boundary_reconciliation_rejected_count += rejected_count;
    m_boundary_reconciliation_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryJointCorrection(bool accepted, double elapsed_milliseconds)
{
    m_boundary_joint_correction_attempt_count++;
    if (accepted)
    {
        m_boundary_joint_correction_accepted_count++;
    }
    else
    {
        m_boundary_joint_correction_fallback_count++;
    }
    m_boundary_joint_correction_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryRescue(bool accepted, bool used_fallback)
{
    m_boundary_rescue_attempt_count++;
    if (accepted)
    {
        m_boundary_rescue_accepted_count++;
    }
    else
    {
        m_boundary_rescue_rejected_count++;
    }
    if (used_fallback) m_boundary_rescue_fallback_count++;
}

void PerformanceCounters::RecordBoundaryRescueExclusions(
    std::size_t hard_failure_count,
    std::size_t invalid_proposal_count,
    std::size_t objective_unavailable_count)
{
    m_boundary_rescue_hard_failure_exclusion_count += hard_failure_count;
    m_boundary_rescue_invalid_proposal_exclusion_count += invalid_proposal_count;
    m_boundary_rescue_objective_unavailable_exclusion_count += objective_unavailable_count;
}

void PerformanceCounters::RecordDependencyPolish(
    std::size_t component_count,
    std::size_t attempt_count,
    std::size_t accepted_count,
    std::size_t fallback_count,
    std::size_t atom_count,
    std::size_t parameter_count,
    std::size_t round_count,
    double elapsed_milliseconds)
{
    m_dependency_polish_component_count += component_count;
    m_dependency_polish_attempt_count += attempt_count;
    m_dependency_polish_accepted_count += accepted_count;
    m_dependency_polish_fallback_count += fallback_count;
    m_dependency_polish_atom_count += atom_count;
    m_dependency_polish_parameter_count += parameter_count;
    m_dependency_polish_round_count += round_count;
    m_dependency_polish_milliseconds += elapsed_milliseconds;
}

double PerformanceCounters::CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start_time).count();
}

std::size_t PerformanceCounters::CountRawSamplingEntries(const SecondStageContext & context)
{
    std::size_t count{ 0 };
    for (const auto & atom_context : context.atom_list)
    {
        count += atom_context.raw_sampling_entries.size();
    }
    return count;
}

std::size_t PerformanceCounters::CountCurrentSolverSymbolicAnalyses() const
{
    std::size_t count{ 0 };
    for (const auto & [key, workspace] : m_solver_workspace_by_key)
    {
        static_cast<void>(key);
        count += workspace.joint_offset.GetSymbolicAnalysisCount();
        count += workspace.joint_polish.GetSymbolicAnalysisCount();
    }
    for (const auto & [key, workspace] :
        m_boundary_joint_correction_workspace_by_key)
    {
        static_cast<void>(key);
        count += workspace.GetSymbolicAnalysisCount();
    }
    return count;
}

static void AppendObjectiveScaleSummary(std::ostringstream & message, const std::vector<double> & scale_list)
{
    if (scale_list.empty())
    {
        message << "unavailable";
        return;
    }
    message
        << array_helper::ComputePercentile(scale_list, 0.5) << "/"
        << array_helper::ComputePercentile(scale_list, 0.99) << "/"
        << std::ranges::max(scale_list);
}

void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset)
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

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode)
{
    if (quiet_mode) return;

    const auto & summary{ topology.summary };
    if (!summary.uses_weighted_graph)
    {
        Logger::Log(LogLevel::Warning,
            "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
        if (Logger::GetLogLevel() >= LogLevel::Debug)
        {
            Logger::Log(LogLevel::Debug,
                "Local-fitting weighted threshold sensitivity is unavailable in binary fallback mode.");
        }
    }
    std::ostringstream message;
    message << "Local-fitting coupling graph mode = "
        << (summary.uses_weighted_graph ? "weighted" : "binary-fallback")
        << std::scientific << std::setprecision(2)
        << ", minimum weight = " << summary.configured_minimum_weight
        << ", candidate/retained/cut edges = "
        << summary.candidate_edge_count << "/"
        << summary.retained_edge_count << "/"
        << summary.cut_edge_count
        << ", weight p50/p95/max = "
        << summary.weight_median << "/"
        << summary.weight_percentile_95 << "/"
        << summary.weight_maximum
        << ", initial components/max atoms/ratio = "
        << summary.component_count << "/"
        << summary.maximum_component_size << "/"
        << std::fixed << std::setprecision(2)
        << summary.maximum_component_ratio << ".";
    Logger::Log(LogLevel::Info, message.str());

    const auto & atom_cutoff_summary{ topology.atom_cutoff_summary };
    std::ostringstream atom_cutoff_message;
    atom_cutoff_message
        << "Local-fitting atom cutoff: atoms="
        << topology.adjacency_list.size()
        << ", limit=" << atom_cutoff_summary.maximum_atom_count_limit
        << ", clusters=" << summary.component_count
        << ", max-atoms=" << summary.maximum_component_size
        << ", cutoff-edges=" << atom_cutoff_summary.cut_edge_count << ".";
    Logger::Log(LogLevel::Info, atom_cutoff_message.str());

    for (const auto & sensitivity : summary.threshold_sensitivity_list)
    {
        std::ostringstream sensitivity_message;
        sensitivity_message
            << std::scientific << std::setprecision(2)
            << "Coupling sensitivity: threshold=" << sensitivity.minimum_weight
            << ", retained/cut="
            << sensitivity.retained_edge_count << "/"
            << sensitivity.cut_edge_count
            << ", components/max-atoms/ratio="
            << sensitivity.component_count << "/"
            << sensitivity.maximum_component_size << "/"
            << std::fixed << std::setprecision(2)
            << sensitivity.maximum_component_ratio << ".";
        Logger::Log(LogLevel::Info, sensitivity_message.str());
    }
}

} // namespace rhbm_gem::core::detail
