#include "core/detail/ClusterHistoryObserver.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <limits>
#include <ranges>
#include <set>
#include <sstream>

namespace rhbm_gem::core::detail {

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
    ClusterHistoryCounters & counters)
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
    counters.RecordObjectiveSampleEvaluation(
        CountObjectiveSamples(samples, domain), domain.unique_sample_count);
    return EvaluateObjectiveContribution(reference, key, samples, domain);
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
static void UpdateClusterHistory(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & key, const std::vector<SampleRef> & samples,
    const ObjectiveDomain & domain, std::string_view source,
    ClusterHistoryCounters & counters,
    ClusterObjectiveState & objective_state,
    const ObjectiveAttemptDiagnostic & diagnostic, ClusterHistoryDiagnostic & history)
{
    history.stored_best_objective = objective_state.best_objective;
    history.best_objective.reset();
    const bool history_complete{ objective_state.best_parameters.atom_index_list == key &&
        objective_state.best_parameters.mdpde_list.size() == key.size() };
    if (!objective_state.best_objective || history_complete)
        history.best_objective = EvaluateBestObjectiveReference(candidate_overlay, key,
            samples, domain, objective_state, counters);
    history.best_reference_unavailable = objective_state.best_objective.has_value() &&
        !history.best_objective.has_value();
    if (history.best_reference_unavailable) return;
    const auto candidate_objective_value{ diagnostic.candidate_objective->GetTotalObjective() };
    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            candidate_overlay.GetState(),
            candidate_overlay.GetBaseline().model_snapshot.node,
            key)
    };
    const auto maximum_transformed_change{ std::ranges::max(transformed_change_summary.maximum_list) };
    auto is_better_than_best{ !history.best_objective.has_value() };
    if (history.best_objective.has_value())
    {
        const auto best_objective_value{ history.best_objective->GetTotalObjective() };
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
        const auto before_step{ objective_state.best_maximum_transformed_change };
        objective_state.best_objective = diagnostic.candidate_objective;
        objective_state.best_parameters = CaptureClusterParameters(candidate_overlay.GetState(), key);
        objective_state.best_maximum_transformed_change = maximum_transformed_change;
        if (candidate_overlay.GetContext().best_trace)
            CaptureBestObjectiveSource(candidate_overlay.GetContext(), key,
                BuildSecondStageModelSnapshot(candidate_overlay.GetContext(), candidate_overlay.GetState()),
                samples, objective_state, history.best_objective, before_step,
                source, !history.best_objective ? "first-best" :
                    (IsBetterAuditObjective(candidate_objective_value, history.best_objective->GetTotalObjective(),
                        kObjectiveStrictTolerance) ? "strict-improvement" : "step-tie-break"),
                diagnostic.trial_count, diagnostic.accepted_factor);
    }
}


void ClusterHistoryObserver::Disable() noexcept
{
    if (m_disabled.exchange(true)) return;
    try { Logger::Log(LogLevel::Debug, "Cluster history observer: status=unavailable, disabled=yes"); }
    catch (...) {}
}

void BeginClusterHistoryObserver(SecondStageContext & context, bool quiet) noexcept
{
    if (quiet || Logger::GetLogLevel() < LogLevel::Debug) return;
    try { context.cluster_history = std::make_shared<ClusterHistoryObserver>(); }
    catch (...) {}
}

void ClusterHistoryObserver::BeginAttempt(SecondStageContext & context,
    const ObjectiveByKey & previous, const FitState & state,
    const CouplingGraphPartition & partition, const ObjectiveDomain & domain,
    std::size_t attempt, std::size_t accepted_iteration) noexcept
{
    if (m_disabled) return;
    try
    {
        ReconcileClusterObjectiveState(previous, state, m_staged);
        BeginBestObjectiveTrace(context, false, domain, attempt, accepted_iteration);
        for (auto & [key, history] : m_staged)
            if (history.best_objective && !history.best_source)
                CaptureBestObjectiveSource(context, key, BuildSecondStageModelSnapshot(context, state),
                    partition.sample_id_list_by_key.at(key), history,
                    history.reset_source ? history.reset_source->objective : std::nullopt,
                    history.reset_source ? history.reset_source->step : 0.0,
                    "iteration-baseline", history.best_reset_reason);
        m_baseline = m_staged;
        m_boundary.clear();
        m_next_observation = 0;
    }
    catch (...) { Disable(); }
}

void ClusterHistoryObserver::ResetPartition(const SecondStageContext & context,
    const CouplingGraphPartition & partition, const ObjectiveDomain & domain,
    const FitState & state) noexcept
{
    if (m_disabled) return;
    try
    {
        auto prior{ std::move(m_staged) };
        m_staged.clear();
        ReconcileClusterObjectiveState(BuildObjectiveByKey(partition, domain, context,
            BuildSecondStageModelSnapshot(context, state)), state, m_staged);
        for (auto & [key, history] : m_staged)
        {
            history.best_reset_reason = "partition-reset";
            const auto it{ prior.find(key) };
            if (it != prior.end()) history.reset_source = it->second.best_source;
        }
    }
    catch (...) { Disable(); }
}

void ClusterHistoryObserver::ResetBackground(const SecondStageContext & context,
    const std::shared_ptr<const FrozenBackground> & previous_background,
    const CouplingGraphPartition & partition, const ObjectiveDomain & domain,
    const FitState & state) noexcept
{
    if (m_disabled) return;
    try
    {
        const auto previous{ BuildObjectiveByKey(partition, domain, context,
            BuildSecondStageModelSnapshot(context, state)) };
        for (const auto & [key, samples] : partition.sample_id_list_by_key)
        {
            const bool changed{ std::ranges::any_of(samples, [&](const SampleRef & sample)
            {
                return !previous_background || previous_background->response_by_atom.at(sample.atom_index) !=
                    context.frozen_background->response_by_atom.at(sample.atom_index);
            }) };
            if (!changed) continue;
            const auto prior{ m_staged.at(key).best_source };
            m_staged.at(key) = ClusterObjectiveState{
                .best_objective = previous.at(key),
                .best_parameters = FitStatePatch::FromState(state, key),
                .best_reset_reason = "background-reset", .reset_source = prior };
        }
    }
    catch (...) { Disable(); }
}

void ClusterHistoryObserver::BeginSearch(const ClusterKey & key) noexcept
{
    if (m_disabled) return;
    try { m_staged.at(key) = m_baseline.at(key); }
    catch (...) { Disable(); }
}

std::shared_ptr<const ClusterHistoryDiagnostic> ClusterHistoryObserver::Local(
    const CandidateEvaluationOverlay & candidate, const ClusterKey & key,
    const std::vector<SampleRef> & samples, const ObjectiveDomain & domain,
    std::string_view source, bool accepted, const ObjectiveAttemptDiagnostic & diagnostic) noexcept
{
    if (m_disabled) return {};
    try
    {
        auto payload{ std::make_shared<ClusterHistoryDiagnostic>() };
        auto & state{ m_staged.at(key) };
        payload->stored_best_objective = state.best_objective;
        if (accepted) UpdateClusterHistory(candidate, key, samples, domain, source, m_counters,
            state, diagnostic, *payload);
        return payload;
    }
    catch (...) { Disable(); return {}; }
}

void ClusterHistoryObserver::BeginBoundary(JointCandidateObjectiveDiagnostic * record) noexcept
{
    if (m_disabled || !record) return;
    try
    {
        const auto id{ ++m_next_observation };
        m_boundary.try_emplace(id);
        record->history_observation = id;
    }
    catch (...) { Disable(); }
}

void ClusterHistoryObserver::BoundaryMember(const CandidateEvaluationOverlay & candidate,
    const ClusterKey & key, const std::vector<SampleRef> & samples, const ObjectiveDomain & domain,
    bool accepted, const ObjectiveAttemptDiagnostic & diagnostic, JointCandidateObjectiveDiagnostic * record) noexcept
{
    if (m_disabled || !record || !record->history_observation) return;
    try
    {
        // Every component trial reads the iteration baseline, not local staged history.
        auto state{ m_baseline.at(key) };
        if (accepted)
        {
            ClusterHistoryDiagnostic payload;
            UpdateClusterHistory(candidate, key, samples, domain, record->source, m_counters,
                state, diagnostic, payload);
            m_boundary.at(record->history_observation).emplace(key, std::move(state));
        }
        else
        {
            record->stored_best = state.best_objective;
            DiagnoseBestObjectiveComparison(record, candidate, key, samples, domain, state);
        }
    }
    catch (...) { Disable(); }
}

void ClusterHistoryObserver::AcceptBoundary(std::size_t observation) noexcept
{
    if (m_disabled || !observation) return;
    try
    {
        for (const auto & [key, history] : m_boundary.at(observation)) m_staged.at(key) = history;
    }
    catch (...) { Disable(); }
}

void ClusterHistoryObserver::Reject(const ClusterKey & key) noexcept
{
    BeginSearch(key);
}

void ClusterHistoryObserver::Publish(const SecondStageContext & context) noexcept
{
    if (m_disabled) return;
    try { LogBestObjectivePublication(context, m_staged); }
    catch (...) { Disable(); }
}

std::optional<ClusterObjectiveState> ClusterHistoryObserver::BaselineSnapshot(const ClusterKey & key) noexcept
{
    if (m_disabled) return std::nullopt;
    try { return m_baseline.at(key); }
    catch (...) { Disable(); return std::nullopt; }
}

std::optional<ClusterObjectiveState> ClusterHistoryObserver::Snapshot(const ClusterKey & key) noexcept
{
    if (m_disabled) return std::nullopt;
    try { return m_staged.at(key); }
    catch (...) { Disable(); return std::nullopt; }
}

} // namespace rhbm_gem::core::detail
