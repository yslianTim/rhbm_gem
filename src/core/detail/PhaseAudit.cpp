#include "core/detail/PhaseAudit.hpp"

#include "core/detail/IterationProcess.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ranges>
#include <sstream>
#include <tuple>

namespace rhbm_gem::core::detail {
namespace {

std::string Quote(std::string_view value)
{
    std::ostringstream out;
    out << '"';
    for (const char character : value)
    {
        const auto c{ static_cast<unsigned char>(character) };
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else out << c;
    }
    out << '"';
    return out.str();
}
std::string Number(double value)
{
    if (!std::isfinite(value)) return "null";
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return out.str();
}
template<class Values> std::string Array(const Values & values)
{
    std::string result{ "[" };
    for (const auto value : values)
    {
        if (result.size() != 1) result += ',';
        result += Number(static_cast<double>(value));
    }
    return result + ']';
}
bool Equal(const FittedGaussianSnapshot & a, const FittedGaussianSnapshot & b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); i++)
        if (a[i].GetAmplitude() != b[i].GetAmplitude() || a[i].GetWidth() != b[i].GetWidth() ||
            a[i].GetOffset() != b[i].GetOffset()) return false;
    return true;
}
std::string Objective(const std::optional<ObjectiveBreakdown> & objective)
{
    if (!objective) return "null";
    return "{\"fit\":" + Number(objective->fit_range_residual_objective) +
        ",\"tail_weighted\":" + Number(objective->GetTailValidationPenalty()) +
        ",\"offset\":" + Number(objective->offset_plausibility_penalty) +
        ",\"total\":" + Number(objective->GetTotalObjective()) + "}";
}
std::string Difference(const std::optional<ObjectiveBreakdown> & a,
    const std::optional<ObjectiveBreakdown> & b)
{
    return a && b ? Number(a->GetTotalObjective() - b->GetTotalObjective()) : "null";
}
struct OperatorAudit
{
    std::string status{ "unavailable" }, reason;
    bool complete{ false }, qualified{ false };
    TransformedChangeSummary summary;
    std::vector<TransformedChange> changes;
    std::string Json() const
    {
        std::string top{ "[" };
        for (std::size_t coordinate = 0; coordinate < 3; coordinate++)
        {
            if (coordinate) top += ',';
            top += '[';
            std::vector<std::size_t> order(changes.size());
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(), [&](auto a, auto b)
                { return changes[a][coordinate] > changes[b][coordinate]; });
            for (std::size_t i = 0; i < std::min<std::size_t>(5, order.size()); i++)
            {
                if (i) top += ',';
                top += "{\"atom_index\":" + std::to_string(order[i]) + ",\"residual\":" +
                    Number(changes[order[i]][coordinate]) + '}';
            }
            top += ']';
        }
        return "{\"status\":" + Quote(status) + ",\"reason\":" + Quote(reason) +
            ",\"complete\":" + (complete ? "true" : "false") +
            ",\"solver_qualified\":" + (qualified ? "true" : "false") +
            ",\"p99\":" + (complete ? Array(summary.percentile_list) : "null") +
            ",\"maximum\":" + (complete ? Array(summary.maximum_list) : "null") +
            ",\"coordinates\":[\"log-peak-height\",\"log-width\",\"offset-to-peak-ratio\"],\"population\":" + Array(summary.population_size_list) + ",\"top_atoms\":" + top + "]}";
    }
};
OperatorAudit Summarize(const IterationProposalResult & proposal, const FitState & state,
    const std::vector<ClusterKey> & keys)
{
    OperatorAudit result;
    ClusterKey population(state.size());
    std::iota(population.begin(), population.end(), 0);
    const auto & evidence{ proposal.fixed_point_operator };
    result.complete = true;
    for (const auto atom : population)
    {
        auto change{ CalculateTransformedChange(evidence.state.at(atom), state.at(atom).mdpde.GetModel()) };
        if (!evidence.shape_available_atom_mask.at(atom))
        {
            change[0] = change[1] = std::numeric_limits<double>::infinity();
            result.complete = false;
        }
        if (!evidence.offset_available_atom_mask.at(atom))
        {
            change[2] = std::numeric_limits<double>::infinity();
            result.complete = false;
        }
        result.changes.push_back(change);
    }
    result.summary = SummarizeActiveDofChanges(result.changes, { population, population });
    const SuspiciousBlockActivity nominal_activity{ std::vector<char>(state.size(), 0),
        std::vector<char>(state.size(), 0), std::vector<char>(state.size(), 0) };
    result.qualified = AreActiveCoordinatesSolverQualified(population, keys, nominal_activity,
        proposal.local_refit_status_by_atom, proposal.health_by_key);
    // A protected offset can trigger a separate unrestricted shape solve whose status is
    // not exposed by BuildIterationProposal. Do not certify it using the protected solve.
    const bool qualification_unavailable{ std::ranges::any_of(proposal.block_activity.offset_fixed_atom_mask,
        [](auto value) { return value != 0; }) };
    if (qualification_unavailable) result.qualified = false;
    result.status = result.complete && result.qualified ? "available" : "unavailable";
    result.reason = !result.complete ? "incomplete-operator" : (qualification_unavailable ? "unrestricted-solver-qualification-unavailable" : (!result.qualified ? "solver-unqualified" : ""));
    return result;
}
bool Reproduced(const OperatorAudit & a, const OperatorAudit & b)
{
    if (a.status != "available" || b.status != "available" || a.changes.size() != b.changes.size()) return false;
    for (std::size_t i = 0; i < a.changes.size(); i++)
        for (std::size_t c = 0; c < 3; c++)
            if (!std::isfinite(a.changes[i][c]) ||
                std::abs(a.changes[i][c] - b.changes[i][c]) > 1.0e-10 + 1.0e-8 * std::abs(b.changes[i][c])) return false;
    return true;
}
} // namespace

PhaseAudit::PhaseAudit(const SecondStageContext & context, const ObjectiveDomain & domain,
    const FitState & baseline, std::vector<ClusterKey> keys, std::size_t attempt, std::size_t domain_id)
    : m_context(context), m_domain(domain), m_baseline(baseline), m_keys(std::move(keys)),
      m_attempt(attempt), m_domain_id(domain_id)
{
    for (const auto & key : m_keys) m_worker_events.try_emplace(key);
    m_context.phase_audit.reset();
    m_context.best_trace.reset();
    Add("baseline", {}, BuildSecondStageModelSnapshot(m_context, baseline).node, {}, 1.0, "baseline", "", false, true);
}
std::string PhaseAudit::Add(std::string stage, ClusterKey key, FittedGaussianSnapshot state,
    FittedGaussianSnapshot parent, double factor, std::string disposition, std::string reason,
    bool probe, bool recertify)
{
    const auto worker{ m_worker_events.find(key) };
    const bool local{ (stage == "local-search" || stage == "local-polish") && worker != m_worker_events.end() };
    // Each cluster worker owns its event buffer until the selection synchronization point.
    std::unique_lock lock(m_mutex, std::defer_lock);
    if (!local) lock.lock();
    auto & events{ local ? worker->second : m_events };
    const auto stem{ std::to_string(m_attempt) + "/" + stage + "/" + Array(key) + "/" };
    std::size_t sequence{ 1 };
    for (const auto & event : events) if (event.stage == stage && event.key == key) sequence++;
    const auto id{ stem + std::to_string(sequence) };
    std::string parent_id;
    if (!m_events.empty())
    {
        parent_id = m_events.front().id;
        if (!parent.empty() && !Equal(m_events.front().state, parent))
        {
            const auto found{ std::ranges::find_if(events | std::views::reverse, [&](const auto & event)
                { return (event.key == key || event.key.empty()) && Equal(event.state, parent); }) };
            if (found != events.rend()) parent_id = found->id;
            else
            {
                parent_id = id + "/parent";
                events.push_back({ parent_id, m_events.front().id, "parent", "reference", "",
                    key, std::move(parent), 1.0, false, recertify });
            }
        }
    }
    events.push_back({ id, parent_id, std::move(stage), std::move(disposition), std::move(reason),
        std::move(key), std::move(state), factor, probe, recertify });
    return id;
}
void PhaseAudit::Capture(std::string_view stage, const ClusterKey & key, const FitStateView & state,
    const FitStateView * parent, double factor, std::string_view disposition, std::string_view reason,
    bool probe, bool recertify) noexcept
{
    try
    {
        Add(std::string(stage), key, BuildSecondStageModelSnapshot(m_context, state).node,
            parent ? BuildSecondStageModelSnapshot(m_context, *parent).node : FittedGaussianSnapshot{},
            factor, std::string(disposition), std::string(reason), probe, recertify);
    }
    catch (...) { m_capture_failures++; }
}
void PhaseAudit::CaptureState(std::string_view stage, const FitState & state, bool probe) noexcept
{
    try
    {
        FittedGaussianSnapshot parent;
        std::string parent_id;
        const auto parent_stage{ stage == "assembly-after-polish" ? "assembly-before-polish" :
            stage == "boundary-final" ? "assembly-after-polish" : stage == "final-selection" ? "boundary-final" : "baseline" };
        for (const auto & event : m_events) if (event.stage == parent_stage)
        { parent = event.state; parent_id = event.id; }
        Add(std::string(stage), {}, BuildSecondStageModelSnapshot(m_context, state).node,
            std::move(parent), 1.0, "observed", "", probe, true);
        if (!parent_id.empty()) m_events.back().parent_id = parent_id;
    }
    catch (...) { m_capture_failures++; }
}
void PhaseAudit::CaptureOperator(const FixedPointOperatorEvidence & evidence) noexcept
{
    try
    {
        const bool complete{ std::ranges::all_of(evidence.shape_available_atom_mask, [](auto v) { return v != 0; }) &&
            std::ranges::all_of(evidence.offset_available_atom_mask, [](auto v) { return v != 0; }) };
        Add("unrestricted-proposal", {}, complete ? evidence.state : FittedGaussianSnapshot{}, {},
            1.0, complete ? "observed" : "unavailable", complete ? "" : "incomplete-production-operator", true, true);
    }
    catch (...) { m_capture_failures++; }
}
void PhaseAudit::Missing(std::string_view stage, const ClusterKey & key, std::string_view reason) noexcept
{
    try { Add(std::string(stage), key, {}, {}, 1.0, "unavailable", std::string(reason), false, false); }
    catch (...) { m_capture_failures++; }
}
void PhaseAudit::MergeWorkers()
{
    for (auto & [key, events] : m_worker_events)
    {
        for (auto & event : events) m_events.push_back(std::move(event));
        events.clear();
    }
}
void PhaseAudit::CaptureSearchAssembly() noexcept
{
    try
    {
        MergeWorkers();
        auto state{ m_baseline };
        // Called after cluster workers have joined; each key has one terminal search event.
        for (const auto & event : m_events)
            if (event.stage == "local-search" && event.disposition == "accepted")
                for (const auto atom : event.key)
                    state.at(atom).mdpde = GaussianModel3DWithUncertainty{ event.state.at(atom), {} };
        CaptureState("assembly-before-polish", state, true);
    }
    catch (...) { m_capture_failures++; }
}

void PhaseAudit::Finish(const FitOptions & options, const std::vector<double> & ridge,
    const SuspiciousBlockActivity & activity, const IterationProposalResult & production,
    const FitState & final_state) noexcept
{
    const auto start{ std::chrono::steady_clock::now() };
    std::size_t objective_count{ 0 }, operator_count{ 0 }, failures{ m_capture_failures.load() };
    try
    {
        MergeWorkers();
        std::map<std::string, std::optional<ObjectiveBreakdown>> objectives;
        std::map<std::string, OperatorAudit> operators;
        const auto evaluate = [&](const FittedGaussianSnapshot & state) -> std::optional<ObjectiveBreakdown>
        {
            if (state.empty()) return std::nullopt;
            objective_count++;
            try
            {
                auto result{ EvaluateAuditObjective(m_domain, m_context, { state, m_context.frozen_background }) };
                if (!result) failures++;
                return result;
            }
            catch (...) { failures++; return std::nullopt; }
        };
        for (const auto & event : m_events)
        {
            objectives[event.id] = evaluate(event.state);
            if (!event.recertify || event.state.empty()) continue;
            operator_count++;
            try
            {
                auto state{ m_baseline };
                for (std::size_t i = 0; i < state.size(); i++)
                    state[i].mdpde = GaussianModel3DWithUncertainty{ event.state.at(i), {} };
                ClusterSolverWorkspaceMap workspaces;
                for (const auto & key : m_keys) workspaces.try_emplace(key);
                auto quiet_options{ options };
                quiet_options.quiet_mode = true;
                quiet_options.thread_size = 1;
                const auto proposal{ BuildIterationProposal(m_context, m_keys, state, quiet_options,
                    ridge, activity, workspaces, "phase-audit") };
                operators[event.id] = Summarize(proposal, state, m_keys);
                if (operators[event.id].status != "available") failures++;
            }
            catch (const std::exception & error) { operators[event.id].reason = error.what(); failures++; }
            catch (...) { operators[event.id].reason = "operator-exception"; failures++; }
        }
        const auto baseline_id{ m_events.front().id };
        const auto production_audit{ Summarize(production, m_baseline, m_keys) };
        const bool reproduced{ Reproduced(operators.at(baseline_id), production_audit) };
        const auto final_snapshot{ BuildSecondStageModelSnapshot(m_context, final_state).node };
        // IDs are per stage/key, independent of worker scheduling.
        std::stable_sort(m_events.begin(), m_events.end(), [](const auto & a, const auto & b)
            { return std::tie(a.key, a.id) < std::tie(b.key, b.id); });
        for (const auto & event : m_events)
        {
            const auto objective{ objectives.at(event.id) };
            const auto parent{ objectives.find(event.parent_id) };
            const auto parent_objective{ parent != objectives.end() ? parent->second : std::nullopt };
            bool retained{ !event.state.empty() && event.disposition != "rejected" && event.disposition != "unavailable" };
            ClusterKey members{ event.key };
            if (members.empty()) { members.resize(final_snapshot.size()); std::iota(members.begin(), members.end(), 0); }
            if (retained)
                for (const auto atom : members)
                    if (!Equal({ event.state.at(atom) }, { final_snapshot.at(atom) })) { retained = false; break; }
            std::string probes{ "[" };
            if (event.probe && !event.state.empty() && parent != objectives.end())
            {
                const auto parent_event{ std::ranges::find(m_events, event.parent_id, &Event::id) };
                if (parent_event != m_events.end() && !parent_event->state.empty())
                    for (const auto alpha : { 1.0, 0.5, 0.125, 0.03125, 0.001 })
                    {
                        std::optional<ObjectiveBreakdown> sampled;
                        if (alpha == 1.0) sampled = objective;
                        else if (const auto models = BuildDampedModelList(parent_event->state, event.state, alpha)) sampled = evaluate(*models);
                        if (probes.size() != 1) probes += ',';
                        probes += "{\"alpha\":" + Number(alpha) + ",\"objective\":" + Objective(sampled) +
                            ",\"delta_parent\":" + Difference(sampled, parent_objective) + '}';
                    }
            }
            const auto op{ operators.find(event.id) };
            const std::string line{ "Second-stage phase audit: schema=1, payload={\"attempt\":" + std::to_string(m_attempt) +
                ",\"domain_id\":" + std::to_string(m_domain_id) + ",\"candidate_id\":" + Quote(event.id) +
                ",\"parent_id\":" + Quote(event.parent_id) + ",\"stage\":" + Quote(event.stage) +
                ",\"key\":" + Array(event.key) + ",\"factor\":" + Number(event.factor) +
                ",\"disposition\":" + Quote(event.disposition) + ",\"reason\":" + Quote(event.reason) +
                ",\"final_retained\":" + (retained ? "true" : "false") +
                ",\"retention_rule\":\"member-parameters-identical\",\"objective\":" + Objective(objective) +
                ",\"delta_baseline\":" + Difference(objective, objectives.at(baseline_id)) +
                ",\"delta_parent\":" + Difference(objective, parent_objective) +
                ",\"operator_reproduced\":" + (reproduced ? "true" : "false") +
                ",\"production_operator\":" + (event.stage == "baseline" ? production_audit.Json() : "null") +
                ",\"operator\":" + (op == operators.end() ? "null" : op->second.Json()) +
                ",\"direction_samples\":" + probes + "]}" };
            Logger::Log(LogLevel::Debug, line);
        }
    }
    catch (const std::exception & error)
    {
        failures++;
        Logger::Log(LogLevel::Debug, "Second-stage phase audit error: " + Quote(error.what()));
    }
    catch (...) { failures++; }
    const auto elapsed{ std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() };
    Logger::Log(LogLevel::Debug, "Second-stage phase audit counters: schema=1, payload={\"attempt\":" +
        std::to_string(m_attempt) + ",\"objective_evaluations\":" + std::to_string(objective_count) +
        ",\"objective_sample_evaluations\":" + std::to_string(objective_count * m_domain.unique_sample_count) +
        ",\"operator_evaluations\":" + std::to_string(operator_count) + ",\"failures\":" +
        std::to_string(failures) + ",\"elapsed_ms\":" + Number(elapsed) + "}");
}

std::shared_ptr<PhaseAudit> BeginPhaseAudit(const SecondStageContext & context, bool quiet,
    const ObjectiveDomain & domain, const FitState & baseline, const std::vector<ClusterKey> & keys,
    std::size_t attempt, std::size_t domain_id) noexcept
{
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE
    if (!quiet && Logger::GetLogLevel() >= LogLevel::Debug)
        try { return std::make_shared<PhaseAudit>(context, domain, baseline, keys, attempt, domain_id); }
        catch (...) { Logger::Log(LogLevel::Debug, "Second-stage phase audit error: capture initialization failed"); }
#else
    (void)context; (void)quiet; (void)domain; (void)baseline; (void)keys; (void)attempt; (void)domain_id;
#endif
    return {};
}
std::string_view PhaseAuditRejectionReason(const ObjectiveAttemptDiagnostic & diagnostic)
{
    if (diagnostic.best_reference_unavailable) return "best-reference-unavailable";
    if (!diagnostic.candidate_objective) return "objective-unavailable";
    if (diagnostic.rejected_by_previous) return diagnostic.rejected_by_best ? "previous+best" : "previous";
    if (diagnostic.rejected_by_best) return "best";
    return "strict-improvement";
}
} // namespace rhbm_gem::core::detail
