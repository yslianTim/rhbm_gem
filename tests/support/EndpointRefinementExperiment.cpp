#include "support/EndpointRefinementExperiment.hpp"
#include "support/SolverFailureCapture.hpp"
#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <stdexcept>

namespace second_stage_test {
namespace {
namespace j = boost::json;
namespace d = rhbm_gem::core::detail;
using rhbm_gem::RHBMEstimationStatus;
struct Run
{
    EndpointPolicy policy{EndpointPolicy::Legacy};
    int budget{};
    bool probe{}, compare{};
    std::filesystem::path directory;
    std::mutex mutex;
    j::array records;
    std::ofstream ledger;
    std::atomic<bool> failed{};
    std::size_t calls{}, triggered{}, accepted{}, equations{}, references{};
};
std::atomic<Run *> active{};
struct Activate
{
    Run * previous;
    explicit Activate(Run * run) : previous{active.exchange(run)} {}
    ~Activate() { active.store(previous); }
};
j::value Number(double x) { return std::isfinite(x) ? j::value(x) : j::value(nullptr); }
template<class V> j::array Vector(const V & values)
{ j::array out; for (const auto x : values) out.push_back(Number(x)); return out; }
template<class M> j::object Matrix(const M & matrix)
{
    j::array values;
    for (Eigen::Index r=0;r<matrix.rows();++r) for (Eigen::Index c=0;c<matrix.cols();++c) values.push_back(Number(matrix(r,c)));
    return {{"rows",matrix.rows()},{"cols",matrix.cols()},{"values",std::move(values)}};
}
template<class M> j::array Model(const M & model)
{ return {Number(model.GetAmplitude()),Number(model.GetWidth()),Number(model.GetOffset())}; }
j::object Estimate(const rhbm_gem::GaussianModel3DWithUncertainty & value)
{ return {{"model",Model(value.GetModel())},{"uncertainty",Model(value.GetStandardDeviationModel())}}; }
void Write(const std::filesystem::path & path, const j::value & value)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream file(path); file.exceptions(std::ios::failbit | std::ios::badbit);
    file << j::serialize(value);
}
j::object Snapshot(const d::SecondStageContext & context, const d::FitState & state,
    const std::vector<d::ClusterKey> & keys, const rhbm_gem::core::FitOptions & options, const std::vector<double> & ridge)
{
    j::array atoms, models, background_models, background_response;
    for (std::size_t i=0;i<state.size();++i)
    {
        const auto & fit{state[i]};
        j::object model{{"alpha",fit.alpha_r},{"ols",Estimate(fit.ols)},{"mdpde",Estimate(fit.mdpde)}};
        if (fit.fit_result)
        {
            const auto & r{*fit.fit_result};
            model["fit"] = j::object{{"status",static_cast<int>(r.status)},{"beta_ols",Vector(r.beta_ols)},
                {"beta_mdpde",Vector(r.beta_mdpde)},{"variance",Number(r.sigma_square)},
                {"weights",Vector(r.data_weight.diagonal())},{"covariance",Vector(r.data_covariance.diagonal())},
                {"iterations",r.diagnostics.iterations},
                {"beta_change",Number(r.diagnostics.squared_beta_change.value_or(NAN))},
                {"variance_change",Number(r.diagnostics.relative_variance_change.value_or(NAN))}};
        }
        models.push_back(std::move(model));
        const auto & atom{context.atom_list[i]};
        j::array samples, neighbors;
        for (const auto & s : atom.raw_sampling_entries) samples.push_back(j::object{
            {"position",j::value_from(s.point.position)},{"distance",Number(s.point.distance)},
            {"selected",s.point.is_selected},{"response",Number(s.response)}});
        for (const auto & n : atom.neighbor_atom_sample_list) neighbors.push_back(j::array{n.atom_index,Number(n.distance)});
        j::object a{{"index",i},{"alpha",atom.alpha_r},{"samples",std::move(samples)},
            {"neighbors",std::move(neighbors)},{"neighbor_offsets",j::value_from(atom.neighbor_atom_sample_offset_list)},
            {"unselected_distances",j::value_from(atom.unselected_distance_list_by_sample)},
            {"refit_design",atom.refit_design.ExperimentSnapshot()}};
        if (atom.atom) a["identity"] = j::object{{"serial_id",atom.atom->GetSerialID()},
            {"chain_id",atom.atom->GetChainID()},{"sequence_id",atom.atom->GetSequenceID()},
            {"component_id",atom.atom->GetComponentID()},{"atom_id",atom.atom->GetAtomID()},
            {"alternate_indicator",atom.atom->GetIndicator()},{"position",j::value_from(atom.atom->GetPosition())}};
        atoms.push_back(std::move(a));
    }
    if (context.frozen_background)
    {
        for (const auto & model : context.frozen_background->model_by_atom) background_models.push_back(Model(model));
        for (const auto & r : context.frozen_background->response_by_atom) background_response.push_back(Vector(r));
    }
    return {{"schema_version",1},{"state",std::move(models)},{"atoms",std::move(atoms)},
        {"background_present",bool(context.frozen_background)},{"background_models",std::move(background_models)},
        {"background_response",std::move(background_response)},{"cluster_keys",j::value_from(keys)},{"ridge",Vector(ridge)},
        {"options",j::object{{"threads",options.thread_size},{"quiet",options.quiet_mode},
            {"exclude_hydrogen",options.exclude_hydrogen},{"halo_depth",options.second_stage_boundary_halo_depth},
            {"polish",options.enable_second_stage_dependency_polish},{"polish_iterations",options.second_stage_dependency_polish_max_iterations}}}};
}
j::object Operator(const d::FixedPointOperatorEvidence & evidence, const d::FitState & state)
{
    const auto assessment{d::AssessNominalOperator(evidence,state)};
    const auto mean{d::QualifiedNominalResidualMeanSquare(evidence,state)};
    j::array models, shape, offsets;
    for (const auto & model : evidence.state) models.push_back(Model(model));
    for (const auto & s : evidence.shape_solves) shape.push_back(j::object{
        {"status",s.status ? j::value(static_cast<int>(*s.status)) : j::value(nullptr)},
        {"variance",Number(s.variance.value_or(NAN))},{"iterations",s.diagnostics.iterations},
        {"beta_change",Number(s.diagnostics.squared_beta_change.value_or(NAN))},
        {"variance_change",Number(s.diagnostics.relative_variance_change.value_or(NAN))}});
    for (const auto & [key,s] : evidence.offset_solves) offsets.push_back(j::object{
        {"key",j::value_from(key)},{"status",static_cast<int>(s.status)},{"offset",Vector(s.offset)},
        {"iterations",s.diagnostics.iterations},{"change",Number(s.diagnostics.normalized_change.value_or(NAN))},
        {"scale",Number(s.diagnostics.robust_scale.value_or(NAN))}});
    return {{"state",std::move(models)},{"shape_solves",std::move(shape)},{"offset_solves",std::move(offsets)},
        {"shape_available",j::value_from(evidence.shape_available_atom_mask)},
        {"offset_available",j::value_from(evidence.offset_available_atom_mask)},
        {"qualified",assessment.certificate.solver_qualified},{"complete",assessment.certificate.operator_complete},
        {"nominal_p99",Vector(assessment.diagnostics.operator_nominal_residual.percentile_list)},
        {"nominal_max",Vector(assessment.diagnostics.operator_nominal_residual.maximum_list)},
        {"recovery_residual_mean_square",mean ? Number(*mean) : j::value(nullptr)}};
}
} // namespace

const char * EndpointPolicyName(EndpointPolicy policy)
{
    switch (policy)
    {
        case EndpointPolicy::Legacy: return "legacy";
        case EndpointPolicy::FailedOnly: return "failed-only";
        case EndpointPolicy::FreshResidual: return "fresh-residual";
    }
    throw std::invalid_argument("Unknown endpoint policy.");
}
bool ShouldRefineEndpoint(EndpointPolicy policy, RHBMEstimationStatus status, bool fresh_pass)
{
    return policy != EndpointPolicy::Legacy &&
        (status != RHBMEstimationStatus::SUCCESS || (policy == EndpointPolicy::FreshResidual && !fresh_pass));
}

struct ScopedSecondStageEndpointExperiment::Impl
{
    Run run;
    std::unique_ptr<Activate> activation;
};
ScopedSecondStageEndpointExperiment::ScopedSecondStageEndpointExperiment()
{
    const char * directory{std::getenv("RHBM_TEST_ENDPOINT_DIR")};
    if (!directory || !*directory) return;
    impl = std::make_unique<Impl>(); auto & run{impl->run}; run.directory = directory;
    const std::string policy{std::getenv("RHBM_TEST_ENDPOINT_POLICY") ? std::getenv("RHBM_TEST_ENDPOINT_POLICY") : "legacy"};
    if (policy == "failed-only") run.policy = EndpointPolicy::FailedOnly;
    else if (policy == "fresh-residual") run.policy = EndpointPolicy::FreshResidual;
    else if (policy != "legacy") throw std::invalid_argument("Unknown test endpoint policy.");
    run.budget = std::getenv("RHBM_TEST_ENDPOINT_BUDGET") ? std::stoi(std::getenv("RHBM_TEST_ENDPOINT_BUDGET")) : 0;
    run.compare = std::getenv("RHBM_TEST_ENDPOINT_COMPARE") && std::string(std::getenv("RHBM_TEST_ENDPOINT_COMPARE")) == "1";
    if ((run.compare || run.policy != EndpointPolicy::Legacy) && run.budget < 11)
        throw std::invalid_argument("A frozen refinement budget is required.");
    std::filesystem::create_directories(run.directory);
    if (std::filesystem::exists(run.directory/"solves.jsonl")) throw std::runtime_error("Use a fresh endpoint experiment directory.");
    run.ledger.open(run.directory/"solves.jsonl"); run.ledger.exceptions(std::ios::failbit | std::ios::badbit);
    impl->activation = std::make_unique<Activate>(&run);
}
ScopedSecondStageEndpointExperiment::~ScopedSecondStageEndpointExperiment()
{
    if (!impl) return;
    impl->activation.reset();
    try
    {
        auto & run{impl->run}; run.ledger.close();
        Write(run.directory/"session.json",j::object{{"schema_version",1},{"policy",EndpointPolicyName(run.policy)},
            {"budget",run.budget},{"complete",!run.failed.load()},{"calls",run.calls},{"triggered",run.triggered},
            {"accepted",run.accepted},{"candidate_and_diagnostic_equations",run.equations},{"reference_updates",run.references}});
    }
    catch (...) { /* The runner requires a complete session marker. */ }
}
bool IsEndpointOperatorProbe() noexcept
{ const auto * run{active.load()}; return run && run->probe; }

rhbm_gem::RHBMBetaEstimateResult ApplyEndpointExperiment(const ShapeFixture & f)
{
    auto * run{active.load()};
    if (!run || (run->policy == EndpointPolicy::Legacy && !run->probe)) return f.expected;
    const auto fresh{EvaluateMDPDEEquations(f.dataset,f.alpha,f.expected.beta_mdpde,f.expected.sigma_square,f.options.data_weight_min)};
    const bool pass{fresh.valid && fresh.scaled.lpNorm<Eigen::Infinity>() <= 1e-8};
    const bool trigger{ShouldRefineEndpoint(run->policy,f.expected.status,pass)};
    auto result{f.expected};
    auto record{CurrentSolverCaptureMember()};
    record["policy"] = EndpointPolicyName(run->policy); record["triggered"] = trigger;
    record["original_status"] = static_cast<int>(f.expected.status); record["original_equations"] = EquationJSON(fresh);
    record["original_iterations"] = f.expected.diagnostics.iterations;
    int equations{1}, references{}; bool accepted{};
    if (trigger)
    {
        auto refined{RefineMDPDEEndpoint(f,run->budget,&fresh)};
        result = std::move(refined.result); accepted = refined.accepted;
        if (auto * root = refined.evidence.if_contains("root")) root->as_object().erase("trace");
        equations = static_cast<int>(refined.evidence.at("candidate_equation_evaluations").as_int64());
        references = static_cast<int>(refined.evidence.at("reference_updates").as_int64());
        record["refinement"] = std::move(refined.evidence);
        if (!accepted) record["rejected_input"] = j::object{{"X",Matrix(f.dataset.X)},
            {"y",Vector(f.dataset.y)},{"alpha",f.alpha},{"weight_floor",f.options.data_weight_min}};
    }
    record["effective_status"] = static_cast<int>(result.status);
    record["effective_beta"] = Vector(result.beta_mdpde); record["effective_variance"] = Number(result.sigma_square);
    std::lock_guard lock(run->mutex);
    ++run->calls; run->triggered += trigger; run->accepted += accepted;
    run->equations += static_cast<std::size_t>(equations); run->references += static_cast<std::size_t>(references);
    if (run->probe) run->records.push_back(std::move(record));
    else
    {
        try { run->ledger << j::serialize(record) << '\n'; }
        catch (...) { run->failed = true; throw; }
    }
    return result;
}

void CaptureEndpointInitialState(const d::SecondStageContext & context, const d::FitState & state,
    const rhbm_gem::core::FitOptions & options)
{
    auto * run{active.load()};
    if (run) Write(run->directory/"initial.json",Snapshot(context,state,{},options,{}));
}

void CompareEndpointOperators(const char * label, const d::SecondStageContext & context,
    const std::vector<d::ClusterKey> & keys, const d::FitState & state,
    const rhbm_gem::core::FitOptions & options, const std::vector<double> & ridge,
    const d::FixedPointOperatorEvidence & baseline)
{
    auto * run{active.load()};
    if (!run || !run->compare || run->probe || run->policy != EndpointPolicy::Legacy) return;
    try
    {
        const auto directory{run->directory/label};
        if (std::filesystem::exists(directory/"comparison.json")) return;
        const auto before{Snapshot(context,state,keys,options,ridge)};
        Write(directory/"input-before.json",before);
        const auto expected{Operator(baseline,state)};
        j::array variants;
        for (const auto policy : {EndpointPolicy::Legacy,EndpointPolicy::FailedOnly,EndpointPolicy::FreshResidual})
        {
            Run probe; probe.policy=policy; probe.budget=run->budget; probe.probe=true;
            d::FixedPointOperatorEvidence evidence;
            { Activate scope(&probe); evidence=d::EvaluateNominalOperator(context,keys,state,options,ridge); }
            auto output{Operator(evidence,state)};
            output["original_solver_qualified"] = expected.at("qualified");
            output["policy"] = EndpointPolicyName(policy);
            output["offset_solves_equal"] = output.at("offset_solves") == expected.at("offset_solves");
            if (policy == EndpointPolicy::Legacy)
            {
                auto numerical{output}; numerical.erase("policy"); numerical.erase("offset_solves_equal");
                numerical.erase("original_solver_qualified");
                output["legacy_exact"] = numerical == expected;
            }
            j::array delta;
            for (std::size_t i=0;i<state.size();++i)
                delta.push_back(Vector(evidence.state[i].ToVector()-baseline.state[i].ToVector()));
            output["delta_abc_from_legacy"] = std::move(delta);
            output["solves"] = std::move(probe.records);
            output["candidate_and_diagnostic_equations"] = probe.equations;
            output["reference_updates"] = probe.references;
            variants.push_back(std::move(output));
        }
        const auto after{Snapshot(context,state,keys,options,ridge)};
        Write(directory/"input-after.json",after);
        Write(directory/"comparison.json",j::object{{"schema_version",1},{"label",label},
            {"input_unchanged",before==after},{"baseline",expected},{"variants",std::move(variants)}});
        if (before != after) throw std::runtime_error("Operator probe changed its input.");
    }
    catch (...) { run->failed=true; throw; }
}

} // namespace second_stage_test
