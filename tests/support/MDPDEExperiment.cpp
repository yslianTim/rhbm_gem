#include "support/MDPDEExperiment.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace second_stage_test {
namespace {
namespace json = boost::json;
using Clock = std::chrono::steady_clock;
thread_local json::array * production_trace{};
json::value Number(double x) { return std::isfinite(x) ? json::value(x) : json::value(nullptr); }
json::array Vector(const Eigen::VectorXd & v)
{
    json::array a; for (double x : v) a.push_back(Number(x)); return a;
}
Eigen::MatrixXd ReadMatrix(std::istream & in)
{
    int rows{}, cols{};
    if (!(in >> rows >> cols) || rows < 0 || rows > 1000000 || cols < 0 || cols > 1000)
        throw std::runtime_error("Invalid shape matrix dimensions.");
    Eigen::MatrixXd m(rows, cols);
    for (int r = 0; r < rows; ++r) for (int c = 0; c < cols; ++c)
    {
        std::string token;
        if (!(in >> token)) throw std::runtime_error("Truncated shape fixture.");
        m(r,c) = std::stod(token);
    }
    return m;
}
Eigen::VectorXd Coordinates(const Eigen::VectorXd & beta, double v)
{
    Eigen::VectorXd u(3); u << beta(0), std::log(beta(1)), std::log(v); return u;
}
json::object Endpoint(const ShapeFixture & f, const Eigen::VectorXd & beta, double v)
{
    auto e{ EvaluateMDPDEEquations(f.dataset, f.alpha, beta, v, f.options.data_weight_min) };
    auto result{ EquationJSON(e) };
    result["beta"] = Vector(beta); result["variance"] = Number(v);
    if (beta.size() == 2 && beta(1) > 0.0)
    {
        result["amplitude"] = Number(std::exp(beta(0)) * std::pow(2.0 * std::acos(-1.0) / beta(1), 1.5));
        result["width"] = Number(1.0 / std::sqrt(beta(1)));
    }
    return result;
}
json::object FixedPoint(const ShapeFixture & f, const std::string & backend, bool own_initialization)
{
    const auto start{ Clock::now() };
    Eigen::VectorXd beta{ own_initialization
        ? MDPDETestBeta(f.dataset, Eigen::VectorXd::Ones(f.dataset.y.size()), backend)
        : f.expected.beta_ols };
    double v{ (f.dataset.y - f.dataset.X * beta).squaredNorm() / (f.dataset.y.size() - 1) };
    json::array trace;
    std::string stop{ "budget-exhausted" };
    int linear_solves{ own_initialization ? 1 : 0 }, evaluations{};
    for (int i = 0; i < 10000; ++i)
    {
        const auto e{ EvaluateMDPDEEquations(f.dataset, f.alpha, beta, v, f.options.data_weight_min) };
        ++evaluations;
        if (e.weights.size() != f.dataset.y.size() || e.rank != f.dataset.X.cols())
        { stop = e.reason; break; }
        const auto previous{ beta }; const double old_v{ v };
        beta = MDPDETestBeta(f.dataset, e.weights, backend); ++linear_solves;
        v = MDPDETestVariance(f.dataset, f.alpha, e.weights, beta);
        auto row{ Endpoint(f, beta, v) }; ++evaluations;
        row["iteration"] = i + 1;
        row["squared_beta_change"] = Number((beta - previous).squaredNorm());
        row["relative_variance_change"] = Number(std::abs(v - old_v) /
            std::max({std::abs(v), std::abs(old_v), f.options.data_weight_min}));
        const bool valid{ row.at("valid").as_bool() };
        const bool reference{ row.at("reference_pass").as_bool() };
        trace.push_back(std::move(row));
        if (reference) { stop = "fresh-residual"; break; }
        if (!valid) { stop = trace.back().at("reason").as_string().c_str(); break; }
        if (v == old_v && (beta.array() == previous.array()).all()) { stop = "stalled"; break; }
    }
    auto out{ Endpoint(f, beta, v) };
    out["method"] = backend + (own_initialization ? "-own-init" : "-production-init");
    out["stop"] = stop; out["trace"] = std::move(trace);
    out["linear_solves"] = linear_solves; out["equation_evaluations"] = evaluations + 1;
    out["seconds"] = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}

json::object RootJSON(const rhbm_gem::mdpde_detail::MDPDERootResult & root,
    const std::vector<rhbm_gem::mdpde_detail::RootEvaluation> & evaluations, const std::string & label)
{
    auto out{ EquationJSON(root.equations) };
    out["beta"] = Vector(root.beta); out["variance"] = Number(root.variance);
    out["amplitude"] = Number(std::exp(root.beta(0)) * std::pow(2.0 * std::acos(-1.0) / root.beta(1), 1.5));
    out["width"] = Number(1.0 / std::sqrt(root.beta(1)));
    out["method"] = "root-" + label; out["native_status"] = root.native_status;
    out["stop"] = root.stop; out["linear_solves"] = 0;
    out["equation_evaluations"] = root.equation_evaluations;
    out["jacobian_evaluations"] = root.jacobian_evaluations;
    out["iterations"] = root.iterations;
    out["initial_u"] = Vector(root.initial_u); out["u"] = Vector(root.u);
    json::array trace;
    for (const auto & e : evaluations)
    {
        json::object row{{"evaluation",trace.size()+1},{"u",Vector(e.u)},
            {"residual_inf",Number(e.residual_inf.value_or(NAN))}};
        if (!e.reason.empty()) row["reason"] = e.reason;
        else { row["denominator"] = Number(e.denominator); row["floor_count"] = e.floor_count; }
        trace.push_back(std::move(row));
    }
    out["trace"] = std::move(trace);
    if (root.jacobian_condition) out["jacobian_condition"] = Number(*root.jacobian_condition);
    if (root.estimated_remaining_u_error) out["estimated_remaining_u_error"] = Number(*root.estimated_remaining_u_error);
    out["verification_equation_evaluations"] = root.verification_equation_evaluations;
    return out;
}

json::object Root(const ShapeFixture & f, const Eigen::VectorXd & initial, const std::string & label,
    int total_budget = 2000)
{
    const auto start{Clock::now()};
    std::vector<rhbm_gem::mdpde_detail::RootEvaluation> trace;
    const auto root{rhbm_gem::mdpde_detail::SolveMDPDERoot(f.dataset,f.alpha,
        f.options.data_weight_min,initial,total_budget,&trace)};
    auto out{RootJSON(root,trace,label)};
    out["seconds"] = std::chrono::duration<double>(Clock::now()-start).count();
    return out;
}

json::object BranchJSON(const rhbm_gem::mdpde_detail::MDPDEBranchComparison & branch)
{
    json::object out{{"pass",branch.pass},{"coordinate_tolerance",1e-6},{"weight_tolerance",1e-6}};
    if (branch.relative_coordinate_difference)
    {
        json::array change;
        for (double x : *branch.relative_coordinate_difference) change.push_back(Number(x));
        out["relative_coordinate_difference"] = std::move(change);
    }
    if (branch.weight_max_difference) out["weight_max_difference"] = Number(*branch.weight_max_difference);
    if (branch.floor_masks_equal) out["floor_masks_equal"] = *branch.floor_masks_equal;
    return out;
}
} // namespace

void RecordMDPDEIteration(const Eigen::VectorXd & beta, double variance, double change, double variance_change)
{
    if (production_trace) production_trace->push_back(boost::json::object{
        {"beta", Vector(beta)}, {"variance", Number(variance)},
        {"squared_beta_change", Number(change)}, {"relative_variance_change", Number(variance_change)}});
}

ShapeFixture ReadShapeFixture(const std::string & path)
{
    std::ifstream in(path); in.imbue(std::locale::classic());
    std::string kind; int version{}, status{};
    ShapeFixture f;
    if (!(in >> kind >> version) || kind != "shape" || version != 1)
        throw std::runtime_error("Expected shape version 1 fixture.");
    in >> f.alpha >> f.options.thread_size >> f.options.max_iterations >> f.options.tolerance >> f.options.data_weight_min;
    f.dataset.X = ReadMatrix(in); f.dataset.y = ReadMatrix(in);
    std::string variance;
    if (!(in >> status >> variance)) throw std::runtime_error("Missing expected shape result.");
    f.expected.status = static_cast<rhbm_gem::RHBMEstimationStatus>(status);
    f.expected.sigma_square = std::stod(variance);
    f.expected.beta_ols = ReadMatrix(in); f.expected.beta_mdpde = ReadMatrix(in);
    std::string extra;
    if (in >> extra) throw std::runtime_error("Unexpected shape fixture content.");
    return f;
}

boost::json::object EquationJSON(const MDPDEEquationEvidence & e)
{
    const double norm{ e.scaled.size() && e.scaled.allFinite()
        ? e.scaled.lpNorm<Eigen::Infinity>() : std::numeric_limits<double>::infinity() };
    return {{"valid", e.valid}, {"reason", e.reason}, {"raw_residual", Vector(e.raw)},
        {"scaled_residual", Vector(e.scaled)}, {"residual_inf", Number(norm)},
        {"denominator", Number(e.denominator)}, {"rank", e.rank},
        {"singular_values", Vector(e.singular_values)}, {"condition", Number(e.condition)},
        {"floor_count", e.floor_count}, {"equation_pass", e.valid && norm <= 1.0e-8},
        {"reference_pass", e.valid && norm <= 1.0e-10}};
}

boost::json::object CompareMDPDEBranches(const Eigen::VectorXd & candidate_beta, double candidate_variance,
    const MDPDEEquationEvidence & candidate, const Eigen::VectorXd & reference_beta, double reference_variance,
    const MDPDEEquationEvidence & reference, double floor)
{
    return BranchJSON(rhbm_gem::mdpde_detail::CompareMDPDEBranches(candidate_beta,candidate_variance,
        candidate,reference_beta,reference_variance,reference,floor));
}

EndpointRefinementResult RefineMDPDEEndpoint(const ShapeFixture & f, int equation_budget,
    const MDPDEEquationEvidence * initial)
{
    rhbm_gem::mdpde_detail::EndpointRefinementEvidence evidence;
    std::vector<rhbm_gem::mdpde_detail::RootEvaluation> trace;
    auto result{rhbm_gem::mdpde_detail::RefineMDPDEEndpoint(f.dataset,f.alpha,f.options,
        f.expected,equation_budget,initial,&evidence,&trace)};
    const auto & d{*result.refinement};
    json::object out{{"schema_version",1},{"accepted",d.accepted},{"reason",d.reason},
        {"original_status",static_cast<int>(f.expected.status)},
        {"original_iterations",f.expected.diagnostics.iterations},
        {"original_squared_beta_change",Number(f.expected.diagnostics.squared_beta_change.value_or(NAN))},
        {"original_relative_variance_change",Number(f.expected.diagnostics.relative_variance_change.value_or(NAN))},
        {"equation_budget",equation_budget},{"candidate_equation_evaluations",d.candidate_equation_evaluations},
        {"candidate_linear_solves",evidence.candidate_linear_solves},{"reference_updates",d.reference_updates},
        {"reference_equation_evaluations",d.reference_updates},{"original_equations",EquationJSON(evidence.original)}};
    if (evidence.root) out["root"] = RootJSON(*evidence.root,trace,"endpoint");
    if (evidence.reference)
    {
        auto reference{EquationJSON(*evidence.reference)};
        reference["beta"] = Vector(evidence.reference_beta);
        reference["variance"] = Number(evidence.reference_variance);
        reference["stop"] = evidence.reference_stop;
        out["reference"] = std::move(reference);
    }
    if (evidence.branch) out["branch"] = BranchJSON(*evidence.branch);
    const bool accepted{d.accepted};
    return {std::move(result),accepted,std::move(out)};
}

boost::json::object CompareMDPDE(const ShapeFixture & f, bool perturb)
{
    json::array trace;
    production_trace = &trace;
    rhbm_gem::RHBMBetaEstimateResult baseline;
    const auto start{ Clock::now() };
    try { baseline = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha, f.dataset, f.options); }
    catch (...) { production_trace = nullptr; throw; }
    production_trace = nullptr;
    const double elapsed{ std::chrono::duration<double>(Clock::now() - start).count() };
    auto production{ Endpoint(f, baseline.beta_mdpde, baseline.sigma_square) };
    production["method"] = "production"; production["native_status"] = static_cast<int>(baseline.status);
    production["seconds"] = elapsed; production["linear_solves"] = baseline.diagnostics.iterations + 1;
    production["iterations"] = baseline.diagnostics.iterations;
    production["equation_evaluations"] = 0;
    production["verification_equation_evaluations"] = baseline.diagnostics.iterations + 1;
    for (auto & value : trace)
    {
        auto & row{ value.as_object() };
        Eigen::VectorXd b(2);
        for (int k = 0; k < 2; ++k)
        {
            const auto & parameter{row.at("beta").at(static_cast<std::size_t>(k))};
            b(k) = parameter.is_double() ? parameter.as_double() : std::numeric_limits<double>::quiet_NaN();
        }
        const double v{ row.at("variance").is_double() ? row.at("variance").as_double() : 0.0 };
        row["equation"] = EquationJSON(EvaluateMDPDEEquations(f.dataset, f.alpha, b, v, f.options.data_weight_min));
    }
    production["trace"] = std::move(trace);
    const bool exact{ baseline.status == f.expected.status && baseline.sigma_square == f.expected.sigma_square &&
        baseline.beta_ols.size() == f.expected.beta_ols.size() &&
        baseline.beta_mdpde.size() == f.expected.beta_mdpde.size() &&
        (baseline.beta_ols.array() == f.expected.beta_ols.array()).all() &&
        (baseline.beta_mdpde.array() == f.expected.beta_mdpde.array()).all() };
    json::array methods; methods.push_back(std::move(production));
    const auto qr_beta{ MDPDETestBeta(f.dataset, Eigen::VectorXd::Ones(f.dataset.y.size()), "qr") };
    const double roundoff_bound{64.0 * std::numeric_limits<double>::epsilon() *
        (f.dataset.y.norm() + f.dataset.X.norm() * qr_beta.norm())};
    if ((f.dataset.y - f.dataset.X * qr_beta).norm() <= roundoff_bound)
        return {{"schema_version", 1}, {"exact_replay", exact}, {"rows", f.dataset.X.rows()},
            {"alpha", f.alpha}, {"classification", "roundoff-exact-fit-boundary"},
            {"roundoff_residual_bound", roundoff_bound}, {"methods", std::move(methods)}};
    for (const auto * backend : {"normal", "qr", "svd"}) methods.push_back(FixedPoint(f, backend, false));
    for (const auto * backend : {"qr", "svd"}) methods.push_back(FixedPoint(f, backend, true));
    const double initial_v{ (f.dataset.y - f.dataset.X * baseline.beta_ols).squaredNorm() / (f.dataset.y.size() - 1) };
    methods.push_back(Root(f, Coordinates(baseline.beta_ols, initial_v), "ols"));
    const auto end{ Coordinates(baseline.beta_mdpde, baseline.sigma_square) };
    methods.push_back(Root(f, end, "endpoint"));
    if (perturb && end.allFinite()) for (int k = 0; k < 3; ++k) for (const int sign : {-1, 1})
    {
        auto u{ end }; u(k) += sign * 1.0e-3;
        methods.push_back(Root(f, u, "perturb-" + std::to_string(k) + "-" + std::to_string(sign)));
    }
    return {{"schema_version", 1}, {"exact_replay", exact}, {"rows", f.dataset.X.rows()},
        {"alpha", f.alpha}, {"methods", std::move(methods)}};
}
} // namespace second_stage_test
