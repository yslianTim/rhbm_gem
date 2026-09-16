#include "support/MDPDEExperiment.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
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

struct RootFunction
{
    const ShapeFixture & fixture;
    mutable json::array trace;
    mutable int evaluations{};
    mutable int jacobians{};
    mutable std::string failure_reason;
    int budget{1992};
    int operator()(const Eigen::VectorXd & u, Eigen::VectorXd & residual) const
    {
        if (evaluations >= budget) return -1;
        ++evaluations;
        Eigen::VectorXd beta(2); beta << u(0), std::exp(u(1));
        const double v{ std::exp(u(2)) };
        auto e{ EvaluateMDPDEEquations(fixture.dataset, fixture.alpha, beta, v,
            fixture.options.data_weight_min) };
        if (e.scaled.size() != 3 || !e.scaled.allFinite() || e.reason == "invalid-model")
        {
            failure_reason = e.reason == "invalid-model" ? "invalid-model-evaluation" : "nonfinite-evaluation";
            trace.push_back(json::object{{"evaluation", evaluations}, {"u", Vector(u)},
                {"residual_inf", nullptr}, {"reason", failure_reason}});
            return -1;
        }
        residual = e.scaled;
        trace.push_back(json::object{{"evaluation", evaluations}, {"u", Vector(u)},
            {"residual_inf", residual.lpNorm<Eigen::Infinity>()},
            {"denominator", Number(e.denominator)}, {"floor_count", e.floor_count}});
        return 0;
    }
    int df(const Eigen::VectorXd & u, Eigen::MatrixXd & jac) const
    {
        ++jacobians; jac.resize(3,3);
        for (int k = 0; k < 3; ++k)
        {
            const double h{std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1.0,std::abs(u(k)))};
            auto plus{u}, minus{u}; plus(k) += h; minus(k) -= h;
            Eigen::VectorXd p(3), m(3);
            if ((*this)(plus,p) < 0 || (*this)(minus,m) < 0) return -1;
            jac.col(k) = (p-m)/(2*h);
        }
        return 0;
    }
};

json::object Root(const ShapeFixture & f, const Eigen::VectorXd & initial, const std::string & label,
    int total_budget = 2000)
{
    if (total_budget < 9) throw std::invalid_argument("Root budget must include endpoint verification.");
    const int root_solve_budget{total_budget - 8};
    const auto start{ Clock::now() };
    auto u{ initial };
    RootFunction function{ f, {}, 0, 0, {}, root_solve_budget };
    Eigen::HybridNonLinearSolver<RootFunction> solver(function);
    solver.parameters.maxfev = root_solve_budget;
    solver.parameters.xtol = 1.0e-12;
    solver.parameters.factor = 1.0;
    const bool valid_start{u.allFinite() && std::isfinite(std::exp(u(1))) && std::exp(u(1)) > 0.0 &&
        std::isfinite(std::exp(u(2))) && std::exp(u(2)) > 0.0};
    int status{};
    if (valid_start) status = static_cast<int>(solver.solve(u));
    Eigen::VectorXd beta(2); beta << u(0), std::exp(u(1));
    auto out{ Endpoint(f, beta, std::exp(u(2))) };
    out["method"] = "root-" + label; out["native_status"] = status;
    out["stop"] = !valid_start ? "invalid-start" : function.evaluations >= root_solve_budget ? "budget-exhausted" :
        !function.failure_reason.empty() ? function.failure_reason :
        (out.at("reference_pass").as_bool() ? "fresh-residual" : "native-stop-without-reference");
    out["linear_solves"] = 0;
    out["equation_evaluations"] = function.evaluations;
    out["jacobian_evaluations"] = function.jacobians;
    out["iterations"] = valid_start ? solver.iter : 0;
    out["initial_u"] = Vector(initial); out["u"] = Vector(u);
    out["trace"] = std::move(function.trace);
    // Independent central-difference check uses the reserved equation calls.
    Eigen::Matrix3d jac;
    bool valid{ u.allFinite() };
    int verification_evaluations{1};
    for (int k = 0; valid && k < 3; ++k)
    {
        const double h{ std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(u(k))) };
        auto plus{ u }, minus{ u }; plus(k) += h; minus(k) -= h;
        auto evaluate = [&](const Eigen::VectorXd & x) {
            Eigen::VectorXd b(2); b << x(0), std::exp(x(1));
            return EvaluateMDPDEEquations(f.dataset, f.alpha, b, std::exp(x(2)), f.options.data_weight_min).scaled;
        };
        const auto p{ evaluate(plus) }, m{ evaluate(minus) };
        verification_evaluations += 2;
        valid = p.size() == 3 && m.size() == 3 && p.allFinite() && m.allFinite();
        if (valid) jac.col(k) = (p - m) / (2.0 * h);
    }
    if (valid)
    {
        const Eigen::JacobiSVD<Eigen::Matrix3d> svd(jac, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const auto e{ EvaluateMDPDEEquations(f.dataset, f.alpha, beta, std::exp(u(2)), f.options.data_weight_min) };
        ++verification_evaluations;
        out["jacobian_condition"] = Number(svd.singularValues()(0) / svd.singularValues()(2));
        out["estimated_remaining_u_error"] = Number(svd.solve(e.scaled).lpNorm<Eigen::Infinity>());
    }
    out["verification_equation_evaluations"] = verification_evaluations;
    out["seconds"] = std::chrono::duration<double>(Clock::now() - start).count();
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
    json::object out{{"pass", false}, {"coordinate_tolerance", 1e-6}, {"weight_tolerance", 1e-6}};
    if (!candidate.valid || !reference.valid || candidate.weights.size() != reference.weights.size()) return out;
    const auto u{Coordinates(candidate_beta, candidate_variance)}, v{Coordinates(reference_beta, reference_variance)};
    Eigen::VectorXd change(3);
    for (int k = 0; k < 3; ++k) change(k) = std::abs(u(k)-v(k)) / std::max({1.0,std::abs(u(k)),std::abs(v(k))});
    const double weights{(candidate.weights-reference.weights).lpNorm<Eigen::Infinity>()};
    const bool floors{((candidate.weights.array() == floor) == (reference.weights.array() == floor)).all()};
    out["relative_coordinate_difference"] = Vector(change);
    out["weight_max_difference"] = Number(weights); out["floor_masks_equal"] = floors;
    out["pass"] = change.allFinite() && change.maxCoeff() <= 1e-6 && weights <= 1e-6 && floors;
    return out;
}

EndpointRefinementResult RefineMDPDEEndpoint(const ShapeFixture & f, int equation_budget,
    const MDPDEEquationEvidence * initial)
{
    using rhbm_gem::RHBMEstimationStatus;
    EndpointRefinementResult output{f.expected, false, {}};
    auto & out{output.evidence};
    out = {{"schema_version", 1}, {"accepted", false}, {"reason", "ineligible"},
        {"original_status", static_cast<int>(f.expected.status)},
        {"original_iterations", f.expected.diagnostics.iterations},
        {"original_squared_beta_change", Number(f.expected.diagnostics.squared_beta_change.value_or(NAN))},
        {"original_relative_variance_change", Number(f.expected.diagnostics.relative_variance_change.value_or(NAN))},
        {"equation_budget", equation_budget}, {"candidate_equation_evaluations", 1},
        {"candidate_linear_solves", 0}, {"reference_updates", 0}, {"reference_equation_evaluations", 0}};
    // Preserve the original numerical endpoint on rejection, including its weights/covariance.
    output.result.status = RHBMEstimationStatus::NUMERICAL_FALLBACK;
    auto reference{initial ? *initial : EvaluateMDPDEEquations(f.dataset, f.alpha, f.expected.beta_mdpde,
        f.expected.sigma_square, f.options.data_weight_min)};
    out["original_equations"] = EquationJSON(reference);
    if (f.dataset.X.rows() <= f.dataset.X.cols() || f.dataset.X.cols() != 2 ||
        !f.dataset.X.allFinite() || !f.dataset.y.allFinite() ||
        f.dataset.X.colPivHouseholderQr().rank() != 2)
    { out["reason"] = reference.reason; return output; }
    const auto qr{MDPDETestBeta(f.dataset, Eigen::VectorXd::Ones(f.dataset.y.size()), "qr")};
    out["candidate_linear_solves"] = 1;
    const double roundoff{64.0 * std::numeric_limits<double>::epsilon() *
        (f.dataset.y.norm() + f.dataset.X.norm() * qr.norm())};
    if ((f.dataset.y-f.dataset.X*qr).norm() <= roundoff)
    { out["reason"] = "roundoff-exact-fit-boundary"; return output; }
    if (!reference.valid) { out["reason"] = reference.reason; return output; }
    if (equation_budget < 11) { out["reason"] = "budget-exhausted"; return output; }
    auto root{Root(f, Coordinates(f.expected.beta_mdpde, f.expected.sigma_square), "endpoint", equation_budget-2)};
    const int work{1 + static_cast<int>(root.at("equation_evaluations").as_int64()) +
        static_cast<int>(root.at("verification_equation_evaluations").as_int64())};
    out["candidate_equation_evaluations"] = work;
    out["root"] = root;
    if (root.at("stop") == "budget-exhausted") { out["reason"] = "budget-exhausted"; return output; }
    if (!root.at("equation_pass").as_bool()) { out["reason"] = "root-equations-unqualified"; return output; }
    Eigen::Vector2d beta;
    beta << root.at("beta").at(0).as_double(), root.at("beta").at(1).as_double();
    const double variance{root.at("variance").as_double()};
    const auto candidate{EvaluateMDPDEEquations(f.dataset,f.alpha,beta,variance,f.options.data_weight_min)};
    out["candidate_equation_evaluations"] = work+1;
    auto reference_beta{f.expected.beta_mdpde};
    double reference_variance{f.expected.sigma_square};
    int updates{};
    std::string reference_stop{"budget-exhausted"};
    for (;;)
    {
        if (!reference.valid) { reference_stop = reference.reason; break; }
        if (reference.scaled.lpNorm<Eigen::Infinity>() <= 1e-10) { reference_stop = "fresh-residual"; break; }
        if (updates + f.expected.diagnostics.iterations >= 10000) break;
        const auto previous{reference_beta}; const double old_v{reference_variance};
        reference_beta = MDPDETestBeta(f.dataset,reference.weights,"normal");
        reference_variance = MDPDETestVariance(f.dataset,f.alpha,reference.weights,reference_beta);
        ++updates;
        reference = EvaluateMDPDEEquations(f.dataset,f.alpha,reference_beta,reference_variance,f.options.data_weight_min);
        if (reference_variance == old_v && (reference_beta.array() == previous.array()).all() &&
            (!reference.valid || reference.scaled.lpNorm<Eigen::Infinity>() > 1e-10))
        { reference_stop = "stalled"; break; }
    }
    auto reference_json{EquationJSON(reference)};
    reference_json["beta"] = Vector(reference_beta); reference_json["variance"] = Number(reference_variance);
    reference_json["stop"] = reference_stop;
    out["reference"] = std::move(reference_json);
    out["reference_updates"] = updates; out["reference_equation_evaluations"] = updates;
    if (reference_stop != "fresh-residual") { out["reason"] = "reference-unqualified"; return output; }
    auto branch{CompareMDPDEBranches(beta,variance,candidate,reference_beta,reference_variance,reference,f.options.data_weight_min)};
    out["branch"] = branch;
    if (!branch.at("pass").as_bool()) { out["reason"] = "branch-mismatch"; return output; }
    output.result.beta_mdpde = beta; output.result.sigma_square = variance;
    output.result.data_weight = candidate.weights.asDiagonal();
    output.result.data_covariance = MDPDETestCovariance(variance,candidate.weights);
    output.result.status = RHBMEstimationStatus::SUCCESS;
    output.accepted = true; out["accepted"] = true; out["reason"] = "accepted";
    return output;
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
