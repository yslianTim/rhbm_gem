#include "utils/hrl/MDPDEEndpointRefinement.hpp"

#include <unsupported/Eigen/NonLinearOptimization>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace rhbm_gem::mdpde_detail {
namespace {

Eigen::VectorXd Coordinates(const Eigen::VectorXd & beta, double variance)
{
    Eigen::VectorXd u(3);
    u << beta(0), std::log(beta(1)), std::log(variance);
    return u;
}

std::optional<double> Residual(const MDPDEEquationEvidence & evidence)
{
    if (evidence.scaled.size() == 0 || !evidence.scaled.allFinite()) return std::nullopt;
    return evidence.scaled.lpNorm<Eigen::Infinity>();
}

struct RootFunction
{
    const RHBMMemberDataset & dataset;
    double alpha, floor;
    int budget;
    std::vector<RootEvaluation> * trace;
    mutable int evaluations{}, jacobians{};
    mutable std::string failure_reason;

    int operator()(const Eigen::VectorXd & u, Eigen::VectorXd & residual) const
    {
        if (evaluations >= budget) return -1;
        ++evaluations;
        Eigen::VectorXd beta(2);
        beta << u(0), std::exp(u(1));
        const auto e{ EvaluateMDPDEEquations(dataset, alpha, beta, std::exp(u(2)), floor) };
        if (e.scaled.size() != 3 || !e.scaled.allFinite() || e.reason == "invalid-model")
        {
            failure_reason = e.reason == "invalid-model" ? "invalid-model-evaluation" : "nonfinite-evaluation";
            if (trace) trace->push_back({u, {}, 0.0, 0, failure_reason});
            return -1;
        }
        residual = e.scaled;
        if (trace) trace->push_back({u, residual.lpNorm<Eigen::Infinity>(), e.denominator, e.floor_count, {}});
        return 0;
    }

    int df(const Eigen::VectorXd & u, Eigen::MatrixXd & jac) const
    {
        ++jacobians;
        jac.resize(3, 3);
        for (int k = 0; k < 3; ++k)
        {
            const double h{ std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(u(k))) };
            auto plus{u}, minus{u};
            plus(k) += h;
            minus(k) -= h;
            Eigen::VectorXd p(3), m(3);
            if ((*this)(plus, p) < 0 || (*this)(minus, m) < 0) return -1;
            jac.col(k) = (p - m) / (2 * h);
        }
        return 0;
    }
};
} // namespace

MDPDERootResult SolveMDPDERoot(const RHBMMemberDataset & dataset, double alpha, double floor,
    const Eigen::VectorXd & initial, int total_budget, std::vector<RootEvaluation> * trace)
{
    if (total_budget < 9) throw std::invalid_argument("Root budget must include endpoint verification.");
    MDPDERootResult result;
    result.initial_u = initial;
    auto & u{ result.u };
    u = initial;
    const int root_budget{ total_budget - 8 };
    RootFunction function{dataset, alpha, floor, root_budget, trace, 0, 0, {}};
    Eigen::HybridNonLinearSolver<RootFunction> solver(function);
    solver.parameters.maxfev = root_budget;
    solver.parameters.xtol = 1.0e-12;
    solver.parameters.factor = 1.0;
    const bool valid_start{ u.allFinite() && std::isfinite(std::exp(u(1))) && std::exp(u(1)) > 0.0 &&
        std::isfinite(std::exp(u(2))) && std::exp(u(2)) > 0.0 };
    if (valid_start) result.native_status = static_cast<int>(solver.solve(u));
    result.beta.resize(2);
    result.beta << u(0), std::exp(u(1));
    result.variance = std::exp(u(2));
    result.equations = EvaluateMDPDEEquations(dataset, alpha, result.beta, result.variance, floor);
    result.stop = !valid_start ? "invalid-start" : function.evaluations >= root_budget ? "budget-exhausted" :
        !function.failure_reason.empty() ? function.failure_reason :
        (result.equations.valid && Residual(result.equations).value_or(INFINITY) <= 1e-10 ?
            "fresh-residual" : "native-stop-without-reference");
    result.equation_evaluations = function.evaluations;
    result.jacobian_evaluations = function.jacobians;
    result.iterations = valid_start ? static_cast<int>(solver.iter) : 0;

    // Preserve the experiment's independent endpoint verification and its reserved work.
    Eigen::Matrix3d jac;
    bool valid{ u.allFinite() };
    result.verification_equation_evaluations = 1;
    for (int k = 0; valid && k < 3; ++k)
    {
        const double h{ std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(u(k))) };
        auto plus{u}, minus{u};
        plus(k) += h;
        minus(k) -= h;
        auto evaluate = [&](const Eigen::VectorXd & x) {
            Eigen::VectorXd beta(2);
            beta << x(0), std::exp(x(1));
            return EvaluateMDPDEEquations(dataset, alpha, beta, std::exp(x(2)), floor).scaled;
        };
        const auto p{ evaluate(plus) }, m{ evaluate(minus) };
        result.verification_equation_evaluations += 2;
        valid = p.size() == 3 && m.size() == 3 && p.allFinite() && m.allFinite();
        if (valid) jac.col(k) = (p - m) / (2.0 * h);
    }
    if (valid)
    {
        const Eigen::JacobiSVD<Eigen::Matrix3d> svd(jac, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const auto e{ EvaluateMDPDEEquations(dataset, alpha, result.beta, result.variance, floor) };
        ++result.verification_equation_evaluations;
        result.jacobian_condition = svd.singularValues()(0) / svd.singularValues()(2);
        result.estimated_remaining_u_error = svd.solve(e.scaled).lpNorm<Eigen::Infinity>();
    }
    return result;
}

MDPDEBranchComparison CompareMDPDEBranches(const Eigen::VectorXd & candidate_beta, double candidate_variance,
    const MDPDEEquationEvidence & candidate, const Eigen::VectorXd & reference_beta, double reference_variance,
    const MDPDEEquationEvidence & reference, double floor)
{
    MDPDEBranchComparison result;
    if (!candidate.valid || !reference.valid || candidate.weights.size() != reference.weights.size()) return result;
    const auto u{ Coordinates(candidate_beta, candidate_variance) }, v{ Coordinates(reference_beta, reference_variance) };
    std::array<double, 3> change{};
    bool coordinates_match{ true };
    for (int k = 0; k < 3; ++k)
    {
        auto & difference{change[static_cast<std::size_t>(k)]};
        difference = std::abs(u(k) - v(k)) / std::max({1.0, std::abs(u(k)), std::abs(v(k))});
        coordinates_match &= std::isfinite(difference) && difference <= 1e-6;
    }
    result.relative_coordinate_difference = change;
    result.weight_max_difference = (candidate.weights - reference.weights).lpNorm<Eigen::Infinity>();
    result.floor_masks_equal = ((candidate.weights.array() == floor) == (reference.weights.array() == floor)).all();
    result.pass = coordinates_match && *result.weight_max_difference <= 1e-6 && *result.floor_masks_equal;
    return result;
}

RHBMBetaEstimateResult RefineMDPDEEndpoint(const RHBMMemberDataset & dataset, double alpha,
    const RHBMExecutionOptions & options, const RHBMBetaEstimateResult & endpoint, int equation_budget,
    const MDPDEEquationEvidence * initial, EndpointRefinementEvidence * evidence, std::vector<RootEvaluation> * trace)
{
    auto result{ endpoint };
    auto & diagnostics{ result.refinement.emplace() };
    diagnostics.reason = "ineligible";
    diagnostics.candidate_equation_evaluations = 1;
    auto reference{ initial ? *initial : EvaluateMDPDEEquations(dataset, alpha, endpoint.beta_mdpde,
        endpoint.sigma_square, options.data_weight_min) };
    diagnostics.original_residual = Residual(reference);
    if (evidence) evidence->original = reference;
    if (dataset.X.rows() <= dataset.X.cols() || dataset.X.cols() != 2 ||
        !dataset.X.allFinite() || !dataset.y.allFinite() || dataset.y.size() != dataset.X.rows() ||
        dataset.X.colPivHouseholderQr().rank() != 2)
    { diagnostics.reason = reference.reason; return result; }
    const Eigen::VectorXd qr{ dataset.X.colPivHouseholderQr().solve(dataset.y) };
    if (evidence) evidence->candidate_linear_solves = 1;
    const double roundoff{ 64.0 * std::numeric_limits<double>::epsilon() *
        (dataset.y.norm() + dataset.X.norm() * qr.norm()) };
    if ((dataset.y - dataset.X * qr).norm() <= roundoff)
    { diagnostics.reason = "roundoff-exact-fit-boundary"; return result; }
    if (!reference.valid) { diagnostics.reason = reference.reason; return result; }
    if (equation_budget < 11) { diagnostics.reason = "budget-exhausted"; return result; }
    const auto root{ SolveMDPDERoot(dataset, alpha, options.data_weight_min,
        Coordinates(endpoint.beta_mdpde, endpoint.sigma_square), equation_budget - 2, trace) };
    diagnostics.candidate_equation_evaluations += root.equation_evaluations + root.verification_equation_evaluations;
    diagnostics.candidate_residual = Residual(root.equations);
    if (evidence) evidence->root = root;
    if (root.stop == "budget-exhausted") { diagnostics.reason = "budget-exhausted"; return result; }
    if (!root.equations.valid || diagnostics.candidate_residual.value_or(INFINITY) > 1e-8)
    { diagnostics.reason = "root-equations-unqualified"; return result; }
    const auto candidate{ EvaluateMDPDEEquations(dataset, alpha, root.beta, root.variance, options.data_weight_min) };
    ++diagnostics.candidate_equation_evaluations;
    auto reference_beta{ endpoint.beta_mdpde };
    double reference_variance{ endpoint.sigma_square };
    std::string reference_stop{ "budget-exhausted" };
    for (;;)
    {
        if (!reference.valid) { reference_stop = reference.reason; break; }
        if (Residual(reference).value_or(INFINITY) <= 1e-10) { reference_stop = "fresh-residual"; break; }
        if (diagnostics.reference_updates + endpoint.diagnostics.iterations >= 10000) break;
        const auto previous{ reference_beta };
        const double old_v{ reference_variance };
        reference_beta = CalculateMDPDEBeta(dataset, reference.weights);
        reference_variance = CalculateMDPDEVariance(dataset, alpha, reference.weights, reference_beta);
        ++diagnostics.reference_updates;
        reference = EvaluateMDPDEEquations(dataset, alpha, reference_beta, reference_variance, options.data_weight_min);
        if (reference_variance == old_v && (reference_beta.array() == previous.array()).all() &&
            (!reference.valid || Residual(reference).value_or(INFINITY) > 1e-10))
        { reference_stop = "stalled"; break; }
    }
    diagnostics.reference_residual = Residual(reference);
    diagnostics.reference_stop = reference_stop;
    if (evidence)
    {
        evidence->reference = reference;
        evidence->reference_beta = reference_beta;
        evidence->reference_variance = reference_variance;
        evidence->reference_stop = reference_stop;
    }
    if (reference_stop != "fresh-residual") { diagnostics.reason = "reference-unqualified"; return result; }
    const auto branch{ CompareMDPDEBranches(root.beta, root.variance, candidate,
        reference_beta, reference_variance, reference, options.data_weight_min) };
    if (evidence) evidence->branch = branch;
    diagnostics.relative_coordinate_difference = branch.relative_coordinate_difference;
    diagnostics.weight_max_difference = branch.weight_max_difference;
    diagnostics.floor_masks_equal = branch.floor_masks_equal;
    if (!branch.pass) { diagnostics.reason = "branch-mismatch"; return result; }
    const auto covariance{ CalculateMDPDECovariance(root.variance, candidate.weights) };
    if (!covariance.diagonal().allFinite() || (covariance.diagonal().array() <= 0.0).any())
    { diagnostics.reason = "invalid-covariance"; return result; }
    result.beta_mdpde = root.beta;
    result.sigma_square = root.variance;
    result.data_weight = candidate.weights.asDiagonal();
    result.data_covariance = covariance;
    diagnostics.accepted = true;
    diagnostics.reason = "accepted";
    return result;
}

RHBMBetaEstimateResult ApplyFailedOnlyRefinement(const RHBMMemberDataset & dataset, double alpha,
    const RHBMExecutionOptions & options, const RHBMBetaEstimateResult & endpoint)
{
    if (endpoint.status == RHBMEstimationStatus::SUCCESS) return endpoint;
    return RefineMDPDEEndpoint(dataset, alpha, options, endpoint);
}

} // namespace rhbm_gem::mdpde_detail
