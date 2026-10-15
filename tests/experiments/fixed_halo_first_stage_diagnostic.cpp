#include "core/detail/FittingWorkset.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/gaussian_fit/FittingRanges.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "utils/hrl/MDPDEEndpointRefinement.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <boost/json.hpp>

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace json = boost::json;
using Clock = std::chrono::steady_clock;

constexpr std::array<int, 4> kFailedHaloIds{16664, 57190, 62463, 63216};

json::value Number(double value)
{
    return std::isfinite(value) ? json::value(value) : json::value(nullptr);
}

json::array Vector(const Eigen::VectorXd & values)
{
    json::array result;
    for (const auto value : values) result.emplace_back(Number(value));
    return result;
}

json::array Matrix(const Eigen::MatrixXd & values)
{
    json::array result;
    for (Eigen::Index row = 0; row < values.rows(); ++row)
    {
        json::array item;
        for (Eigen::Index column = 0; column < values.cols(); ++column)
            item.emplace_back(Number(values(row, column)));
        result.emplace_back(std::move(item));
    }
    return result;
}

std::string StatusText(rhbm_gem::RHBMEstimationStatus status)
{
    using rhbm_gem::RHBMEstimationStatus;
    switch (status)
    {
        case RHBMEstimationStatus::SUCCESS: return "SUCCESS";
        case RHBMEstimationStatus::SINGLE_MEMBER: return "SINGLE_MEMBER";
        case RHBMEstimationStatus::INSUFFICIENT_DATA: return "INSUFFICIENT_DATA";
        case RHBMEstimationStatus::MAX_ITERATIONS_REACHED: return "MAX_ITERATIONS_REACHED";
        case RHBMEstimationStatus::NUMERICAL_FALLBACK: return "NUMERICAL_FALLBACK";
    }
    return "UNKNOWN";
}

std::string QualificationText(rhbm_gem::RHBMSolveQualification qualification)
{
    using rhbm_gem::RHBMSolveQualification;
    switch (qualification)
    {
        case RHBMSolveQualification::NativeSuccess: return "NativeSuccess";
        case RHBMSolveQualification::RefinedSuccess: return "RefinedSuccess";
        case RHBMSolveQualification::Unqualified: return "Unqualified";
    }
    return "Unqualified";
}

json::object FitRecord(const rhbm_gem::RHBMMemberDataset & dataset, double alpha,
    const rhbm_gem::RHBMExecutionOptions & options,
    const rhbm_gem::RHBMBetaEstimateResult & fit, double seconds)
{
    using namespace rhbm_gem;
    const auto equations = mdpde_detail::EvaluateMDPDEEquations(
        dataset, alpha, fit.beta_mdpde, fit.sigma_square, options.data_weight_min);
    const auto model = linearization_service::DecodeParameterVector(fit.beta_mdpde);
    Eigen::JacobiSVD<Eigen::MatrixXd> design_svd(dataset.X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const double design_condition = design_svd.singularValues().size() >= 2 &&
        design_svd.singularValues()(1) > 0.0 ?
        design_svd.singularValues()(0) / design_svd.singularValues()(1) :
        std::numeric_limits<double>::infinity();
    const double residual = equations.valid && equations.scaled.size() > 0 && equations.scaled.allFinite() ?
        equations.scaled.lpNorm<Eigen::Infinity>() : std::numeric_limits<double>::infinity();
    json::object result{
        {"max_iterations", options.max_iterations},
        {"native_status", StatusText(fit.status)},
        {"qualification", QualificationText(fit.Qualification())},
        {"iterations", fit.diagnostics.iterations},
        {"final_squared_beta_change", fit.diagnostics.squared_beta_change ?
            Number(*fit.diagnostics.squared_beta_change) : json::value(nullptr)},
        {"final_relative_variance_change", fit.diagnostics.relative_variance_change ?
            Number(*fit.diagnostics.relative_variance_change) : json::value(nullptr)},
        {"beta_ols", Vector(fit.beta_ols)},
        {"beta_mdpde", Vector(fit.beta_mdpde)},
        {"A", Number(model.GetAmplitude())},
        {"B", Number(model.GetWidth())},
        {"variance", Number(fit.sigma_square)},
        {"variance_above_weight_floor", fit.sigma_square > options.data_weight_min},
        {"estimating_equation_residual_inf", Number(residual)},
        {"estimating_equations_valid", equations.valid},
        {"estimating_equations_reason", equations.reason},
        {"equation_effective_rank", equations.rank},
        {"equation_condition", Number(equations.condition)},
        {"weight_floor_count", equations.floor_count},
        {"design_rank", dataset.X.colPivHouseholderQr().rank()},
        {"design_condition", Number(design_condition)},
        {"wall_seconds", seconds}
    };
    if (fit.refinement)
    {
        const auto & refinement = *fit.refinement;
        result["recovery"] = json::object{
            {"attempted", true}, {"accepted", refinement.accepted}, {"reason", refinement.reason},
            {"candidate_equation_evaluations", refinement.candidate_equation_evaluations},
            {"reference_updates", refinement.reference_updates},
            {"reference_stop", refinement.reference_stop},
            {"original_residual", refinement.original_residual ? Number(*refinement.original_residual) : json::value(nullptr)},
            {"candidate_residual", refinement.candidate_residual ? Number(*refinement.candidate_residual) : json::value(nullptr)},
            {"reference_residual", refinement.reference_residual ? Number(*refinement.reference_residual) : json::value(nullptr)},
            {"weight_max_difference", refinement.weight_max_difference ? Number(*refinement.weight_max_difference) : json::value(nullptr)},
            {"floor_masks_equal", refinement.floor_masks_equal ? json::value(*refinement.floor_masks_equal) : json::value(nullptr)}
        };
    }
    else result["recovery"] = json::object{{"attempted", false}};
    return result;
}

json::object FixedPointTrajectory(const rhbm_gem::RHBMMemberDataset & dataset, double alpha,
    const rhbm_gem::RHBMExecutionOptions & options)
{
    using namespace rhbm_gem;
    Eigen::VectorXd beta = mdpde_detail::CalculateMDPDEBeta(
        dataset, Eigen::VectorXd::Ones(dataset.y.size()));
    double variance = (dataset.y - dataset.X * beta).squaredNorm() /
        static_cast<double>(dataset.y.size() - 1);
    double minimum_variance = variance;
    double previous_variance_delta{};
    std::array<double, 2> previous_beta_delta{};
    std::array<int, 2> beta_direction_reversals{};
    int variance_direction_reversals{};
    int nonfinite_intermediates{};
    bool converged{};
    int iterations{};
    json::array trajectory;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration)
    {
        const auto equations = mdpde_detail::EvaluateMDPDEEquations(
            dataset, alpha, beta, variance, options.data_weight_min);
        if (!equations.valid)
        {
            ++nonfinite_intermediates;
            break;
        }
        const auto next_beta = mdpde_detail::CalculateMDPDEBeta(dataset, equations.weights);
        const double next_variance = mdpde_detail::CalculateMDPDEVariance(
            dataset, alpha, equations.weights, next_beta);
        const double beta_change = (next_beta - beta).squaredNorm();
        const double variance_delta = next_variance - variance;
        const double variance_change = std::abs(variance_delta) /
            std::max({std::abs(next_variance), std::abs(variance), options.data_weight_min});
        const auto endpoint_equations = mdpde_detail::EvaluateMDPDEEquations(
            dataset, alpha, next_beta, next_variance, options.data_weight_min);
        const double equation_residual = endpoint_equations.valid && endpoint_equations.scaled.allFinite() ?
            endpoint_equations.scaled.lpNorm<Eigen::Infinity>() :
            std::numeric_limits<double>::infinity();
        if (!next_beta.allFinite() || !std::isfinite(next_variance) || !std::isfinite(beta_change) ||
            !std::isfinite(variance_change))
            ++nonfinite_intermediates;
        if (std::isfinite(previous_variance_delta) && previous_variance_delta * variance_delta < 0.0)
            ++variance_direction_reversals;
        for (int k = 0; k < 2; ++k)
        {
            const double delta = next_beta(k) - beta(k);
            if (previous_beta_delta[static_cast<std::size_t>(k)] * delta < 0.0)
                ++beta_direction_reversals[static_cast<std::size_t>(k)];
            previous_beta_delta[static_cast<std::size_t>(k)] = delta;
        }
        trajectory.emplace_back(json::object{
            {"iteration", iteration}, {"beta_squared_change", Number(beta_change)},
            {"relative_variance_change", Number(variance_change)}, {"variance", Number(next_variance)},
            {"equation_residual_inf", Number(equation_residual)},
            {"weight_floor_count", endpoint_equations.floor_count}
        });
        minimum_variance = std::min(minimum_variance, next_variance);
        iterations = iteration;
        if (beta_change < options.tolerance && variance_change < options.tolerance)
        {
            converged = true;
            break;
        }
        beta = next_beta;
        variance = next_variance;
        previous_variance_delta = variance_delta;
    }
    return {
        {"converged", converged}, {"iterations", iterations},
        {"final_squared_beta_change", trajectory.empty() ? json::value(nullptr) :
            trajectory.back().as_object().at("beta_squared_change")},
        {"final_relative_variance_change", trajectory.empty() ? json::value(nullptr) :
            trajectory.back().as_object().at("relative_variance_change")},
        {"minimum_variance", Number(minimum_variance)},
        {"variance_direction_reversals", variance_direction_reversals},
        {"beta_direction_reversals", json::array{beta_direction_reversals[0], beta_direction_reversals[1]}},
        {"nonfinite_intermediates", nonfinite_intermediates}, {"iterations_trace", std::move(trajectory)}
    };
}

json::object DiagnoseAtom(const rhbm_gem::AtomObject & atom,
    const rhbm_gem::JointInitializationAtom & record,
    std::size_t joint_support_rows, const std::string & model_path, const std::string & map_path,
    bool normalization_applied, double normalization_divisor)
{
    using namespace rhbm_gem;
    const auto view = AtomLocalPotentialView::For(atom);
    const auto raw_samples = view.GetRawSamplingEntries(false);
    const auto samples = view.GetSamplingEntries(FittingStage::First);
    std::vector<double> responses;
    responses.reserve(samples.size());
    for (const auto & sample : samples) responses.push_back(sample.response);
    const core::detail::PreparedLocalGaussianDesign design(samples, 0.0, core::detail::kSignalDistanceMax);
    const auto dataset = design.BuildDataset(responses, GaussianModel3D{0.0, 1.0, 0.0});
    const double alpha = view.GetAlphaR(FittingStage::First);

    const auto baseline_start = Clock::now();
    RHBMExecutionOptions baseline_options{.thread_size = 1};
    const auto baseline = rhbm_helper::EstimateBetaMDPDE(alpha, dataset, baseline_options);
    const double baseline_seconds = std::chrono::duration<double>(Clock::now() - baseline_start).count();

    json::array selected;
    for (const auto & sample : raw_samples)
        selected.emplace_back(json::object{
            {"distance", sample.point.distance}, {"response", sample.response},
            {"selected", sample.point.is_selected}
        });

    json::object runs;
    for (const int budget : {100, 300, 1000})
    {
        RHBMExecutionOptions options{.thread_size = 1};
        options.max_iterations = budget;
        const auto start = Clock::now();
        const auto fit = budget == 100 ? baseline : rhbm_helper::EstimateBetaMDPDE(alpha, dataset, options);
        const double seconds = budget == 100 ? baseline_seconds :
            std::chrono::duration<double>(Clock::now() - start).count();
        runs["native_" + std::to_string(budget)] = FitRecord(dataset, alpha, options, fit, seconds);
    }

    const auto recovery_start = Clock::now();
    mdpde_detail::EndpointRefinementEvidence refinement_evidence;
    const auto recovered = mdpde_detail::RefineMDPDEEndpoint(dataset, alpha, baseline_options, baseline,
        mdpde_detail::kEndpointEquationBudget, nullptr, &refinement_evidence);
    const double recovery_seconds = std::chrono::duration<double>(Clock::now() - recovery_start).count();
    auto recovery_record = FitRecord(dataset, alpha, baseline_options, recovered, recovery_seconds);
    if (recovered.refinement)
    {
        auto & recovery = recovery_record["recovery"].as_object();
        recovery["method"] = "existing-failed-only-endpoint-refinement";
        recovery["root"] = refinement_evidence.root ? json::value(json::object{
            {"stop", refinement_evidence.root->stop},
            {"iterations", refinement_evidence.root->iterations},
            {"equation_evaluations", refinement_evidence.root->equation_evaluations},
            {"residual_inf", refinement_evidence.root->equations.valid &&
                refinement_evidence.root->equations.scaled.allFinite() ?
                Number(refinement_evidence.root->equations.scaled.lpNorm<Eigen::Infinity>()) : json::value(nullptr)},
            {"effective_rank", refinement_evidence.root->equations.rank},
            {"condition", Number(refinement_evidence.root->equations.condition)},
            {"jacobian_condition", refinement_evidence.root->jacobian_condition ?
                Number(*refinement_evidence.root->jacobian_condition) : json::value(nullptr)},
            {"estimated_remaining_u_error", refinement_evidence.root->estimated_remaining_u_error ?
                Number(*refinement_evidence.root->estimated_remaining_u_error) : json::value(nullptr)}
        }) : json::value(nullptr);
        recovery["branch_check_passed"] = refinement_evidence.branch ?
            json::value(refinement_evidence.branch->pass) : json::value(nullptr);
        recovery["candidate_linear_solves"] = refinement_evidence.candidate_linear_solves;
        recovery["covariance_valid"] = recovered.data_covariance.size() > 0 &&
            recovered.data_covariance.diagonal().allFinite() &&
            (recovered.data_covariance.diagonal().array() > 0.0).all();
    }
    runs["existing_failed_only_refinement"] = std::move(recovery_record);

    const auto ols_model = linearization_service::DecodeParameterVector(baseline.beta_ols);
    const Eigen::VectorXd qr = dataset.X.colPivHouseholderQr().solve(dataset.y);
    const double exact_fit_roundoff = 64.0 * std::numeric_limits<double>::epsilon() *
        (dataset.y.norm() + dataset.X.norm() * qr.norm());
    RHBMExecutionOptions trajectory_options{.thread_size = 1};
    trajectory_options.max_iterations = 1000;
    return {
        {"atom_id", record.id},
        {"source", json::object{
            {"model", model_path},
            {"map", map_path},
            {"map_normalization_applied", normalization_applied != 0},
            {"map_normalization_divisor", Number(normalization_divisor)},
            {"selection", "chain A non-hydrogen contributors; FirstStage selected samples"}
        }},
        {"native_status", StatusText(baseline.status)},
        {"joint_support_rows", joint_support_rows},
        {"first_stage_sample_count", record.sample_count},
        {"selected_raw_sample_count", samples.size()},
        {"local_alpha", Number(alpha)},
        {"native_solver_options", json::object{
            {"max_iterations", 100}, {"tolerance", baseline_options.tolerance},
            {"data_weight_min", baseline_options.data_weight_min}, {"thread_size", baseline_options.thread_size}
        }},
        {"first_stage_native_diagnostics", json::object{
            {"iterations", baseline.diagnostics.iterations},
            {"final_squared_beta_change", baseline.diagnostics.squared_beta_change ?
                Number(*baseline.diagnostics.squared_beta_change) : json::value(nullptr)},
            {"final_relative_variance_change", baseline.diagnostics.relative_variance_change ?
                Number(*baseline.diagnostics.relative_variance_change) : json::value(nullptr)}
        }},
        {"initial_ols", json::object{
            {"beta", Vector(baseline.beta_ols)}, {"A", Number(ols_model.GetAmplitude())},
            {"B", Number(ols_model.GetWidth())},
            {"residual_norm", Number((dataset.y - dataset.X * qr).norm())},
            {"roundoff_exact_fit_threshold", Number(exact_fit_roundoff)},
            {"roundoff_exact_fit", (dataset.y - dataset.X * qr).norm() <= exact_fit_roundoff}
        }},
        {"samples", std::move(selected)},
        {"dataset", json::object{
            {"selection_rule", "selected raw samples with distance in [0, signal-distance-max] and response > 0"},
            {"design_matrix", Matrix(dataset.X)}, {"log_response", Vector(dataset.y)}
        }},
        {"native_fixed_point_trajectory_to_1000", FixedPointTrajectory(dataset, alpha, trajectory_options)},
        {"runs", std::move(runs)}
    };
}
}

int main(int argc, char ** argv)
{
    try
    {
        if (argc != 4) throw std::invalid_argument("Usage: fixed_halo_first_stage_diagnostic MODEL MAP OUTPUT.json");
        auto model = rhbm_gem::ReadModel(argv[1]);
        auto map = rhbm_gem::ReadMap(argv[2]);
        model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
            return atom.GetElement() != Element::HYDROGEN && atom.GetChainID() == "A";
        });
        const auto normalization_divisor = map->MapValueArrayNormalization();
        model->EditAnalysis().InitializeFromSelection();
        const auto problem = rhbm_gem::core::BuildJointProblem(*map, *model);
        const auto workset = rhbm_gem::core::detail::MakeJointFittingWorkset(*model, problem);
        rhbm_gem::core::FitOptions fit_options;
        fit_options.thread_size = 1;
        fit_options.quiet_mode = true;
        const auto initialization = rhbm_gem::core::detail::PrepareJointWidthInitialization(
            *map, *model, workset, fit_options);

        json::array failures;
        for (std::size_t index = 0; index < workset.contributors.size(); ++index)
        {
            const auto & record = initialization.atoms.at(index);
            if (workset.target_mask[index] || !record.native_status ||
                *record.native_status == static_cast<int>(rhbm_gem::RHBMEstimationStatus::SUCCESS)) continue;
            const auto atom_id = workset.contributors[index]->GetSerialID();
            if (std::find(kFailedHaloIds.begin(), kFailedHaloIds.end(), atom_id) == kFailedHaloIds.end()) continue;
            failures.emplace_back(DiagnoseAtom(*workset.contributors[index], record,
                problem.Input().support[index].size(), argv[1], argv[2], normalization_divisor.has_value(),
                normalization_divisor.value_or(1.0)));
        }
        if (failures.size() != kFailedHaloIds.size())
            throw std::runtime_error("Expected all four reported failed halo IDs; captured " + std::to_string(failures.size()));

        json::object output{
            {"fixture_contract", "6z6u-first-stage-mdpde-diagnostic-v1"},
            {"source", json::object{{"model", argv[1]}, {"map", argv[2]},
                {"normalization", "MapValueArrayNormalization() before sampling"},
                {"sampling", "FibonacciDeterministic"}, {"joint_chain", "A"},
                {"supported_halos", static_cast<std::uint64_t>(problem.Input().atom_ids.size() -
                    problem.Input().selection_domain->target_indices.size())},
                {"observations", static_cast<std::uint64_t>(problem.Input().observations.size())},
                {"targets", static_cast<std::uint64_t>(problem.Input().selection_domain->target_indices.size())},
                {"map_normalization_applied", normalization_divisor.has_value()},
                {"map_normalization_divisor", Number(normalization_divisor.value_or(1.0))}}},
            {"native_defaults", json::object{{"max_iterations", 100}, {"tolerance", 1.0e-5},
                {"data_weight_min", 1.0e-8}}},
            {"comparison_cases", json::array{"native_100", "native_300", "native_1000",
                "existing_failed_only_refinement"}},
            {"failed_halos", std::move(failures)}
        };
        std::ofstream file(argv[3]);
        file.exceptions(std::ios::failbit | std::ios::badbit);
        file << json::serialize(output) << '\n';
        std::cout << json::serialize(json::object{{"output", argv[3]},
            {"failed_halos", output.at("failed_halos").as_array().size()},
            {"normalization_divisor", output.at("source").as_object().at("map_normalization_divisor")}}) << '\n';
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
