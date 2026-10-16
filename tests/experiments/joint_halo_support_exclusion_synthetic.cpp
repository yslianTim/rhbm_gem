#include "support/JointHaloObservationExclusion.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace json = boost::json;
namespace core = rhbm_gem::core;
namespace joint = core::joint_component;
namespace audit = joint_test::halo_exclusion;
using Input = core::JointProblemInput;
using Position = audit::Position;
constexpr std::array<double, 3> SimulationRadii{2.5, 3.5, 5.0};
constexpr std::array<double, 4> NoiseSigmas{0.0, 0.005, 0.02, 0.05};
constexpr std::array<std::uint64_t, 10> Seeds{1001,1002,1003,1004,1005,1006,1007,1008,1009,1010};
constexpr std::array<double, 3> TargetTruth{1.4, 0.72, 0.12};
constexpr std::array<double, 3> HaloTruth{0.7, 0.9, 0.18};
constexpr std::array<double, 2> FirstStageHalo{0.58, 0.82};
constexpr double TargetInitialWidth = 0.81;
constexpr Position HaloPosition{0.0, 0.0, 2.0};

struct Scenario
{
    Input input;
    std::vector<Position> positions;
    std::vector<double> target_truth, halo_truth, clean_observations;
};

struct MethodRun
{
    std::uint64_t seed{};
    bool candidate_numerically_passed{}, production_formally_qualified{};
    std::array<double, 3> error{};
    std::array<double, 3> estimate{};
    std::string point_quality{"Unavailable"}, parameter_stability{"NotAssessed"}, local_conditioning{"NotAssessed"};
    double target_prediction_rmse{}, full_prediction_rmse{}, retained_halo_rmse{};
    double objective{};
    std::vector<core::JointCheck> checks;
};

json::value Number(double value)
{ return std::isfinite(value) ? json::value(value) : json::value(nullptr); }

void Write(const std::filesystem::path & path, const json::value & value)
{
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    const auto temporary = std::filesystem::path(path.string() + ".tmp");
    std::ofstream output(temporary);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << json::serialize(value) << '\n';
    output.close();
    std::filesystem::rename(temporary, path);
}

joint::BasisValues SimulationBasis(double square, double width, double support_radius)
{
    if (square > support_radius * support_radius) return {};
    const double distance = std::sqrt(square);
    auto value = joint::EvaluateKernel(square, width, std::max(2.5, support_radius));
    value.charge = distance < 1e-5 ? std::sqrt(2.0 / std::numbers::pi) / width :
        std::erf(distance / width / std::sqrt(2.0)) / distance;
    return value;
}

Scenario BuildScenario(double simulation_radius)
{
    Scenario scenario;
    scenario.input.atom_ids = {"target", "halo"};
    scenario.input.support.resize(2);
    scenario.input.selection_domain.emplace();
    scenario.input.selection_domain->target_indices = {0};
    scenario.input.selection_domain->contributor_policy = std::string(rhbm_gem::FixedHaloContributorPolicy);
    scenario.input.selection_domain->fixed_halo_selection = rhbm_gem::JointFixedHaloSelectionProvenance{1,1,1,0};

    for (int z = -5; z <= 5; ++z)
        for (int y = -5; y <= 5; ++y)
            for (int x = -5; x <= 5; ++x)
            {
                const Position position{0.5 * x, 0.5 * y, 0.5 * z};
                const double target_square = position[0] * position[0] + position[1] * position[1] + position[2] * position[2];
                if (target_square > 6.25) continue;
                const auto row = scenario.positions.size();
                scenario.positions.push_back(position);
                scenario.input.row_ids.push_back(std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z));
                scenario.input.support[0].push_back({row, target_square});
                const double dx = position[0] - HaloPosition[0];
                const double dy = position[1] - HaloPosition[1];
                const double dz = position[2] - HaloPosition[2];
                const double halo_square = dx * dx + dy * dy + dz * dz;
                if (halo_square <= 6.25) scenario.input.support[1].push_back({row, halo_square});

                const auto target = SimulationBasis(target_square, TargetTruth[1], 2.5);
                const auto halo = SimulationBasis(halo_square, HaloTruth[1], simulation_radius);
                scenario.target_truth.push_back(TargetTruth[0] * target.gaussian + TargetTruth[2] * target.charge);
                scenario.halo_truth.push_back(HaloTruth[0] * halo.gaussian + HaloTruth[2] * halo.charge);
                scenario.input.observations.push_back(scenario.target_truth.back() + scenario.halo_truth.back());
            }
    if (scenario.input.support[1].empty()) throw std::logic_error("Synthetic Halo has no Joint support.");
    scenario.clean_observations = scenario.input.observations;
    return scenario;
}

double Rmse(const std::vector<double> & errors)
{
    if (errors.empty()) return std::numeric_limits<double>::quiet_NaN();
    long double sum{};
    for (const double value : errors) sum += static_cast<long double>(value) * value;
    return std::sqrt(static_cast<double>(sum / static_cast<long double>(errors.size())));
}

bool HasPassed(const std::vector<core::JointCheck> & checks, std::string_view name)
{
    return std::any_of(checks.begin(), checks.end(), [&](const auto & check) {
        return check.name == name && check.status == rhbm_gem::JointCheckStatus::Passed;
    });
}

MethodRun Fit(const Scenario & scenario, std::span<const std::uint8_t> retained,
    std::string_view method, std::uint64_t seed)
{
    const std::size_t row_count = scenario.input.observations.size();
    Input noisy = scenario.input;
    core::JointProblem original(noisy);
    core::JointProblem fit_problem = original;
    if (method == "A")
    {
        rhbm_gem::JointFixedHaloParameter fixed;
        fixed.contributor_index = 1;
        fixed.atom_id = "halo";
        fixed.amplitude = FirstStageHalo[0];
        fixed.width = FirstStageHalo[1];
        fixed.charge = 0.0;
        fixed.first_stage_sample_count = 20;
        fixed.first_stage_status = static_cast<int>(rhbm_gem::RHBMEstimationStatus::SUCCESS);
        fixed.first_stage_fit_qualified = true;
        fit_problem = original.WithFixedHalos({fixed});
    }
    else
    {
        auto masked = audit::BuildMaskedInput(noisy, retained, false);
        if (method == "D")
            for (std::size_t row = 0; row < row_count; ++row)
                masked.input.observations.at(row) -= scenario.halo_truth[row];
        fit_problem = core::JointProblem(std::move(masked.input));
        fit_problem = core::JointProblemAccess::WithReferenceNormalization(fit_problem, original);
    }
    const auto fit = core::FitJointComponents(fit_problem, {TargetInitialWidth});
    MethodRun result;
    result.seed = seed;
    result.checks = fit.evidence;
    result.objective = fit.objective.value_or(std::numeric_limits<double>::quiet_NaN());
    if (!fit.assembled_state || fit.assembled_state->ac.size() != 2 || fit.assembled_state->b.size() != 1)
        return result;
    result.estimate = {fit.assembled_state->ac[0], fit.assembled_state->b[0], fit.assembled_state->ac[1]};
    for (std::size_t parameter = 0; parameter < result.error.size(); ++parameter)
        result.error[parameter] = result.estimate[parameter] - TargetTruth[parameter];
    bool parameter_valid{};
    if (fit.target_evidence && !fit.target_evidence->point_atoms.empty())
    {
        const auto & point = fit.target_evidence->point_atoms.front();
        result.point_quality = point.point_quality;
        result.parameter_stability = point.parameter_stability;
        result.local_conditioning = point.local_conditioning;
        parameter_valid = point.parameter_validity == rhbm_gem::JointCheckStatus::Passed &&
            point.local_conditioning == "Adequate";
    }

    const bool replay = HasPassed(fit.evidence, "independent-replay");
    result.candidate_numerically_passed = fit.RuntimeConvergence() == rhbm_gem::JointCheckStatus::Passed &&
        fit.TargetRuntimeConvergence() == rhbm_gem::JointCheckStatus::Passed && replay && parameter_valid;
    result.production_formally_qualified = result.candidate_numerically_passed &&
        result.point_quality == "Available" && result.parameter_stability == "Stable";

    std::vector<double> target_errors, full_errors, retained_halo_errors;
    const auto & original_rows = scenario.input.support[0];
    for (std::size_t row = 0; row < row_count; ++row)
    {
        double target_prediction{};
        for (const auto & support : original_rows)
            if (support.row == row)
            {
                const auto basis = joint::EvaluateKernel(support.squared_distance, result.estimate[1], 2.5);
                target_prediction += result.estimate[0] * basis.gaussian + result.estimate[2] * basis.charge;
            }
        target_errors.push_back(target_prediction - scenario.target_truth[row]);
        double halo_prediction{};
        if (method == "A")
        {
            for (const auto & support : scenario.input.support[1])
                if (support.row == row)
                    halo_prediction += FirstStageHalo[0] *
                        joint::EvaluateKernel(support.squared_distance, FirstStageHalo[1], 2.5).gaussian;
        }
        else if (method == "D") halo_prediction = scenario.halo_truth[row];
        full_errors.push_back(target_prediction + halo_prediction - scenario.clean_observations[row]);
        if (retained[row]) retained_halo_errors.push_back(scenario.halo_truth[row]);
    }
    result.target_prediction_rmse = Rmse(target_errors);
    result.full_prediction_rmse = Rmse(full_errors);
    result.retained_halo_rmse = Rmse(retained_halo_errors);
    return result;
}

struct ParameterSummary
{
    std::size_t count{};
    double bias{}, variance_population{}, variance_sample{}, mse{}, rmse{}, median_absolute{}, p95_absolute{}, maximum_absolute{};
};

double Percentile(std::vector<double> values, double probability)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const auto index = static_cast<std::size_t>(std::ceil(probability * static_cast<double>(values.size() - 1)));
    return values.at(index);
}

ParameterSummary Summarize(const std::vector<double> & errors)
{
    ParameterSummary result;
    result.count = errors.size();
    if (errors.empty())
    {
        result.bias = result.variance_population = result.variance_sample = result.mse = result.rmse =
            result.median_absolute = result.p95_absolute = result.maximum_absolute = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    const long double count = static_cast<long double>(errors.size());
    long double sum{}, square{};
    std::vector<double> absolute;
    absolute.reserve(errors.size());
    for (const double error : errors)
    { sum += error; square += static_cast<long double>(error) * error; absolute.push_back(std::abs(error)); }
    result.bias = static_cast<double>(sum / count);
    result.mse = static_cast<double>(square / count);
    long double centered{};
    for (const double error : errors) centered += (error - result.bias) * (error - result.bias);
    result.variance_population = static_cast<double>(centered / count);
    result.variance_sample = errors.size() > 1 ? static_cast<double>(centered / static_cast<long double>(errors.size() - 1)) :
        std::numeric_limits<double>::quiet_NaN();
    result.rmse = std::sqrt(result.mse);
    result.median_absolute = Percentile(absolute, 0.5);
    result.p95_absolute = Percentile(absolute, 0.95);
    result.maximum_absolute = *std::max_element(absolute.begin(), absolute.end());
    return result;
}

json::object StatsJson(const ParameterSummary & value)
{
    return {{"sample_count", value.count}, {"bias", Number(value.bias)},
        {"variance_population", Number(value.variance_population)}, {"variance_sample", Number(value.variance_sample)},
        {"mse", Number(value.mse)}, {"rmse", Number(value.rmse)},
        {"median_absolute_error", Number(value.median_absolute)}, {"p95_absolute_error", Number(value.p95_absolute)},
        {"maximum_absolute_error", Number(value.maximum_absolute)}};
}

json::object RunJson(const MethodRun & run)
{
    json::array checks;
    const auto status_name = [](rhbm_gem::JointCheckStatus status) -> const char * {
        switch (status)
        {
        case rhbm_gem::JointCheckStatus::Passed: return "Passed";
        case rhbm_gem::JointCheckStatus::Failed: return "Failed";
        case rhbm_gem::JointCheckStatus::Unavailable: return "Unavailable";
        case rhbm_gem::JointCheckStatus::NotRun: return "NotRun";
        }
        return "Unknown";
    };
    for (const auto & check : run.checks)
        if (check.name == "independent-replay" || check.name == "parameter-validity" ||
            check.name == "conditional-ac-optimality" || check.name == "width-stationarity")
            checks.emplace_back(json::object{{"name", check.name}, {"status", status_name(check.status)},
                {"value", check.value ? Number(*check.value) : json::value(nullptr)},
                {"threshold", check.threshold ? Number(*check.threshold) : json::value(nullptr)}, {"reason", check.reason}});
    return {{"seed", run.seed}, {"candidate_numerically_passed", run.candidate_numerically_passed},
        {"production_formally_qualified", run.production_formally_qualified},
        {"A_estimate", Number(run.estimate[0])}, {"B_estimate_angstrom", Number(run.estimate[1])},
        {"C_estimate", Number(run.estimate[2])}, {"A_error", Number(run.error[0])},
        {"B_error_angstrom", Number(run.error[1])}, {"C_error", Number(run.error[2])},
        {"target_prediction_rmse", Number(run.target_prediction_rmse)},
        {"full_prediction_rmse", Number(run.full_prediction_rmse)},
        {"retained_true_halo_rmse", Number(run.retained_halo_rmse)}, {"objective", Number(run.objective)},
        {"point_quality", run.point_quality}, {"parameter_stability", run.parameter_stability},
        {"local_conditioning", run.local_conditioning}, {"checks", std::move(checks)}};
}

json::object GroupJson(const std::vector<MethodRun> & runs)
{
    std::array<std::vector<double>, 3> candidate_errors, formal_errors;
    std::vector<double> target_prediction, full_prediction, retained_halo, formal_target_prediction;
    std::size_t candidate_passed{}, candidate_failed{}, formal_passed{};
    std::map<std::string, std::size_t> stability_counts, conditioning_counts, quality_counts;
    json::array by_seed;
    for (const auto & run : runs)
    {
        by_seed.push_back(RunJson(run));
        if (run.candidate_numerically_passed)
        {
            ++candidate_passed;
            for (std::size_t parameter = 0; parameter < candidate_errors.size(); ++parameter)
                candidate_errors[parameter].push_back(run.error[parameter]);
            target_prediction.push_back(run.target_prediction_rmse);
            full_prediction.push_back(run.full_prediction_rmse);
            retained_halo.push_back(run.retained_halo_rmse);
        }
        else ++candidate_failed;
        if (run.production_formally_qualified)
        {
            ++formal_passed;
            for (std::size_t parameter = 0; parameter < formal_errors.size(); ++parameter)
                formal_errors[parameter].push_back(run.error[parameter]);
            formal_target_prediction.push_back(run.target_prediction_rmse);
        }
        ++stability_counts[run.parameter_stability];
        ++conditioning_counts[run.local_conditioning];
        ++quality_counts[run.point_quality];
    }
    json::array candidate_quality, formal_quality;
    for (const auto & values : candidate_errors) candidate_quality.push_back(StatsJson(Summarize(values)));
    for (const auto & values : formal_errors) formal_quality.push_back(StatsJson(Summarize(values)));
    json::object stability, conditioning, point_quality;
    for (const auto & [name, count] : stability_counts) stability[name] = count;
    for (const auto & [name, count] : conditioning_counts) conditioning[name] = count;
    for (const auto & [name, count] : quality_counts) point_quality[name] = count;
    return {{"attempted_repetitions", runs.size()}, {"candidate_numerically_passed", candidate_passed},
        {"candidate_numerically_failed", candidate_failed}, {"production_formally_qualified", formal_passed},
        {"parameter_order", json::array{"A", "B_angstrom", "C"}},
        {"candidate_diagnostic_parameter_bias_variance_mse", std::move(candidate_quality)},
        {"production_qualified_parameter_bias_variance_mse", std::move(formal_quality)},
        {"candidate_diagnostic_target_prediction_rmse", StatsJson(Summarize(target_prediction))},
        {"candidate_diagnostic_full_prediction_rmse", StatsJson(Summarize(full_prediction))},
        {"candidate_diagnostic_retained_true_halo_rmse", StatsJson(Summarize(retained_halo))},
        {"production_qualified_target_prediction_rmse", StatsJson(Summarize(formal_target_prediction))},
        {"parameter_stability_counts", std::move(stability)}, {"local_conditioning_counts", std::move(conditioning)},
        {"point_quality_counts", std::move(point_quality)}, {"repetitions", std::move(by_seed)}};
}

json::object RunCondition(const Scenario & scenario, double sigma)
{
    std::array<std::vector<MethodRun>, 4> runs;
    const std::array<std::string, 4> methods{"A", "B", "C", "D"};
    const auto exact = audit::HaloSupportExclusionRows(scenario.input);
    const std::array<Position, 1> halos{HaloPosition};
    const auto distances = audit::NearestHaloDistances(scenario.positions, halos);
    const auto geometry = audit::RetainedRows(distances, 2.0);
    const audit::RowMask all(scenario.input.observations.size(), 1);
    json::array noise_realizations;
    for (const auto seed : Seeds)
    {
        Input noisy = scenario.input;
        std::vector<double> observations(noisy.observations.size());
        const auto noise = audit::FillGaussianNoise(scenario.input.observations, observations, sigma, seed);
        noisy.observations = std::move(observations);
        Scenario paired = scenario;
        paired.input.observations = noisy.observations;
        noise_realizations.emplace_back(json::object{{"seed", seed}, {"realized_mean", noise.mean},
            {"realized_standard_deviation", noise.standard_deviation}});
        const std::array<std::span<const std::uint8_t>, 4> masks{all, exact, geometry, all};
        for (std::size_t method = 0; method < methods.size(); ++method)
        {
            const auto result = Fit(paired, masks[method], methods[method], seed);
            runs[method].push_back(result);
        }
        (void)noise;
    }
    json::object groups;
    for (std::size_t method = 0; method < methods.size(); ++method)
        groups[methods[method]] = GroupJson(runs[method]);
    const auto exact_summary = audit::AnalyzeCoverage(scenario.input, exact,
        std::array<double,1>{2.0}, distances);
    const auto geometry_summary = audit::AnalyzeCoverage(scenario.input, geometry,
        std::array<double,1>{2.0}, distances);
    const double original_row_count = static_cast<double>(scenario.input.observations.size());
    return {{"noise_sigma_normalized_map_units", sigma}, {"paired_seeds", json::value_from(Seeds)},
        {"noise_realizations", std::move(noise_realizations)},
        {"coverage", json::object{{"original_rows", scenario.input.observations.size()},
            {"exact_retained_rows", exact_summary.retained_rows},
            {"exact_retained_fraction", static_cast<double>(exact_summary.retained_rows) / original_row_count},
            {"geometry_r2_retained_rows", geometry_summary.retained_rows},
            {"geometry_r2_retained_fraction", static_cast<double>(geometry_summary.retained_rows) / original_row_count},
            {"target_components_after_exact_mask", exact_summary.target_components}}},
        {"groups", std::move(groups)}};
}

std::size_t PeakRssBytes()
{
#ifdef _WIN32
    return 0;
#else
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss < 0) return 0;
#ifdef __APPLE__
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
#endif
#endif
}

void Run(const std::filesystem::path & output_path)
{
    json::array conditions;
    for (const double simulation_radius : SimulationRadii)
    {
        const auto scenario = BuildScenario(simulation_radius);
        for (const double sigma : NoiseSigmas)
            conditions.emplace_back(json::object{{"simulation_halo_support_radius_angstrom", simulation_radius},
                {"joint_support_radius_angstrom", 2.5}, {"condition", RunCondition(scenario, sigma)}});
    }
    json::object threads{{"Eigen", Eigen::nbThreads()},
        {"OpenMP", std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset"},
        {"OpenBLAS", std::getenv("OPENBLAS_NUM_THREADS") ? std::getenv("OPENBLAS_NUM_THREADS") : "unset"},
        {"Accelerate", std::getenv("VECLIB_MAXIMUM_THREADS") ? std::getenv("VECLIB_MAXIMUM_THREADS") : "unset"}};
    Write(output_path, json::object{{"contract", "joint-halo-support-exclusion-synthetic-robustness-v1"},
        {"target_count", 1}, {"halo_count", 1}, {"joint_support_radius_angstrom", 2.5},
        {"simulation_support_regimes_angstrom", json::value_from(SimulationRadii)},
        {"noise_sigmas_normalized_map_units", json::value_from(NoiseSigmas)},
        {"target_truth", json::array{TargetTruth[0], TargetTruth[1], TargetTruth[2]}},
        {"halo_truth", json::array{HaloTruth[0], HaloTruth[1], HaloTruth[2]}},
        {"first_stage_halo_prior", json::object{{"A", FirstStageHalo[0]}, {"B_angstrom", FirstStageHalo[1]},
            {"C", 0.0}, {"qualification", "deterministic synthetic qualified prior; not Oracle truth"}}},
        {"target_initial_width", TargetInitialWidth},
        {"noise_layer", "Joint observation rows; identical realization supplied to A/B/C/D in each paired seed"},
        {"support_mismatch_contract", "Target source support remains 2.5 A; Halo Gaussian and charge sources use the stated radius; Joint Halo support remains 2.5 A"},
        {"groups", json::array{"A-production-fixed-halo", "B-exact-halo-support-exclusion", "C-geometry-exclusion-r2", "D-oracle-fullabc-control"}},
        {"threads", std::move(threads)}, {"peak_rss_bytes", PeakRssBytes()}, {"conditions", std::move(conditions)}});
}
}

int main(int argc, char ** argv)
{
    try
    {
        Eigen::setNbThreads(1);
        if (argc != 2) throw std::invalid_argument("Usage: joint_halo_support_exclusion_synthetic OUTPUT.json");
        Run(argv[1]);
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
