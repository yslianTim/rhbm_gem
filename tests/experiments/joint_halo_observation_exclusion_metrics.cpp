#include "support/JointHaloObservationExclusion.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {
namespace json = boost::json;
namespace core = rhbm_gem::core;
namespace audit = joint_test::halo_exclusion;
namespace joint = core::joint_component;
using Path = std::filesystem::path;

struct TruthAtom
{
    int element{};
    double charge{}, width{};
};

struct Statistics
{
    std::size_t count{};
    double bias{}, mae{}, rmse{}, median_absolute_error{}, p95_absolute_error{}, maximum_absolute_error{};
};

json::value Number(double value)
{ return std::isfinite(value) ? json::value(value) : json::value(nullptr); }

json::value ReadJson(const Path & path)
{
    std::ifstream input(path);
    input.exceptions(std::ios::failbit | std::ios::badbit);
    std::ostringstream contents;
    contents << input.rdbuf();
    return json::parse(contents.str());
}

void Write(const Path & path, const json::value & value)
{
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    const Path temporary(path.string() + ".tmp");
    std::ofstream output(temporary);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << json::serialize(value) << '\n';
    output.close();
    std::filesystem::rename(temporary, path);
}

double Percentile(std::vector<double> values, double probability)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const auto index = static_cast<std::size_t>(std::ceil(probability * static_cast<double>(values.size() - 1)));
    return values.at(index);
}

Statistics Summarize(const std::vector<double> & errors)
{
    Statistics result;
    result.count = errors.size();
    if (errors.empty())
    {
        result.bias = result.mae = result.rmse = result.median_absolute_error =
            result.p95_absolute_error = result.maximum_absolute_error = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    std::vector<double> absolute;
    absolute.reserve(errors.size());
    long double sum{}, sum_square{}, absolute_sum{};
    for (const double error : errors)
    {
        const double magnitude = std::abs(error);
        sum += error;
        sum_square += static_cast<long double>(error) * error;
        absolute_sum += magnitude;
        absolute.push_back(magnitude);
    }
    const long double denominator = static_cast<long double>(errors.size());
    result.bias = static_cast<double>(sum / denominator);
    result.mae = static_cast<double>(absolute_sum / denominator);
    result.rmse = std::sqrt(static_cast<double>(sum_square / denominator));
    result.median_absolute_error = Percentile(absolute, 0.5);
    result.p95_absolute_error = Percentile(absolute, 0.95);
    result.maximum_absolute_error = *std::max_element(absolute.begin(), absolute.end());
    return result;
}

json::object StatisticsJson(const Statistics & statistics)
{
    return {{"count", statistics.count}, {"bias", Number(statistics.bias)}, {"mae", Number(statistics.mae)},
        {"rmse", Number(statistics.rmse)}, {"median_absolute_error", Number(statistics.median_absolute_error)},
        {"p95_absolute_error", Number(statistics.p95_absolute_error)},
        {"maximum_absolute_error", Number(statistics.maximum_absolute_error)}};
}

std::uint64_t RowIdentityHash(const std::vector<std::string> & row_ids)
{
    std::uint64_t hash = 14695981039346656037ull;
    for (const auto & id : row_ids)
    {
        for (const char character : id)
        { hash ^= static_cast<unsigned char>(character); hash *= 1099511628211ull; }
        hash ^= '\n'; hash *= 1099511628211ull;
    }
    return hash;
}

std::string Hex(std::uint64_t value)
{
    std::ostringstream result;
    result << std::hex << std::setw(16) << std::setfill('0') << value;
    return result.str();
}

const char * DistanceStratum(double distance)
{
    if (distance < 2.5) return "<2.5A";
    if (distance < 5.0) return "2.5-5A";
    return ">=5A";
}

json::object ErrorStatistics(const std::vector<double> & amplitude, const std::vector<double> & width,
    const std::vector<double> & charge)
{
    return {{"A_raw_units", StatisticsJson(Summarize(amplitude))}, {"B_angstrom", StatisticsJson(Summarize(width))},
        {"C_raw_units", StatisticsJson(Summarize(charge))}};
}

void Run(const Path & model_path, const Path & map_path, const Path & sidecar_path,
    const Path & fits_path, const Path & output_path)
{
    const auto sidecar = ReadJson(sidecar_path).as_object();
    const auto fits = ReadJson(fits_path).as_object();
    auto model = rhbm_gem::ReadModel(model_path);
    auto map = rhbm_gem::ReadMap(map_path);
    const auto map_divisor = map->MapValueArrayNormalization().value_or(1.0);
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement() != Element::HYDROGEN && atom.GetChainID() == "A";
    });
    const auto problem = core::BuildJointProblem(*map, *model);
    const auto & input = problem.Input();
    const auto initialization = joint::InitializeContributors(*map, *model, problem);
    const auto & target_indices = input.selection_domain->target_indices;
    std::unordered_map<std::string, TruthAtom> truth;
    for (const auto & value : sidecar.at("atoms").as_array())
    {
        const auto & atom = value.as_object();
        if (atom.at("element").as_int64() == static_cast<int>(Element::HYDROGEN)) continue;
        const auto & width = atom.at("effective_gaussian_width");
        if (width.is_null()) continue;
        truth.emplace(std::to_string(atom.at("serial_id").as_int64()), TruthAtom{
            static_cast<int>(atom.at("element").as_int64()), atom.at("charge_used").as_double(), width.as_double()});
    }

    const std::size_t row_count = input.observations.size();
    std::vector<double> true_target(row_count), true_halo(row_count);
    std::vector<audit::Position> row_positions, target_positions, halo_positions;
    row_positions.reserve(input.row_ids.size());
    for (const auto & id : input.row_ids) row_positions.push_back(map->GetGridPosition(std::stoull(id)));
    std::vector<double> target_halo_distances;
    target_halo_distances.reserve(target_indices.size());
    for (std::size_t contributor = 0; contributor < input.atom_ids.size(); ++contributor)
    {
        const auto * atom = model->FindAtomPtr(std::stoi(input.atom_ids[contributor]));
        if (!atom) throw std::runtime_error("Contributor is missing from the simulation model.");
        const auto position = atom->GetPosition();
        const bool target = std::binary_search(target_indices.begin(), target_indices.end(), contributor);
        if (target) target_positions.push_back(position);
        else halo_positions.push_back(position);

        const auto found = truth.find(input.atom_ids[contributor]);
        if (found == truth.end()) throw std::runtime_error("Contributor is missing from the truth sidecar.");
        const auto & parameters = found->second;
        const double amplitude = static_cast<double>(parameters.element) / map_divisor;
        const double charge = parameters.charge / map_divisor;
        for (const auto & support : input.support[contributor])
        {
            const auto basis = joint::EvaluateKernel(support.squared_distance, parameters.width, 2.5);
            const double value = amplitude * basis.gaussian + charge * basis.charge;
            (target ? true_target : true_halo)[support.row] += value;
        }
    }
    target_halo_distances = audit::NearestHaloDistances(target_positions, halo_positions);
    const auto nearest_rows = audit::NearestHaloDistances(row_positions, halo_positions);
    std::vector<double> truth_reconstruction_error(row_count);
    for (std::size_t row = 0; row < row_count; ++row)
        truth_reconstruction_error[row] = true_target[row] + true_halo[row] - input.observations[row];

    json::array runs;
    for (const auto & run_value : fits.at("runs").as_array())
    {
        const auto & run = run_value.as_object();
        const std::string method = json::value_to<std::string>(run.at("method"));
        const double radius = run.at("radius_angstrom").as_double();
        const auto retained = method == "A" ? audit::RowMask(row_count, 1) :
            method == "B" && radius == 2.5 ? audit::HaloSupportExclusionRows(input) :
            audit::RetainedRows(nearest_rows, radius);
        std::vector<std::string> retained_ids;
        retained_ids.reserve(row_count);
        for (std::size_t row = 0; row < row_count; ++row) if (retained[row]) retained_ids.push_back(input.row_ids[row]);
        if (static_cast<std::size_t>(run.at("observation_count").as_int64()) != retained_ids.size() ||
            json::value_to<std::string>(run.at("row_ids_fnv1a64")) != Hex(RowIdentityHash(retained_ids)))
            throw std::runtime_error("Controlled-fit row identities do not match the reconstructed geometry mask.");

        std::vector<double> estimated_target(row_count);
        std::vector<double> run_first_stage_halo(row_count);
        std::size_t fixed_halo_count{};
        std::map<std::string, std::array<double, 2>> run_fixed_halos;
        if (const auto * fixed_values = run.if_contains("fixed_halos"))
            for (const auto & fixed_value : fixed_values->as_array())
            {
                const auto & fixed = fixed_value.as_object();
                if (!fixed.at("first_stage_fit_qualified").as_bool() ||
                    fixed.at("charge").as_double() != 0.0)
                    throw std::runtime_error("Controlled fit contains an unqualified fixed Halo.");
                run_fixed_halos.emplace(json::value_to<std::string>(fixed.at("atom_id")),
                    std::array<double, 2>{fixed.at("amplitude").as_double(), fixed.at("width").as_double()});
            }
        if (method == "A" || method == "C")
            for (std::size_t contributor = 0; contributor < input.atom_ids.size(); ++contributor)
            {
                if (std::binary_search(target_indices.begin(), target_indices.end(), contributor)) continue;
                const auto & support = input.support[contributor];
                const bool retained_support = method == "A" ||
                    std::any_of(support.begin(), support.end(), [&](const auto & entry) { return retained[entry.row]; });
                if (!retained_support) continue;
                ++fixed_halo_count;
                const auto found = run_fixed_halos.find(input.atom_ids[contributor]);
                const double amplitude = found != run_fixed_halos.end() ? found->second[0] :
                    initialization.atoms.at(contributor).mdpde[0];
                const double width = found != run_fixed_halos.end() ? found->second[1] :
                    initialization.atoms.at(contributor).mdpde[1];
                for (const auto & entry : support)
                    run_first_stage_halo[entry.row] += amplitude *
                        joint::EvaluateKernel(entry.squared_distance, width, 2.5).gaussian;
            }
        std::map<std::string, std::array<double, 3>> estimates;
        std::map<std::string, std::string> point_quality, parameter_validity, local_conditioning, parameter_stability;
        for (const auto & target_value : run.at("targets").as_array())
        {
            const auto & target = target_value.as_object();
            const std::string id = json::value_to<std::string>(target.at("atom_id"));
            if (!target.at("amplitude").is_null() && !target.at("charge").is_null() && !target.at("width").is_null())
                estimates[id] = {target.at("amplitude").as_double(), target.at("charge").as_double(), target.at("width").as_double()};
            point_quality[id] = json::value_to<std::string>(target.at("point_quality"));
            parameter_validity[id] = json::value_to<std::string>(target.at("parameter_validity"));
            local_conditioning[id] = json::value_to<std::string>(target.at("local_conditioning"));
            parameter_stability[id] = json::value_to<std::string>(target.at("parameter_stability"));
        }
        const auto count_statuses = [](const auto & statuses) {
            std::map<std::string, std::size_t> counts;
            for (const auto & [id, status] : statuses) { (void)id; ++counts[status]; }
            json::object result;
            for (const auto & [status, count] : counts) result[status] = count;
            return result;
        };

        std::vector<double> amplitude_errors, width_errors, charge_errors;
        std::vector<double> qualified_amplitude_errors, qualified_width_errors, qualified_charge_errors;
        std::map<std::string, std::vector<double>> strata_a, strata_b, strata_c;
        json::array target_records;
        std::size_t targets_without_estimate{}, qualified_target_estimates{};
        bool complete_target_estimate_set = true;
        for (std::size_t target = 0; target < target_indices.size(); ++target)
        {
            const auto contributor = target_indices[target];
            const auto & id = input.atom_ids.at(contributor);
            const auto & reference = truth.at(id);
            const double distance = target_halo_distances.at(target);
            const std::string stratum = DistanceStratum(distance);
            const auto found = estimates.find(id);
            complete_target_estimate_set &= found != estimates.end();
            if (found == estimates.end()) ++targets_without_estimate;
            const bool qualified = found != estimates.end() && point_quality.at(id) == "Available" &&
                parameter_validity.at(id) == "Passed" && local_conditioning.at(id) == "Adequate" &&
                parameter_stability.at(id) == "Stable";
            if (qualified) ++qualified_target_estimates;
            json::object record{{"atom_id", id}, {"nearest_halo_distance", Number(distance)},
                {"distance_stratum", stratum}, {"point_quality", point_quality.at(id)},
                {"parameter_validity", parameter_validity.at(id)},
                {"local_conditioning", local_conditioning.at(id)},
                {"parameter_stability", parameter_stability.at(id)},
                {"qualified_estimate", qualified},
                {"A_truth_raw_units", reference.element},
                {"B_truth_angstrom", reference.width},
                {"C_truth_raw_units", reference.charge},
                {"A_estimate_raw_units", nullptr}, {"B_estimate_angstrom", nullptr},
                {"C_estimate_raw_units", nullptr}, {"A_error_raw_units", nullptr},
                {"B_error_angstrom", nullptr}, {"C_error_raw_units", nullptr}};
            if (found != estimates.end())
            {
                const auto & estimate = found->second;
                const double error_a = estimate[0] * map_divisor - reference.element;
                const double error_c = estimate[1] * map_divisor - reference.charge;
                const double error_b = estimate[2] - reference.width;
                amplitude_errors.push_back(error_a); width_errors.push_back(error_b); charge_errors.push_back(error_c);
                strata_a[stratum].push_back(error_a); strata_b[stratum].push_back(error_b); strata_c[stratum].push_back(error_c);
                if (qualified)
                {
                    qualified_amplitude_errors.push_back(error_a);
                    qualified_width_errors.push_back(error_b);
                    qualified_charge_errors.push_back(error_c);
                }
                record["A_estimate_raw_units"] = estimate[0] * map_divisor;
                record["B_estimate_angstrom"] = estimate[2];
                record["C_estimate_raw_units"] = estimate[1] * map_divisor;
                record["A_error_raw_units"] = error_a;
                record["B_error_angstrom"] = error_b;
                record["C_error_raw_units"] = error_c;
            }
            target_records.emplace_back(std::move(record));
        }
        json::object strata;
        for (const auto & [name, ignored] : strata_a)
        {
            (void)ignored;
            strata[name] = ErrorStatistics(strata_a[name], strata_b[name], strata_c[name]);
        }

        std::vector<double> actual_prediction_error, oracle_prediction_error,
            target_prediction_error, retained_prediction_error, excluded_prediction_error,
            retained_target_error, excluded_target_error;
        for (std::size_t row = 0; row < row_count; ++row)
        {
            if (!complete_target_estimate_set) break;
            const double halo = run_first_stage_halo[row];
            const double actual_residual = estimated_target[row] + halo - input.observations[row];
            const double oracle_residual = estimated_target[row] + true_halo[row] - input.observations[row];
            const double target_error = estimated_target[row] - true_target[row];
            actual_prediction_error.push_back(actual_residual);
            oracle_prediction_error.push_back(oracle_residual);
            target_prediction_error.push_back(target_error);
            if (retained[row])
            {
                retained_prediction_error.push_back(actual_residual);
                retained_target_error.push_back(target_error);
            }
            else
            {
                excluded_prediction_error.push_back(actual_residual);
                excluded_target_error.push_back(target_error);
            }
        }
        const auto quality = run.at("runtime_convergence");
        const auto target_quality = run.at("target_runtime_convergence");
        runs.emplace_back(json::object{{"method", method}, {"radius_angstrom", radius},
            {"runtime_convergence", quality}, {"target_runtime_convergence", target_quality},
            {"elapsed_seconds", run.at("elapsed_seconds")}, {"peak_rss_bytes", run.at("peak_rss_bytes")},
            {"observation_count", retained_ids.size()}, {"excluded_count", row_count - retained_ids.size()},
            {"target_count", target_indices.size()}, {"fixed_halo_count", fixed_halo_count},
            {"targets_without_finite_estimate", targets_without_estimate},
            {"qualified_target_estimate_count", qualified_target_estimates},
            {"prediction_diagnostics_available", complete_target_estimate_set},
            {"point_quality_counts", count_statuses(point_quality)},
            {"parameter_validity_counts", count_statuses(parameter_validity)},
            {"local_conditioning_counts", count_statuses(local_conditioning)},
            {"parameter_stability_counts", count_statuses(parameter_stability)},
            {"finite_candidate_parameter_errors", ErrorStatistics(amplitude_errors, width_errors, charge_errors)},
            {"qualified_parameter_errors", ErrorStatistics(qualified_amplitude_errors,
                qualified_width_errors, qualified_charge_errors)},
            {"target_halo_distance_strata", std::move(strata)}, {"per_target_errors", std::move(target_records)},
            {"prediction_errors_normalized_map_units", json::object{
                {"target_only_original", StatisticsJson(Summarize(target_prediction_error))},
                {"actual_model_original", StatisticsJson(Summarize(actual_prediction_error))},
                {"actual_model_retained", StatisticsJson(Summarize(retained_prediction_error))},
                {"actual_model_excluded", StatisticsJson(Summarize(excluded_prediction_error))},
                {"oracle_halo_assisted_original", StatisticsJson(Summarize(oracle_prediction_error))},
                {"target_only_retained", StatisticsJson(Summarize(retained_target_error))},
                {"target_only_excluded", StatisticsJson(Summarize(excluded_target_error))}}}});
    }

    const double truth_rmse = Summarize(truth_reconstruction_error).rmse;
    Write(output_path, json::object{{"contract", "joint-halo-observation-exclusion-scientific-validation-v1"},
        {"model", model_path.string()}, {"map", map_path.string()}, {"simulation_truth", sidecar_path.string()},
        {"model_sha256", json::value_to<std::string>(sidecar.at("source").as_object().at("model_sha256"))},
        {"map_sha256", json::value_to<std::string>(sidecar.at("output").as_object().at("map_sha256"))},
        {"map_normalization_divisor", map_divisor}, {"original_observation_count", row_count},
        {"original_observation_scale", problem.ObservationScale()}, {"fit_noise", fits.at("noise")},
        {"target_count", target_indices.size()},
        {"supported_halo_count", input.atom_ids.size() - target_indices.size()},
        {"sidecar_truth_reconstruction_rmse_normalized_map_units", truth_rmse},
        {"distance_strata_angstrom", json::array{"<2.5", "2.5-5.0", ">=5.0"}},
        {"runs", std::move(runs)}});
}
}

int main(int argc, char ** argv)
{
    try
    {
        Eigen::setNbThreads(1);
        if (argc != 6) throw std::invalid_argument(
            "Usage: joint_halo_observation_exclusion_metrics MODEL MAP SIDECAR FITS.json OUTPUT.json");
        Run(argv[1], argv[2], argv[3], argv[4], argv[5]);
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
