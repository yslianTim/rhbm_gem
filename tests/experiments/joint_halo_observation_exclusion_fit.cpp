#include "support/JointHaloObservationExclusion.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <boost/hash2/sha2.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace json = boost::json;
namespace core = rhbm_gem::core;
namespace audit = joint_test::halo_exclusion;
namespace joint = core::joint_component;
using Path = std::filesystem::path;
using Clock = std::chrono::steady_clock;
constexpr std::array<double, 3> ExperimentRadii{1.0, 2.0, 2.5};

struct NoiseSettings
{
    double sigma{};
    std::uint64_t seed{};
};

json::value Number(double value)
{ return std::isfinite(value) ? json::value(value) : json::value(nullptr); }

const char * StatusName(rhbm_gem::JointCheckStatus status)
{
    switch (status)
    {
    case rhbm_gem::JointCheckStatus::Passed: return "Passed";
    case rhbm_gem::JointCheckStatus::Failed: return "Failed";
    case rhbm_gem::JointCheckStatus::Unavailable: return "Unavailable";
    case rhbm_gem::JointCheckStatus::NotRun: return "NotRun";
    }
    return "Unknown";
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

struct AppliedNoise
{
    std::string digest;
    audit::NoiseSummary summary;
};

AppliedNoise AddNormalizedMapNoise(rhbm_gem::MapObject & map, const NoiseSettings & settings)
{
    if (!std::isfinite(settings.sigma) || settings.sigma < 0.0)
        throw std::invalid_argument("Noise sigma must be finite and nonnegative.");
    boost::hash2::sha2_256 hash;
    constexpr std::string_view tag = "normalized-full-map-voxels-gaussian-noise-v1";
    hash.update(tag.data(), tag.size());
    for (unsigned shift = 0; shift < 64; shift += 8)
    {
        const auto byte = static_cast<unsigned char>((settings.seed >> shift) & 0xff);
        hash.update(&byte, 1);
    }
    const std::size_t voxel_count = map.GetMapValueArraySize();
    auto array = std::make_unique<double[]>(voxel_count);
    const auto summary = audit::FillGaussianNoise(
        std::span<const double>(map.GetMapValueArray(), voxel_count),
        std::span<double>(array.get(), voxel_count), settings.sigma, settings.seed);
    for (std::size_t index = 0; index < voxel_count; ++index)
    {
        const auto bits = std::bit_cast<std::uint64_t>(array[index]);
        std::array<unsigned char, 8> bytes{};
        for (unsigned byte = 0; byte < bytes.size(); ++byte)
            bytes[byte] = static_cast<unsigned char>((bits >> (8 * byte)) & 0xff);
        hash.update(bytes.data(), bytes.size());
    }
    map.SetMapValueArray(std::move(array));
    return {boost::hash2::to_string(hash.result()), summary};
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

std::vector<double> TargetWidths(const core::JointProblemInput & input,
    const rhbm_gem::JointInitialization & initialization)
{
    std::vector<double> widths;
    widths.reserve(input.selection_domain->target_indices.size());
    for (const auto contributor : input.selection_domain->target_indices)
    {
        const double width = initialization.b.at(contributor);
        if (!std::isfinite(width) || !(width > 0.0))
            throw std::runtime_error("A selected target has no finite positive FirstStage width.");
        widths.push_back(width);
    }
    return widths;
}

json::object FitJson(const std::string & method, double radius, const core::JointProblem & problem,
    const rhbm_gem::JointInitialization & initialization,
    const std::vector<std::size_t> & source_contributors, const core::JointFitResult & fit,
    double elapsed_seconds)
{
    const auto & input = problem.Input();
    std::size_t sweeps{};
    json::array components, evidence, targets, fixed_halos;
    for (const auto & component : fit.components)
    {
        sweeps += component.search_sweeps;
        json::array conditioning, stability;
        std::map<std::string, std::size_t> conditioning_counts, stability_counts;
        for (const auto & block : component.local_conditioning)
        {
            ++conditioning_counts[block.conditioning_classification];
            conditioning.emplace_back(json::object{{"block_id", block.block_id},
                {"classification", block.conditioning_classification},
                {"free_ac_status", block.free_ac_status}, {"free_ac_rank", block.free_ac_rank ?
                    json::value(*block.free_ac_rank) : json::value(nullptr)},
                {"free_ac_condition", block.free_ac_condition ? Number(*block.free_ac_condition) : json::value(nullptr)},
                {"projected_width_status", block.projected_width_status},
                {"projected_width_rank", block.projected_width_rank ?
                    json::value(*block.projected_width_rank) : json::value(nullptr)},
                {"projected_width_condition", block.projected_width_condition ?
                    Number(*block.projected_width_condition) : json::value(nullptr)}});
        }
        json::value stability_summary = nullptr;
        if (component.parameter_stability)
        {
            const auto & result = *component.parameter_stability;
            for (const auto & block : result.blocks)
            {
                ++stability_counts[block.status];
                stability.emplace_back(json::object{{"block_id", block.block_id}, {"status", block.status},
                    {"prediction_status", block.prediction_status},
                    {"objective_consistency_status", block.objective_consistency_status},
                    {"log_width_parameter_change_inf", block.log_width_parameter_change_inf ?
                        Number(*block.log_width_parameter_change_inf) : json::value(nullptr)},
                    {"scaled_ac_parameter_change_inf", block.scaled_ac_parameter_change_inf ?
                        Number(*block.scaled_ac_parameter_change_inf) : json::value(nullptr)},
                    {"local_prediction_change_relative", block.local_prediction_change_relative ?
                        Number(*block.local_prediction_change_relative) : json::value(nullptr)}});
            }
            json::object counts;
            for (const auto & [name, count] : stability_counts) counts[name] = count;
            stability_summary = json::object{{"trajectory_stability_status", result.trajectory_stability_status},
                {"perturbation_status", result.perturbation_status},
                {"maximum_eta_change", result.maximum_eta_change ? Number(*result.maximum_eta_change) : json::value(nullptr)},
                {"maximum_scaled_ac_change", result.maximum_scaled_ac_change ?
                    Number(*result.maximum_scaled_ac_change) : json::value(nullptr)},
                {"recent_sweep_count", result.recent_sweep_count}, {"block_status_counts", std::move(counts)},
                {"perturbation_directions", result.perturbation_directions},
                {"perturbation_seeds", result.perturbation_seeds}, {"perturbation_solves", result.perturbation_solves},
                {"perturbation_evaluations", result.perturbation_evaluations}, {"seconds", result.seconds}};
        }
        json::object conditioning_summary;
        for (const auto & [name, count] : conditioning_counts) conditioning_summary[name] = count;
        components.emplace_back(json::object{{"id", component.id}, {"search_completed", component.search_completed},
            {"search_sweeps", component.search_sweeps}, {"stop_reason", component.stop_reason},
            {"profile_evaluations", component.profile_evaluations}, {"reference_evaluations", component.reference_evaluations},
            {"accepted_updates", component.accepted_updates}, {"native_status", component.native_status},
            {"runtime_convergence", StatusName(component.RuntimeConvergence())},
            {"target_runtime_convergence", StatusName(component.TargetRuntimeConvergence())},
            {"local_conditioning_seconds", component.local_conditioning_seconds},
            {"local_conditioning_status_counts", std::move(conditioning_summary)},
            {"local_conditioning_blocks", std::move(conditioning)},
            {"parameter_stability", std::move(stability_summary)},
            {"parameter_stability_blocks", std::move(stability)}});
    }
    for (const auto & check : fit.evidence)
        evidence.emplace_back(json::object{{"name", check.name}, {"status", StatusName(check.status)},
            {"scope", check.scope == rhbm_gem::JointEvidenceScope::AssembledGlobal ? "assembled-global" : "component-local"},
            {"value", check.value ? Number(*check.value) : json::value(nullptr)},
            {"threshold", check.threshold ? Number(*check.threshold) : json::value(nullptr)},
            {"reason", check.reason}});
    for (const auto & fixed : problem.FixedHalos())
        fixed_halos.emplace_back(json::object{{"atom_id", fixed.atom_id},
            {"contributor_index", fixed.contributor_index}, {"amplitude", fixed.amplitude},
            {"width", fixed.width}, {"charge", fixed.charge},
            {"first_stage_sample_count", fixed.first_stage_sample_count},
            {"first_stage_fit_qualified", fixed.first_stage_fit_qualified}});

    const auto & target_indices = input.selection_domain->target_indices;
    for (std::size_t target = 0; target < target_indices.size(); ++target)
    {
        const auto contributor = target_indices[target];
        const auto source_contributor = source_contributors.at(contributor);
        const auto & seed = initialization.atoms.at(source_contributor);
        json::value amplitude = nullptr, charge = nullptr, width = nullptr;
        if (fit.assembled_state && 2 * target + 1 < fit.assembled_state->ac.size() &&
            target < fit.assembled_state->b.size())
        {
            amplitude = Number(fit.assembled_state->ac[2 * target]);
            charge = Number(fit.assembled_state->ac[2 * target + 1]);
            width = Number(fit.assembled_state->b[target]);
        }
        json::value point_quality = "Unavailable", point_reason = "no-target-evidence";
        json::value local_conditioning = "NotAssessed", parameter_stability = "NotAssessed";
        json::value parameter_validity = "Unavailable";
        if (fit.target_evidence)
        {
            const auto point = std::find_if(fit.target_evidence->point_atoms.begin(), fit.target_evidence->point_atoms.end(),
                [&](const auto & candidate) { return candidate.atom == target; });
            if (point != fit.target_evidence->point_atoms.end())
            {
                point_quality = point->point_quality;
                point_reason = point->reason;
                local_conditioning = point->local_conditioning;
                parameter_stability = point->parameter_stability;
                parameter_validity = StatusName(point->parameter_validity);
            }
        }
        targets.emplace_back(json::object{{"atom_id", input.atom_ids.at(contributor)},
            {"seed_first_stage_qualified", seed.first_stage_fit_qualified}, {"seed_reason", seed.reason},
            {"seed_source", seed.seed_source}, {"initial_width", Number(initialization.b.at(source_contributor))},
            {"amplitude", amplitude}, {"charge", charge}, {"width", width},
            {"point_quality", std::move(point_quality)}, {"point_reason", std::move(point_reason)},
            {"parameter_validity", std::move(parameter_validity)},
            {"local_conditioning", std::move(local_conditioning)},
            {"parameter_stability", std::move(parameter_stability)}});
    }

    std::size_t retained_halos{};
    if (input.selection_domain->contributor_policy == rhbm_gem::FixedHaloContributorPolicy)
        retained_halos = input.atom_ids.size() - target_indices.size();
    std::ostringstream row_hash;
    row_hash << std::hex << std::setw(16) << std::setfill('0') << RowIdentityHash(input.row_ids);
    return {{"method", method}, {"radius_angstrom", radius}, {"observation_count", input.observations.size()},
        {"target_count", target_indices.size()}, {"halo_count", retained_halos},
        {"contributor_policy", input.selection_domain->contributor_policy},
        {"fixed_halos", std::move(fixed_halos)},
        {"row_ids_fnv1a64", row_hash.str()}, {"observation_scale", problem.ObservationScale()},
        {"elapsed_seconds", elapsed_seconds}, {"peak_rss_bytes", PeakRssBytes()},
        {"search_completed", fit.search_completed}, {"search_sweeps", sweeps},
        {"runtime_convergence", StatusName(fit.RuntimeConvergence())},
        {"target_runtime_convergence", StatusName(fit.TargetRuntimeConvergence())},
        {"objective", fit.objective ? Number(*fit.objective) : json::value(nullptr)},
        {"solver_method", fit.solver_provenance.search_method},
        {"fixed_neighbor_core_atoms", fit.solver_provenance.fixed_neighbor_core_atoms},
        {"initialization_valid", fit.initialization.valid}, {"initialization_reason", fit.initialization.reason},
        {"costs", json::object{{"initialization_seconds", fit.costs.initialization_seconds},
            {"search_seconds", fit.costs.search_seconds}, {"assessment_seconds", fit.costs.assessment_seconds},
            {"assembly_seconds", fit.costs.assembly_seconds}}},
        {"components", std::move(components)}, {"evidence", std::move(evidence)}, {"targets", std::move(targets)}};
}

json::object ExceptionJson(const std::string & method, double radius, const std::string & message)
{ return {{"method", method}, {"radius_angstrom", radius}, {"status", "exception"}, {"error", message}}; }

void Run(const Path & model_path, const Path & map_path, const Path & output_path,
    std::optional<NoiseSettings> noise = std::nullopt)
{
    auto model = rhbm_gem::ReadModel(model_path);
    auto map = rhbm_gem::ReadMap(map_path);
    const auto normalization = map->MapValueArrayNormalization();
    std::optional<std::string> noise_digest;
    std::optional<audit::NoiseSummary> noise_summary;
    if (noise)
    {
        const auto applied = AddNormalizedMapNoise(*map, *noise);
        noise_digest = applied.digest;
        noise_summary = applied.summary;
    }
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement() != Element::HYDROGEN && atom.GetChainID() == "A";
    });
    const auto original_problem = core::BuildJointProblem(*map, *model);
    const auto & original_input = original_problem.Input();
    const auto initialization = joint::InitializeContributors(*map, *model, original_problem);
    const auto widths = TargetWidths(original_input, initialization);
    std::vector<audit::Position> observations, halos;
    observations.reserve(original_input.row_ids.size());
    for (const auto & row_id : original_input.row_ids)
        observations.push_back(map->GetGridPosition(std::stoull(row_id)));
    const auto & targets = original_input.selection_domain->target_indices;
    for (std::size_t contributor = 0; contributor < original_input.atom_ids.size(); ++contributor)
    {
        if (std::binary_search(targets.begin(), targets.end(), contributor)) continue;
        const auto * atom = model->FindAtomPtr(std::stoi(original_input.atom_ids[contributor]));
        if (!atom) throw std::runtime_error("A supported Halo is missing from the input model.");
        halos.push_back(atom->GetPosition());
    }
    const auto halo_distances = audit::NearestHaloDistances(observations, halos);

    json::array runs;
    json::object root{{"contract", "joint-halo-observation-exclusion-controlled-fit-v1"},
        {"model", model_path.string()}, {"map", map_path.string()}, {"selection", "non-hydrogen chain A"},
        {"mask_selection_basis", "nearest supported Halo distance; fixed before truth evaluation"},
        {"experiment_radii_angstrom", json::value_from(ExperimentRadii)},
        {"map_normalization_divisor", normalization ? Number(*normalization) : json::value(1.0)},
        {"target_count", targets.size()}, {"halo_count", halos.size()},
        {"original_observation_count", original_input.observations.size()},
        {"original_observation_scale", original_problem.ObservationScale()},
        {"noise", noise ? json::value(json::object{{"layer", "full normalized map voxels before FirstStage and Joint construction"},
            {"distribution", "independent N(0,sigma^2)"}, {"sigma_normalized_map_units", noise->sigma},
            {"seed", noise->seed}, {"realization_sha256", *noise_digest},
            {"realized_mean", noise_summary->mean},
            {"realized_standard_deviation", noise_summary->standard_deviation},
            {"first_stage_uses_same_noisy_map", true}}) : json::value(nullptr)},
        {"first_stage_initialization_valid", initialization.valid},
        {"first_stage_initialization_reason", initialization.reason},
        {"first_stage_halos_qualified", static_cast<std::size_t>(std::count_if(initialization.atoms.begin(), initialization.atoms.end(),
            [](const auto & atom) { return !atom.target && atom.first_stage_fit_qualified; }))},
        {"thread_settings", json::object{{"Eigen", Eigen::nbThreads()},
            {"OpenMP", std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset"},
            {"OpenBLAS", std::getenv("OPENBLAS_NUM_THREADS") ? std::getenv("OPENBLAS_NUM_THREADS") : "unset"},
            {"Accelerate", std::getenv("VECLIB_MAXIMUM_THREADS") ? std::getenv("VECLIB_MAXIMUM_THREADS") : "unset"}}},
        {"runs", runs}};
    Write(output_path, root);

    auto record = [&](json::object outcome) {
        runs.emplace_back(std::move(outcome));
        root["runs"] = runs;
        Write(output_path, root);
    };
    auto execute = [&](const std::string & method, double radius,
        const std::vector<std::size_t> & source_contributors, auto && fit_factory) {
        std::cout << "Starting " << method << " R=" << radius << " A\n" << std::flush;
        const auto started = Clock::now();
        try
        {
            auto result = fit_factory();
            const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
            auto outcome = FitJson(method, radius, *result.problem, initialization, source_contributors, result, seconds);
            const std::string status = json::value_to<std::string>(outcome.at("runtime_convergence"));
            const std::string target_status = json::value_to<std::string>(outcome.at("target_runtime_convergence"));
            std::cout << "Finished " << method << " R=" << radius << " A: " << status << ", target " <<
                target_status << ", " << seconds << " s\n" << std::flush;
            record(std::move(outcome));
        }
        catch (const std::exception & error)
        {
            std::cout << "Failed " << method << " R=" << radius << " A: " << error.what() << '\n' << std::flush;
            record(ExceptionJson(method, radius, error.what()));
        }
    };

    std::vector<std::size_t> original_contributors(original_input.atom_ids.size());
    std::iota(original_contributors.begin(), original_contributors.end(), 0);
    execute("A", 0.0, original_contributors,
        [&] { return core::FitJointComponentsInitialized(original_problem, initialization); });
    if (noise)
    {
        const std::array<std::pair<std::string, audit::RowMask>, 2> comparisons{{
            {"B", audit::HaloSupportExclusionRows(original_input)},
            {"R2", audit::RetainedRows(halo_distances, 2.0)}}};
        for (const auto & [method, retained] : comparisons)
        {
            const double radius = method == "B" ? 2.5 : 2.0;
            auto targets_only = audit::BuildMaskedInput(original_input, retained, false);
            auto b_problem = core::JointProblem(std::move(targets_only.input));
            b_problem = core::JointProblemAccess::WithReferenceNormalization(b_problem, original_problem);
            execute(method, radius, targets_only.source_contributors,
                [&] { return core::FitJointComponents(b_problem, widths); });
        }
    }
    else for (const double radius : ExperimentRadii)
    {
        const auto retained = radius == 2.5 ? audit::HaloSupportExclusionRows(original_input) :
            audit::RetainedRows(halo_distances, radius);
        auto targets_only = audit::BuildMaskedInput(original_input, retained, false);
        auto b_problem = core::JointProblem(std::move(targets_only.input));
        b_problem = core::JointProblemAccess::WithReferenceNormalization(b_problem, original_problem);
        execute("B", radius, targets_only.source_contributors,
            [&] { return core::FitJointComponents(b_problem, widths); });

        auto with_halos = audit::BuildMaskedInput(original_input, retained, true);
        auto c_problem = core::JointProblem(std::move(with_halos.input));
        std::vector<rhbm_gem::JointFixedHaloParameter> fixed_halos;
        const auto & masked_input = c_problem.Input();
        for (std::size_t contributor = 0; contributor < masked_input.atom_ids.size(); ++contributor)
        {
            if (std::binary_search(masked_input.selection_domain->target_indices.begin(),
                masked_input.selection_domain->target_indices.end(), contributor)) continue;
            const auto source = with_halos.source_contributors.at(contributor);
            std::string failure;
            auto fixed = core::detail::MakeFixedHaloParameter(initialization.atoms.at(source), contributor,
                masked_input.atom_ids[contributor], failure);
            if (!fixed) throw std::runtime_error("A retained Halo lacks a qualified FirstStage snapshot: " + failure);
            fixed_halos.push_back(std::move(*fixed));
        }
        c_problem = c_problem.WithFixedHalos(std::move(fixed_halos));
        c_problem = core::JointProblemAccess::WithReferenceNormalization(c_problem, original_problem);
        execute("C", radius, with_halos.source_contributors,
            [&] { return core::FitJointComponents(c_problem, widths); });
    }
}
}

int main(int argc, char ** argv)
{
    try
    {
        Eigen::setNbThreads(1);
        if (argc != 4 && argc != 6)
            throw std::invalid_argument("Usage: joint_halo_observation_exclusion_fit MODEL MAP OUTPUT.json [NOISE_SIGMA SEED]");
        std::optional<NoiseSettings> noise;
        if (argc == 6) noise = NoiseSettings{std::stod(argv[4]), std::stoull(argv[5])};
        Run(argv[1], argv[2], argv[3], noise);
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
