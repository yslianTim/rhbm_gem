#include "support/JointHaloObservationExclusion.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace json = boost::json;
namespace core = rhbm_gem::core;
namespace audit = joint_test::halo_exclusion;
using Path = std::filesystem::path;
using Clock = std::chrono::steady_clock;

json::value Number(double value)
{ return std::isfinite(value) ? json::value(value) : json::value(nullptr); }

void Write(const Path & path, const json::value & value)
{
    const Path temporary(path.string() + ".tmp");
    std::ofstream output(temporary);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << json::serialize(value) << '\n';
    output.close();
    std::filesystem::rename(temporary, path);
}

json::object CoverageJson(double radius, const audit::CoverageSummary & summary)
{
    json::array targets, radial, proximity, observation_proximity;
    for (const auto & target : summary.targets)
    {
        json::array shells;
        for (const auto & shell : target.radial)
            shells.emplace_back(json::object{{"original", shell.original}, {"retained", shell.retained}});
        targets.emplace_back(json::object{{"contributor_index", target.contributor_index}, {"atom_id", target.atom_id},
            {"nearest_halo_distance", Number(target.nearest_halo_distance)},
            {"original_rows", target.original_rows}, {"retained_rows", target.retained_rows},
            {"retained_fraction", target.original_rows ? Number(static_cast<double>(target.retained_rows) /
                static_cast<double>(target.original_rows)) : json::value(nullptr)}, {"radial", std::move(shells)}});
    }
    for (const auto & shell : summary.radial)
        radial.emplace_back(json::object{{"original_target_memberships", shell.original},
            {"retained_target_memberships", shell.retained}});
    for (const auto count : summary.halo_proximity_bins) proximity.emplace_back(count);
    for (const auto count : summary.observation_halo_distance_bins) observation_proximity.emplace_back(count);
    return {{"radius_angstrom", radius}, {"original_rows", summary.original_rows},
        {"excluded_rows", summary.excluded_rows}, {"retained_rows", summary.retained_rows},
        {"retained_fraction", summary.original_rows ? static_cast<double>(summary.retained_rows) /
            static_cast<double>(summary.original_rows) : 0.0},
        {"mixed_target_halo_rows", summary.mixed_target_halo_rows},
        {"excluded_mixed_target_halo_rows", summary.excluded_mixed_target_halo_rows},
        {"rows_shared_by_multiple_targets", summary.rows_shared_by_targets},
        {"retained_rows_shared_by_multiple_targets", summary.retained_rows_shared_by_targets},
        {"max_targets_supporting_one_row", summary.max_targets_per_row},
        {"target_components", summary.target_components},
        {"targets_without_retained_rows", summary.targets_without_retained_rows},
        {"target_radial_diversity_counts", json::value_from(summary.target_radial_diversity)},
        {"target_halo_distance_bins_lt1_1to1.5_1.5to2_2to2.5_ge2.5", std::move(proximity)},
        {"observation_halo_distance_bins_lt0.5_0.5to1_1to1.5_1.5to2_2to2.5_ge2.5", std::move(observation_proximity)},
        {"radial_target_memberships", std::move(radial)}, {"targets", std::move(targets)}};
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

void Run(const Path & model_path, const Path & map_path, const Path & output_path)
{
    const auto started = Clock::now();
    auto model = rhbm_gem::ReadModel(model_path);
    auto map = rhbm_gem::ReadMap(map_path);
    const auto normalization = map->MapValueArrayNormalization();
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement() != Element::HYDROGEN && atom.GetChainID() == "A";
    });
    const auto problem = core::BuildJointProblem(*map, *model);
    const auto & input = problem.Input();

    std::vector<audit::Position> observations;
    observations.reserve(input.row_ids.size());
    for (const auto & row_id : input.row_ids)
        observations.push_back(map->GetGridPosition(std::stoull(row_id)));
    std::vector<audit::Position> halos, target_positions;
    const auto & target_indices = input.selection_domain->target_indices;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        const auto * object = model->FindAtomPtr(std::stoi(input.atom_ids[atom]));
        if (!object) throw std::runtime_error("Joint contributor is missing from the model.");
        const audit::Position position = object->GetPosition();
        if (std::binary_search(target_indices.begin(), target_indices.end(), atom)) target_positions.push_back(position);
        else halos.push_back(position);
    }
    const auto nearest_rows = audit::NearestHaloDistances(observations, halos);
    const auto nearest_targets = audit::NearestHaloDistances(target_positions, halos);
    json::array radius_results;
    for (const double radius : audit::DefaultRadii)
    {
        const auto retained = audit::RetainedRows(nearest_rows, radius);
        const auto summary = audit::AnalyzeCoverage(input, retained, nearest_targets, nearest_rows);
        radius_results.emplace_back(CoverageJson(radius, summary));
    }
    if (!output_path.parent_path().empty()) std::filesystem::create_directories(output_path.parent_path());
    Write(output_path, json::object{{"contract", "joint-halo-observation-exclusion-coverage-v1"},
        {"model", model_path.string()}, {"map", map_path.string()},
        {"selection", "non-hydrogen chain A"},
        {"observation_domain", "production BuildJointProblem union of selected-target 2.5-A supports"},
        {"map_normalization_divisor", normalization ? Number(*normalization) : json::value(1.0)},
        {"observation_scale", problem.ObservationScale()}, {"target_count", target_indices.size()},
        {"halo_count", halos.size()}, {"observation_count", input.observations.size()},
        {"thread_settings", json::object{{"Eigen", Eigen::nbThreads()},
            {"OpenMP", std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset"},
            {"OpenBLAS", std::getenv("OPENBLAS_NUM_THREADS") ? std::getenv("OPENBLAS_NUM_THREADS") : "unset"},
            {"Accelerate", std::getenv("VECLIB_MAXIMUM_THREADS") ? std::getenv("VECLIB_MAXIMUM_THREADS") : "unset"}}},
        {"elapsed_seconds", std::chrono::duration<double>(Clock::now() - started).count()},
        {"peak_rss_bytes", PeakRssBytes()}, {"radii_angstrom", json::value_from(audit::DefaultRadii)},
        {"radial_shell_edges_angstrom", json::value_from(audit::RadialShellEdges)},
        {"results", std::move(radius_results)}});
}
}

int main(int argc, char ** argv)
{
    try
    {
        Eigen::setNbThreads(1);
        if (argc != 4) throw std::invalid_argument("Usage: joint_halo_observation_exclusion MODEL MAP OUTPUT.json");
        Run(argv[1], argv[2], argv[3]);
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
