#include "support/JointHaloObservationExclusion.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <string>
#include <unordered_set>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace json = boost::json;
namespace core = rhbm_gem::core;
namespace audit = joint_test::halo_exclusion;
using Path = std::filesystem::path;
namespace joint = core::joint_component;
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

json::value SpectrumJson(const audit::DesignSpectrum & spectrum)
{
    return json::object{{"rows", spectrum.rows}, {"columns", spectrum.columns}, {"rank", spectrum.rank},
        {"available", spectrum.available}, {"full_rank", spectrum.full_rank},
        {"minimum_singular_value", spectrum.available ? Number(spectrum.minimum_singular_value) : json::value(nullptr)},
        {"maximum_singular_value", spectrum.available ? Number(spectrum.maximum_singular_value) : json::value(nullptr)},
        {"condition", spectrum.full_rank ? Number(spectrum.condition) : json::value(nullptr)},
        {"ill_conditioned", spectrum.ill_conditioned}};
}

bool ValidTargetSeed(const rhbm_gem::JointInitializationAtom & seed)
{
    return seed.target && seed.first_stage_fit_qualified && seed.seed_source == "fitted" &&
        std::isfinite(seed.mdpde[0]) && seed.mdpde[0] >= 0.0 &&
        std::isfinite(seed.mdpde[1]) && seed.mdpde[1] > 0.0 && std::isfinite(seed.mdpde[2]);
}

struct ComponentBlocks
{
    joint::ComponentProblemView view;
    joint::StructuralBlockPartition partition;
};

std::vector<ComponentBlocks> ExistingBlocks(const joint::ProblemData & data)
{
    std::vector<ComponentBlocks> result;
    result.reserve(data.partition.components.size());
    const joint::FixedNeighborSearchPolicy production_policy;
    for (const auto & component : data.partition.components)
    {
        auto view = joint::MakeComponentProblemView(data.layout, component);
        joint::StructuralTopology topology(view.domain, view.atom_ids, view.local_layout);
        auto partition = topology.BuildCores(production_policy.core_atoms);
        result.push_back({std::move(view), std::move(partition)});
    }
    return result;
}

audit::LocalDesignMetrics TargetDesign(const core::JointProblemInput & input,
    const rhbm_gem::JointInitialization & initialization, std::size_t contributor,
    const audit::RowMask & retained)
{
    const auto & seed = initialization.atoms.at(contributor);
    const auto & support = input.support.at(contributor);
    std::size_t rows{};
    for (const auto & entry : support) if (retained.at(entry.row)) ++rows;
    Eigen::MatrixXd ac = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rows), 2);
    Eigen::MatrixXd width = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rows), 1);
    Eigen::Index row{};
    for (const auto & entry : support)
    {
        if (!retained.at(entry.row)) continue;
        const auto basis = joint::EvaluateKernel(entry.squared_distance, seed.mdpde[1], 2.5);
        ac(row, 0) = basis.gaussian;
        ac(row, 1) = basis.charge;
        width(row, 0) = seed.mdpde[0] * basis.gaussian_log_width + seed.mdpde[2] * basis.charge_log_width;
        ++row;
    }
    return audit::AnalyzeLocalDesign(ac, width);
}

struct BlockDesign
{
    bool available{};
    std::string reason;
    std::size_t retained_rows{};
    std::vector<std::size_t> contributors;
    audit::LocalDesignMetrics metrics;
};

BlockDesign BlockDesignAtRadius(const rhbm_gem::JointInitialization & initialization, const joint::ProblemData & data,
    const ComponentBlocks & component, const joint::StructuralCore & block, const audit::RowMask & retained)
{
    BlockDesign result;
    std::vector<Eigen::Index> rows;
    rows.reserve(block.affected_rows.size());
    for (const auto local_row : block.affected_rows)
    {
        const auto parent_row = component.view.local_to_parent_rows.at(static_cast<std::size_t>(local_row));
        if (retained.at(static_cast<std::size_t>(parent_row))) rows.push_back(local_row);
    }
    result.retained_rows = rows.size();
    std::vector<std::size_t> target_contributors;
    target_contributors.reserve(block.atoms.size());
    for (const auto local_atom : block.atoms)
    {
        const auto component_atom = static_cast<std::size_t>(local_atom);
        const auto parent_target = static_cast<std::size_t>(component.view.local_to_parent_atoms.at(component_atom));
        const auto contributor = data.target_to_contributor.at(parent_target);
        if (!ValidTargetSeed(initialization.atoms.at(contributor)))
        {
            result.reason = "unqualified-first-stage-target";
            return result;
        }
        target_contributors.push_back(contributor);
    }
    Eigen::MatrixXd ac = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rows.size()),
        static_cast<Eigen::Index>(2 * block.atoms.size()));
    Eigen::MatrixXd width = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rows.size()),
        static_cast<Eigen::Index>(block.atoms.size()));
    for (std::size_t local_atom = 0; local_atom < block.atoms.size(); ++local_atom)
    {
        const auto atom_index = static_cast<std::size_t>(block.atoms[local_atom]);
        const auto & seed = initialization.atoms.at(target_contributors[local_atom]);
        for (const auto support : component.view.domain.atoms[atom_index])
        {
            const auto position = std::lower_bound(rows.begin(), rows.end(), support.row);
            if (position == rows.end() || *position != support.row) continue;
            const auto row = static_cast<Eigen::Index>(position - rows.begin());
            const auto basis = joint::EvaluateKernel(support.square, seed.mdpde[1], 2.5);
            ac(row, static_cast<Eigen::Index>(2 * local_atom)) = basis.gaussian;
            ac(row, static_cast<Eigen::Index>(2 * local_atom + 1)) = basis.charge;
            width(row, static_cast<Eigen::Index>(local_atom)) =
                seed.mdpde[0] * basis.gaussian_log_width + seed.mdpde[2] * basis.charge_log_width;
        }
    }
    result.available = true;
    result.contributors = std::move(target_contributors);
    result.metrics = audit::AnalyzeLocalDesign(ac, width);
    return result;
}

json::object TargetDesignJson(const std::string & atom_id, double halo_distance,
    std::size_t original_rows, std::size_t retained_rows, const audit::LocalDesignMetrics & metrics, bool valid)
{
    return {{"atom_id", atom_id}, {"nearest_halo_distance", Number(halo_distance)},
        {"original_rows", original_rows}, {"retained_rows", retained_rows}, {"available", valid},
        {"ac", SpectrumJson(metrics.ac)}, {"abc", SpectrumJson(metrics.abc)},
        {"projected_width", SpectrumJson(metrics.projected_width)},
        {"minimum_projected_width_fraction", metrics.ac.available && metrics.projected_width.available ?
            Number(metrics.minimum_projected_width_fraction) : json::value(nullptr)}};
}

json::object BlockDesignJson(const core::JointProblemInput & input, const ComponentBlocks & component,
    const joint::StructuralCore & block, const BlockDesign & design)
{
    json::array atom_ids;
    for (const auto contributor : design.contributors) atom_ids.emplace_back(input.atom_ids.at(contributor));
    if (!design.available)
        return {{"component_seed_atom_id", component.view.atom_ids.size() == 0 ? "" : component.view.atom_ids[0]},
            {"block_id", block.id}, {"core_size", block.atoms.size()},
            {"original_affected_rows", block.affected_rows.size()}, {"retained_affected_rows", design.retained_rows},
            {"available", false}, {"reason", design.reason}};
    return {{"component_seed_atom_id", component.view.atom_ids.size() == 0 ? "" : component.view.atom_ids[0]},
        {"block_id", block.id}, {"core_size", design.contributors.size()},
        {"original_affected_rows", block.affected_rows.size()},
        {"retained_affected_rows", design.metrics.ac.rows}, {"available", true},
        {"ac", SpectrumJson(design.metrics.ac)}, {"abc", SpectrumJson(design.metrics.abc)},
        {"projected_width", SpectrumJson(design.metrics.projected_width)},
        {"minimum_projected_width_fraction", Number(design.metrics.minimum_projected_width_fraction)},
        {"target_contributor_indices", json::value_from(design.contributors)}, {"target_atom_ids", std::move(atom_ids)}};
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

json::object IdentifiabilityJson(double radius, const audit::CoverageSummary & coverage,
    const core::JointProblem & problem, const rhbm_gem::JointInitialization & initialization,
    const std::vector<ComponentBlocks> & components, const audit::RowMask & retained)
{
    const auto & input = problem.Input();
    const auto & data = core::JointProblemAccess::Get(problem);
    json::array targets, blocks, unqualified_targets;
    std::size_t valid_targets{}, invalid_targets{}, target_ac_deficient{}, target_abc_deficient{},
        target_width_deficient{}, target_ill_conditioned{}, block_count{}, block_unavailable{},
        block_ac_deficient{}, block_abc_deficient{}, block_width_deficient{}, block_ill_conditioned{};
    std::unordered_set<std::size_t> targets_in_rank_deficient_blocks;
    for (std::size_t index = 0; index < input.selection_domain->target_indices.size(); ++index)
    {
        const auto contributor = input.selection_domain->target_indices[index];
        const auto & seed = initialization.atoms.at(contributor);
        const bool valid = ValidTargetSeed(seed);
        if (valid) ++valid_targets; else ++invalid_targets;
        if (!valid) unqualified_targets.emplace_back(seed.id);
        const auto metrics = valid ? TargetDesign(input, initialization, contributor, retained) : audit::LocalDesignMetrics{};
        const auto & coverage_target = coverage.targets[index];
        if (valid)
        {
            target_ac_deficient += metrics.ac.rank < metrics.ac.columns;
            target_abc_deficient += metrics.abc.rank < metrics.abc.columns;
            target_width_deficient += metrics.projected_width.rank < metrics.projected_width.columns;
            target_ill_conditioned += metrics.ac.ill_conditioned || metrics.abc.ill_conditioned ||
                metrics.projected_width.ill_conditioned;
        }
        targets.emplace_back(TargetDesignJson(seed.id, coverage_target.nearest_halo_distance,
            coverage_target.original_rows, coverage_target.retained_rows, metrics, valid));
    }
    for (const auto & component : components)
    {
        for (const auto & block : component.partition.cores)
        {
            ++block_count;
            const auto design = BlockDesignAtRadius(initialization, data, component, block, retained);
            blocks.emplace_back(BlockDesignJson(input, component, block, design));
            if (!design.available) { ++block_unavailable; continue; }
            const auto & metrics = design.metrics;
            const bool ac_deficient = metrics.ac.rank < metrics.ac.columns;
            const bool abc_deficient = metrics.abc.rank < metrics.abc.columns;
            const bool width_deficient = metrics.projected_width.rank < metrics.projected_width.columns;
            block_ac_deficient += ac_deficient;
            block_abc_deficient += abc_deficient;
            block_width_deficient += width_deficient;
            block_ill_conditioned += metrics.ac.ill_conditioned || metrics.abc.ill_conditioned ||
                metrics.projected_width.ill_conditioned;
            if (ac_deficient || abc_deficient || width_deficient)
                targets_in_rank_deficient_blocks.insert(design.contributors.begin(), design.contributors.end());
        }
    }
    return {{"radius_angstrom", radius}, {"first_stage_target_seeds_qualified", valid_targets},
        {"first_stage_target_seeds_unqualified", invalid_targets}, {"target_singleton_ac_rank_deficient", target_ac_deficient},
        {"unqualified_target_atom_ids", std::move(unqualified_targets)},
        {"target_singleton_abc_rank_deficient", target_abc_deficient},
        {"target_singleton_projected_width_rank_deficient", target_width_deficient},
        {"target_singleton_ill_conditioned", target_ill_conditioned},
        {"fixed_neighbor_core_size", joint::FixedNeighborSearchPolicy{}.core_atoms},
        {"fixed_neighbor_blocks", block_count}, {"blocks_unavailable", block_unavailable},
        {"blocks_ac_rank_deficient", block_ac_deficient}, {"blocks_abc_rank_deficient", block_abc_deficient},
        {"blocks_projected_width_rank_deficient", block_width_deficient},
        {"blocks_ill_conditioned", block_ill_conditioned},
        {"targets_in_rank_deficient_blocks", targets_in_rank_deficient_blocks.size()},
        {"targets", std::move(targets)}, {"blocks", std::move(blocks)}};
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
    const auto initialization = joint::InitializeContributors(*map, *model, problem);
    const auto blocks = ExistingBlocks(core::JointProblemAccess::Get(problem));
    std::size_t qualified_halos{}, unqualified_halos{};
    json::array unqualified_halo_ids;
    for (const auto & seed : initialization.atoms)
    {
        if (seed.target) continue;
        if (seed.first_stage_fit_qualified) ++qualified_halos;
        else { ++unqualified_halos; unqualified_halo_ids.emplace_back(seed.id); }
    }

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
        radius_results.emplace_back(json::object{{"coverage", CoverageJson(radius, summary)},
            {"identifiability", IdentifiabilityJson(radius, summary, problem, initialization, blocks, retained)}});
    }
    std::filesystem::create_directories(output_path.parent_path());
    Write(output_path, json::object{{"contract", "joint-halo-observation-exclusion-feasibility-v1"},
        {"model", model_path.string()}, {"map", map_path.string()},
        {"selection", "non-hydrogen chain A"}, {"observation_domain", "production BuildJointProblem union of selected-target 2.5-A supports"},
        {"map_normalization_divisor", normalization ? Number(*normalization) : json::value(1.0)},
        {"observation_scale", problem.ObservationScale()}, {"target_count", target_indices.size()},
        {"halo_count", halos.size()}, {"observation_count", input.observations.size()},
        {"first_stage_initialization_valid", initialization.valid}, {"first_stage_initialization_reason", initialization.reason},
        {"first_stage_halos_qualified", qualified_halos}, {"first_stage_halos_unqualified", unqualified_halos},
        {"unqualified_halo_atom_ids", std::move(unqualified_halo_ids)},
        {"original_target_components", core::JointProblemAccess::Get(problem).partition.components.size()},
        {"thread_settings", json::object{{"Eigen", Eigen::nbThreads()}, {"OpenMP", std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset"},
            {"OpenBLAS", std::getenv("OPENBLAS_NUM_THREADS") ? std::getenv("OPENBLAS_NUM_THREADS") : "unset"},
            {"Accelerate", std::getenv("VECLIB_MAXIMUM_THREADS") ? std::getenv("VECLIB_MAXIMUM_THREADS") : "unset"}}},
        {"elapsed_seconds", std::chrono::duration<double>(Clock::now() - started).count()},
        {"peak_rss_bytes", PeakRssBytes()},
        {"radii_angstrom", json::value_from(audit::DefaultRadii)},
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
