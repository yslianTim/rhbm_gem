#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/GroupPotentialFitting.hpp"
#include "core/detail/HaloMedianReconstruction.hpp"
#include "core/detail/MapSampler.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "data/detail/JointStageAdapter.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <boost/json.hpp>
#include <Eigen/SparseCore>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace rg = rhbm_gem;
namespace core = rg::core;
namespace detail = core::detail;
namespace joint = core::joint_component;
namespace json = boost::json;
using Position = std::array<double, 3>;
using Profile = std::array<double, 3>; // A, B, C
using Clock = std::chrono::steady_clock;

struct Scenario
{
    std::string name;
    double halo_distance{};
    bool contrast{};
    bool wide_charge_contrast{};
    double noise{};
    bool map_edge_halo{};
    std::size_t extra_free_halo_observations{};
    bool real_geometry{}, real_full_structure{};
    bool require_short_bonded_halo_pair{};
    int sequence_start{}, sequence_end{};
    std::uint64_t noise_seed{};
};

struct Halo
{
    int id{};
    Position position{};
    Profile truth{};
    std::vector<std::size_t> rows;
};

struct MethodResult
{
    bool available{};
    std::string reason;
    std::vector<std::optional<Profile>> profiles;
    std::optional<double> a_rmse, b_rmse, c_rmse;
    std::optional<double> a_bias, b_bias, c_bias;
    std::optional<double> a_median_abs, b_median_abs, c_median_abs;
    std::optional<double> a_p95_abs, b_p95_abs, c_p95_abs;
    std::optional<double> a_max_abs, b_max_abs, c_max_abs;
    std::optional<double> halo_signal_rmse, shared_signal_rmse, domain_signal_rmse, max_signal_error;
    std::optional<double> peeling_rmse, tail_ratio_mae, local_a_rmse, local_b_rmse;
    std::optional<double> group_a_error, group_b_error;
    std::optional<double> target_endpoint_drift_max;
    std::optional<double> local_effective_sample_mean, local_information_trace_mean;
    std::optional<std::size_t> local_effective_sample_min;
    std::optional<double> group_alpha_g, group_mdpde_a, group_mdpde_b;
    std::optional<double> group_posterior_a_rmse, group_posterior_b_rmse;
    std::optional<double> group_prior_uncertainty_a, group_prior_uncertainty_b;
    std::optional<double> group_member_uncertainty_a_mean, group_member_uncertainty_b_mean;
    std::size_t peeling_available{}, peeling_total{}, local_available{}, local_qualified_count{};
    std::size_t local_sample_count_total{}, local_effective_sample_count_total{}, local_information_count{};
    std::size_t group_members{}, group_outlier_count{};
    double seconds{};
};

double Distance(const Position & lhs, const Position & rhs)
{
    double square{};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const double delta = lhs[axis] - rhs[axis];
        square += delta * delta;
    }
    return std::sqrt(square);
}

double SupportSquare(const Position & lhs, const Position & rhs)
{
    double square{};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const double delta = lhs[axis] - rhs[axis];
        square += delta * delta;
    }
    return square;
}

std::size_t Grid(const std::array<int, 3> & size, int x, int y, int z)
{
    return static_cast<std::size_t>(x) + static_cast<std::size_t>(size[0]) *
        (static_cast<std::size_t>(y) + static_cast<std::size_t>(size[1]) * static_cast<std::size_t>(z));
}

rg::AtomObject MakeAtom(int id, const Position & position, bool target)
{
    rg::AtomObject atom;
    atom.SetSerialID(id);
    atom.SetSequenceID(id);
    atom.SetComponentID("ALA");
    atom.SetAtomID(target ? "CA" : "CB");
    atom.SetChainID("A");
    atom.SetIndicator(".");
    atom.SetComponentKey(target ? static_cast<ComponentKey>(1) : static_cast<ComponentKey>(id));
    atom.SetAtomKey(static_cast<AtomKey>(target ? static_cast<int>(Spot::CA) : id));
    atom.SetResidue(Residue::ALA);
    atom.SetElement(Element::CARBON);
    atom.SetSpot(target ? Spot::CA : Spot::CB);
    atom.SetStructure(Structure::FREE);
    atom.SetPosition(position);
    return atom;
}

struct Synthetic
{
    rg::ModelObject model;
    std::unique_ptr<rg::MapObject> map;
    std::unique_ptr<rg::MapObject> target_only_map;
    std::vector<Position> target_positions;
    std::vector<Profile> target_truth;
    std::unordered_map<int,Profile> target_truth_by_id;
    std::vector<Halo> halos;
    std::set<std::pair<int,int>> covalent_pairs;
    std::optional<double> map_normalization_divisor;
    std::array<int, 3> grid_size{51, 63, 16};
    std::array<double, 3> spacing{0.5, 0.5, 0.5};
    std::array<double, 3> origin{-3.0, -3.0, -3.0};
};

struct RealGeometrySource
{
    std::filesystem::path model, map, truth;
};

std::optional<RealGeometrySource> real_geometry_source;

Synthetic MakeSynthetic(const Scenario & scenario)
{
    Synthetic out;
    out.target_positions.reserve(20);
    out.target_truth.assign(20, Profile{1.0, 1.10, 0.10});
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    atoms.reserve(27);
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 4; ++x)
        {
            const Position position{6.2 * x, 6.2 * y, 0.0};
            out.target_positions.push_back(position);
            atoms.push_back(std::make_unique<rg::AtomObject>(MakeAtom(
                static_cast<int>(out.target_positions.size()), position, true)));
            out.target_truth_by_id.emplace(static_cast<int>(out.target_positions.size()),
                out.target_truth.back());
        }

    const Position target0 = out.target_positions.front();
    const Position halo0{target0[0], target0[1], 3.0};
    Position halo1{};
    if (scenario.halo_distance > 5.0)
        halo1 = {out.target_positions[1][0], out.target_positions[1][1], 3.0};
    else
    {
        const double angle = 2.0 * std::asin(scenario.halo_distance / 6.0);
        halo1 = {target0[0] + 3.0 * std::sin(angle), target0[1], 3.0 * std::cos(angle)};
    }
    const Profile similar0{0.96, 1.08, 0.08}, similar1{1.04, 1.12, 0.12};
    const Profile contrast0{1.55, scenario.wide_charge_contrast ? 1.20 : 1.12, 0.38};
    const Profile contrast1{0.55, scenario.wide_charge_contrast ? 0.65 : 0.72, -0.28};
    const auto add_halo = [&](int id, Position position, Profile truth) {
        out.halos.push_back({id, position, truth, {}});
        atoms.push_back(std::make_unique<rg::AtomObject>(MakeAtom(id, position, false)));
    };
    add_halo(101, halo0, scenario.contrast ? contrast0 : similar0);
    add_halo(102, halo1, scenario.contrast ? contrast1 : similar1);
    // One halo in each corner has exactly 3, 2, 1, and 0 distinct OriginalInput rows.
    if (scenario.extra_free_halo_observations)
        add_halo(103, {20.6, scenario.extra_free_halo_observations == 3 ? 27.9 : 27.8, -3.0},
            {1.0, 1.10, 0.10});
    add_halo(104, {-2.5, 28.10, -2.5}, {1.0, 1.10, 0.10});
    add_halo(105, {-1.6, -3.40, -3.0}, {1.0, 1.10, 0.10});
    add_halo(106, {21.1, -3.20, -2.5}, {1.0, 1.10, 0.10});
    if (scenario.map_edge_halo)
        add_halo(107, {22.0, 27.5, 4.2}, {0.70, 1.05, -0.15});

    out.model = rg::ModelObject(std::move(atoms));
    out.model.SelectAtoms([](const rg::AtomObject & atom) { return atom.GetSerialID() <= 20; });
    out.model.EditAnalysis().InitializeFromSelection();

    const std::size_t voxel_count = static_cast<std::size_t>(out.grid_size[0]) *
        static_cast<std::size_t>(out.grid_size[1]) * static_cast<std::size_t>(out.grid_size[2]);
    auto raw_values = std::make_unique<double[]>(voxel_count);
    auto target_values = std::make_unique<double[]>(voxel_count);
    for (int z = 0; z < out.grid_size[2]; ++z)
        for (int y = 0; y < out.grid_size[1]; ++y)
            for (int x = 0; x < out.grid_size[0]; ++x)
            {
                const Position position{out.origin[0] + out.spacing[0] * x,
                    out.origin[1] + out.spacing[1] * y, out.origin[2] + out.spacing[2] * z};
                double target_value{}, halo_value{};
                for (std::size_t atom = 0; atom < out.target_positions.size(); ++atom)
                {
                    const double square = SupportSquare(position, out.target_positions[atom]);
                    if (square <= 6.25)
                    {
                        const auto basis = joint::EvaluateKernel(square, out.target_truth[atom][1], 2.5);
                        target_value += out.target_truth[atom][0] * basis.gaussian + out.target_truth[atom][2] * basis.charge;
                    }
                }
                for (const auto & halo : out.halos)
                {
                    const double square = SupportSquare(position, halo.position);
                    if (square <= 6.25)
                    {
                        const auto basis = joint::EvaluateKernel(square, halo.truth[1], 2.5);
                        halo_value += halo.truth[0] * basis.gaussian + halo.truth[2] * basis.charge;
                    }
                }
                const double deterministic_noise = scenario.noise *
                    (std::sin(0.017 * static_cast<double>(Grid(out.grid_size, x, y, z))) +
                     0.5 * std::sin(0.113 * static_cast<double>(Grid(out.grid_size, x, y, z))));
                const auto row = Grid(out.grid_size, x, y, z);
                target_values[row] = target_value;
                raw_values[row] = target_value + halo_value + deterministic_noise;
            }
    out.map = std::make_unique<rg::MapObject>(out.grid_size, out.spacing, out.origin, std::move(raw_values));
    out.target_only_map = std::make_unique<rg::MapObject>(out.grid_size, out.spacing, out.origin, std::move(target_values));
    return out;
}

std::unordered_map<int,Profile> ReadSimulationTruth(const std::filesystem::path & path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read 6Z6U simulation truth: " + path.string());
    json::parse_options options;
    options.numbers = json::number_precision::precise;
    const auto metadata = json::parse(std::string(std::istreambuf_iterator<char>(input), {}), {}, options);
    std::unordered_map<int,Profile> truth;
    for (const auto & atom : metadata.at("atoms").as_array())
    {
        const int id = json::value_to<int>(atom.at("serial_id"));
        const double gaussian_width = json::value_to<double>(atom.at("effective_gaussian_width"));
        const double charge_width = json::value_to<double>(atom.at("effective_charge_width"));
        if (gaussian_width != charge_width)
            throw std::runtime_error("The real-geometry runner requires one common Gaussian/charge width per atom.");
        truth.emplace(id, Profile{static_cast<double>(json::value_to<int>(atom.at("element"))),
            gaussian_width, json::value_to<double>(atom.at("charge_used"))});
    }
    return truth;
}

void AddTruthAtomToGrid(const rg::AtomObject & atom, const Profile & profile,
    const std::array<int,3> & dimensions, const std::array<double,3> & spacing,
    const std::array<double,3> & origin, std::vector<double> & values)
{
    std::array<int,3> lower{}, upper{};
    const auto position = atom.GetPosition();
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        lower[axis] = std::max(0, static_cast<int>(std::ceil(
            (position[axis] - 2.5 - origin[axis]) / spacing[axis])));
        upper[axis] = std::min(dimensions[axis] - 1, static_cast<int>(std::floor(
            (position[axis] + 2.5 - origin[axis]) / spacing[axis])));
    }
    for (int z = lower[2]; z <= upper[2]; ++z)
        for (int y = lower[1]; y <= upper[1]; ++y)
            for (int x = lower[0]; x <= upper[0]; ++x)
            {
                const Position grid_position{origin[0] + spacing[0] * x,
                    origin[1] + spacing[1] * y, origin[2] + spacing[2] * z};
                const double square = SupportSquare(grid_position, position);
                if (square > 6.25) continue;
                const auto basis = joint::EvaluateKernel(square, profile[1], 2.5);
                const auto index = static_cast<std::size_t>(x) + static_cast<std::size_t>(dimensions[0]) *
                    (static_cast<std::size_t>(y) + static_cast<std::size_t>(dimensions[1]) *
                    static_cast<std::size_t>(z));
                values.at(index) += profile[0] * basis.gaussian + profile[2] * basis.charge;
            }
}

std::pair<int,int> OrderedPair(int left, int right)
{ return left < right ? std::pair{left,right} : std::pair{right,left}; }

Synthetic MakeRealGeometrySynthetic(const Scenario & scenario)
{
    if (!real_geometry_source) throw std::logic_error("Real-geometry inputs were not configured.");
    auto source_model = rg::ReadModel(real_geometry_source->model.string());
    if (!source_model) throw std::runtime_error("Could not load the real-geometry model.");
    const auto truths = ReadSimulationTruth(real_geometry_source->truth);
    Synthetic out;
    std::unordered_set<int> target_ids;
    std::array<double,3> minimum{std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
    std::array<double,3> maximum{-std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};

    if (scenario.real_full_structure)
    {
        auto source_map = rg::ReadMap(real_geometry_source->map.string());
        if (!source_map) throw std::runtime_error("Could not load the 6Z6U map.");
        out.model = std::move(*source_model);
        out.model.SelectAtoms([](const rg::AtomObject & atom) {
            return atom.GetElement() != Element::HYDROGEN && atom.GetChainID() == "A";
        });
        if (out.model.GetSelectedAtomCount() != 1418)
            throw std::runtime_error("The 6Z6U chain-A target count differs from the recorded 1,418 baseline.");
        out.model.EditAnalysis().InitializeFromSelection();
        out.map = std::move(source_map);
        out.map_normalization_divisor = out.map->MapValueArrayNormalization();
        const double truth_scale = out.map_normalization_divisor ? 1.0 / *out.map_normalization_divisor : 1.0;
        out.grid_size = out.map->GetGridSize();
        out.spacing = out.map->GetGridSpacing();
        out.origin = out.map->GetOrigin();
        for (const auto & bond : out.model.GetBondList())
            out.covalent_pairs.insert(OrderedPair(bond->GetAtomSerialID1(), bond->GetAtomSerialID2()));
        std::size_t maximum_id{};
        for (const auto & owned : out.model.GetAtomList())
            maximum_id = std::max(maximum_id, static_cast<std::size_t>(owned->GetSerialID()));
        out.target_truth.assign(maximum_id, Profile{});
        std::vector<double> target_values(out.map->GetMapValueArraySize(), 0.0);
        for (const auto * atom : out.model.GetSelectedAtoms())
        {
            const auto found = truths.find(atom->GetSerialID());
            if (found == truths.end()) throw std::runtime_error("Missing 6Z6U truth for a selected Target.");
            target_ids.insert(atom->GetSerialID());
            auto normalized_truth = found->second;
            normalized_truth[0] *= truth_scale;
            normalized_truth[2] *= truth_scale;
            out.target_truth_by_id.emplace(atom->GetSerialID(), normalized_truth);
            out.target_truth.at(static_cast<std::size_t>(atom->GetSerialID() - 1)) = normalized_truth;
            out.target_positions.push_back(atom->GetPosition());
            AddTruthAtomToGrid(*atom, normalized_truth, out.grid_size, out.spacing, out.origin, target_values);
        }
        for (const auto & owned : out.model.GetAtomList())
        {
            const int id = owned->GetSerialID();
            if (owned->GetElement() == Element::HYDROGEN || target_ids.contains(id)) continue;
            const auto found = truths.find(id);
            if (found == truths.end()) throw std::runtime_error("Missing 6Z6U truth for a Halo candidate.");
            auto normalized_truth = found->second;
            normalized_truth[0] *= truth_scale;
            normalized_truth[2] *= truth_scale;
            out.halos.push_back({id, owned->GetPosition(), normalized_truth, {}});
        }
        auto target_array = std::make_unique<double[]>(target_values.size());
        std::copy(target_values.begin(), target_values.end(), target_array.get());
        out.target_only_map = std::make_unique<rg::MapObject>(out.grid_size, out.spacing, out.origin,
            std::move(target_array));
        return out;
    }

    std::vector<const rg::AtomObject *> seed_targets;
    for (const auto & atom : source_model->GetAtomList())
    {
        if (atom->GetElement() == Element::HYDROGEN || atom->GetChainID() != "A" ||
            atom->GetSequenceID() < scenario.sequence_start || atom->GetSequenceID() > scenario.sequence_end)
            continue;
        seed_targets.push_back(atom.get());
        const auto position = atom->GetPosition();
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
            minimum[axis] = std::min(minimum[axis], position[axis]);
            maximum[axis] = std::max(maximum[axis], position[axis]);
        }
    }
    if (seed_targets.empty()) throw std::runtime_error("The requested real-geometry sequence window has no targets.");
    for (const auto * atom : seed_targets) target_ids.insert(atom->GetSerialID());
    const double padding = 5.0;
    std::vector<std::unique_ptr<rg::AtomObject>> local_atoms;
    std::size_t maximum_id{};
    for (const auto & owned : source_model->GetAtomList())
    {
        if (owned->GetElement() == Element::HYDROGEN) continue;
        if (owned->GetChainID() == "A" && !target_ids.contains(owned->GetSerialID())) continue;
        const auto position = owned->GetPosition();
        bool inside = true;
        for (std::size_t axis = 0; axis < 3; ++axis)
            inside &= position[axis] >= minimum[axis] - padding && position[axis] <= maximum[axis] + padding;
        if (!inside) continue;
        maximum_id = std::max(maximum_id, static_cast<std::size_t>(owned->GetSerialID()));
        local_atoms.push_back(std::make_unique<rg::AtomObject>(*owned));
    }
    out.model = rg::ModelObject(std::move(local_atoms));
    out.model.SelectAtoms([&target_ids](const rg::AtomObject & atom) {return target_ids.contains(atom.GetSerialID());});
    out.model.EditAnalysis().InitializeFromSelection();
    for (const auto & bond : source_model->GetBondList())
        out.covalent_pairs.insert(OrderedPair(bond->GetAtomSerialID1(), bond->GetAtomSerialID2()));

    constexpr double grid_spacing = 0.5;
    const double grid_padding = scenario.map_edge_halo ? 2.5 : 3.0;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        out.origin[axis] = std::floor((minimum[axis] - grid_padding) / grid_spacing) * grid_spacing;
        out.grid_size[axis] = static_cast<int>(std::ceil((maximum[axis] + grid_padding - out.origin[axis]) /
            grid_spacing)) + 1;
        out.spacing[axis] = grid_spacing;
    }
    const auto voxel_count = static_cast<std::size_t>(out.grid_size[0]) *
        static_cast<std::size_t>(out.grid_size[1]) * static_cast<std::size_t>(out.grid_size[2]);
    std::vector<double> raw_values(voxel_count, 0.0), target_values(voxel_count, 0.0);
    for (const auto * atom : out.model.GetSelectedAtoms())
    {
        const auto profile = truths.find(atom->GetSerialID());
        if (profile == truths.end()) throw std::runtime_error("Missing 6Z6U truth for a selected local Target.");
        target_ids.insert(atom->GetSerialID());
        out.target_truth_by_id.emplace(atom->GetSerialID(), profile->second);
        out.target_positions.push_back(atom->GetPosition());
        AddTruthAtomToGrid(*atom, profile->second, out.grid_size, out.spacing, out.origin, target_values);
        AddTruthAtomToGrid(*atom, profile->second, out.grid_size, out.spacing, out.origin, raw_values);
    }
    for (const auto & owned : out.model.GetAtomList())
    {
        const int id = owned->GetSerialID();
        if (owned->GetElement() == Element::HYDROGEN || target_ids.contains(id)) continue;
        const auto profile = truths.find(id);
        if (profile == truths.end()) throw std::runtime_error("Missing 6Z6U truth for a local Halo.");
        out.halos.push_back({id, owned->GetPosition(), profile->second, {}});
        AddTruthAtomToGrid(*owned, profile->second, out.grid_size, out.spacing, out.origin, raw_values);
    }
    if (scenario.noise > 0.0)
    {
        std::mt19937_64 generator(scenario.noise_seed);
        std::normal_distribution<double> normal(0.0, scenario.noise);
        for (double & value : raw_values) value += normal(generator);
    }
    auto raw_array = std::make_unique<double[]>(voxel_count);
    auto target_array = std::make_unique<double[]>(voxel_count);
    std::copy(raw_values.begin(), raw_values.end(), raw_array.get());
    std::copy(target_values.begin(), target_values.end(), target_array.get());
    out.map = std::make_unique<rg::MapObject>(out.grid_size, out.spacing, out.origin, std::move(raw_array));
    out.target_only_map = std::make_unique<rg::MapObject>(out.grid_size, out.spacing, out.origin,
        std::move(target_array));
    out.target_truth.assign(maximum_id, Profile{});
    for (const auto & [id, profile] : out.target_truth_by_id)
        out.target_truth.at(static_cast<std::size_t>(id - 1)) = profile;
    return out;
}

std::optional<Profile> TargetPoint(const rg::JointAnalysisResult & result, std::size_t target_index)
{
    for (const auto & component : result.components)
    {
        if (!component.state) continue;
        for (std::size_t local = 0; local < component.atoms.size(); ++local)
            if (component.atoms[local] == target_index &&
                component.state->ac.size() >= 2 * local + 2 && component.state->b.size() > local)
                return Profile{component.state->ac[2 * local], component.state->b[local],
                    component.state->ac[2 * local + 1]};
    }
    return std::nullopt;
}

using ProfileMap = std::unordered_map<int, Profile>;

std::optional<Profile> TargetTruth(const Synthetic & synthetic, int atom_id)
{
    if (const auto found = synthetic.target_truth_by_id.find(atom_id);
        found != synthetic.target_truth_by_id.end()) return found->second;
    if (atom_id > 0 && static_cast<std::size_t>(atom_id) <= synthetic.target_truth.size())
        return synthetic.target_truth[static_cast<std::size_t>(atom_id - 1)];
    return std::nullopt;
}

ProfileMap GeometryProfiles(const detail::HaloMedianBuildResult & built, const rg::ModelObject & model)
{
    ProfileMap out;
    for (const auto & halo : built.snapshot.halos)
        out.emplace(std::stoi(halo.atom_id), Profile{halo.amplitude, halo.width, halo.charge});
    detail::GeometryLocalTargetDonorIndex donor_index(built.target_donors);
    std::set<int> targets;
    for (const auto & donor : built.target_donors) targets.insert(std::stoi(donor.atom_id));
    for (const auto & atom : model.GetAtomList())
    {
        if (atom->GetElement() == Element::HYDROGEN || targets.contains(atom->GetSerialID()) ||
            out.contains(atom->GetSerialID()))
            continue;
        const auto donors = donor_index.Query(atom->GetPosition());
        if (donors.empty()) continue;
        out.emplace(atom->GetSerialID(), detail::GeometryLocalTargetMedian(donors));
    }
    return out;
}

ProfileMap FullABCProfiles(const detail::HaloFullABCReconstructionResult & result)
{
    ProfileMap out;
    for (const auto & halo : result.halos)
        if (halo.reconstructed_parameters)
            out.emplace(std::stoi(halo.classification.atom_id), *halo.reconstructed_parameters);
    return out;
}

std::string HaloParameterSourceName(detail::HaloFullABCParameterSource source)
{
    switch(source)
    {
    case detail::HaloFullABCParameterSource::Unavailable: return "Unavailable";
    case detail::HaloFullABCParameterSource::PreFixedGlobalMedian: return "PreFixedGlobalMedian";
    case detail::HaloFullABCParameterSource::FullABC: return "FullABC";
    case detail::HaloFullABCParameterSource::ObservationACFixedWidthRecovery:
        return "ObservationACFixedWidthRecovery";
    case detail::HaloFullABCParameterSource::PostFailureGlobalMedianRecovery:
        return "PostFailureGlobalMedianRecovery";
    }
    return "Unknown";
}

double Median(std::vector<double> values)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const auto middle = values.size() / 2;
    return values.size() % 2 ? values[middle] : 0.5 * (values[middle - 1] + values[middle]);
}

double Quantile(std::vector<double> values, double probability)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double location = probability * static_cast<double>(values.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(location));
    const auto upper = static_cast<std::size_t>(std::ceil(location));
    return values[lower] + (location - static_cast<double>(lower)) * (values[upper] - values[lower]);
}

std::optional<double> Rmse(const std::vector<double> & errors)
{
    if (errors.empty()) return std::nullopt;
    double sum{};
    for (const double error : errors) sum += error * error;
    return std::sqrt(sum / static_cast<double>(errors.size()));
}

void ScoreParameters(MethodResult & result, const std::vector<Halo> & halos)
{
    std::array<std::vector<double>, 3> errors, absolute;
    for (std::size_t index = 0; index < halos.size(); ++index)
    {
        if (index >= result.profiles.size() || !result.profiles[index]) return;
        for (std::size_t parameter = 0; parameter < 3; ++parameter)
        {
            const auto error = (*result.profiles[index])[parameter] - halos[index].truth[parameter];
            errors[parameter].push_back(error);
            absolute[parameter].push_back(std::abs(error));
        }
    }
    result.a_rmse = Rmse(errors[0]); result.b_rmse = Rmse(errors[1]); result.c_rmse = Rmse(errors[2]);
    std::array<std::optional<double> *, 3> biases{&result.a_bias, &result.b_bias, &result.c_bias};
    std::array<std::optional<double> *, 3> medians{&result.a_median_abs, &result.b_median_abs, &result.c_median_abs};
    std::array<std::optional<double> *, 3> p95s{&result.a_p95_abs, &result.b_p95_abs, &result.c_p95_abs};
    std::array<std::optional<double> *, 3> maxima{&result.a_max_abs, &result.b_max_abs, &result.c_max_abs};
    for (std::size_t parameter = 0; parameter < 3; ++parameter)
    {
        double sum{};
        for (const double error : errors[parameter]) sum += error;
        *biases[parameter] = sum / static_cast<double>(errors[parameter].size());
        *medians[parameter] = Median(absolute[parameter]);
        *p95s[parameter] = Quantile(absolute[parameter], 0.95);
        *maxima[parameter] = *std::max_element(absolute[parameter].begin(), absolute[parameter].end());
    }
}

std::optional<ProfileMap> ObservationACProfiles(const core::JointProblem & problem,
    const rg::JointAnalysisResult & target_result, const detail::HaloMedianBuildResult & geometry_medians,
    const detail::HaloFullABCPreclassification & preclassification,
    const Synthetic & synthetic, std::string & reason)
{
    const auto & input = problem.OriginalInput();
    const auto global = geometry_medians.snapshot.global_median;
    const auto geometry_profiles = GeometryProfiles(geometry_medians, synthetic.model);
    std::unordered_map<int, std::size_t> original_index;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
        original_index.emplace(std::stoi(input.atom_ids[atom]), atom);
    std::vector<const detail::HaloFullABCClassificationRecord *> free_halos;
    std::vector<double> response = input.observations;
    for (std::size_t atom = 0; atom < input.selection_domain->target_indices.size(); ++atom)
    {
        const auto contributor = input.selection_domain->target_indices[atom];
        const auto point = TargetPoint(target_result, contributor);
        if (!point) { reason = "target-endpoint-missing"; return std::nullopt; }
        for (const auto & support : input.support[contributor])
        {
            const auto basis = joint::EvaluateKernel(support.squared_distance, (*point)[1], 2.5);
            response[support.row] -= (*point)[0] * basis.gaussian + (*point)[2] * basis.charge;
        }
    }
    for (const auto & halo : preclassification.halos)
    {
        if (halo.role == detail::HaloFullABCRole::FreeFullABC)
        {
            free_halos.push_back(&halo);
            continue;
        }
        if (!halo.original_contributor_index) continue;
        const auto & profile = global;
        const auto contributor = *halo.original_contributor_index;
        for (const auto row : halo.original_support_rows)
        {
            const auto support = std::find_if(input.support[contributor].begin(), input.support[contributor].end(),
                [row](const auto & entry) { return entry.row == row; });
            if (support == input.support[contributor].end()) continue;
            const auto basis = joint::EvaluateKernel(support->squared_distance, profile[1], 2.5);
            response[row] -= profile[0] * basis.gaussian + profile[2] * basis.charge;
        }
    }
    if (free_halos.empty()) { reason = "no-free-halo-columns"; return ProfileMap{}; }
    const auto columns = static_cast<Eigen::Index>(2 * free_halos.size());
    const auto rows = static_cast<Eigen::Index>(input.observations.size());
    std::vector<Eigen::Triplet<double>> triplets;
    std::size_t support_count{};
    for (const auto * halo : free_halos) support_count += halo->original_support_rows.size();
    triplets.reserve(2 * support_count);
    for (std::size_t local = 0; local < free_halos.size(); ++local)
    {
        const int id = std::stoi(free_halos[local]->atom_id);
        const auto found = original_index.find(id);
        if (found == original_index.end()) { reason = "free-halo-original-mapping-missing"; return std::nullopt; }
        const auto profile = geometry_profiles.find(id);
        if (profile == geometry_profiles.end()) { reason = "geometry-width-unavailable"; return std::nullopt; }
        const auto contributor = found->second;
        for (const auto row : free_halos[local]->original_support_rows)
        {
            const auto support = std::find_if(input.support[contributor].begin(), input.support[contributor].end(),
                [row](const auto & entry) { return entry.row == row; });
            if (support == input.support[contributor].end()) continue;
            const auto basis = joint::EvaluateKernel(support->squared_distance, profile->second[1], 2.5);
            triplets.emplace_back(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(2 * local), basis.gaussian);
            triplets.emplace_back(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(2 * local + 1), basis.charge);
        }
    }
    joint::Sparse design(rows, columns);
    design.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::VectorXd y(rows), weights = Eigen::VectorXd::Ones(rows);
    for (Eigen::Index row = 0; row < rows; ++row) y(row) = response[static_cast<std::size_t>(row)];
    const auto solved = joint::SolveLinear(design, y, weights, false, true);
    if (!solved.valid)
    {
        reason = "observation-ac-solver:" + solved.reason;
        return std::nullopt;
    }
    ProfileMap result;
    for (const auto & halo : preclassification.halos)
    {
        const int id = std::stoi(halo.atom_id);
        if (halo.role == detail::HaloFullABCRole::FixedGlobalTargetMedian)
            result.emplace(id, global);
    }
    for (std::size_t local = 0; local < free_halos.size(); ++local)
    {
        const int id = std::stoi(free_halos[local]->atom_id);
        const auto width_source = geometry_profiles.find(id);
        if (width_source == geometry_profiles.end()) continue;
        result[id] = Profile{solved.beta(static_cast<Eigen::Index>(2 * local)), width_source->second[1],
            solved.beta(static_cast<Eigen::Index>(2 * local + 1))};
    }
    (void)synthetic;
    return result;
}

struct SignalScore
{
    std::optional<double> all, shared, domain, maximum;
};

SignalScore ScoreSignal(const core::JointProblem & problem, const Synthetic & synthetic,
    const ProfileMap & profiles, const rg::JointAnalysisResult & target_result)
{
    const auto & input = problem.OriginalInput();
    std::unordered_map<int, std::size_t> original_index;
    for (std::size_t i = 0; i < input.atom_ids.size(); ++i) original_index.emplace(std::stoi(input.atom_ids[i]), i);
    std::vector<double> true_signal(input.observations.size(), 0.0), estimate(input.observations.size(), 0.0);
    std::vector<double> estimated_total(input.observations.size(), 0.0);
    std::vector<std::size_t> shared_count(input.observations.size(), 0);
    for (const auto & halo : synthetic.halos)
    {
        const auto found = original_index.find(halo.id);
        if (found == original_index.end()) continue;
        const auto & support = input.support[found->second];
        std::set<std::size_t> visited;
        for (const auto & entry : support)
        {
            if (!visited.insert(entry.row).second) continue;
            const auto true_basis = joint::EvaluateKernel(entry.squared_distance, halo.truth[1], 2.5);
            true_signal[entry.row] += halo.truth[0] * true_basis.gaussian + halo.truth[2] * true_basis.charge;
            ++shared_count[entry.row];
            const auto profile = profiles.find(halo.id);
            if (profile == profiles.end()) continue;
            const auto estimate_basis = joint::EvaluateKernel(entry.squared_distance, profile->second[1], 2.5);
            estimate[entry.row] += profile->second[0] * estimate_basis.gaussian + profile->second[2] * estimate_basis.charge;
        }
    }
    for (const auto target : input.selection_domain->target_indices)
    {
        const auto fitted = TargetPoint(target_result, target);
        if (!fitted) continue;
        for (const auto & entry : input.support[target])
        {
            const auto fitted_basis = joint::EvaluateKernel(entry.squared_distance, (*fitted)[1], 2.5);
            estimated_total[entry.row] += (*fitted)[0] * fitted_basis.gaussian + (*fitted)[2] * fitted_basis.charge;
        }
    }
    for (std::size_t row = 0; row < estimated_total.size(); ++row) estimated_total[row] += estimate[row];
    std::vector<double> domain_errors;
    std::vector<double> all_errors, shared_errors;
    double maximum{};
    for (std::size_t row = 0; row < input.observations.size(); ++row)
    {
        const double error = estimate[row] - true_signal[row];
        all_errors.push_back(error);
        if (shared_count[row] >= 2) shared_errors.push_back(error);
        maximum = std::max(maximum, std::abs(error));
        domain_errors.push_back(estimated_total[row] - input.observations[row]);
    }
    return {Rmse(all_errors), Rmse(shared_errors), Rmse(domain_errors), maximum};
}

std::size_t PeakRssBytes()
{
#ifdef _WIN32
    return 0;
#else
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
#endif
#endif
}

json::value Number(const std::optional<double> & value)
{
    return value && std::isfinite(*value) ? json::value(*value) : json::value(nullptr);
}

json::object MethodJson(const MethodResult & result)
{
    return {
        {"available", result.available}, {"reason", result.reason}, {"seconds", result.seconds},
        {"halo_A_rmse", Number(result.a_rmse)}, {"halo_B_rmse", Number(result.b_rmse)},
        {"halo_C_rmse", Number(result.c_rmse)}, {"halo_A_bias", Number(result.a_bias)},
        {"halo_B_bias", Number(result.b_bias)}, {"halo_C_bias", Number(result.c_bias)},
        {"halo_A_median_abs", Number(result.a_median_abs)}, {"halo_B_median_abs", Number(result.b_median_abs)},
        {"halo_C_median_abs", Number(result.c_median_abs)}, {"halo_A_p95_abs", Number(result.a_p95_abs)},
        {"halo_B_p95_abs", Number(result.b_p95_abs)}, {"halo_C_p95_abs", Number(result.c_p95_abs)},
        {"halo_A_max_abs", Number(result.a_max_abs)}, {"halo_B_max_abs", Number(result.b_max_abs)},
        {"halo_C_max_abs", Number(result.c_max_abs)}, {"halo_signal_rmse", Number(result.halo_signal_rmse)},
        {"shared_support_rmse", Number(result.shared_signal_rmse)},
        {"original_domain_prediction_rmse", Number(result.domain_signal_rmse)},
        {"maximum_signal_error", Number(result.max_signal_error)},
        {"peeling_rmse", Number(result.peeling_rmse)}, {"tail_ratio_mae", Number(result.tail_ratio_mae)},
        {"local_A_rmse", Number(result.local_a_rmse)}, {"local_B_rmse", Number(result.local_b_rmse)},
        {"local_qualified_count", result.local_qualified_count},
        {"local_sample_count_total", result.local_sample_count_total},
        {"local_effective_sample_count_total", result.local_effective_sample_count_total},
        {"local_effective_sample_count_mean", Number(result.local_effective_sample_mean)},
        {"local_effective_sample_count_min", result.local_effective_sample_min ?
            json::value(*result.local_effective_sample_min) : json::value(nullptr)},
        {"local_information_count", result.local_information_count},
        {"local_information_trace_mean", Number(result.local_information_trace_mean)},
        {"group_posterior_A_error", Number(result.group_a_error)},
        {"group_posterior_B_error", Number(result.group_b_error)},
        {"group_alpha_g", Number(result.group_alpha_g)},
        {"group_mdpde_A", Number(result.group_mdpde_a)},
        {"group_mdpde_B", Number(result.group_mdpde_b)},
        {"group_posterior_A_rmse", Number(result.group_posterior_a_rmse)},
        {"group_posterior_B_rmse", Number(result.group_posterior_b_rmse)},
        {"group_prior_uncertainty_A", Number(result.group_prior_uncertainty_a)},
        {"group_prior_uncertainty_B", Number(result.group_prior_uncertainty_b)},
        {"group_member_uncertainty_A_mean", Number(result.group_member_uncertainty_a_mean)},
        {"group_member_uncertainty_B_mean", Number(result.group_member_uncertainty_b_mean)},
        {"group_outlier_count", result.group_outlier_count},
        {"target_endpoint_drift_max_abs", Number(result.target_endpoint_drift_max)},
        {"peeling_available", result.peeling_available}, {"peeling_total", result.peeling_total},
        {"local_available", result.local_available}, {"group_members", result.group_members}
    };
}

json::value FailureBlockJson(const std::optional<detail::HaloFullABCBlockFailureDiagnostic> & failure)
{
    if(!failure) return nullptr;
    json::array core_ids,active_amplitudes,column_norms,design_singulars,widths,width_derivative_norms,
        projected_width_singulars,projected_width_norms,eta_steps;
    for(const auto & id:failure->core_atom_ids) core_ids.emplace_back(id);
    for(const auto index:failure->active_amplitude_indices) active_amplitudes.emplace_back(index);
    for(const auto value:failure->local_column_norms) column_norms.emplace_back(value);
    for(const auto value:failure->local_design_singular_values) design_singulars.emplace_back(value);
    for(const auto value:failure->initial_widths) widths.emplace_back(value);
    for(const auto value:failure->width_derivative_column_norms) width_derivative_norms.emplace_back(value);
    for(const auto value:failure->projected_width_singular_values) projected_width_singulars.emplace_back(value);
    for(const auto value:failure->projected_width_column_norms) projected_width_norms.emplace_back(value);
    for(const auto value:failure->proposed_eta_step) eta_steps.emplace_back(value);
    const auto optional_size=[](const std::optional<std::size_t> & value) {
        return value ? json::value(*value) : json::value(nullptr);
    };
    return json::object{
        {"component_index",failure->component_index},{"component_id",failure->component_id},
        {"sweep",failure->sweep_index},{"block",failure->block_index},{"status",failure->status},
        {"reason",failure->reason},{"core_atom_ids",std::move(core_ids)},
        {"affected_rows",failure->affected_rows},{"local_design_rows",failure->local_design_rows},
        {"local_design_columns",failure->local_design_columns},{"local_design_rank",optional_size(failure->local_design_rank)},
        {"local_design_rank_threshold",Number(failure->local_design_rank_threshold)},
        {"local_design_singular_values",std::move(design_singulars)},
        {"condition_indicator",Number(failure->condition_indicator)},
        {"local_column_norms",std::move(column_norms)},
        {"initial_profile_valid",failure->initial_profile_valid},
        {"initial_profile_trusted",failure->initial_profile_trusted},
        {"initial_profile_reason",failure->initial_profile_reason},
        {"local_search_stop_reason",failure->local_search_stop_reason},
        {"active_amplitude_indices",std::move(active_amplitudes)},
        {"first_invalid_coefficient",optional_size(failure->first_invalid_coefficient)},
        {"first_invalid_atom_id",failure->first_invalid_atom_id},
        {"first_invalid_parameter",failure->first_invalid_parameter},
        {"solver_valid",failure->solver_valid},{"solver_reason",failure->solver_reason},
        {"primal_feasibility",Number(failure->primal_feasibility)},
        {"scaled_projected_kkt",Number(failure->scaled_projected_kkt)},
        {"local_response_norm",Number(failure->local_response_norm)},
        {"observation_scale",Number(failure->observation_scale)},
        {"initial_widths",std::move(widths)},{"width_upper_bound",Number(failure->width_upper_bound)},
        {"width_derivative_column_norms",std::move(width_derivative_norms)},
        {"width_gradient_inf",Number(failure->width_gradient_inf)},
        {"projected_width_column_norms",std::move(projected_width_norms)},
        {"projected_width_singular_values",std::move(projected_width_singulars)},
        {"projected_width_rank_threshold",Number(failure->projected_width_rank_threshold)},
        {"width_step_inf",Number(failure->width_step_inf)},
        {"proposed_eta_step",std::move(eta_steps)},
        {"lm_subproblem_status",failure->lm_subproblem_status},
        {"lm_actual_decrease",Number(failure->lm_actual_decrease)},
        {"lm_predicted_decrease",Number(failure->lm_predicted_decrease)},
        {"lm_ratio",Number(failure->lm_ratio)},
        {"weak_width_step_fraction",Number(failure->weak_width_step_fraction)},
        {"width_bound_active",failure->width_bound_active},
        {"local_objective_before",Number(failure->local_objective_before)},
        {"local_objective_after",Number(failure->local_objective_after)},
        {"global_objective_before",Number(failure->global_objective_before)},
        {"global_objective_after",Number(failure->global_objective_after)},
        {"global_replay_delta",Number(failure->global_replay_delta)},
        {"local_global_delta_error",Number(failure->local_global_delta_error)},
        {"local_final_ac_kkt",Number(failure->local_final_ac_kkt)},
        {"global_ac_kkt",Number(failure->global_ac_kkt)},
        {"global_width_gradient_inf_norm",Number(failure->global_width_gradient_inf_norm)},
        {"replay_prediction_discrepancy",Number(failure->replay_prediction_discrepancy)},
        {"replay_gradient_discrepancy",Number(failure->replay_gradient_discrepancy)},
        {"replay_kkt_discrepancy",Number(failure->replay_kkt_discrepancy)}
    };
}

json::array ComponentSearchesJson(const std::vector<detail::HaloFullABCComponentSearchDiagnostic> & searches)
{
    json::array output;
    for(const auto & search:searches)
    {
        json::array blocks,candidates;
        json::array local_conditioning;
        for(const auto & candidate:search.candidate_profiles)
            candidates.emplace_back(json::object{{"atom_id",candidate.atom_id},
                {"parameters",json::array{candidate.parameters[0],candidate.parameters[1],candidate.parameters[2]}}});
        for(const auto & block:search.blocks)
            blocks.emplace_back(json::object{{"sweep",block.sweep},{"block",block.block},
                {"affected_rows",block.affected_rows},{"status",block.status},{"reason",block.reason},
                {"local_search_stop_reason",block.local_search_stop_reason},
                {"profile_evaluations",block.profile_evaluations},{"accepted_updates",block.accepted_updates},
                {"local_objective_before",Number(block.local_objective_before)},
                {"local_objective_after",Number(block.local_objective_after)},
                {"global_objective_before",Number(block.global_objective_before)},
                {"global_objective_after",Number(block.global_objective_after)},
                {"global_replay_delta",Number(block.global_replay_delta)},
                {"local_global_delta_error",Number(block.local_global_delta_error)},
                {"local_final_ac_kkt",Number(block.local_final_ac_kkt)}});
        for(const auto & condition:search.local_conditioning)
            local_conditioning.emplace_back(json::object{{"block_id",condition.block_id},
                {"status",condition.status},{"reason",condition.reason},
                {"core_atom_count",condition.core_atom_count},{"affected_row_count",condition.affected_row_count},
                {"active_amplitude_count",condition.active_amplitude_count},
                {"free_ac_columns",condition.free_ac_columns},{"free_ac_rank",condition.free_ac_rank ?
                    json::value(*condition.free_ac_rank) : json::value(nullptr)},
                {"free_ac_rank_threshold",Number(condition.free_ac_rank_threshold)},
                {"free_ac_sigma_min",Number(condition.free_ac_sigma_min)},
                {"free_ac_sigma_max",Number(condition.free_ac_sigma_max)},
                {"free_ac_condition",Number(condition.free_ac_condition)},
                {"free_ac_conditioning",condition.free_ac_conditioning},
                {"projected_width_rank",condition.projected_width_rank ?
                    json::value(*condition.projected_width_rank) : json::value(nullptr)},
                {"projected_width_rank_threshold",Number(condition.projected_width_rank_threshold)},
                {"projected_width_sigma_min",Number(condition.projected_width_sigma_min)},
                {"projected_width_sigma_max",Number(condition.projected_width_sigma_max)},
                {"projected_width_condition",Number(condition.projected_width_condition)},
                {"projected_width_conditioning",condition.projected_width_conditioning},
                {"local_widths",condition.local_widths},{"projected_width_column_norms",condition.projected_width_column_norms},
                {"local_ac_projected_kkt",Number(condition.local_ac_projected_kkt)},
                {"local_projected_width_kkt",Number(condition.local_projected_width_kkt)}});
        output.emplace_back(json::object{{"component_id",search.component_id},{"reason",search.reason},
            {"search_converged",search.search_converged},{"sweep_count",search.sweep_count},
            {"objective",Number(search.objective)},{"global_ac_kkt",Number(search.global_ac_kkt)},
            {"global_width_gradient_inf_norm",Number(search.global_width_gradient_inf_norm)},
            {"eta_change_inf",Number(search.eta_change_inf)},
            {"objective_replay_error",Number(search.objective_replay_error)},
            {"cache_replay_error",Number(search.cache_replay_error)},
            {"endpoint_verified",search.endpoint_verified},{"endpoint_qualified",search.endpoint_qualified},
            {"endpoint_failure",search.endpoint_failure},{"trust_failure",search.trust_failure},
            {"local_conditioning",std::move(local_conditioning)},
            {"candidate_profiles",std::move(candidates)},{"blocks",std::move(blocks)}});
    }
    return output;
}

json::array HaloObservationGeometryJson(const detail::HaloFullABCReconstructionResult & result,
    const core::JointProblem & problem,const rg::ModelObject & model)
{
    const auto & input=problem.OriginalInput();
    const auto count=result.halos.size();
    std::vector<std::size_t> parent(count);
    std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](std::size_t index) {
        while(parent[index]!=index)
        {
            parent[index]=parent[parent[index]];
            index=parent[index];
        }
        return index;
    };
    std::map<std::size_t,std::vector<std::size_t>> row_users;
    for(std::size_t halo=0;halo<count;++halo)
        if(result.halos[halo].classification.role==detail::HaloFullABCRole::FreeFullABC)
            for(const auto row:result.halos[halo].classification.original_support_rows)
                row_users[row].push_back(halo);
    std::map<std::pair<std::size_t,std::size_t>,std::vector<std::size_t>> shared_by_pair;
    for(const auto & [row,users]:row_users)
        for(std::size_t left=0;left<users.size();++left)
            for(std::size_t right=left+1;right<users.size();++right)
            {
                const auto a=std::min(users[left],users[right]);
                const auto b=std::max(users[left],users[right]);
                shared_by_pair[{a,b}].push_back(row);
            }
    for(const auto & [pair,rows]:shared_by_pair)
    {
        (void)rows;
        const auto left=root(pair.first),right=root(pair.second);
        if(left!=right) parent[std::max(left,right)]=std::min(left,right);
    }
    std::map<std::size_t,std::set<std::size_t>> component_rows;
    std::map<std::size_t,std::vector<std::size_t>> component_halos;
    for(std::size_t halo=0;halo<count;++halo)
    {
        if(result.halos[halo].classification.role!=detail::HaloFullABCRole::FreeFullABC) continue;
        const auto component=root(halo);
        component_halos[component].push_back(halo);
        component_rows[component].insert(result.halos[halo].classification.original_support_rows.begin(),
            result.halos[halo].classification.original_support_rows.end());
    }

    json::array output;
    for(std::size_t halo=0;halo<count;++halo)
    {
        const auto & classification=result.halos[halo].classification;
        const int id=std::stoi(classification.atom_id);
        const auto * atom=model.FindAtomPtr(id);
        const auto position=atom ? atom->GetPosition() : Position{};
        json::array rows,shared_rows,component_ids;
        if(classification.original_contributor_index)
        {
            std::map<std::size_t,double> radial_by_row;
            for(const auto & support:input.support.at(*classification.original_contributor_index))
                radial_by_row.emplace(support.row,support.squared_distance);
            for(const auto & [row,square]:radial_by_row)
                rows.emplace_back(json::object{{"row",row},{"row_id",input.row_ids.at(row)},
                    {"squared_radial_distance",square}});
        }
        for(const auto & [row,users]:row_users)
            if(users.size()>1 && std::find(users.begin(),users.end(),halo)!=users.end())
                shared_rows.emplace_back(json::object{{"row",row},{"row_id",input.row_ids.at(row)}});
        std::size_t local_effective_observations{};
        if(classification.role==detail::HaloFullABCRole::FreeFullABC)
        {
            const auto component=root(halo);
            local_effective_observations=component_rows.at(component).size();
            for(const auto member:component_halos.at(component))
                component_ids.emplace_back(result.halos[member].classification.atom_id);
        }
        output.emplace_back(json::object{{"halo_atom_id",classification.atom_id},
            {"position",json::array{position[0],position[1],position[2]}},
            {"distinct_support_rows",classification.distinct_original_observations},
            {"role",classification.role==detail::HaloFullABCRole::FreeFullABC ? "FreeFullABC" :
                "PreFixedGlobalMedian"},{"support",std::move(rows)},
            {"free_component_atom_ids",std::move(component_ids)},
            {"local_effective_observation_count",local_effective_observations},
            {"shared_rows",std::move(shared_rows)}});
    }
    return output;
}

void RunPostprocessing(const Synthetic & synthetic, const core::JointProblem & problem,
    const rg::JointAnalysisResult & snapshot, const detail::HaloMedianBuildResult * geometry,
    const detail::HaloFullABCReconstructionResult * fullabc, const std::string & run_id,
    MethodResult & result)
{
    rg::ModelObject model = synthetic.model;
    rg::data_internal::ApplyJointStageEstimates(model, snapshot, run_id);
    double target_drift{};
    for (const auto target : problem.OriginalInput().selection_domain->target_indices)
    {
        const int id = std::stoi(problem.OriginalInput().atom_ids.at(target));
        const auto expected = TargetPoint(snapshot, target);
        const auto actual = rg::AtomLocalPotentialView::For(*model.FindAtomPtr(id)).GetStageEstimate(
            rg::FittingStage::Second).point;
        if (!expected || !actual) throw std::runtime_error("Frozen Target endpoint is unavailable during postprocessing.");
        target_drift = std::max(target_drift, std::abs(actual->GetAmplitude() - (*expected)[0]));
        target_drift = std::max(target_drift, std::abs(actual->GetWidth() - (*expected)[1]));
        target_drift = std::max(target_drift, std::abs(actual->GetOffset() - (*expected)[2]));
    }
    result.target_endpoint_drift_max = target_drift;
    if (target_drift != 0.0) throw std::runtime_error("A reconstruction method changed the frozen Target endpoint.");
    detail::JointPostprocessingMetrics metrics;
    detail::JointPostprocessingContext context(*synthetic.map, model, problem, &snapshot, &metrics,
        geometry, fullabc);
    auto peelings = context.BuildPostFitPeelingSamples();
    std::vector<std::pair<int, rg::PostFitPeelingResult>> peeling_updates;
    for (const auto & [id, peeling] : peelings) peeling_updates.emplace_back(id, peeling);
    model.EditAnalysis().ApplyPostFitPeelingsBatch(std::move(peeling_updates));

    std::vector<double> peeling_errors, tail_errors, local_a_errors, local_b_errors;
    std::vector<double> local_information_traces, group_posterior_a_errors, group_posterior_b_errors;
    std::vector<double> group_uncertainty_a, group_uncertainty_b;
    std::vector<std::pair<int, rg::PostJointLocalResult>> local_updates;
    core::FitOptions options;
    options.thread_size = 1;
    options.quiet_mode = true;
    options.exclude_hydrogen = true;
    for (const auto * atom : model.GetSelectedAtoms())
    {
        const int id = atom->GetSerialID();
        const auto peeled = peelings.find(id);
        if (peeled == peelings.end()) continue;
        const auto raw = rg::AtomLocalPotentialView::For(*atom).GetRawSamplingEntriesView();
        result.peeling_total += peeled->second.samples.size();
        for (std::size_t sample = 0; sample < peeled->second.samples.size() && sample < raw.size(); ++sample)
        {
            if (!peeled->second.samples[sample].response) continue;
            ++result.peeling_available;
            const SamplingPointList point{raw[sample].point};
            const auto target_truth = detail::BuildLocalPotentialSampleList(*synthetic.target_only_map, point);
            if (!target_truth.empty())
                peeling_errors.push_back(*peeled->second.samples[sample].response - target_truth.front().response);
        }
        const auto assessment = rg::AtomLocalPotentialView::For(*atom).AssessLocalFittingPeelingRatio(1.0, 2.0);
        if (assessment.ratio)
        {
            double raw_target{}, raw_sum{};
            std::size_t count{};
            for (std::size_t sample = 0; sample < raw.size(); ++sample)
            {
                if (raw[sample].point.distance < 1.0 || raw[sample].point.distance > 2.0 ||
                    !peeled->second.samples[sample].response) continue;
                const SamplingPointList point{raw[sample].point};
                const auto target_truth = detail::BuildLocalPotentialSampleList(*synthetic.target_only_map, point);
                if (target_truth.empty()) continue;
                raw_target += target_truth.front().response;
                raw_sum += raw[sample].response;
                ++count;
            }
            if (count && std::abs(raw_sum) > 1e-12)
                tail_errors.push_back(*assessment.ratio - (1.0 - raw_target / raw_sum));
        }
        auto local = detail::FitPostJointLocalTargetForResearch(context, model, id,
            peeled->second, options, metrics);
        result.local_sample_count_total += local.sample_count;
        if (local.status == rg::EvidenceStatus::Available)
        {
            ++result.local_qualified_count;
            result.local_effective_sample_count_total += local.effective_sample_count;
            result.local_effective_sample_min = result.local_effective_sample_min ?
                std::min(*result.local_effective_sample_min, local.effective_sample_count) :
                std::optional<std::size_t>{local.effective_sample_count};
        }
        if (local.information)
        {
            ++result.local_information_count;
            local_information_traces.push_back(local.information->precision.trace());
        }
        if (local.status == rg::EvidenceStatus::Available)
        {
            const auto truth = TargetTruth(synthetic,id);
            if (truth)
            {
                local_a_errors.push_back(local.fit.mdpde.GetModel().GetAmplitude() - (*truth)[0]);
                local_b_errors.push_back(local.fit.mdpde.GetModel().GetWidth() - (*truth)[1]);
            }
        }
        local_updates.emplace_back(id, std::move(local));
    }
    model.EditAnalysis().ApplyPostJointLocalResultsBatch(std::move(local_updates));
    detail::RunJointGroupPotentialFitting(model, options, &metrics);
    result.peeling_rmse = Rmse(peeling_errors);
    result.tail_ratio_mae = tail_errors.empty() ? std::nullopt :
        std::optional<double>{std::accumulate(tail_errors.begin(), tail_errors.end(), 0.0,
            [](double sum, double value) { return sum + std::abs(value); }) / static_cast<double>(tail_errors.size())};
    result.local_a_rmse = Rmse(local_a_errors);
    result.local_b_rmse = Rmse(local_b_errors);
    result.local_available = local_a_errors.size();
    if (result.local_qualified_count)
        result.local_effective_sample_mean = static_cast<double>(result.local_effective_sample_count_total) /
            static_cast<double>(result.local_qualified_count);
    if (!local_information_traces.empty())
        result.local_information_trace_mean = std::accumulate(local_information_traces.begin(),
            local_information_traces.end(), 0.0) / static_cast<double>(local_information_traces.size());
    if (const auto keys = model.GetAnalysisView().CollectAtomGroupKeys(); !keys.empty())
    {
        const auto summary = model.GetAnalysisView().GetGroupParameterSummary(keys.front());
        if (summary && summary->inference)
        {
            result.group_members = summary->member_ids.size();
            result.group_alpha_g = summary->inference->alpha_g;
            result.group_mdpde_a = summary->inference->mdpde.GetAmplitude();
            result.group_mdpde_b = summary->inference->mdpde.GetWidth();
            result.group_prior_uncertainty_a =
                summary->inference->prior.GetStandardDeviationModel().GetAmplitude();
            result.group_prior_uncertainty_b =
                summary->inference->prior.GetStandardDeviationModel().GetWidth();
            for (std::size_t index = 0; index < summary->member_ids.size() &&
                index < summary->inference->member_results.size(); ++index)
            {
                const int member_id = summary->member_ids[index];
                const auto & member = summary->inference->member_results[index];
                result.group_outlier_count += member.is_outlier;
                const auto truth = TargetTruth(synthetic,member_id);
                if (!truth) continue;
                group_posterior_a_errors.push_back(member.posterior.GetModel().GetAmplitude() - (*truth)[0]);
                group_posterior_b_errors.push_back(member.posterior.GetModel().GetWidth() - (*truth)[1]);
                group_uncertainty_a.push_back(member.posterior.GetStandardDeviationModel().GetAmplitude());
                group_uncertainty_b.push_back(member.posterior.GetStandardDeviationModel().GetWidth());
            }
            std::vector<double> group_truth_a,group_truth_b;
            for (const int member_id : summary->member_ids)
                if (const auto truth = TargetTruth(synthetic,member_id))
                {
                    group_truth_a.push_back((*truth)[0]);
                    group_truth_b.push_back((*truth)[1]);
                }
            if (!group_truth_a.empty())
            {
                result.group_a_error = summary->inference->prior.GetModel().GetAmplitude() - Median(group_truth_a);
                result.group_b_error = summary->inference->prior.GetModel().GetWidth() - Median(group_truth_b);
            }
        }
    }
    result.group_posterior_a_rmse = Rmse(group_posterior_a_errors);
    result.group_posterior_b_rmse = Rmse(group_posterior_b_errors);
    if (!group_uncertainty_a.empty())
        result.group_member_uncertainty_a_mean = std::accumulate(group_uncertainty_a.begin(), group_uncertainty_a.end(), 0.0) /
            static_cast<double>(group_uncertainty_a.size());
    if (!group_uncertainty_b.empty())
        result.group_member_uncertainty_b_mean = std::accumulate(group_uncertainty_b.begin(), group_uncertainty_b.end(), 0.0) /
            static_cast<double>(group_uncertainty_b.size());
}

json::object RunScenario(const Scenario & scenario, std::ofstream & summary_csv, std::ofstream & halo_csv)
{
    const auto started = Clock::now();
    auto synthetic = scenario.real_geometry ? MakeRealGeometrySynthetic(scenario) : MakeSynthetic(scenario);
    const auto problem = core::BuildPureHaloSupportExclusionProblem(*synthetic.map, synthetic.model);
    if (scenario.real_full_structure)
    {
        const auto & input = problem.OriginalInput();
        const auto target_count = input.selection_domain->target_indices.size();
        const auto supported_halo_count = input.atom_ids.size() - target_count;
        if (target_count != 1418 || supported_halo_count != 659 || input.observations.size() != 229783)
            throw std::runtime_error("The 6Z6U Partial Selection no longer matches the recorded chain-A baseline.");
    }
    const auto workset = core::detail::MakeJointFittingWorkset(synthetic.model, problem);
    core::RunPotentialSamplingWorkflow(*synthetic.map, synthetic.model, workset.contributors,
        SphereSamplingMethod::FibonacciDeterministic, 1, "P3 FirstStage", true);
    core::FitOptions options;
    options.thread_size = 1;
    options.quiet_mode = true;
    options.exclude_hydrogen = true;
    auto initialization = core::detail::RunJointFirstStageInitializationFromPreparedSamples(
        synthetic.model, workset, options);
    const auto fit_start = Clock::now();
    const auto fit = core::FitJointComponentsInitialized(problem, initialization);
    const double target_joint_seconds = std::chrono::duration<double>(Clock::now() - fit_start).count();
    const auto snapshot = core::CaptureJointAnalysisResult(fit);
    json::object case_json{
        {"name", scenario.name},
        {"halo_center_distance_requested", scenario.real_geometry ? json::value(nullptr) : json::value(scenario.halo_distance)},
        {"contrast", scenario.contrast}, {"wide_charge_contrast", scenario.wide_charge_contrast},
        {"noise_amplitude", scenario.noise}, {"map_edge_halo", scenario.map_edge_halo},
        {"noise_seed", scenario.real_geometry ? json::value(scenario.noise_seed) : json::value(nullptr)},
        {"real_geometry", scenario.real_geometry},
        {"requires_short_bonded_C_O_or_C_N_halo_pair",scenario.require_short_bonded_halo_pair},
        {"map_normalization_divisor", synthetic.map_normalization_divisor ?
            json::value(*synthetic.map_normalization_divisor) : json::value(nullptr)},
        {"target_sequence_range", scenario.real_geometry ? json::value(json::array{scenario.sequence_start,scenario.sequence_end}) : json::value(nullptr)},
        {"extra_free_halo_observations", scenario.extra_free_halo_observations},
        {"target_joint_qualified", rg::IsQualifiedJointAnalysisResult(snapshot)},
        {"target_joint_seconds", target_joint_seconds}, {"first_stage_qualified_width_donors", 0},
        {"first_stage_donor_coverage", nullptr}, {"target_endpoint_drift_max_abs", nullptr},
        {"geometry_kNN4", nullptr}, {"observation_AC_fixed_kNN_width", nullptr},
        {"halo_fullabc", nullptr}, {"oracle_halo_truth", nullptr}, {"geometry_peeling", nullptr},
        {"fullabc_peeling", nullptr}, {"solver_failure_reasons", json::array{}}
    };
    std::size_t qualified_width_donors{};
    for (const auto & atom : initialization.atoms) qualified_width_donors +=
        rg::IsQualifiedFirstStageEstimate(atom.native_status, atom.first_stage_fit_qualified,
            atom.first_stage_qualification, atom.first_stage_recovery);
    case_json["first_stage_qualified_width_donors"] = qualified_width_donors;
    case_json["initialization_valid"] = initialization.valid;
    case_json["initialization_reason"] = initialization.reason;
    case_json["first_stage_donor_coverage"] = initialization.atoms.empty() ? json::value(nullptr) :
        json::value(static_cast<double>(qualified_width_donors) / static_cast<double>(initialization.atoms.size()));
    case_json["fit_search_completed"] = fit.search_completed;
    case_json["fit_runtime_convergence"] = static_cast<int>(fit.RuntimeConvergence());
    case_json["fit_target_runtime_convergence"] = static_cast<int>(fit.TargetRuntimeConvergence());
    case_json["fitting_domain_rows"] = problem.Input().observations.size();
    case_json["fit_components"] = json::array{};
    for (const auto & component : fit.components)
    {
        json::object component_json{{"id", component.id}, {"stop_reason", component.stop_reason},
            {"search_completed", component.search_completed}, {"search_sweeps", component.search_sweeps},
            {"atom_count", component.atoms.size()}, {"row_count", component.rows.size()},
            {"native_status", component.native_status}, {"has_state", component.state.has_value()}};
        component_json["atom_ids"] = json::array{};
        for (const auto atom : component.atoms)
            if (atom < problem.Input().atom_ids.size())
                component_json["atom_ids"].as_array().emplace_back(problem.Input().atom_ids[atom]);
        component_json["evidence"] = json::array{};
        for (const auto & evidence : component.evidence)
            component_json["evidence"].as_array().emplace_back(json::object{
                {"name", evidence.name}, {"status", static_cast<int>(evidence.status)},
                {"scope", static_cast<int>(evidence.scope)}, {"reason", evidence.reason},
                {"value", evidence.value ? json::value(*evidence.value) : json::value(nullptr)},
                {"threshold", evidence.threshold ? json::value(*evidence.threshold) : json::value(nullptr)}});
        case_json["fit_components"].as_array().emplace_back(std::move(component_json));
    }
    case_json["fit_evidence"] = json::array{};
    for (const auto & evidence : fit.evidence)
        case_json["fit_evidence"].as_array().emplace_back(json::object{
            {"name", evidence.name}, {"status", static_cast<int>(evidence.status)},
            {"scope", static_cast<int>(evidence.scope)}, {"reason", evidence.reason},
            {"value", evidence.value ? json::value(*evidence.value) : json::value(nullptr)},
            {"threshold", evidence.threshold ? json::value(*evidence.threshold) : json::value(nullptr)}});
    if (!rg::IsQualifiedJointAnalysisResult(snapshot))
    {
        case_json["target_joint_failure"] = "runtime-or-target-qualification-not-passed";
        return case_json;
    }

    const std::string run_id = "p3-" + scenario.name;
    const auto targets = problem.OriginalInput().selection_domain->target_indices;
    std::vector<Profile> target_endpoint;
    for (const auto target : targets)
    {
        const auto point = TargetPoint(snapshot, target);
        if (!point) throw std::runtime_error("Qualified Target Joint endpoint is missing a selected point.");
        target_endpoint.push_back(*point);
    }
    const auto geometry_started = Clock::now();
    const auto geometry = detail::BuildHaloMedianReconstruction(synthetic.model, problem, snapshot, run_id);
    const double geometry_estimation_seconds = std::chrono::duration<double>(Clock::now() - geometry_started).count();
    const auto fullabc_started = Clock::now();
    const bool capture_solver_diagnostics=scenario.name=="three_observation_free_stress" ||
        scenario.name=="four_observation_free_stress";
    const auto fullabc = detail::EstimateHaloFullABCJoint(
        synthetic.model,problem,snapshot,run_id,capture_solver_diagnostics);
    const double fullabc_estimation_seconds = std::chrono::duration<double>(Clock::now() - fullabc_started).count();
    const auto preclassification = detail::BuildHaloFullABCPreclassification(synthetic.model, problem, snapshot);
    case_json["free_halo_count"] = fullabc.free_halo_count;
    case_json["fixed_halo_count"] = fullabc.fixed_halo_count;
    case_json["fullabc_solver_invoked"] = fullabc.diagnostics.solver_invoked;
    case_json["fullabc_solver_qualified"] = fullabc.diagnostics.solver_qualified;
    case_json["fullabc_solver_status"] = static_cast<int>(fullabc.diagnostics.solver_status);
    case_json["fullabc_unavailable_reason"] = fullabc.unavailable_reason;
    case_json["fullabc_seconds"] = fullabc_estimation_seconds;
    case_json["fullabc_peak_workspace_estimate_bytes"] = fullabc.diagnostics.estimated_peak_memory_bytes ?
        json::value(*fullabc.diagnostics.estimated_peak_memory_bytes) : json::value(nullptr);
    case_json["fullabc_free_support_entries"] = fullabc.diagnostics.free_support_entries;
    case_json["fullabc_zero_sensitivity_rows"] = fullabc.diagnostics.zero_sensitivity_rows;
    case_json["fullabc_original_observation_rows"] = fullabc.diagnostics.original_observation_rows;
    case_json["fullabc_observation_scale"] = Number(fullabc.diagnostics.observation_scale);
    case_json["fullabc_objective"] = Number(fullabc.diagnostics.objective);
    case_json["fullabc_qualified_first_stage_width_donors"] =
        fullabc.diagnostics.qualified_first_stage_width_donors;
    case_json["fullabc_first_stage_width_p95"] = Number(fullabc.diagnostics.first_stage_width_p95);
    case_json["fullabc_width_upper_bound"] = Number(fullabc.diagnostics.width_upper_bound);
    json::array solver_failures;
    for (const auto & evidence : fullabc.diagnostics.evidence)
        if (evidence.status == rg::JointCheckStatus::Failed || evidence.status == rg::JointCheckStatus::Unavailable)
            solver_failures.emplace_back(json::object{{"name", evidence.name}, {"status", static_cast<int>(evidence.status)},
                {"reason", evidence.reason},
                {"value", evidence.value ? json::value(*evidence.value) : json::value(nullptr)},
                {"threshold", evidence.threshold ? json::value(*evidence.threshold) : json::value(nullptr)}});
    case_json["solver_failure_reasons"] = std::move(solver_failures);
    if(capture_solver_diagnostics || scenario.real_geometry)
    {
        case_json["halo_observation_geometry"]=HaloObservationGeometryJson(fullabc,problem,synthetic.model);
        case_json["fullabc_first_failing_block"]=FailureBlockJson(fullabc.diagnostics.first_failing_block);
        case_json["fullabc_component_searches"]=ComponentSearchesJson(fullabc.diagnostics.component_searches);
    }

    std::unordered_map<int, std::size_t> original_index;
    const auto & input = problem.OriginalInput();
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
        original_index.emplace(std::stoi(input.atom_ids[atom]), atom);
    for (auto & halo : synthetic.halos)
        if (const auto found = original_index.find(halo.id); found != original_index.end())
            for (const auto & support : input.support[found->second]) halo.rows.push_back(support.row);
    for (auto & halo : synthetic.halos)
    {
        std::sort(halo.rows.begin(), halo.rows.end());
        halo.rows.erase(std::unique(halo.rows.begin(), halo.rows.end()), halo.rows.end());
    }
    std::vector<const Halo *> geometry_halos;
    for (const auto & halo : synthetic.halos)
        if (!halo.rows.empty() || !scenario.real_full_structure) geometry_halos.push_back(&halo);
    double minimum_halo_distance = std::numeric_limits<double>::infinity();
    std::size_t center_pairs_below_1_5{}, support_overlap_pairs{}, shared_voxel_pairs{};
    std::size_t exact_center_pairs{}, nonbonded_pairs_below_1_0{}, close_bonded_pairs_below_1_5{},
        short_bonded_co_cn_pairs{};
    json::array shared_voxel_counts;
    std::map<std::string, std::size_t> distance_bins{
        {"below_1_5", 0}, {"1_5_to_2", 0}, {"2_to_3", 0}, {"3_to_4", 0}, {"4_to_5", 0}, {"at_least_5", 0}};
    for (std::size_t left = 0; left < geometry_halos.size(); ++left)
        for (std::size_t right = left + 1; right < geometry_halos.size(); ++right)
        {
            const double distance = Distance(geometry_halos[left]->position, geometry_halos[right]->position);
            minimum_halo_distance = std::min(minimum_halo_distance, distance);
            std::string bin = distance < 1.5 ? "below_1_5" : distance < 2.0 ? "1_5_to_2" :
                distance < 3.0 ? "2_to_3" : distance < 4.0 ? "3_to_4" : distance < 5.0 ? "4_to_5" : "at_least_5";
            ++distance_bins[bin];
            center_pairs_below_1_5 += distance < 1.5;
            exact_center_pairs += distance <= 1e-9;
            const bool covalent = synthetic.covalent_pairs.contains(
                OrderedPair(geometry_halos[left]->id,geometry_halos[right]->id));
            if (covalent && distance < 1.5)
            {
                ++close_bonded_pairs_below_1_5;
                const auto * left_atom = synthetic.model.FindAtomPtr(geometry_halos[left]->id);
                const auto * right_atom = synthetic.model.FindAtomPtr(geometry_halos[right]->id);
                if (left_atom && right_atom)
                {
                    const auto left_element = left_atom->GetElement();
                    const auto right_element = right_atom->GetElement();
                    const bool c_o = (left_element == Element::CARBON && right_element == Element::OXYGEN) ||
                        (right_element == Element::CARBON && left_element == Element::OXYGEN);
                    const bool c_n = (left_element == Element::CARBON && right_element == Element::NITROGEN) ||
                        (right_element == Element::CARBON && left_element == Element::NITROGEN);
                    short_bonded_co_cn_pairs += c_o || c_n;
                }
            }
            if (!covalent && distance < 1.0) ++nonbonded_pairs_below_1_0;
            std::vector<std::size_t> intersection;
            std::set_intersection(geometry_halos[left]->rows.begin(), geometry_halos[left]->rows.end(),
                geometry_halos[right]->rows.begin(), geometry_halos[right]->rows.end(),
                std::back_inserter(intersection));
            if (!intersection.empty())
            {
                ++support_overlap_pairs;
                shared_voxel_pairs += intersection.size();
                shared_voxel_counts.emplace_back(intersection.size());
            }
        }
    json::array halo_target_distances;
    for (const auto * halo : geometry_halos)
    {
        double nearest = std::numeric_limits<double>::infinity();
        for (const auto & target : synthetic.target_positions)
            nearest = std::min(nearest, Distance(halo->position, target));
        halo_target_distances.emplace_back(json::object{
            {"atom_id", halo->id}, {"minimum_center_distance", nearest}});
    }
    case_json["halo_halo_minimum_center_distance"] = minimum_halo_distance;
    case_json["halo_halo_center_pairs_below_1_5_angstrom"] = center_pairs_below_1_5;
    case_json["close_bonded_halo_pairs_below_1_5_angstrom"] = close_bonded_pairs_below_1_5;
    case_json["short_bonded_C_O_or_C_N_pairs_below_1_5_angstrom"] = short_bonded_co_cn_pairs;
    case_json["nonbonded_halo_pairs_below_1_0_angstrom"] = nonbonded_pairs_below_1_0;
    case_json["exact_halo_center_pairs"] = exact_center_pairs;
    case_json["shared_voxels_per_overlap_pair"] = std::move(shared_voxel_counts);
    if (exact_center_pairs || nonbonded_pairs_below_1_0)
        throw std::runtime_error("A scene contains coincident or unreasonably close nonbonded Halo centers.");
    if (scenario.require_short_bonded_halo_pair && short_bonded_co_cn_pairs == 0)
        throw std::runtime_error("The required physical short C-O/C-N Halo bond pair was not retained.");
    case_json["halo_target_minimum_center_distances"] = std::move(halo_target_distances);
    case_json["halo_halo_center_distance_bins"] = json::object{};
    for (const auto & [bin, count] : distance_bins) case_json["halo_halo_center_distance_bins"].as_object()[bin] = count;
    case_json["support_overlap_pairs"] = support_overlap_pairs;
    case_json["pairwise_shared_voxel_count"] = shared_voxel_pairs;
    case_json["n_h_distribution"] = json::object{};
    for (const auto & halo : fullabc.halos)
    {
        const auto key = std::to_string(halo.classification.distinct_original_observations);
        auto & value = case_json["n_h_distribution"].as_object()[key];
        value = value.is_int64() ? json::value(value.as_int64() + 1) : json::value(1);
    }
    std::unordered_map<int, std::size_t> observed_by_id;
    for (const auto & halo : fullabc.halos)
        observed_by_id.emplace(std::stoi(halo.classification.atom_id),
            halo.classification.distinct_original_observations);
    const std::size_t expected_halo_103 = scenario.extra_free_halo_observations;
    if (!scenario.real_geometry && ((expected_halo_103 && observed_by_id[103] != expected_halo_103) || observed_by_id[104] != 2 ||
        observed_by_id[105] != 1 || observed_by_id[106] != 0)
       )
        throw std::runtime_error("Synthetic distinct OriginalInput observation counts differ from the scenario contract: " +
            std::to_string(observed_by_id[103]) + "," + std::to_string(observed_by_id[104]) + "," +
            std::to_string(observed_by_id[105]) + "," + std::to_string(observed_by_id[106]));
    case_json["original_domain_rows"] = input.observations.size();
    case_json["original_domain_scale"] = problem.ObservationScale();
    case_json["target_endpoint_drift_max_abs"] = 0.0;

    MethodResult geometry_result;
    geometry_result.available = true;
    geometry_result.profiles.resize(synthetic.halos.size());
    const auto geometry_profiles = GeometryProfiles(geometry, synthetic.model);
    for (std::size_t index = 0; index < synthetic.halos.size(); ++index)
        if (const auto found = geometry_profiles.find(synthetic.halos[index].id); found != geometry_profiles.end())
            geometry_result.profiles[index] = found->second;
        else if (geometry_profiles.empty()) geometry_result.available = false;
    ScoreParameters(geometry_result, synthetic.halos);
    const auto geometry_signal = ScoreSignal(problem, synthetic, geometry_profiles, snapshot);
    geometry_result.halo_signal_rmse = geometry_signal.all;
    geometry_result.shared_signal_rmse = geometry_signal.shared;
    geometry_result.domain_signal_rmse = geometry_signal.domain;
    geometry_result.max_signal_error = geometry_signal.maximum;
    geometry_result.seconds = geometry_estimation_seconds;
    case_json["geometry_kNN4"] = MethodJson(geometry_result);

    std::string ac_reason;
    std::optional<ProfileMap> ac_profiles;
    const auto ac_started = Clock::now();
    if (!scenario.real_full_structure)
        ac_profiles = ObservationACProfiles(problem, snapshot, geometry, preclassification, synthetic, ac_reason);
    else ac_reason = "not-run-for-full-chain-partial-selection";
    MethodResult ac_result;
    ac_result.seconds = std::chrono::duration<double>(Clock::now() - ac_started).count();
    ac_result.available = ac_profiles.has_value();
    ac_result.reason = ac_reason;
    ac_result.profiles.resize(synthetic.halos.size());
    if (ac_profiles)
    {
        for (std::size_t index = 0; index < synthetic.halos.size(); ++index)
            if (const auto found = ac_profiles->find(synthetic.halos[index].id); found != ac_profiles->end())
                ac_result.profiles[index] = found->second;
        ScoreParameters(ac_result, synthetic.halos);
        const auto signal = ScoreSignal(problem, synthetic, *ac_profiles, snapshot);
        ac_result.halo_signal_rmse = signal.all;
        ac_result.shared_signal_rmse = signal.shared;
        ac_result.domain_signal_rmse = signal.domain;
        ac_result.max_signal_error = signal.maximum;
    }
    case_json["observation_AC_fixed_kNN_width"] = MethodJson(ac_result);

    MethodResult full_result;
    full_result.available = fullabc.available;
    full_result.reason = fullabc.unavailable_reason;
    full_result.profiles.resize(synthetic.halos.size());
    const auto full_profiles = FullABCProfiles(fullabc);
    if (fullabc.available)
    {
        for (std::size_t index = 0; index < synthetic.halos.size(); ++index)
            if (const auto found = full_profiles.find(synthetic.halos[index].id); found != full_profiles.end())
                full_result.profiles[index] = found->second;
        ScoreParameters(full_result, synthetic.halos);
        const auto signal = ScoreSignal(problem, synthetic, full_profiles, snapshot);
        full_result.halo_signal_rmse = signal.all;
        full_result.shared_signal_rmse = signal.shared;
        full_result.domain_signal_rmse = signal.domain;
        full_result.max_signal_error = signal.maximum;
        const auto post_start = Clock::now();
        RunPostprocessing(synthetic, problem, snapshot, nullptr, &fullabc, run_id, full_result);
        full_result.seconds = fullabc_estimation_seconds +
            std::chrono::duration<double>(Clock::now() - post_start).count();
        case_json["fullabc_peeling"] = MethodJson(full_result);
    }
    case_json["halo_fullabc"] = MethodJson(full_result);

    if (capture_solver_diagnostics || scenario.real_geometry)
    {
        json::array recovery_results;
        const std::array recovery_methods{
            detail::HaloFullABCRecoveryMethod::FixedWidthACPerHalo,
            detail::HaloFullABCRecoveryMethod::PostFailureGlobalMedian,
            detail::HaloFullABCRecoveryMethod::FixedWidthACPerComponent,
            detail::HaloFullABCRecoveryMethod::IsolateFailedComponent};
        std::vector<detail::HaloFullABCRecoveryMethod> methods;
        if (scenario.real_full_structure) methods.push_back(detail::HaloFullABCRecoveryMethod::FixedWidthACPerComponent);
        else if (scenario.real_geometry)
            methods.assign(recovery_methods.begin(),recovery_methods.begin()+3);
        else methods.assign(recovery_methods.begin(),recovery_methods.end());
        for (const auto method:methods)
        {
            const auto recovered_started=Clock::now();
            const auto recovered=detail::EstimateHaloFullABCJointWithRecovery(
                synthetic.model,problem,snapshot,run_id,method);
            if (capture_solver_diagnostics && (!recovered.available || !recovered.diagnostics.reconstruction_qualified ||
                !recovered.diagnostics.assembly_replay_qualified || recovered.diagnostics.solver_qualified))
                throw std::runtime_error("P1 recovery qualification contract failed for " +
                    recovered.diagnostics.recovery_method + ": " + recovered.unavailable_reason);
            if (capture_solver_diagnostics && (!recovered.diagnostics.assembled_prediction_replay_error ||
                !recovered.diagnostics.assembled_objective_replay_error ||
                !recovered.diagnostics.original_domain_objective ||
                !recovered.diagnostics.original_domain_normalization_scale))
                throw std::runtime_error("P1 mixed-source OriginalInput replay evidence is incomplete.");
            const auto recovered_free_halo=std::find_if(recovered.halos.begin(),recovered.halos.end(),
                [](const auto & halo) {return halo.classification.atom_id=="103";});
            if (capture_solver_diagnostics && (recovered_free_halo==recovered.halos.end() || !recovered_free_halo->qualified ||
                !recovered_free_halo->reconstructed_parameters)
               )
                throw std::runtime_error("P1 recovery omitted the failed Free Halo.");
            if (capture_solver_diagnostics && (method==detail::HaloFullABCRecoveryMethod::FixedWidthACPerHalo ||
                method==detail::HaloFullABCRecoveryMethod::FixedWidthACPerComponent)
               )
            {
                if (recovered_free_halo->parameter_source!=
                    detail::HaloFullABCParameterSource::ObservationACFixedWidthRecovery)
                    throw std::runtime_error("P1 fixed-width AC recovery source was not recorded.");
                const auto check=std::find_if(recovered.diagnostics.recovery.begin(),
                    recovered.diagnostics.recovery.end(),[](const auto & value) {
                        return std::find(value.atom_ids.begin(),value.atom_ids.end(),"103")!=value.atom_ids.end();
                    });
                if (check==recovered.diagnostics.recovery.end() || !check->qualified || !check->full_ac_rank ||
                    !check->independent_replay_passed)
                    throw std::runtime_error("P1 fixed-width AC recovery failed rank or replay qualification.");
            }
            else if (capture_solver_diagnostics && recovered_free_halo->parameter_source!=
                detail::HaloFullABCParameterSource::PostFailureGlobalMedianRecovery)
                throw std::runtime_error("P1 median recovery source was not recorded.");
            if (capture_solver_diagnostics && (method==detail::HaloFullABCRecoveryMethod::PostFailureGlobalMedian ||
                method==detail::HaloFullABCRecoveryMethod::IsolateFailedComponent)
               )
                for (const std::string_view atom_id:{"101","102"})
                {
                    const auto retained=std::find_if(recovered.halos.begin(),recovered.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==atom_id;
                    });
                    if (retained==recovered.halos.end() || retained->parameter_source!=
                        detail::HaloFullABCParameterSource::FullABC)
                        throw std::runtime_error("P1 local fallback discarded an independent FullABC component.");
                }
            if (capture_solver_diagnostics && method==detail::HaloFullABCRecoveryMethod::IsolateFailedComponent)
                for (const std::string_view atom_id:{"101","102"})
                {
                    const auto retained=std::find_if(recovered.halos.begin(),recovered.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==atom_id;
                    });
                    if (retained==recovered.halos.end() || retained->parameter_source!=
                        detail::HaloFullABCParameterSource::FullABC)
                        throw std::runtime_error("P1 failed-component isolation changed a qualified component.");
                }
            MethodResult metrics;
            metrics.available=recovered.available;
            metrics.reason=recovered.unavailable_reason;
            metrics.profiles.resize(synthetic.halos.size());
            const auto recovered_profiles=FullABCProfiles(recovered);
            for (std::size_t index=0;index<synthetic.halos.size();++index)
                if (const auto found=recovered_profiles.find(synthetic.halos[index].id);
                    found!=recovered_profiles.end()) metrics.profiles[index]=found->second;
            if (recovered.available)
            {
                ScoreParameters(metrics,synthetic.halos);
                const auto signal=ScoreSignal(problem,synthetic,recovered_profiles,snapshot);
                metrics.halo_signal_rmse=signal.all;
                metrics.shared_signal_rmse=signal.shared;
                metrics.domain_signal_rmse=signal.domain;
                metrics.max_signal_error=signal.maximum;
                const auto post_started=Clock::now();
                RunPostprocessing(synthetic,problem,snapshot,nullptr,&recovered,run_id,metrics);
                metrics.seconds=std::chrono::duration<double>(Clock::now()-post_started).count();
            }
            json::array sources;
            std::size_t fullabc_source_count{},ac_recovery_count{},median_recovery_count{},
                pre_fixed_count{},unavailable_count{};
            for (const auto & halo:recovered.halos)
            {
                switch(halo.parameter_source)
                {
                case detail::HaloFullABCParameterSource::FullABC: ++fullabc_source_count; break;
                case detail::HaloFullABCParameterSource::ObservationACFixedWidthRecovery: ++ac_recovery_count; break;
                case detail::HaloFullABCParameterSource::PostFailureGlobalMedianRecovery:
                    ++median_recovery_count; break;
                case detail::HaloFullABCParameterSource::PreFixedGlobalMedian: ++pre_fixed_count; break;
                case detail::HaloFullABCParameterSource::Unavailable: ++unavailable_count; break;
                }
                sources.emplace_back(json::object{{"atom_id",halo.classification.atom_id},
                    {"source",HaloParameterSourceName(halo.parameter_source)},
                    {"qualified",halo.qualified},{"reason",halo.qualification_reason},
                    {"parameters",halo.reconstructed_parameters ? json::value(json::array{
                        (*halo.reconstructed_parameters)[0],(*halo.reconstructed_parameters)[1],
                        (*halo.reconstructed_parameters)[2]}) : json::value(nullptr)}});
            }
            json::array recovery_checks;
            for (const auto & check:recovered.diagnostics.recovery)
            {
                json::array atom_ids;
                for (const auto & atom_id:check.atom_ids) atom_ids.emplace_back(atom_id);
                recovery_checks.emplace_back(json::object{{"component_id",check.component_id},
                    {"method",check.method},{"status",check.status},{"reason",check.reason},
                    {"atom_ids",std::move(atom_ids)},
                    {"attempted_halos",check.attempted_halos},{"qualified_halos",check.qualified_halos},
                    {"ac_rank",check.ac_rank ? json::value(*check.ac_rank) : json::value(nullptr)},
                    {"ac_columns",check.ac_columns ? json::value(*check.ac_columns) : json::value(nullptr)},
                    {"ac_condition",Number(check.ac_condition)},
                    {"fixed_width",Number(check.fixed_width)},
                    {"primary_ac_kkt",Number(check.primary_ac_kkt)},
                    {"independent_ac_kkt",Number(check.independent_ac_kkt)},
                    {"objective_replay_error",Number(check.objective_replay_error)},
                    {"full_ac_rank",check.full_ac_rank},
                    {"independent_replay_passed",check.independent_replay_passed},
                    {"fullabc_refit_attempted",check.fullabc_refit_attempted},
                    {"fullabc_refit_qualified",check.fullabc_refit_qualified},
                    {"qualified",check.qualified}});
            }
            recovery_results.emplace_back(json::object{{"method",recovered.diagnostics.recovery_method},
                {"available",recovered.available},{"reason",recovered.unavailable_reason},
                {"solver_qualified",recovered.diagnostics.solver_qualified},
                {"reconstruction_qualified",recovered.diagnostics.reconstruction_qualified},
                {"assembly_replay_qualified",recovered.diagnostics.assembly_replay_qualified},
                {"assembled_prediction_replay_error",Number(
                    recovered.diagnostics.assembled_prediction_replay_error)},
                {"assembled_objective_replay_error",Number(
                    recovered.diagnostics.assembled_objective_replay_error)},
                {"original_domain_objective",Number(recovered.diagnostics.original_domain_objective)},
                {"original_domain_normalization_scale",Number(
                    recovered.diagnostics.original_domain_normalization_scale)},
                {"source_counts",json::object{{"fullabc",fullabc_source_count},
                    {"observation_ac_recovery",ac_recovery_count},{"median_recovery",median_recovery_count},
                    {"pre_fixed_median",pre_fixed_count},{"unavailable",unavailable_count}}},
                {"seconds",std::chrono::duration<double>(Clock::now()-recovered_started).count()},
                {"sources",std::move(sources)},{"recovery_checks",std::move(recovery_checks)},
                {"metrics",MethodJson(metrics)}});
        }
        case_json["p1_recovery_methods"]=std::move(recovery_results);
    }

    MethodResult post_geometry = geometry_result;
    auto geometry_snapshot = snapshot;
    geometry_snapshot.halo_median_reconstruction = geometry.snapshot;
    const auto post_geometry_start = Clock::now();
    RunPostprocessing(synthetic, problem, geometry_snapshot, &geometry, nullptr, run_id, post_geometry);
    post_geometry.seconds = std::chrono::duration<double>(Clock::now() - post_geometry_start).count();
    case_json["geometry_peeling"] = MethodJson(post_geometry);

    MethodResult oracle;
    oracle.available = true;
    oracle.profiles.resize(synthetic.halos.size());
    ProfileMap oracle_profiles;
    for (std::size_t index = 0; index < synthetic.halos.size(); ++index)
    {
        oracle.profiles[index] = synthetic.halos[index].truth;
        oracle_profiles.emplace(synthetic.halos[index].id, synthetic.halos[index].truth);
    }
    ScoreParameters(oracle, synthetic.halos);
    const auto oracle_signal = ScoreSignal(problem, synthetic, oracle_profiles, snapshot);
    oracle.halo_signal_rmse = 0.0;
    oracle.shared_signal_rmse = 0.0;
    oracle.domain_signal_rmse = oracle_signal.domain;
    oracle.max_signal_error = 0.0;
    case_json["oracle_halo_truth"] = MethodJson(oracle);

    case_json["rss_process_peak_bytes"] = PeakRssBytes();
    case_json["total_case_seconds"] = std::chrono::duration<double>(Clock::now() - started).count();
    case_json["halos"] = json::array{};
    for (std::size_t index = 0; index < synthetic.halos.size(); ++index)
    {
        const auto & halo = synthetic.halos[index];
        json::object record{{"atom_id", halo.id}, {"observation_count", halo.rows.size()},
            {"truth_A", halo.truth[0]}, {"truth_B", halo.truth[1]}, {"truth_C", halo.truth[2]},
            {"position", json::array{halo.position[0], halo.position[1], halo.position[2]}},
            {"geometry", json::array{}}, {"observation_AC", json::array{}},
            {"fullabc", json::array{}}};
        const auto geom = geometry_profiles.find(halo.id);
        if (geom != geometry_profiles.end())
            record["geometry"] = json::array{geom->second[0], geom->second[1], geom->second[2]};
        if (ac_profiles && ac_profiles->contains(halo.id))
        {
            const auto & profile = ac_profiles->at(halo.id);
            record["observation_AC"] = json::array{profile[0], profile[1], profile[2]};
        }
        if (full_profiles.contains(halo.id))
        {
            const auto & profile = full_profiles.at(halo.id);
            record["fullabc"] = json::array{profile[0], profile[1], profile[2]};
        }
        case_json["halos"].as_array().emplace_back(record);
        double nearest_target_distance=std::numeric_limits<double>::infinity();
        for(const auto & target:synthetic.target_positions)
            nearest_target_distance=std::min(nearest_target_distance,Distance(halo.position,target));
        halo_csv << scenario.name << ',' << halo.id << ',' << halo.rows.size() << ',' << halo.truth[0] << ','
            << halo.truth[1] << ',' << halo.truth[2] << ',' << nearest_target_distance << ','
            << record.at("geometry").as_array().size() << ',' << record.at("observation_AC").as_array().size() << ','
            << record.at("fullabc").as_array().size() << '\n';
    }
    const auto write_summary = [&](const char * method_name, const MethodResult & method) {
        summary_csv << scenario.name << ',' << method_name << ',' << method.available << ',' << method.reason << ','
            << synthetic.halos.size() << ',' << fullabc.fixed_halo_count << ',' << fullabc.free_halo_count << ','
            << (method.a_rmse ? std::to_string(*method.a_rmse) : "") << ','
            << (method.b_rmse ? std::to_string(*method.b_rmse) : "") << ','
            << (method.c_rmse ? std::to_string(*method.c_rmse) : "") << ','
            << (method.halo_signal_rmse ? std::to_string(*method.halo_signal_rmse) : "") << ','
            << (method.peeling_rmse ? std::to_string(*method.peeling_rmse) : "") << ','
            << (method.tail_ratio_mae ? std::to_string(*method.tail_ratio_mae) : "") << ','
            << (method.local_a_rmse ? std::to_string(*method.local_a_rmse) : "") << ','
            << (method.local_b_rmse ? std::to_string(*method.local_b_rmse) : "") << ','
            << method.seconds << ',' << PeakRssBytes() << '\n';
    };
    write_summary("geometry-kNN4", post_geometry);
    write_summary("observation-AC-fixed-kNN-width", ac_result);
    write_summary("observation-driven-fullabc", full_result);
    write_summary("oracle-truth", oracle);
    return case_json;
}

void RunP2Suite(const std::filesystem::path & output, std::vector<Scenario> scenarios,
    std::string contract)
{
    if (!real_geometry_source) throw std::logic_error("P2 source paths were not configured.");
    const bool production_map_normalization = scenarios.size() == 1 && scenarios.front().real_full_structure;
    std::filesystem::create_directories(output);
    std::ofstream summary(output / "summary.csv");
    std::ofstream halos(output / "halos.csv");
    if (!summary || !halos) throw std::runtime_error("Cannot create P2 result files.");
    summary << "scenario,method,available,reason,halo_count,fixed_count,free_count,A_RMSE,B_RMSE,C_RMSE,halo_signal_RMSE,peeling_RMSE,tail_ratio_MAE,local_A_RMSE,local_B_RMSE,seconds,process_peak_RSS_bytes\n";
    halos << "scenario,atom_id,n_h,truth_A,truth_B,truth_C,nearest_target_distance_angstrom,geometry_fields,observation_AC_fields,fullabc_fields\n";
    json::array cases;
    for (const auto & scenario : scenarios)
    {
        std::cerr << "[P2] " << scenario.name << '\n';
        cases.emplace_back(RunScenario(scenario, summary, halos));
        summary.flush();
        halos.flush();
    }
    std::ofstream output_json(output / "results.json");
    if (!output_json) throw std::runtime_error("Cannot create P2 JSON output.");
    output_json << json::serialize(json::object{
        {"contract",std::move(contract)},
        {"model",real_geometry_source->model.string()},
        {"map",real_geometry_source->map.string()},
        {"simulation_truth",real_geometry_source->truth.string()},
        {"map_normalization_applied",production_map_normalization},
        {"noise_seeds",json::array{4101,4102,4103}},
        {"scenarios",std::move(cases)},
        {"process_peak_rss_bytes",PeakRssBytes()},
        {"method_notes",json::array{
            "The chain-A structural coordinates come from 6Z6U and all map values are generated from the matching simulation truth sidecar.",
            "Each noise seed re-estimates one Target endpoint, then all Halo methods use that same frozen endpoint and observation realization.",
            "Short bonded Halo contacts are retained; exact coincidences and nonbonded center distances below 1.0 Angstrom are rejected as geometry errors.",
            "Group Alpha is reported per case; alpha=0 is not counted as group-alpha training evidence."}}}) << '\n';
}

}

int main(int argc, char ** argv)
try
{
    if (argc > 1 && std::string_view(argv[1]) == "--real-geometry")
    {
        if (argc < 6) throw std::runtime_error("usage: --real-geometry model.cif map truth.json output [--first]");
        real_geometry_source = RealGeometrySource{argv[2],argv[3],argv[4]};
        std::vector<Scenario> scenarios;
        for (const auto seed : {4101ULL,4102ULL,4103ULL})
        {
            Scenario scenario;
            scenario.name = "6z6u_real_geometry_seq145_147_seed_" + std::to_string(seed);
            scenario.noise = 0.003;
            scenario.real_geometry = true;
            scenario.sequence_start = 145;
            scenario.sequence_end = 147;
            scenario.noise_seed = seed;
            scenario.require_short_bonded_halo_pair = true;
            scenarios.push_back(scenario);
        }
        Scenario second_region;
        second_region.name = "6z6u_real_geometry_seq80_82_seed_4101";
        second_region.noise = 0.003;
        second_region.real_geometry = true;
        second_region.sequence_start = 80;
        second_region.sequence_end = 82;
        second_region.noise_seed = 4101;
        scenarios.push_back(second_region);
        Scenario boundary;
        boundary.name = "6z6u_real_geometry_map_boundary_seq5_8_seed_4102";
        boundary.noise = 0.003;
        boundary.map_edge_halo = true;
        boundary.real_geometry = true;
        boundary.sequence_start = 5;
        boundary.sequence_end = 8;
        boundary.noise_seed = 4102;
        scenarios.push_back(boundary);
        if (argc > 6 && std::string_view(argv[6]) == "--first") scenarios.resize(1);
        RunP2Suite(argv[5],std::move(scenarios),
            "fixed-target-halo-fullabc-p2-real-geometry-multiseed-v1");
        std::cout << "P2 real-geometry results written to " << argv[5] << '\n';
        return 0;
    }
    if (argc > 1 && std::string_view(argv[1]) == "--6z6u-partial-selection")
    {
        if (argc != 6) throw std::runtime_error("usage: --6z6u-partial-selection model.cif map truth.json output");
        real_geometry_source = RealGeometrySource{argv[2],argv[3],argv[4]};
        Scenario scenario;
        scenario.name = "6z6u_chain_a_partial_selection";
        scenario.real_geometry = true;
        scenario.real_full_structure = true;
        RunP2Suite(argv[5],{scenario},"fixed-target-halo-fullabc-p2-6z6u-chain-a-paired-v1");
        std::cout << "P2 6Z6U Partial Selection results written to " << argv[5] << '\n';
        return 0;
    }
    const std::filesystem::path output = argc > 1 ? argv[1] : "docs/research/results/halo_fullabc_p3";
    std::filesystem::create_directories(output);
    std::ofstream summary(output / "summary.csv");
    std::ofstream halos(output / "halos.csv");
    if (!summary || !halos) throw std::runtime_error("Cannot create P3 result files.");
    summary << "scenario,method,available,reason,halo_count,fixed_count,free_count,A_RMSE,B_RMSE,C_RMSE,halo_signal_RMSE,peeling_RMSE,tail_ratio_MAE,local_A_RMSE,local_B_RMSE,seconds,process_peak_RSS_bytes\n";
    halos << "scenario,atom_id,n_h,truth_A,truth_B,truth_C,nearest_target_distance_angstrom,geometry_fields,observation_AC_fields,fullabc_fields\n";
    std::vector<Scenario> scenarios{
        {"separated_6p2_similar_low_noise", 6.2, false, false, 0.003, false},
        {"separated_6p2_similar_high_noise", 6.2, false, false, 0.005, false},
        {"overlap_4_similar_low_noise", 4.0, false, false, 0.003, false},
        {"overlap_3_contrast_low_noise", 3.0, true, false, 0.003, false},
        {"overlap_2_width_charge_low_noise", 2.0, true, true, 0.003, false},
        {"overlap_1p6_contrast_low_noise", 1.6, true, false, 0.003, false},
        {"overlap_4_width_charge_contrast", 4.0, true, true, 0.003, false},
        {"boundary_halo_reconstruction_only", 3.0, false, false, 0.003, true},
        {"three_observation_free_stress", 3.0, true, false, 0.003, false, 3},
        {"four_observation_free_stress", 4.0, true, false, 0.003, false, 4}
    };
    if (argc > 2)
    {
        const std::string_view requested = argv[2];
        if (requested == "--first") scenarios.resize(1);
        else
        {
            const auto found = std::find_if(scenarios.begin(), scenarios.end(), [&](const auto & scenario) {
                return scenario.name == requested;
            });
            if (found == scenarios.end()) throw std::runtime_error("Unknown P3 scenario: " + std::string(requested));
            const Scenario selected = *found;
            scenarios.assign(1, selected);
        }
    }
    json::array cases;
    for (const auto & scenario : scenarios)
    {
        std::cerr << "[P3] " << scenario.name << '\n';
        cases.emplace_back(RunScenario(scenario, summary, halos));
        summary.flush(); halos.flush();
    }
    json::object result{
        {"contract", "fixed-target-observation-driven-halo-fullabc-p3-v2-nearest-target-distance"},
        {"map_grid", json::array{51, 63, 16}}, {"spacing_angstrom", 0.5},
        {"target_count", 20}, {"synthetic_scenarios", std::move(cases)},
        {"process_peak_rss_bytes", PeakRssBytes()},
        {"method_notes", json::array{
            "Every scenario uses the actual first-stage local initializer and qualified target-only Pure Exclusion Joint endpoint.",
            "Observation A/C uses the same n_h<3 global-median fixed set and a sparse constrained A/C solve at geometry kNN-4 widths.",
            "The oracle row is a signal and parameter diagnostic only; it does not enter estimation.",
            "The original isolated synthetic scenarios remain separate from the real-geometry P2 suite; no 1.5-Angstrom center-distance gate is applied."
        }}
    };
    std::ofstream json_file(output / "results.json");
    if (!json_file) throw std::runtime_error("Cannot create P3 JSON output.");
    json_file << json::serialize(result) << '\n';
    std::cout << "P3 results written to " << output << '\n';
    return 0;
}
catch (const std::exception & error)
{
    std::cerr << "joint_halo_fullabc_validation: " << error.what() << '\n';
    return 1;
}
