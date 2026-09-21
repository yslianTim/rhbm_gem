#include "utils/domain/FileFingerprint.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "support/MDPDEExperiment.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>

namespace second_stage_test {
namespace {
namespace j = boost::json;
namespace fs = std::filesystem;
namespace sim = rhbm_gem::core::simulation;
using namespace rhbm_gem;
using Position = std::array<double, 3>;
j::value ReadJSON(const fs::path & path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Missing experiment input: " + path.string());
    j::parse_options options; options.numbers = j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in), {}), {}, options);
}
double Distance(const Position & a, const Position & b)
{
    double s{}; for (std::size_t k = 0; k < 3; ++k) s += (a[k]-b[k])*(a[k]-b[k]); return std::sqrt(s);
}
struct Target { std::size_t atom{}; double alpha{}; SamplingPointList points; };
struct Layers
{
    std::vector<double> estimator, generator, grid, roundtrip, original, neighbors;
};

// Contributor selection for prediction uses the entire atom list. The generator,
// not a query at the sampling-point center, determines voxel contributions.
Layers Evaluate(const sim::SimulationAtomPreparationResult & atoms, double width, double cutoff,
    const Target & target, const MapObject & grid, const MapObject & roundtrip, const MapObject * original)
{
    Layers out;
    const auto gs{SampleExperimentPoints(grid, target.points)};
    const auto rs{SampleExperimentPoints(roundtrip, target.points)};
    const auto os{original ? SampleExperimentPoints(*original, target.points) : rs};
    ElectricPotential potential; potential.SetBlurringWidth(width);
    for (std::size_t i = 0; i < target.points.size(); ++i)
    {
        const auto & p{ target.points[i].position };
        double estimator{}, generator{}, neighbors{};
        for (std::size_t a = 0; a < atoms.atom_list.size(); ++a)
        {
            const auto & atom{ atoms.atom_list[a] };
            const double r{ Distance(p, atom.position) };
            const GaussianModel3D model{static_cast<double>(atom.element), width, atom.charge_used};
            if (a == target.atom || (Distance(atom.position, atoms.atom_list[target.atom].position) <= 5.0 && r <= 2.5))
            {
                const double value{ model.ResponseAtDistance(r) };
                estimator += value; if (a != target.atom) neighbors += value;
            }
            if (r <= cutoff) generator += potential.GetPotentialValue(atom.element, r, atom.charge_used);
        }
        out.estimator.push_back(estimator); out.generator.push_back(generator); out.neighbors.push_back(neighbors);
        out.grid.push_back(gs[i].response); out.roundtrip.push_back(rs[i].response); out.original.push_back(os[i].response);
    }
    return out;
}

std::pair<bool, bool> StencilFlags(const MapObject & map, const Position & p,
    const sim::SimulationAtomPreparationResult & atoms, double cutoff)
{
    const auto index{ map.GetIndexFromPosition(p) }, size{ map.GetGridSize() };
    const auto origin{map.GetOrigin()}, h{map.GetGridSpacing()};
    bool boundary{}, crosses{};
    for (std::size_t k = 0; k < 3; ++k) boundary |= index[k] - 1 < 0 || index[k] + 2 >= size[k];
    for (const auto & atom : atoms.atom_list)
    {
        if (std::abs(Distance(p,atom.position)-cutoff) > 4.0 * *std::max_element(h.begin(),h.end())) continue;
        bool inside{}, outside{};
        for (int x = -1; x <= 2; ++x) for (int y = -1; y <= 2; ++y) for (int z = -1; z <= 2; ++z)
        {
            Position node; const std::array<int,3> delta{x,y,z};
            for (std::size_t k = 0; k < 3; ++k) node[k] = origin[k] + h[k] * std::clamp(index[k]+delta[k], 0, size[k]-1);
            if (Distance(node, atom.position) <= cutoff) inside = true; else outside = true;
        }
        crosses |= inside && outside;
    }
    return {crosses, boundary};
}

void WriteComparison(std::ofstream & samples, std::ofstream & fits, const std::string & label,
    const sim::SimulationAtomPreparationResult & atoms, double width, double cutoff,
    const Target & target, const MapObject & grid, const MapObject & roundtrip, const MapObject * original)
{
    const auto values{ Evaluate(atoms, width, cutoff, target, grid, roundtrip, original) };
    const auto & atom{ atoms.atom_list.at(target.atom) };
    const GaussianModel3D model{static_cast<double>(atom.element), width, atom.charge_used};
    const double peak{ model.EvaluateAtDistance(0.0).signal };
    const double h{ grid.GetGridSpacing()[0] };
    LocalPotentialSampleList entries;
    for (std::size_t i = 0; i < target.points.size(); ++i)
    {
        const auto & p{target.points[i]};
        const auto flags{StencilFlags(grid, p.position, atoms, cutoff)};
        samples << label << ',' << atom.serial_id << ',' << h << ',' << i << ',' << p.distance << ','
            << p.is_selected << ',' << flags.first << ',' << flags.second << ',' << peak << ','
            << values.estimator[i] << ',' << values.generator[i] << ',' << values.grid[i] << ','
            << values.roundtrip[i] << ',' << (original ? values.original[i] : std::numeric_limits<double>::quiet_NaN()) << '\n';
        entries.push_back({values.estimator[i], p});
    }
    const core::detail::PreparedLocalGaussianDesign design(entries, 0.0, 1.0);
    const std::array<const std::vector<double> *,5> layers{&values.estimator, &values.generator,
        &values.grid, &values.roundtrip, &values.original};
    const std::array<const char *,5> names{"estimator", "generator", "grid", "roundtrip", "original"};
    for (std::size_t layer = 0; layer < (original ? 5u : 4u); ++layer)
    {
        auto adjusted{ *layers[layer] };
        j::array membership;
        double log_sum{}, log_square{}, log_max{}; std::size_t common{};
        for (std::size_t i = 0; i < adjusted.size(); ++i)
        {
            adjusted[i] -= values.neighbors[i];
            if (target.points[i].distance > 1.0) continue;
            const auto evaluation{model.EvaluateAtDistance(target.points[i].distance)};
            const double corrected{adjusted[i] - (evaluation.response - evaluation.signal)};
            if (corrected <= 0.0) continue;
            membership.push_back(i);
            const double truth_corrected{values.estimator[i] - values.neighbors[i] - (evaluation.response - evaluation.signal)};
            if (truth_corrected <= 0.0) continue;
            const double delta{std::log(corrected) - std::log(truth_corrected)};
            log_sum += delta; log_square += delta*delta; log_max = std::max(log_max,std::abs(delta)); ++common;
        }
        ShapeFixture f; f.alpha = target.alpha; f.dataset = design.BuildDataset(adjusted, model);
        f.expected = rhbm_helper::EstimateBetaMDPDE(f.alpha, f.dataset, f.options);
        auto comparison{CompareMDPDE(f)};
        comparison["case"] = label; comparison["serial_id"] = atom.serial_id;
        comparison["h"] = h; comparison["layer"] = names[layer]; comparison["truth_A"] = static_cast<double>(atom.element);
        comparison["truth_B"] = width; comparison["fixed_C"] = atom.charge_used;
        comparison["membership"] = std::move(membership); comparison["common_count"] = common;
        comparison["log_bias"] = common ? j::value(log_sum/static_cast<double>(common)) : j::value(nullptr);
        comparison["log_rmse"] = common ? j::value(std::sqrt(log_square/static_cast<double>(common))) : j::value(nullptr);
        comparison["log_max"] = log_max;
        // Forward comparisons keep endpoint diagnostics; full trajectories live in solver replay results.
        for (auto & method : comparison.at("methods").as_array()) method.as_object().erase("trace");
        fits << j::serialize(comparison) << '\n';
    }
}

void GridCase(std::ofstream & samples, std::ofstream & fits, const fs::path & output,
    const std::string & label, const sim::SimulationAtomPreparationResult & atoms,
    const std::vector<Target> & targets, double width, double cutoff, double h,
    const Position & origin, const std::array<int,3> & size, const MapObject * original)
{
    MapObject grid(size, {h,h,h}, origin);
    core::MapSimulationRequest request; request.job_count = 4; request.cutoff_distance = cutoff;
    request.potential_model_choice = core::PotentialModel::SINGLE_GAUS;
    sim::PopulateMapValueArray(grid, atoms, request, width);
    const auto file{ output / (label + "-" + std::to_string(h) + ".map") };
    WriteMap(file, grid); const auto roundtrip{ReadMap(file)};
    std::size_t voxel_mismatches{};
    for (std::size_t i = 0; i < grid.GetMapValueArraySize(); ++i)
        if (roundtrip->GetMapValue(i) != static_cast<double>(static_cast<float>(grid.GetMapValue(i)))) ++voxel_mismatches;
    if (voxel_mismatches) throw std::runtime_error("CCP4 voxel round-trip mismatch.");
    for (const auto & target : targets) WriteComparison(samples, fits, label, atoms, width, cutoff, target, grid, *roundtrip, original);
}
} // namespace

void RunForwardExperiment(const std::string & manifest_path, const std::string & map_path,
    const std::string & capture_directory, const std::string & output_directory)
{
    const auto manifest{ReadJSON(manifest_path)};
    if (manifest.at("schema_version")!=1) throw std::runtime_error("Historical experiment only supports frozen manifest v1.");
    if (rhbm_gem::FileSha256(map_path) != j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Map hash does not match the simulation manifest.");
    const auto & settings{manifest.at("settings")};
    const double width{j::value_to<double>(settings.at("blurring_width"))};
    const double cutoff{j::value_to<double>(settings.at("cutoff_distance"))};
    sim::SimulationAtomPreparationResult atoms;
    for (const auto & a : manifest.at("atoms").as_array()) atoms.atom_list.push_back(sim::SimulationAtom{
        .serial_id=j::value_to<int>(a.at("serial_id")), .element=static_cast<Element>(j::value_to<int>(a.at("element"))),
        .position=j::value_to<Position>(a.at("position")), .charge_used=j::value_to<double>(a.at("charge_used"))});
    j::value final_context; std::set<std::size_t> failed_indices;
    for (const auto & file : fs::directory_iterator(capture_directory))
    {
        const auto name{file.path().filename().string()};
        if (!name.starts_with("shape-") || !name.ends_with(".txt.json")) continue;
        const auto metadata{ReadJSON(file.path())};
        if (metadata.at("phase") == "final") final_context = ReadJSON(fs::path(capture_directory) /
            j::value_to<std::string>(metadata.at("context_file")));
        if (!name.starts_with("shape-success-") &&
            (metadata.at("phase") == "final" || metadata.at("phase") == "recovery-current"))
            failed_indices.insert(j::value_to<std::size_t>(metadata.at("indices").at(0)));
    }
    if (final_context.is_null() || failed_indices.empty()) throw std::runtime_error("Missing final/recovery experiment captures.");
    const auto original{ReadMap(map_path)};
    std::vector<Target> targets;
    for (const auto & a : final_context.at("atoms").as_array())
    {
        const int serial{j::value_to<int>(a.at("identity").at("serial_id"))};
        const auto it{std::find_if(atoms.atom_list.begin(), atoms.atom_list.end(), [&](const auto & atom){return atom.serial_id == serial;})};
        if (it == atoms.atom_list.end()) throw std::runtime_error("Unmatched captured atom.");
        const auto & truth_record{manifest.at("atoms").at(static_cast<std::size_t>(it-atoms.atom_list.begin()))};
        for (const auto * field : {"serial_id","chain_id","sequence_id","component_id","atom_id","alternate_indicator","position"})
            if (truth_record.at(field) != a.at("identity").at(field)) throw std::runtime_error("Captured atom identity mismatch.");
        Target target{static_cast<std::size_t>(it-atoms.atom_list.begin()), j::value_to<double>(a.at("alpha")), {}};
        for (const auto & s : a.at("samples").as_array()) target.points.push_back({j::value_to<double>(s.at("distance")),
            j::value_to<Position>(s.at("position")), s.at("selected").as_bool()});
        const auto replay{SampleExperimentPoints(*original, target.points)};
        for (std::size_t i = 0; i < replay.size(); ++i)
            if (replay[i].response != j::value_to<double>(a.at("samples").at(i).at("response")))
                throw std::runtime_error("Captured map sample did not replay exactly.");
        targets.push_back(std::move(target));
    }
    if (targets.size() != atoms.atom_list.size()) throw std::runtime_error("Incomplete captured atom population.");
    const fs::path output{output_directory}; fs::create_directories(output);
    std::ofstream samples(output / "samples.csv"), fits(output / "fits.jsonl");
    samples.exceptions(std::ios::failbit | std::ios::badbit); fits.exceptions(std::ios::failbit | std::ios::badbit);
    samples << std::setprecision(17) << "case,serial_id,h,sample,distance,selected,crosses_cutoff,map_boundary,peak,estimator,generator,grid,roundtrip,original\n";
    const auto origin{j::value_to<Position>(settings.at("origin"))};
    GridCase(samples, fits, output, "fold", atoms, targets, width, cutoff, 0.1, origin,
        j::value_to<std::array<int,3>>(settings.at("grid_size")), original.get());
    for (const auto index : failed_indices) for (double h : {0.1,0.05,0.025})
    {
        const auto & target{targets.at(index)}; Position local_origin; std::array<int,3> size;
        for (std::size_t k = 0; k < 3; ++k)
        {
            double lo{std::numeric_limits<double>::max()}, hi{std::numeric_limits<double>::lowest()};
            for (const auto & p : target.points) {lo = std::min(lo,p.position[k]); hi = std::max(hi,p.position[k]);}
            const double first{std::floor((lo-origin[k])/h)-2};
            local_origin[k] = origin[k]+first*h; size[k] = static_cast<int>(std::ceil((hi-origin[k])/h)-first)+4;
        }
        GridCase(samples, fits, output, "failure-"+std::to_string(atoms.atom_list[target.atom].serial_id),
            atoms, {target}, width, cutoff, h, local_origin, size, nullptr);
    }
    for (int mode = 0; mode < 4; ++mode) for (double h : {0.1,0.05,0.025})
    {
        sim::SimulationAtomPreparationResult small;
        small.atom_list.push_back(sim::SimulationAtom{.serial_id=1,.element=Element::CARBON,
            .position={0.013,0.021,0.037},.charge_used=mode == 0 ? 0.0 : (mode == 1 ? 0.3 : -0.3)});
        if (mode == 3)
        {
            small.atom_list.push_back(sim::SimulationAtom{.serial_id=2,.element=Element::OXYGEN,.position={1.1,0.4,-0.2},.charge_used=0.3});
            small.atom_list.push_back(sim::SimulationAtom{.serial_id=3,.element=Element::NITROGEN,.position={-1.2,0.3,0.2},.charge_used=0.0});
        }
        Target target{0,0.1,{}};
        for (double r : {0.0,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0,2.49,2.5,2.51,2.8})
            for (std::size_t axis = 0; axis < 3; ++axis) for (int sign : {-1,1})
            {auto p{small.atom_list[0].position}; p[axis] += sign*r; target.points.push_back({r,p,true});}
        target.points.push_back({Distance({-3.0,-3.0,-3.0}, small.atom_list[0].position), {-3.0,-3.0,-3.0},true});
        const int n{static_cast<int>(std::round(6.0/h))+1};
        GridCase(samples,fits,output,"small-"+std::to_string(mode),small,{target},0.5,2.5,h,{-3.0,-3.0,-3.0},{n,n,n},nullptr);
    }
}
} // namespace second_stage_test
