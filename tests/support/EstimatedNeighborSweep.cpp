#include "support/EstimatedNeighborSweep.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace second_stage_test::matched {
namespace {
namespace j = boost::json;
namespace fs = std::filesystem;
namespace sim = rhbm_gem::core::simulation;
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("Missing sweep input: " + path.string());
    j::parse_options options; options.numbers = j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in), {}), {}, options);
}
void Write(const fs::path & path, const j::value & value)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit | std::ios::badbit); out << j::serialize(value) << '\n';
}
double Analytic(const Position & position, const std::vector<Atom> & atoms, double cutoff)
{
    double value{};
    for (const auto & atom : atoms)
    {
        const auto b{EvaluateBasis(SquareDistance(position, atom.position), atom.width, cutoff)};
        value += atom.amplitude*b.gaussian + atom.charge*b.charge;
    }
    return value;
}
bool Crosses(const Stencil & stencil, const std::vector<Atom> & atoms, double cutoff)
{
    for (const auto & atom : atoms) for (double radius : {cutoff, std::min(cutoff, 2.5)})
    {
        bool inside{}, outside{};
        for (const auto & slot : stencil.slots)
        {
            if (SquareDistance(slot.position, atom.position) <= radius*radius) inside = true; else outside = true;
        }
        if (inside && outside) return true;
    }
    return false;
}
struct Target
{
    j::value identity;
    std::vector<Position> positions;
    std::vector<double> distances;
    std::vector<bool> selected, crosses;
    std::vector<Stencil> stencils;
    Eigen::VectorXd observations;
    std::array<Eigen::VectorXd, 2> truth_neighbors;
};
} // namespace

FrozenTarget PrepareFrozenTarget(const std::vector<Atom> & state, std::size_t target,
    const std::vector<Position> & positions, const std::vector<Stencil> & stencils, double cutoff)
{
    if (target >= state.size() || positions.size() != stencils.size())
        throw std::invalid_argument("Invalid frozen target membership.");
    for (const auto & atom : state)
        if (!std::isfinite(atom.amplitude) || atom.amplitude < 0.0 || !std::isfinite(atom.width) ||
            atom.width <= 0.0 || !std::isfinite(atom.charge) ||
            !std::all_of(atom.position.begin(), atom.position.end(), [](double x){ return std::isfinite(x); }))
            throw std::invalid_argument("Invalid checkpoint model.");
    auto neighbors{state}; neighbors.erase(neighbors.begin() + static_cast<std::ptrdiff_t>(target));
    FrozenTarget out;
    for (auto & y : out.neighbors) y.resize(static_cast<Eigen::Index>(positions.size()));
    for (std::size_t p=0; p<positions.size(); ++p)
    {
        out.design[0].push_back({{SquareDistance(positions[p], state[target].position), 1.0}});
        out.neighbors[0](static_cast<Eigen::Index>(p)) = Analytic(positions[p], neighbors, cutoff);
        out.neighbors[1](static_cast<Eigen::Index>(p)) = Predict(stencils[p], neighbors, cutoff);
    }
    out.design[1] = MakeDesign(stencils, state[target].position);
    return out;
}

j::array FitFrozenTarget(const FrozenTarget & prepared, const std::vector<double> & distances,
    const Eigen::VectorXd & observations, const Atom & input, double cutoff)
{
    if (distances.size() != static_cast<std::size_t>(observations.size()) || !observations.allFinite())
        throw std::invalid_argument("Invalid frozen target observations.");
    j::array results;
    for (std::size_t mode=0; mode<2; ++mode) for (bool free_charge : {false, true})
    {
        Design design; std::vector<double> response; j::array membership;
        for (std::size_t p=0; p<distances.size(); ++p)
        {
            if (!free_charge && distances[p] > 1.0) continue;
            design.push_back(prepared.design[mode].at(p));
            response.push_back(observations(static_cast<Eigen::Index>(p))-prepared.neighbors[mode](static_cast<Eigen::Index>(p)));
            membership.push_back(p);
        }
        const Eigen::Map<const Eigen::VectorXd> y(response.data(), static_cast<Eigen::Index>(response.size()));
        const auto start{std::chrono::steady_clock::now()};
        auto fit{Fit(design, y, free_charge ? 0.0 : input.charge, free_charge, cutoff)};
        fit["seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const auto basis{EvaluateDesign(design, input.width, cutoff)};
        fit["input_loss"] = (y-input.amplitude*basis.col(0)-input.charge*basis.col(1)).squaredNorm()/static_cast<double>(y.size());
        fit["input"] = j::array{input.amplitude, input.width, input.charge};
        fit["prediction"] = mode == 0 ? "analytic" : "matched";
        fit["parameters"] = free_charge ? "ABC" : "AB";
        fit["membership"] = std::move(membership);
        results.push_back(std::move(fit));
    }
    return results;
}

void RunEstimatedNeighborSweep(const std::string & manifest_path, const std::string & map_path,
    const std::string & index_path, const std::string & output_path)
{
    const fs::path output{output_path};
    if (fs::exists(output/"fits") || fs::exists(output/"samples")) throw std::runtime_error("Sweep output already exists.");
    const auto manifest{Read(manifest_path)}, index{Read(index_path)};
    const auto & settings{manifest.at("settings")};
    if (settings.at("potential_model") != "single_gaus" || index.at("schema_version") != 1 ||
        sim::FileSha256(map_path) != j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Unsupported sweep manifest or mismatched map hash.");
    const double cutoff{j::value_to<double>(settings.at("cutoff_distance"))};
    const auto original{rhbm_gem::ReadMap(map_path)};
    rhbm_gem::MapObject generation(j::value_to<std::array<int, 3>>(settings.at("grid_size")),
        j::value_to<Position>(settings.at("grid_spacing")), j::value_to<Position>(settings.at("origin")));
    const int jobs{std::getenv("OMP_NUM_THREADS") ? std::stoi(std::getenv("OMP_NUM_THREADS")) : 4};
    if (jobs != 1 && jobs != 4) throw std::runtime_error("Sweep requires OMP_NUM_THREADS=1 or 4.");
    std::vector<j::value> contexts;
    std::set<std::string> ids;
    for (const auto & entry : index.at("states").as_array())
    {
        const auto path{j::value_to<std::string>(entry.at("context"))};
        if (!ids.insert(j::value_to<std::string>(entry.at("id"))).second ||
            sim::FileSha256(path) != j::value_to<std::string>(entry.at("sha256")))
            throw std::runtime_error("Duplicate state or mismatched context hash: " + path);
        contexts.push_back(Read(path));
        const auto & c{contexts.back()};
        if (c.at("schema_version") != 1 || c.at("phase") != entry.at("phase") ||
            c.at("attempt") != entry.at("attempt") || c.at("state").as_array().size() != 168 ||
            c.at("atoms").as_array().size() != 168)
            throw std::runtime_error("Invalid sweep context: " + path);
    }
    if (contexts.size() != 8 || manifest.at("atoms").as_array().size() != 168)
        throw std::runtime_error("Sweep requires eight states of the full fold-168 population.");

    std::vector<Target> targets;
    std::vector<Atom> truth;
    std::set<int> serials;
    // Truth is used only for independent forward checks and later subtraction scoring.
    for (std::size_t i=0; i<168; ++i)
    {
        const auto & captured{contexts[0].at("atoms").at(i)};
        Target target; target.identity = captured.at("identity");
        const int serial{j::value_to<int>(target.identity.at("serial_id"))};
        const auto & atoms{manifest.at("atoms").as_array()};
        const auto it{std::find_if(atoms.begin(), atoms.end(), [&](const auto & a){return a.at("serial_id") == serial;})};
        if (it == atoms.end() || !serials.insert(serial).second) throw std::runtime_error("Missing or duplicate sweep identity.");
        for (const char * key : {"serial_id", "chain_id", "sequence_id", "component_id", "atom_id", "alternate_indicator", "position"})
            if (it->at(key) != target.identity.at(key)) throw std::runtime_error("Sweep identity mismatch.");
        for (const auto & c : contexts)
            if (c.at("atoms").at(i).at("index") != i || c.at("atoms").at(i).at("identity") != target.identity ||
                c.at("atoms").at(i).at("samples") != captured.at("samples"))
                throw std::runtime_error("State identity or sample membership differs.");
        truth.push_back({j::value_to<Position>(target.identity.at("position")),
            static_cast<double>(ChemicalDataHelper::GetAtomicNumber(static_cast<Element>(j::value_to<int>(it->at("element"))))),
            j::value_to<double>(settings.at("blurring_width")), j::value_to<double>(it->at("charge_used"))});
        SamplingPointList points;
        for (const auto & s : captured.at("samples").as_array())
        {
            target.positions.push_back(j::value_to<Position>(s.at("position")));
            target.distances.push_back(j::value_to<double>(s.at("distance")));
            target.selected.push_back(s.at("selected").as_bool());
            target.stencils.push_back(MakeStencil(generation, *original, target.positions.back()));
            points.push_back({target.distances.back(), target.positions.back(), target.selected.back()});
        }
        if (points.size() != 200 || std::count_if(target.distances.begin(), target.distances.end(), [](double r){return r<=1.0;}) != 100)
            throw std::runtime_error("Incomplete sweep sample population.");
        const auto replay{SampleExperimentPoints(*original, points)};
        target.observations.resize(200);
        for (std::size_t p=0; p<200; ++p)
        {
            const double value{j::value_to<double>(captured.at("samples").at(p).at("response"))};
            if (!std::isfinite(value) || replay[p].response != value) throw std::runtime_error("Sweep response failed exact map replay.");
            target.observations(static_cast<Eigen::Index>(p)) = value;
        }
        targets.push_back(std::move(target));
    }
    double maximum_float_difference{}, maximum_bound_excess{};
    for (std::size_t i=0; i<targets.size(); ++i)
    {
        auto & target{targets[i]};
        target.truth_neighbors = PrepareFrozenTarget(truth, i, target.positions, target.stencils, cutoff).neighbors;
        for (std::size_t p=0; p<200; ++p)
        {
            double absolute{}, bound{};
            const double prediction{Predict(target.stencils[p], truth, cutoff, false, &absolute, &bound)};
            const double quantized{Predict(target.stencils[p], truth, cutoff, true)};
            const double tau{512*std::numeric_limits<double>::epsilon()*std::max(1.0, absolute)};
            const double observed{target.observations(static_cast<Eigen::Index>(p))};
            const double difference{std::abs(quantized-observed)}, excess{std::abs(prediction-observed)-bound-tau};
            if (difference>tau || excess>0.0) throw std::runtime_error("Sweep forward validation failed.");
            maximum_float_difference = std::max(maximum_float_difference, difference);
            maximum_bound_excess = std::max(maximum_bound_excess, excess);
            target.crosses.push_back(Crosses(target.stencils[p], truth, cutoff));
        }
    }
    fs::create_directories(output/"fits"); fs::create_directories(output/"samples");
    Write(output/"forward-status.json", j::object{{"passed", true}, {"capture_replay", true}, {"samples", 33600},
        {"maximum_float_difference", maximum_float_difference}, {"maximum_bound_excess", maximum_bound_excess}});
    for (std::size_t s=0; s<contexts.size(); ++s)
    {
        const auto id{j::value_to<std::string>(index.at("states").at(s).at("id"))};
        std::vector<Atom> models;
        for (std::size_t i=0; i<168; ++i)
        {
            const auto abc{j::value_to<std::array<double, 3>>(contexts[s].at("state").at(i))};
            models.push_back({j::value_to<Position>(targets[i].identity.at("position")), abc[0], abc[1], abc[2]});
        }
        const auto state{std::move(models)};
        std::vector<std::string> sample_rows(168);
        std::vector<std::exception_ptr> errors(168);
#ifdef USE_OPENMP
        #pragma omp parallel for schedule(static) num_threads(jobs)
#endif
        for (std::size_t i=0; i<168; ++i)
        {
            try
            {
                const auto & t{targets[i]};
                const auto prepared{PrepareFrozenTarget(state, i, t.positions, t.stencils, cutoff)};
                auto fits = FitFrozenTarget(prepared, t.distances, t.observations, state[i], cutoff);
                for (auto & fit : fits)
                {
                    fit.as_object()["state_id"] = id;
                    fit.as_object()["identity"] = t.identity;
                    fit.as_object()["serial_id"] = t.identity.at("serial_id");
                }
                Write(output/"fits"/(std::to_string(s*168+i)+".json"), fits);
                std::ostringstream rows; rows << std::setprecision(17);
                for (std::size_t p=0; p<200; ++p)
                {
                    const auto k{static_cast<Eigen::Index>(p)};
                    rows << id << ',' << j::value_to<int>(t.identity.at("serial_id")) << ',' << p << ','
                        << t.distances[p] << ',' << t.selected[p] << ',' << t.crosses[p] << ',' << t.stencils[p].boundary
                        << ',' << t.observations(k);
                    for (std::size_t mode=0; mode<2; ++mode)
                    {
                        const double estimated{prepared.neighbors[mode](k)}, reference{t.truth_neighbors[mode](k)},
                            common{t.truth_neighbors[1](k)};
                        rows << ',' << estimated << ',' << reference << ',' << estimated-reference
                            << ',' << reference-common << ',' << estimated-common;
                    }
                    rows << '\n';
                }
                sample_rows[i] = rows.str();
            }
            catch (...) { errors[i] = std::current_exception(); }
        }
        for (const auto & error : errors) if (error) std::rethrow_exception(error);
        std::ofstream samples(output/"samples"/(std::to_string(s)+".csv"));
        samples.exceptions(std::ios::failbit | std::ios::badbit);
        samples << "state_id,serial_id,sample,distance,selected,crosses_cutoff,boundary,observation";
        for (const char * mode : {"analytic", "matched"}) for (const char * field : {"estimate", "truth", "estimation", "operator", "total"})
            samples << ',' << mode << '_' << field;
        samples << '\n'; for (const auto & rows : sample_rows) samples << rows;
    }
    Write(output/"fit-index.json", j::object{{"schema_version", 1}, {"experiment", "estimated-neighbor-sweep"},
        {"states", index.at("states")}, {"targets", 1344}, {"fits", 5376}, {"jobs", jobs},
        {"width_bounds", j::array{0.1, 2.0}}, {"profile_points", 129}, {"reference_points", 257},
        {"brent_budget", 128}, {"stationarity_tolerance", 1e-8}, {"reference_tolerance", 1e-6}});
}
} // namespace second_stage_test::matched
