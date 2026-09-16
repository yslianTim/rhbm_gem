#include "support/ObservationMatchedExperiment.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>

namespace second_stage_test::matched {
namespace {
namespace j = boost::json;
namespace fs = std::filesystem;
namespace sim = rhbm_gem::core::simulation;
using rhbm_gem::MapObject;
struct Target {std::size_t atom; SamplingPointList points;};
struct Prepared
{
    std::string label; int serial; double h, cutoff; Atom truth;
    SamplingPointList points;
    std::array<Design,3> design;
    std::array<Eigen::VectorXd,3> neighbors;
    std::array<Eigen::VectorXd,4> observations;
};
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("Missing experiment input: "+path.string());
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},options);
}
void Write(const fs::path & path, const j::value & value)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out << j::serialize(value) << '\n';
}
bool Crosses(const Stencil & stencil, const std::vector<Atom> & atoms, double cutoff)
{
    for (const auto & a : atoms)
    {
        bool inside{},outside{};
        for (const auto & slot : stencil.slots)
        {
            if (SquareDistance(slot.position,a.position)<=cutoff*cutoff) inside=true; else outside=true;
        }
        if (inside && outside) return true;
    }
    return false;
}
double Analytic(const Position & p, const std::vector<Atom> & atoms, double cutoff)
{
    double value{};
    for (const auto & a : atoms)
    {
        const auto b{EvaluateBasis(SquareDistance(p,a.position),a.width,cutoff)};
        value += a.amplitude*b.gaussian+a.charge*b.charge;
    }
    return value;
}
std::vector<Atom> Models(const sim::SimulationAtomPreparationResult & atoms, double width)
{
    std::vector<Atom> result;
    for (const auto & a : atoms.atom_list) result.push_back({a.position,
        static_cast<double>(ChemicalDataHelper::GetAtomicNumber(a.element)),width,a.charge_used});
    return result;
}

void PrepareCase(const fs::path & output, const std::string & label,
    const sim::SimulationAtomPreparationResult & atoms, const std::vector<Target> & targets,
    double width, double cutoff, double h, const Position & origin, const std::array<int,3> & size,
    const MapObject * original, int jobs, std::ofstream & samples, std::vector<Prepared> & prepared,
    j::array & checks)
{
    MapObject grid(size,{h,h,h},origin);
    rhbm_gem::core::MapSimulationRequest request;
    request.job_count=jobs; request.cutoff_distance=cutoff;
    request.potential_model_choice=rhbm_gem::core::PotentialModel::SINGLE_GAUS;
    sim::PopulateMapValueArray(grid,atoms,request,width);
    const auto map_path{output/(label+"-"+std::to_string(h)+".map")};
    rhbm_gem::WriteMap(map_path,grid); const auto roundtrip{rhbm_gem::ReadMap(map_path)};
    auto header_values{std::make_unique<double[]>(grid.GetMapValueArraySize())};
    std::size_t voxel_mismatches{}, original_mismatches{};
    for (std::size_t i=0;i<grid.GetMapValueArraySize();++i)
    {
        header_values[i]=grid.GetMapValue(i);
        voxel_mismatches += roundtrip->GetMapValue(i)!=static_cast<double>(static_cast<float>(header_values[i]));
        if (original) original_mismatches += original->GetMapValue(i)!=roundtrip->GetMapValue(i);
    }
    MapObject header(size,roundtrip->GetGridSpacing(),roundtrip->GetOrigin(),std::move(header_values));
    bool geometry_equal{!original || (original->GetGridSize()==size && original->GetOrigin()==roundtrip->GetOrigin()
        && original->GetGridSpacing()==roundtrip->GetGridSpacing())};
    const auto models{Models(atoms,width)};
    std::size_t count{}, failures{};
    double max_double{}, max_header{}, max_float{}, max_bound_excess{};
    std::set<std::size_t> nodes;
    for (const auto & target : targets)
    {
        Prepared item; item.label=label; item.serial=atoms.atom_list.at(target.atom).serial_id;
        item.h=h; item.cutoff=cutoff; item.truth=models.at(target.atom); item.points=target.points;
        const auto native_grid{SampleExperimentPoints(grid,target.points)};
        const auto native_header{SampleExperimentPoints(header,target.points)};
        const auto native_roundtrip{SampleExperimentPoints(*roundtrip,target.points)};
        const auto native_original{original ? SampleExperimentPoints(*original,target.points) : native_roundtrip};
        std::vector<Stencil> stencils, file_stencils;
        const auto n{static_cast<Eigen::Index>(target.points.size())};
        for (auto & y : item.observations) y.resize(n);
        for (auto & y : item.neighbors) y.resize(n);
        auto neighbors{models}; neighbors.erase(neighbors.begin()+static_cast<std::ptrdiff_t>(target.atom));
        for (std::size_t p=0;p<target.points.size();++p)
        {
            const auto & point{target.points[p]}; const auto row{static_cast<Eigen::Index>(p)};
            stencils.push_back(MakeStencil(grid,grid,point.position));
            file_stencils.push_back(MakeStencil(grid,*roundtrip,point.position));
            for (const auto & s : stencils.back().slots) nodes.insert(s.index);
            for (const auto & s : file_stencils.back().slots) nodes.insert(s.index);
            double absolute{}, file_absolute{}, bound{};
            const double prediction{Predict(stencils.back(),models,cutoff,false,&absolute)};
            const double file_prediction{Predict(file_stencils.back(),models,cutoff,false,&file_absolute,&bound)};
            const double quantized{Predict(file_stencils.back(),models,cutoff,true)};
            const double tau{512*std::numeric_limits<double>::epsilon()*std::max(1.0,absolute)};
            const double file_tau{512*std::numeric_limits<double>::epsilon()*std::max(1.0,file_absolute)};
            const double dg{std::abs(prediction-native_grid[p].response)};
            const double dh{std::abs(file_prediction-native_header[p].response)};
            const double dq{std::abs(quantized-native_roundtrip[p].response)};
            const double excess{std::abs(file_prediction-native_roundtrip[p].response)-bound-file_tau};
            const bool passed{dg<=tau && dh<=file_tau && dq<=file_tau && excess<=0.0 &&
                native_original[p].response==native_roundtrip[p].response};
            ++count; failures += !passed;
            max_double=std::max(max_double,dg); max_header=std::max(max_header,dh); max_float=std::max(max_float,dq);
            max_bound_excess=std::max(max_bound_excess,excess);
            item.observations[0](row)=Analytic(point.position,models,cutoff);
            item.observations[1](row)=native_grid[p].response;
            item.observations[2](row)=native_header[p].response;
            item.observations[3](row)=native_original[p].response;
            item.neighbors[0](row)=Analytic(point.position,neighbors,cutoff);
            item.neighbors[1](row)=Predict(stencils.back(),neighbors,cutoff);
            item.neighbors[2](row)=Predict(file_stencils.back(),neighbors,cutoff);
            item.design[0].push_back({{SquareDistance(point.position,item.truth.position),1.0}});
            samples << label << ',' << item.serial << ',' << h << ',' << p << ',' << point.distance << ','
                << point.is_selected << ',' << Crosses(stencils.back(),models,cutoff) << ','
                << Crosses(file_stencils.back(),models,cutoff) << ',' << stencils.back().boundary << ','
                << file_stencils.back().boundary << ',' << item.observations[0](row) << ',' << native_grid[p].response << ','
                << native_header[p].response << ',' << native_roundtrip[p].response << ',' << native_original[p].response << ','
                << prediction << ',' << file_prediction << ',' << quantized << ',' << tau << ',' << file_tau << ','
                << bound << ',' << passed << '\n';
        }
        item.design[1]=MakeDesign(stencils,item.truth.position);
        item.design[2]=MakeDesign(file_stencils,item.truth.position);
        prepared.push_back(std::move(item));
    }
    const bool passed{failures==0 && voxel_mismatches==0 && original_mismatches==0 && geometry_equal};
    checks.emplace_back(j::object{{"case",label},{"h",h},{"samples",count},{"failures",failures},{"passed",passed},
        {"max_double_difference",max_double},{"max_header_difference",max_header},{"max_float_difference",max_float},
        {"max_quantization_bound_excess",max_bound_excess},{"voxel_mismatches",voxel_mismatches},
        {"original_voxel_mismatches",original_mismatches},{"original_geometry_equal",geometry_equal},
        {"unique_voxels",nodes.size()},{"generation_origin",j::value_from(origin)},
        {"sampling_origin",j::value_from(roundtrip->GetOrigin())},{"sampling_spacing",j::value_from(roundtrip->GetGridSpacing())}});
    Write(output/"forward-checks.json",checks);
    samples.flush();
    if (!passed) throw std::runtime_error("Forward validation failed: "+label+" h="+std::to_string(h));
}

void FitPrepared(const fs::path & output, const std::vector<Prepared> & prepared, int jobs)
{
    fs::create_directories(output/"fits");
    std::vector<std::exception_ptr> errors(prepared.size());
    // Independent targets, fixed per-target arithmetic and deterministic filenames.
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) num_threads(jobs)
#endif
    for (std::size_t i=0;i<prepared.size();++i)
    {
        try
        {
            const auto & p{prepared[i]}; j::array results;
            for (std::size_t layer=0;layer<4;++layer) for (bool matched : {false,true}) for (bool free_charge : {false,true})
            {
                const std::size_t mode{!matched ? 0u : layer<=1 ? 1u : 2u};
                Design design; std::vector<double> response; j::array membership;
                for (std::size_t s=0;s<p.points.size();++s)
                {
                    if (!free_charge && p.points[s].distance>1.0) continue;
                    design.push_back(p.design[mode][s]);
                    response.push_back(p.observations[layer](static_cast<Eigen::Index>(s))-p.neighbors[mode](static_cast<Eigen::Index>(s)));
                    membership.push_back(s);
                }
                const Eigen::Map<const Eigen::VectorXd> y(response.data(),static_cast<Eigen::Index>(response.size()));
                const auto start{std::chrono::steady_clock::now()};
                auto result{Fit(design,y,free_charge ? 0.0 : p.truth.charge,free_charge,p.cutoff)};
                result["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                result["case"]=p.label; result["serial_id"]=p.serial; result["h"]=p.h;
                result["layer"]=std::array{"analytic","grid","header","file"}[layer];
                result["prediction"]=matched ? "matched" : "analytic";
                result["parameters"]=free_charge ? "ABC" : "AB"; result["membership"]=std::move(membership);
                result["truth"]=j::array{p.truth.amplitude,p.truth.width,p.truth.charge};
                results.push_back(std::move(result));
            }
            Write(output/"fits"/(std::to_string(i)+".json"),results);
        }
        catch (...) {errors[i]=std::current_exception();}
    }
    for (const auto & error : errors) if (error) std::rethrow_exception(error);
    Write(output/"fit-index.json",j::object{{"schema_version",1},{"targets",prepared.size()},
        {"fits",prepared.size()*16},{"jobs",jobs},{"width_bounds",j::array{0.1,2.0}},
        {"profile_points",129},{"reference_points",257},{"brent_budget",128},
        {"stationarity_tolerance",1e-8},{"reference_tolerance",1e-6}});
}
} // namespace

void Run(const std::string & manifest_path, const std::string & map_path,
    const std::string & captures, const std::string & output_path)
{
    const fs::path output{output_path}; fs::create_directories(output);
    if (fs::exists(output/"samples.csv") || fs::exists(output/"fits")) throw std::runtime_error("Experiment output already exists.");
    const auto manifest{Read(manifest_path)}; const auto & settings{manifest.at("settings")};
    if (settings.at("potential_model")!="single_gaus" ||
        sim::FileSha256(map_path)!=j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Unsupported generation model or mismatched map hash.");
    const double width{j::value_to<double>(settings.at("blurring_width"))}, cutoff{j::value_to<double>(settings.at("cutoff_distance"))};
    const auto origin{j::value_to<Position>(settings.at("origin"))};
    const auto original{rhbm_gem::ReadMap(map_path)};
    sim::SimulationAtomPreparationResult atoms;
    for (const auto & a : manifest.at("atoms").as_array()) atoms.atom_list.push_back(sim::SimulationAtom{
        .serial_id=j::value_to<int>(a.at("serial_id")),.element=static_cast<Element>(j::value_to<int>(a.at("element"))),
        .position=j::value_to<Position>(a.at("position")),.charge_used=j::value_to<double>(a.at("charge_used"))});
    std::set<fs::path> contexts;
    for (const auto & file : fs::directory_iterator(captures))
    {
        const auto name{file.path().filename().string()};
        if (!name.starts_with("shape-") || !name.ends_with(".txt.json")) continue;
        const auto metadata{Read(file.path())};
        if (metadata.at("phase")=="final") contexts.insert(fs::path(captures)/j::value_to<std::string>(metadata.at("context_file")));
    }
    if (contexts.size()!=1) throw std::runtime_error("Expected exactly one final capture context.");
    const auto context{Read(*contexts.begin())}; std::vector<Target> targets; std::set<std::size_t> identities;
    for (const auto & a : context.at("atoms").as_array())
    {
        const auto & identity{a.at("identity")};
        const auto & truth{manifest.at("atoms").as_array()};
        const auto it{std::find_if(truth.begin(),truth.end(),[&](const auto & v){return v.at("serial_id")==identity.at("serial_id");})};
        if (it==truth.end()) throw std::runtime_error("Missing atom identity.");
        for (const char * key : {"serial_id","chain_id","sequence_id","component_id","atom_id","alternate_indicator","position"})
            if (it->at(key)!=identity.at(key)) throw std::runtime_error("Atom identity mismatch.");
        Target target{static_cast<std::size_t>(it-truth.begin()),{}};
        if (!identities.insert(target.atom).second) throw std::runtime_error("Duplicate captured identity.");
        for (const auto & s : a.at("samples").as_array()) target.points.push_back({j::value_to<double>(s.at("distance")),
            j::value_to<Position>(s.at("position")),s.at("selected").as_bool()});
        const auto replay{SampleExperimentPoints(*original,target.points)};
        for (std::size_t s=0;s<replay.size();++s)
            if (replay[s].response!=j::value_to<double>(a.at("samples").at(s).at("response")))
                throw std::runtime_error("Capture response failed exact replay.");
        targets.push_back(std::move(target));
    }
    if (targets.size()!=168 || identities.size()!=atoms.atom_list.size() ||
        std::any_of(targets.begin(),targets.end(),[](const auto & t){return t.points.size()!=200;}))
        throw std::runtime_error("Incomplete fold-168 population.");
    const int jobs{std::getenv("OMP_NUM_THREADS") ? std::stoi(std::getenv("OMP_NUM_THREADS")) : 4};
    if (jobs!=1 && jobs!=4) throw std::runtime_error("Experiment requires OMP_NUM_THREADS=1 or 4.");
    std::ofstream samples(output/"samples.csv"); samples.exceptions(std::ios::failbit|std::ios::badbit);
    samples << std::setprecision(17) << "case,serial_id,h,sample,distance,selected,crosses_grid,crosses_file,boundary_grid,boundary_file,analytic,grid,header,roundtrip,original,matched_grid,matched_header,matched_float,tau,file_tau,quantization_bound,passed\n";
    std::vector<Prepared> prepared; j::array checks;
    PrepareCase(output,"fold",atoms,targets,width,cutoff,0.1,origin,
        j::value_to<std::array<int,3>>(settings.at("grid_size")),original.get(),jobs,samples,prepared,checks);
    const auto failed{std::find_if(targets.begin(),targets.end(),[&](const auto & t){return atoms.atom_list[t.atom].serial_id==100;})};
    if (failed==targets.end()) throw std::runtime_error("Missing serial 100.");
    for (double h : {0.1,0.05,0.025})
    {
        Position local_origin; std::array<int,3> size;
        for (std::size_t k=0;k<3;++k)
        {
            double lo{std::numeric_limits<double>::max()},hi{std::numeric_limits<double>::lowest()};
            for (const auto & p : failed->points) {lo=std::min(lo,p.position[k]);hi=std::max(hi,p.position[k]);}
            const double first{std::floor((lo-origin[k])/h)-2};
            local_origin[k]=origin[k]+first*h; size[k]=static_cast<int>(std::ceil((hi-origin[k])/h)-first)+4;
        }
        PrepareCase(output,"failure-100",atoms,{*failed},width,cutoff,h,local_origin,size,nullptr,jobs,samples,prepared,checks);
    }
    for (int mode=0;mode<4;++mode) for (double h : {0.1,0.05,0.025})
    {
        sim::SimulationAtomPreparationResult small;
        small.atom_list.push_back(sim::SimulationAtom{.serial_id=1,.element=Element::CARBON,
            .position={0.013,0.021,0.037},.charge_used=mode==0 ? 0.0 : mode==1 ? 0.3 : -0.3});
        if (mode==3)
        {
            small.atom_list.push_back(sim::SimulationAtom{.serial_id=2,.element=Element::OXYGEN,.position={1.1,0.4,-0.2},.charge_used=0.3});
            small.atom_list.push_back(sim::SimulationAtom{.serial_id=3,.element=Element::NITROGEN,.position={-1.2,0.3,0.2},.charge_used=0.0});
        }
        Target target{0,{}};
        for (double r : {0.0,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0,2.49,2.5,2.51,2.8})
            for (std::size_t axis=0;axis<3;++axis) for (int sign : {-1,1})
            {auto p{small.atom_list[0].position};p[axis]+=sign*r;target.points.push_back({r,p,true});}
        target.points.push_back({std::sqrt(SquareDistance({-3,-3,-3},small.atom_list[0].position)),{-3,-3,-3},true});
        const int n{static_cast<int>(std::round(6/h))+1};
        PrepareCase(output,"small-"+std::to_string(mode),small,{target},0.5,2.5,h,{-3,-3,-3},{n,n,n},nullptr,jobs,samples,prepared,checks);
    }
    // No fitting begins until every forward case has passed.
    Write(output/"forward-status.json",j::object{{"passed",true},{"capture_replay",true},{"cases",checks.size()}});
    FitPrepared(output,prepared,jobs);
}
} // namespace second_stage_test::matched
