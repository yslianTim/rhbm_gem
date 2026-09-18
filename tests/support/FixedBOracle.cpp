#include "support/FixedBOracle.hpp"
#include "support/AtomCenteredVoxelUnion.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>

namespace second_stage_test::matched::fixed_b {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
using Sparse=Eigen::SparseMatrix<double>;
constexpr double eps=std::numeric_limits<double>::epsilon();
j::value Number(double v) {return std::isfinite(v) ? j::value(v) : j::value(nullptr);}
j::array Vector(const Eigen::VectorXd & v)
{
    j::array out; for (double n:v) out.push_back(Number(n)); return out;
}
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("Missing oracle input: "+path.string());
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},options);
}
void Write(const fs::path & path,const j::value & value)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';
}
std::ofstream CSV(const fs::path & path,const std::string & header)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<header<<'\n'<<std::setprecision(17); return out;
}
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
double Difference(const Eigen::VectorXd & a,const Eigen::VectorXd & b)
{return ((a-b).array().abs()/(1+a.array().abs().max(b.array().abs()))).maxCoeff();}
j::object Endpoint(const Sparse & x,const Eigen::VectorXd & y,const joint_ac::LinearResult & result,double seconds)
{
    auto out=Certificate(x,y,result.beta);
    out["beta"]=Vector(result.beta); out["valid"]=result.valid; out["reason"]=result.reason;
    out["free_rank"]=result.rank; out["linear_solves"]=result.solves;
    out["constraint_releases"]=result.releases; out["seconds"]=seconds; return out;
}
}

j::object Certificate(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & beta)
{
    j::object out{{"feasible",false},{"kkt_passed",false},{"projected_kkt",nullptr}};
    if (x.rows()!=y.size() || x.cols()!=beta.size() || !beta.allFinite() || !y.allFinite()) return out;
    Eigen::VectorXd norms(x.cols());
    for (Eigen::Index k=0;k<x.cols();++k) norms(k)=x.col(k).norm();
    if (!norms.allFinite() || (norms.array()<=0).any()) return out;
    const Eigen::VectorXd residual=x*beta-y;
    if (!residual.allFinite()) return out;
    const double s=std::max(1.0,y.norm());
    const Eigen::VectorXd u=norms.array()*beta.array()/s;
    const Eigen::VectorXd gradient=(x.transpose()*residual).array()/norms.array()/s;
    Eigen::VectorXd projected=u-gradient; j::array active; bool feasible=true;
    for (Eigen::Index k=0;k<beta.size();k+=2)
    {
        feasible &= beta(k)>=0; projected(k)=std::max(0.0,projected(k));
        if (beta(k)==0) active.push_back(k/2);
    }
    const double kkt=(u-projected).lpNorm<Eigen::Infinity>();
    out["feasible"]=feasible; out["projected_kkt"]=Number(kkt); out["kkt_passed"]=feasible && kkt<=1e-10;
    out["rss"]=Number(residual.squaredNorm()); out["objective"]=Number(.5*residual.squaredNorm());
    out["residual_scale"]=Number(residual.squaredNorm()/static_cast<double>(y.size()));
    out["residual_rmse"]=Number(residual.norm()/std::sqrt(static_cast<double>(y.size())));
    out["residual_max"]=Number(residual.cwiseAbs().maxCoeff()); out["relative_residual"]=Number(residual.norm()/s);
    out["active_atoms"]=active; return out;
}

j::object Fit(const Sparse & x,const Eigen::VectorXd & y,const j::object & spectrum)
{
    const auto weights=Eigen::VectorXd::Ones(y.size());
    auto start=std::chrono::steady_clock::now();
    const auto primary=joint_ac::WeightedSolve(x,y,weights);
    const double primary_seconds=Seconds(start); start=std::chrono::steady_clock::now();
    const auto reference=joint_ac::WeightedSolve(x,y,weights,true);
    auto p=Endpoint(x,y,primary,primary_seconds),r=Endpoint(x,y,reference,Seconds(start));
    const double difference=Difference(primary.beta,reference.beta);
    const bool qualified=primary.valid && reference.valid && p.at("kkt_passed").as_bool() &&
        r.at("kkt_passed").as_bool() && difference<=1e-10 && j::value_to<int>(spectrum.at("rank"))==x.cols();
    return {{"schema_version",1},{"experiment","fixed-b-oracle"},{"alpha",0},{"qualified",qualified},
        {"primary",p},{"reference",r},{"scaled_reference_difference",Number(difference)},
        {"design_spectrum",spectrum},{"design_nonzeros",x.nonZeros()},
        {"linear_solver","sparse-qr-householder-1024"},{"reference_solver","independent-tsqr-8192-svd"},
        {"variance_semantics","descriptive RSS/N; zero permitted"},{"kkt_tolerance",1e-10},{"reference_tolerance",1e-10}};
}

Data Prepare(const std::string & manifest_path,const std::string & map_path,
    const std::string & checkpoint_path,const std::string & output_path,const std::string & experiment)
{
    Eigen::setNbThreads(1);
    const fs::path output(output_path);
    if (fs::exists(output/"dataset.json") || fs::exists(output/"fits") || fs::exists(output/"residuals"))
        throw std::runtime_error("Oracle output already exists.");
    const auto manifest=Read(manifest_path),checkpoint=Read(checkpoint_path);
    const auto & settings=manifest.at("settings");
    if (settings.at("potential_model")!="single_gaus" || settings.at("blurring_width")!=.5 ||
        settings.at("cutoff_distance")!=2.5 || manifest.at("atoms").as_array().size()!=168 ||
        checkpoint.at("schema_version")!=1 || checkpoint.at("atoms").as_array().size()!=168 ||
        checkpoint.at("state").as_array().size()!=168 ||
        sim::FileSha256(map_path)!=j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Invalid oracle manifest/map/checkpoint.");
    const auto original=rhbm_gem::ReadMap(map_path);
    auto generation_ptr=std::make_unique<rhbm_gem::MapObject>(j::value_to<std::array<int,3>>(settings.at("grid_size")),
        j::value_to<Position>(settings.at("grid_spacing")),j::value_to<Position>(settings.at("origin")));
    auto & generation=*generation_ptr;
    std::vector<Atom> truth,states; sim::SimulationAtomPreparationResult generator;
    j::array identities; Eigen::VectorXd truth_beta(336),checkpoint_b(168);
    std::set<int> serials;
    // Manifest preparation order controls generator accumulation, independently of X.
    for (std::size_t a=0;a<168;++a)
    {
        const auto & atom=manifest.at("atoms").at(a);
        if (j::value_to<std::size_t>(atom.at("preparation_index"))!=a) throw std::runtime_error("Changed preparation order.");
        const int serial=j::value_to<int>(atom.at("serial_id"));
        if (!serials.insert(serial).second) throw std::runtime_error("Duplicate atom identity.");
        const auto & captured=checkpoint.at("atoms").as_array();
        const auto found=std::find_if(captured.begin(),captured.end(),[&](const auto & c){return c.at("identity").at("serial_id")==atom.at("serial_id");});
        if (found==captured.end()) throw std::runtime_error("Missing checkpoint identity.");
        const auto & identity=found->at("identity");
        for (const char * key:{"serial_id","chain_id","sequence_id","component_id","atom_id","alternate_indicator","position"})
            if (identity.at(key)!=atom.at(key)) throw std::runtime_error("Checkpoint identity mismatch.");
        const auto index=static_cast<std::size_t>(found-captured.begin());
        if (j::value_to<std::size_t>(found->at("index"))!=index) throw std::runtime_error("Checkpoint index mismatch.");
        const double width=j::value_to<double>(checkpoint.at("state").at(index).at(1));
        if (!(width>0) || !std::isfinite(width)) throw std::runtime_error("Invalid checkpoint B.");
        const auto position=j::value_to<Position>(atom.at("position"));
        const int element=j::value_to<int>(atom.at("element"));
        const double charge=j::value_to<double>(atom.at("charge_used"));
        truth.push_back({position,static_cast<double>(element),.5,charge});
        states.push_back({position,0,width,0}); identities.push_back(identity);
        truth_beta(static_cast<Eigen::Index>(2*a))=element; truth_beta(static_cast<Eigen::Index>(2*a+1))=charge;
        checkpoint_b(static_cast<Eigen::Index>(a))=width;
        generator.atom_list.push_back(sim::SimulationAtom{.serial_id=serial,.element=static_cast<Element>(element),.position=position,.charge_used=charge});
    }
    fs::create_directories(output/"fits"); fs::create_directories(output/"residuals");
    j::array completed; Write(output/"completion.json",j::object{{"complete",false},{"cases",completed}});
    const auto start=std::chrono::steady_clock::now();
    auto grid=atom_union::BuildGrid(states,{},generation,*original);
    rhbm_gem::core::MapSimulationRequest request; request.job_count=1; request.cutoff_distance=2.5;
    request.potential_model_choice=rhbm_gem::core::PotentialModel::SINGLE_GAUS;
    const int jobs=sim::PopulateMapValueArray(generation,generator,request,.5);
    Eigen::VectorXd y64(static_cast<Eigen::Index>(grid.voxels.size())),y32(y64.size());
    auto table=CSV(output/"voxels.csv","row,index,x,y,z,multiplicity,nearest_distance,reference_double,observed,quantization_delta");
    double maximum_forward{},maximum_quantization{}; std::size_t memberships{};
    for (std::size_t p=0;p<grid.voxels.size();++p)
    {
        const auto & v=grid.voxels[p]; double absolute{};
        const double direct=unique_grid::Direct(v.position,truth,2.5,&absolute),reference=generation.GetMapValue(v.index);
        const double error=std::abs(reference-direct),q=v.observed-reference;
        if (error>512*eps*std::max(1.0,absolute) || static_cast<double>(static_cast<float>(reference))!=v.observed)
            throw std::runtime_error("Oracle forward/quantization mismatch at voxel "+std::to_string(v.index));
        maximum_forward=std::max(maximum_forward,error); maximum_quantization=std::max(maximum_quantization,std::abs(q));
        memberships+=v.multiplicity; y64(static_cast<Eigen::Index>(p))=reference; y32(static_cast<Eigen::Index>(p))=v.observed;
        table<<p<<','<<v.index<<','<<v.position[0]<<','<<v.position[1]<<','<<v.position[2]<<','<<v.multiplicity<<','<<v.nearest_distance<<','<<reference<<','<<v.observed<<','<<q<<'\n';
    }
    table.close();
    Write(output/"dataset.json",j::object{{"schema_version",1},{"experiment",experiment},{"atoms",identities},
        {"row_count",grid.voxels.size()},{"memberships",memberships},{"radius",2.5},{"membership_geometry","generation"},
        {"grid_size",settings.at("grid_size")},{"generation_origin",settings.at("origin")},{"generation_spacing",settings.at("grid_spacing")},
        {"header_origin",j::value_from(original->GetOrigin())},{"header_spacing",j::value_from(original->GetGridSpacing())},
        {"checkpoint_b",Vector(checkpoint_b)},{"voxel_table_sha256",sim::FileSha256((output/"voxels.csv").string())}});
    Write(output/"scoring-truth.json",j::object{{"beta",Vector(truth_beta)},{"width",.5}});
    Write(output/"forward-status.json",j::object{{"passed",true},{"maximum_double_difference",maximum_forward},
        {"maximum_quantized_difference",0},{"maximum_quantization_delta",maximum_quantization},{"generator_jobs",jobs},
        {"voxel_count",grid.voxels.size()},{"seconds",Seconds(start)}});
    return {std::move(grid),std::move(states),std::move(generation_ptr),std::move(y64),std::move(y32),std::move(checkpoint_b),std::move(identities)};
}

void Run(const std::string & manifest_path,const std::string & map_path,
    const std::string & checkpoint_path,const std::string & output_path)
{
    const auto data=Prepare(manifest_path,map_path,checkpoint_path,output_path);
    const fs::path output(output_path); j::array completed;
    const auto & grid=data.grid; const auto & states=data.atoms;
    const auto & generation=*data.generation;
    const auto & y64=data.y64; const auto & y32=data.y32; const auto & checkpoint_b=data.checkpoint_b;
    for (bool oracle:{true,false})
    {
        auto model=states; if (oracle) for (auto & atom:model) atom.width=.5;
        auto clock=std::chrono::steady_clock::now();
        const auto x=atom_union::BuildDesign(grid,model,generation);
        const double build_seconds=Seconds(clock); clock=std::chrono::steady_clock::now();
        const auto spectrum=joint_ac::SparseSpectrum(x,Eigen::VectorXd::Ones(x.rows()));
        const double spectrum_seconds=Seconds(clock);
        for (bool quantized:{false,true})
        {
            const std::string name=std::string(oracle ? "true-b-" : "checkpoint-b-")+(quantized ? "float32" : "double");
            const auto & y=quantized ? y32 : y64; auto fit=Fit(x,y,spectrum);
            fit["case"]=name; fit["row_count"]=x.rows(); fit["design_build_seconds"]=build_seconds;
            fit["spectrum_seconds"]=spectrum_seconds; fit["b"]=oracle ? Vector(Eigen::VectorXd::Constant(168,.5)) : Vector(checkpoint_b);
            Eigen::VectorXd beta(336);
            for (std::size_t k=0;k<336;++k) beta(static_cast<Eigen::Index>(k))=j::value_to<double>(fit.at("primary").at("beta").at(k));
            const Eigen::VectorXd prediction=x*beta;
            auto residuals=CSV(output/"residuals"/(name+".csv"),"row,prediction,residual");
            for (Eigen::Index p=0;p<y.size();++p) residuals<<p<<','<<prediction(p)<<','<<y(p)-prediction(p)<<'\n';
            residuals.close(); fit["residual_sha256"]=sim::FileSha256((output/"residuals"/(name+".csv")).string());
            Write(output/"fits"/(name+".json"),fit); completed.emplace_back(name);
            Write(output/"completion.json",j::object{{"complete",completed.size()==4},{"cases",completed}});
        }
    }
}
} // namespace second_stage_test::matched::fixed_b
