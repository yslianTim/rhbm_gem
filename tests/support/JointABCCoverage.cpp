#include "support/JointABCCoverage.hpp"
#include "support/AtomCenteredVoxelUnion.hpp"
#include "support/FixedBOracle.hpp"
#include "support/JointABCCertification.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sys/resource.h>

namespace second_stage_test::matched::coverage {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
using Vector=Eigen::VectorXd;
using Clock=std::chrono::steady_clock;
constexpr std::array<const char *,8> synthetic_names{
    "baseline","weak-1e-2","weak-1e-4","active-a","zero-signal","near-0.10","near-0.02","duplicate"};
j::value Number(double x) {return std::isfinite(x) ? j::value(x) : j::value(nullptr);}
j::array Values(const Vector & v) {j::array out; for (double x:v) out.push_back(Number(x)); return out;}
Vector Parse(const j::value & v)
{
    Vector out(static_cast<Eigen::Index>(v.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k)));
    return out;
}
j::object Resources(Clock::time_point start)
{
    struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    const auto bytes=usage.ru_maxrss;
#else
    const auto bytes=usage.ru_maxrss*1024;
#endif
    return {{"seconds",std::chrono::duration<double>(Clock::now()-start).count()},
        {"process_peak_rss_bytes",bytes}};
}
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("Missing coverage input: "+path.string());
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
// Independent analytic data generator: deliberately does not call EvaluateBasis.
double Forward(const Position & point,const std::vector<Atom> & atoms,double * absolute=nullptr)
{
    double value{},total{};
    for (const auto & atom:atoms)
    {
        const double square=SquareDistance(point,atom.position);
        if (square>6.25) continue;
        const double b=atom.width,r=std::sqrt(square);
        const double g=atom.amplitude*std::pow(2*M_PI*b*b,-1.5)*std::exp(-square/(2*b*b));
        const double c=atom.charge*(r<1e-5 ? std::sqrt(2/M_PI)/b : std::erf(r/b/std::sqrt(2.0))/r);
        value+=g+c; total+=std::abs(g)+std::abs(c);
    }
    if (absolute) *absolute=total;
    return value;
}
struct Data
{
    std::unique_ptr<rhbm_gem::ModelObject> model;
    std::unique_ptr<rhbm_gem::MapObject> generation,observation;
    std::vector<Atom> truth;
    j::array identities;
};
Data MainData(const std::string & model_path,const std::string & map_path,const j::value & manifest)
{
    const auto & settings=manifest.at("settings");
    if (manifest.at("atom_count")!=168 || settings.at("potential_model")!="single_gaus" ||
        settings.at("blurring_width")!=.5 || settings.at("cutoff_distance")!=2.5 ||
        sim::FileSha256(model_path)!=j::value_to<std::string>(manifest.at("source").at("model_sha256")) ||
        sim::FileSha256(map_path)!=j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Coverage model/map/manifest contract mismatch.");
    Data data; data.model=rhbm_gem::ReadModel(model_path); data.observation=rhbm_gem::ReadMap(map_path);
    data.generation=std::make_unique<rhbm_gem::MapObject>(j::value_to<std::array<int,3>>(settings.at("grid_size")),
        j::value_to<Position>(settings.at("grid_spacing")),j::value_to<Position>(settings.at("origin")));
    sim::SimulationAtomPreparationResult generator;
    std::set<int> serials;
    for (const auto & entry:manifest.at("atoms").as_array())
    {
        const int serial=j::value_to<int>(entry.at("serial_id")),element=j::value_to<int>(entry.at("element"));
        if (!serials.insert(serial).second || j::value_to<std::size_t>(entry.at("preparation_index"))!=data.truth.size())
            throw std::runtime_error("Coverage duplicate identity or changed preparation order.");
        const auto * atom=data.model->FindAtomPtr(serial);
        if (!atom) throw std::runtime_error("Coverage model identity missing.");
        const auto identity=Identity(*atom);
        for (const auto & field:identity) if (entry.at(field.key())!=field.value())
            throw std::runtime_error("Coverage model identity mismatch: "+std::to_string(serial));
        const auto position=j::value_to<Position>(entry.at("position"));
        const double charge=j::value_to<double>(entry.at("charge_used"));
        data.truth.push_back({position,static_cast<double>(element),ElementWidth(element),charge});
        data.identities.push_back(identity);
        generator.atom_list.push_back(sim::SimulationAtom{.serial_id=serial,.element=static_cast<Element>(element),
            .position=position,.charge_used=charge});
    }
    if (data.truth.size()!=168) throw std::runtime_error("Coverage atom population mismatch.");
    rhbm_gem::core::MapSimulationRequest request; request.job_count=1; request.cutoff_distance=2.5;
    request.potential_model_choice=rhbm_gem::core::PotentialModel::SINGLE_GAUS;
    if (sim::PopulateMapValueArray(*data.generation,generator,request,.5)!=1)
        throw std::runtime_error("Coverage generator must use one thread.");
    return data;
}
Data SyntheticData(const std::string & name)
{
    Data data; data.truth=SyntheticAtoms(name);
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for (std::size_t k=0;k<data.truth.size();++k)
    {
        const int element=6+static_cast<int>(k%3);
        auto atom=std::make_unique<rhbm_gem::AtomObject>();
        atom->SetSerialID(static_cast<int>(k+1)); atom->SetSequenceID(static_cast<int>(k+1));
        atom->SetChainID("T"); atom->SetComponentID("ALA"); atom->SetIndicator(".");
        atom->SetAtomID(element==6 ? "C" : element==7 ? "N" : "O");
        atom->SetElement(static_cast<Element>(element)); atom->SetPosition(data.truth[k].position);
        data.identities.push_back(Identity(*atom)); atoms.push_back(std::move(atom));
    }
    data.model=std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
    const std::array<int,3> size{33,29,29}; const Position spacing{.3,.3,.3},origin{-3.6,-3.6,-3.6};
    data.generation=std::make_unique<rhbm_gem::MapObject>(size,spacing,origin);
    auto y64=std::make_unique<double[]>(data.generation->GetMapValueArraySize());
    auto y32=std::make_unique<double[]>(data.generation->GetMapValueArraySize());
    for (std::size_t k=0;k<data.generation->GetMapValueArraySize();++k)
    {
        y64[k]=Forward(data.generation->GetGridPosition(k),data.truth);
        y32[k]=static_cast<double>(static_cast<float>(y64[k]));
    }
    data.generation->SetMapValueArray(std::move(y64));
    data.observation=std::make_unique<rhbm_gem::MapObject>(size,spacing,origin,std::move(y32));
    return data;
}
void RunDataset(const std::string & name,Data & data,const fs::path & root,bool certify)
{
    fs::path output=root;
    const auto preparation_start=Clock::now();
    fs::create_directories(output/"fits"); fs::create_directories(output/"residuals");
    fs::create_directories(output/"weak-directions"); fs::create_directories(output/"controls");
    Write(output/"completion.json",j::object{{"complete",false},{"cases",j::array{}}});
    auto grid=atom_union::BuildGrid(data.truth,{},*data.generation,*data.observation);
    if(certify)
    {
        const auto size=data.generation->GetGridSize();
        for(auto & voxel:grid.voxels)
        {
            const auto i=voxel.index; const auto nx=static_cast<std::size_t>(size[0]),ny=static_cast<std::size_t>(size[1]);
            const auto position=sim::GridPosition({static_cast<int>(i%nx),static_cast<int>(i/nx%ny),
                static_cast<int>(i/(nx*ny))},data.generation->GetGridSpacing(),data.generation->GetOrigin());
            if(position!=voxel.position) throw std::runtime_error("Legacy coordinates disagree with explicit FMA snapshot.");
            voxel.position=position;
        }
    }
    Vector y64(static_cast<Eigen::Index>(grid.voxels.size())),y32(y64.size());
    auto table=CSV(output/"voxels.csv","row,index,x,y,z,multiplicity,nearest_distance,reference_double,observed,quantization_delta");
    std::size_t memberships{}; double max_forward{},max_quantization{};
    // Validate all stored map values, including voxels outside the fitted ROI.
    if (data.generation->GetGridSize()!=data.observation->GetGridSize()) throw std::runtime_error("Map dimensions differ.");
    for (std::size_t k=0;k<data.generation->GetMapValueArraySize();++k)
    {
        double absolute{};
        const double direct=Forward(data.generation->GetGridPosition(k),data.truth,&absolute);
        const double generated=data.generation->GetMapValue(k),observed=data.observation->GetMapValue(k);
        const double difference=std::abs(direct-generated);
        if (!std::isfinite(generated) || !std::isfinite(observed) ||
            difference>512*std::numeric_limits<double>::epsilon()*std::max(1.0,absolute) ||
            static_cast<double>(static_cast<float>(generated))!=observed)
            throw std::runtime_error("Coverage forward/quantization mismatch at voxel "+std::to_string(k));
        max_forward=std::max(max_forward,difference); max_quantization=std::max(max_quantization,std::abs(observed-generated));
    }
    for (std::size_t k=0;k<grid.voxels.size();++k)
    {
        const auto & v=grid.voxels[k]; const double generated=data.generation->GetMapValue(v.index);
        y64(static_cast<Eigen::Index>(k))=generated; y32(static_cast<Eigen::Index>(k))=v.observed;
        memberships+=v.multiplicity;
        table<<k<<','<<v.index<<','<<v.position[0]<<','<<v.position[1]<<','<<v.position[2]<<','<<v.multiplicity<<','
            <<v.nearest_distance<<','<<generated<<','<<v.observed<<','<<v.observed-generated<<'\n';
    }
    table.close();
    Vector beta(static_cast<Eigen::Index>(2*data.truth.size())),width(static_cast<Eigen::Index>(data.truth.size()));
    for (std::size_t k=0;k<data.truth.size();++k)
    {
        beta(static_cast<Eigen::Index>(2*k))=data.truth[k].amplitude;
        beta(static_cast<Eigen::Index>(2*k+1))=data.truth[k].charge;
        width(static_cast<Eigen::Index>(k))=data.truth[k].width;
    }
    Write(output/"scoring-truth.json",j::object{{"beta",Values(beta)},{"b",Values(width)},
        {"known_unidentified",name=="zero-signal" || name=="duplicate"}});
    Write(output/"dataset.json",j::object{{"schema_version",1},{"experiment","joint-abc-coverage"},{"name",name},
        {"atoms",data.identities},{"row_count",grid.voxels.size()},{"memberships",memberships},{"radius",2.5},
        {"membership_geometry","generation"},{"grid_size",j::value_from(data.generation->GetGridSize())},
        {"generation_origin",j::value_from(data.generation->GetOrigin())},{"generation_spacing",j::value_from(data.generation->GetGridSpacing())},
        {"header_origin",j::value_from(data.observation->GetOrigin())},{"header_spacing",j::value_from(data.observation->GetGridSpacing())},
        {"voxel_table_sha256",sim::FileSha256((output/"voxels.csv").string())}});
    Write(output/"forward-status.json",j::object{{"passed",true},{"maximum_double_difference",max_forward},
        {"maximum_quantized_difference",0},{"maximum_quantization_delta",max_quantization},
        {"checked_map_voxels",data.generation->GetMapValueArraySize()},
        {"generator",name=="heterogeneous-168" ? "production-single-gaus" : "independent-analytic"},{"generator_jobs",1}});
    Write(output/"preparation-resources.json",Resources(preparation_start));
    const auto initialization=Initialize(*data.model,*data.observation,data.identities);
    Write(output/"initialization.json",initialization.evidence);
    // Optimizer domain contains geometry/support only; truth never supplies a start.
    const auto domain=certify ? certification::Snapshot(grid,data.truth,root) : joint_abc::Domain(grid,data.truth);
    if (name=="heterogeneous-168")
    {
        for (const std::string kind:{"true-b","first-stage-b"})
        {
            if (kind=="first-stage-b" && !initialization.valid) continue;
            auto atoms=data.truth; const Vector b=kind=="true-b" ? width : initialization.b;
            for (std::size_t k=0;k<atoms.size();++k) atoms[k].width=b(static_cast<Eigen::Index>(k));
            const auto x=atom_union::BuildDesign(grid,atoms,*data.generation);
            const auto spectrum=joint_ac::SparseSpectrum(x,Vector::Ones(y64.size()));
            for (bool quantized:{false,true})
            {
                const std::string label=kind+(quantized ? "-float32" : "-double");
                auto fit=fixed_b::Fit(x,quantized ? y32 : y64,spectrum); fit["b"]=Values(b); fit["case"]=label;
                Write(output/"controls"/(label+".json"),fit);
            }
        }
    }
    for (const std::string & variant:certify ? std::vector<std::string>{"legacy","guarded","guarded-log"} : std::vector<std::string>{""})
    {
    if (certify)
    {
        output=root/variant;
        for (const char * sub:{"fits","residuals","weak-directions","controls","audits"}) fs::create_directories(output/sub);
        for (const char * file:{"voxels.csv","dataset.json","scoring-truth.json","forward-status.json","initialization.json","snapshot.json","contributors.csv"})
            fs::copy_file(root/file,output/file);
        for (const auto & file:fs::directory_iterator(root/"controls")) fs::copy_file(file.path(),output/"controls"/file.path().filename());
    }
    j::array completed;
    for (bool quantized:{false,true}) for (const std::string start:{"first-stage","narrower","wider","mixed"})
    {
        const std::string label=start+(quantized ? "-float32" : "-double");
        std::cout<<"Starting "<<variant<<'/'<<name<<'/'<<label<<std::endl;
        j::object fit;
        if (initialization.valid)
        {
            const auto initial=InitialWidths(initialization.b,data.identities,start);
            j::object resources;
            fit=joint_abc::Fit(domain,quantized ? y32 : y64,initial,&resources,variant);
            fit["initial_b"]=Values(initial); fit["resources"]=resources;
            fit["initialization_status"]="valid";
            if (fit.at("primary").at("valid").as_bool())
            {
                const auto & endpoint=fit.at("primary");
                const Vector coefficients=Parse(endpoint.at("beta")),b=Parse(endpoint.at("b"));
                auto atoms=data.truth;
                for (std::size_t k=0;k<atoms.size();++k) atoms[k].width=b(static_cast<Eigen::Index>(k));
                const auto x=certify ? joint_abc::Evaluate(domain,quantized ? y32 : y64,Parse(endpoint.at("eta"))).x :
                    atom_union::BuildDesign(grid,atoms,*data.generation);
                const Vector prediction=x*coefficients; const auto & y=quantized ? y32 : y64;
                auto residual=CSV(output/"residuals"/(label+".csv"),"row,prediction,residual");
                for (Eigen::Index k=0;k<y.size();++k) residual<<k<<','<<prediction(k)<<','<<prediction(k)-y(k)<<'\n';
                residual.close(); fit["residual_sha256"]=sim::FileSha256((output/"residuals"/(label+".csv")).string());
            }
            if (fit.contains("width_spectrum"))
            {
                auto weak=CSV(output/"weak-directions"/(label+".csv"),"serial_id,weakest,second,third");
                for (std::size_t k=0;k<data.truth.size();++k)
                {
                    weak<<j::value_to<int>(data.identities[k].at("serial_id"));
                    for (std::size_t d=0;d<3;++d) weak<<','<<j::value_to<double>(fit.at("width_spectrum").at("weak_directions").at(d).at(k));
                    weak<<'\n';
                }
            }
        }
        else fit={{"execution_complete",false},{"joint_qualified",false},
            {"initialization_status","invalid"},{"qualification_failure","initialization-failed"}};
        fit["case"]=label; fit["dataset"]=name;
        if(certify) {fit["variant"]=variant; fit["observation_snapshot_sha256"]=sim::FileSha256(root/"snapshot.json");}
        Write(output/"fits"/(label+".json"),fit); completed.emplace_back(label);
        Write(output/"completion.json",j::object{{"complete",completed.size()==8},{"cases",completed}});
        std::cout<<name<<'/'<<label<<" qualified="<<fit.at("joint_qualified")<<std::endl;
    }
    }
    if (certify) Write(root/"completion.json",j::object{{"complete",true},{"variants",j::array{"legacy","guarded","guarded-log"}}});
}
} // namespace

double ElementWidth(int element) {return element==8 ? .4 : element==7 ? .45 : .5;}

std::vector<Atom> SyntheticAtoms(const std::string & name)
{
    if (std::find(synthetic_names.begin(),synthetic_names.end(),name)==synthetic_names.end())
        throw std::invalid_argument("Unknown coverage synthetic dataset.");
    std::vector<Atom> atoms;
    for (int z=0;z<2;++z) for (int y=0;y<2;++y) for (int x=0;x<3;++x)
    {
        const int serial=1+x+3*y+6*z,element=6+(serial-1)%3;
        atoms.push_back({{1.2*x,1.2*y,1.2*z},static_cast<double>(element),ElementWidth(element),serial%2 ? .2 : -.2});
    }
    for (const std::size_t k:{1u,5u,9u})
    {
        if (name=="weak-1e-2" || name=="weak-1e-4")
        {
            const double factor=name=="weak-1e-2" ? .01 : .0001;
            atoms[k].amplitude*=factor; atoms[k].charge*=factor;
        }
        if (name=="active-a") {atoms[k].amplitude=0; atoms[k].charge=.2;}
    }
    if (name=="zero-signal") {atoms[5].amplitude=0; atoms[5].charge=0;}
    if (name=="near-0.10" || name=="near-0.02" || name=="duplicate")
    {atoms[3].position=atoms[0].position; atoms[3].position[0]+=name=="near-0.10" ? .1 : name=="near-0.02" ? .02 : 0;}
    return atoms;
}

j::object Identity(const rhbm_gem::AtomObject & atom)
{
    return {{"serial_id",atom.GetSerialID()},{"chain_id",atom.GetChainID()},{"sequence_id",atom.GetSequenceID()},
        {"component_id",atom.GetComponentID()},{"atom_id",atom.GetAtomID()},{"alternate_indicator",atom.GetIndicator()},
        {"position",j::value_from(atom.GetPosition())},{"element",static_cast<int>(atom.GetElement())}};
}

Initialization Initialize(rhbm_gem::ModelObject & model,rhbm_gem::MapObject & map,const j::array & identities)
{
    const auto start=Clock::now();
    Initialization result; result.b=Vector::Constant(static_cast<Eigen::Index>(identities.size()),std::numeric_limits<double>::quiet_NaN());
    result.evidence={{"source","first-stage-float32-map"},{"sampling_method","FibonacciDeterministic"},{"jobs",1},
        {"seed_abc",j::array{0,1,0}},{"uses_peeling",false},{"uses_truth",false},{"valid",false}};
    model.SelectAllAtoms(); model.ApplyElementSelection(Element::HYDROGEN,true);
    if (model.GetSelectedAtomCount()!=identities.size()) throw std::runtime_error("Initialization atom population mismatch.");
    std::set<int> serials;
    for (const auto & identity:identities)
    {
        const int serial=j::value_to<int>(identity.at("serial_id"));
        const auto * atom=model.FindAtomPtr(serial);
        if (!serials.insert(serial).second || !atom || Identity(*atom)!=identity)
            throw std::runtime_error("Initialization identity mismatch.");
    }
    model.EditAnalysis().InitializeFromSelection();
    rhbm_gem::core::RunPotentialSamplingWorkflow(map,model,SphereSamplingMethod::FibonacciDeterministic,1);
    model.EditAnalysis().InitializeLocalFittingSeedModels();
    rhbm_gem::core::FitOptions options; options.thread_size=1; options.quiet_mode=true; options.exclude_hydrogen=true;
    try
    {
        rhbm_gem::core::RunLocalAlphaTraining(model,options,rhbm_gem::FittingStage::First);
        rhbm_gem::core::RunFixedOffsetLocalFitting(model,options,rhbm_gem::FittingStage::First);
    }
    catch (const std::exception & error)
    {
        result.evidence["reason"]=std::string(error.what()); result.evidence["resources"]=Resources(start);
        result.evidence["atoms"]=j::array{}; result.evidence["b0"]=Values(result.b); return result;
    }
    j::array records; std::size_t warnings{};
    for (std::size_t k=0;k<identities.size();++k)
    {
        const auto * atom=model.FindAtomPtr(j::value_to<int>(identities[k].at("serial_id")));
        const auto view=rhbm_gem::AtomLocalPotentialView::For(*atom);
        const auto & local=view.GetGaussianResult(rhbm_gem::FittingStage::First);
        const auto samples=view.GetSamplingEntries(rhbm_gem::FittingStage::First);
        const auto abc=[](const auto & model_value) {return Values(model_value.GetModel().ToVector());};
        j::object row{{"identity",identities[k]},{"alpha",Number(local.alpha_r)},
            {"ols",abc(local.ols)},{"mdpde",abc(local.mdpde)},
            {"raw_sample_count",view.GetRawSamplingEntries(false).size()},{"sample_count",samples.size()}};
        if (local.fit_result)
        {
            const auto & fit=*local.fit_result;
            row["diagnostics"]=j::object{{"native_status",static_cast<int>(fit.status)},
                {"qualification",static_cast<int>(fit.Qualification())},{"iterations",fit.diagnostics.iterations},
                {"squared_beta_change",fit.diagnostics.squared_beta_change ? Number(*fit.diagnostics.squared_beta_change) : j::value(nullptr)},
                {"relative_variance_change",fit.diagnostics.relative_variance_change ? Number(*fit.diagnostics.relative_variance_change) : j::value(nullptr)},
                {"sigma_square",Number(fit.sigma_square)},{"beta_ols",Values(fit.beta_ols)},{"beta_mdpde",Values(fit.beta_mdpde)}};
            warnings+=fit.Qualification()==rhbm_gem::RHBMSolveQualification::Unqualified;
        }
        else {row["diagnostics"]=nullptr; ++warnings;}
        result.b(static_cast<Eigen::Index>(k))=local.mdpde.GetModel().GetWidth(); records.push_back(row);
    }
    result.valid=result.b.allFinite() && (result.b.array()>0).all();
    result.evidence["valid"]=result.valid; result.evidence["reason"]=result.valid ? "valid-widths" : "invalid-widths";
    result.evidence["unqualified_local_count"]=warnings; result.evidence["atoms"]=records;
    result.evidence["b0"]=Values(result.b); result.evidence["resources"]=Resources(start); return result;
}

Vector InitialWidths(const Vector & b,const j::array & identities,const std::string & start)
{
    if (b.size()!=static_cast<Eigen::Index>(identities.size()) || !b.allFinite() || (b.array()<=0).any())
        throw std::invalid_argument("Invalid first-stage initial widths.");
    if (start!="first-stage" && start!="narrower" && start!="wider" && start!="mixed")
        throw std::invalid_argument("Unknown coverage initialization.");
    Vector result=b;
    for (Eigen::Index k=0;k<result.size();++k)
    {
        const int serial=j::value_to<int>(identities.at(static_cast<std::size_t>(k)).at("serial_id"));
        if (start=="narrower" || (start=="mixed" && serial%2)) result(k)*=.8;
        else if (start=="wider" || start=="mixed") result(k)*=1.2;
    }
    return result;
}

void Run(const std::string & model,const std::string & map,const std::string & manifest,const std::string & output_path,bool certify)
{
    Eigen::setNbThreads(1); const fs::path output(output_path);
    if (fs::exists(output/"completion.json") || fs::exists(output/"datasets"))
        throw std::runtime_error("Coverage output already exists.");
    fs::create_directories(output); j::array completed;
    Write(output/"completion.json",j::object{{"complete",false},{"datasets",completed}});
    const auto metadata=Read(manifest);
    if (metadata.at("schema_version")==1)
    {
        const auto fixture=Read(fs::path(__FILE__).parent_path().parent_path()/"benchmarks/joint_abc_coverage.json");
        if (sim::FileSha256(manifest)!=j::value_to<std::string>(fixture.at("input_hashes").at("manifest")))
            throw std::runtime_error("Unknown v1 coverage width contract.");
    }
    else if (!certify || metadata.at("schema_version")!=2)
        throw std::runtime_error("Unsupported coverage manifest version.");
    else
    {
        const auto & kernel=metadata.at("kernel"); const auto & support=metadata.at("support");
        const j::object policy{{"id","element-scaled-v1"},{"oxygen",.8},{"nitrogen",.9},{"other",1.0}};
        if(kernel.at("version")!="single_gaus-v1" || kernel.at("width_policy")!=policy ||
            kernel.as_object().contains("effective_charge_width") || kernel.at("near_zero_distance")!=1e-5 ||
            kernel.at("charge_term_cutoff")!=2.5 || !kernel.at("minimum_charge_width").is_null() ||
            support.at("version")!="sphere-fma-v1" || support.at("outer_cutoff")!=2.5 ||
            support.at("coordinates")!="fma(index,spacing,origin)" ||
            support.at("squared_distance")!="fma(dz,dz,fma(dy,dy,dx*dx))" ||
            support.at("comparison")!="squared_distance<=cutoff*cutoff")
            throw std::runtime_error("Changed manifest v2 kernel/support policy.");
        for(const auto & atom:metadata.at("atoms").as_array())
            if(atom.at("effective_gaussian_width")!=ElementWidth(j::value_to<int>(atom.at("element"))) ||
                atom.at("effective_charge_width")!=atom.at("effective_gaussian_width"))
                throw std::runtime_error("Wrong per-atom effective widths.");
    }
    auto main=MainData(model,map,metadata);
    RunDataset("heterogeneous-168",main,output/"datasets"/"heterogeneous-168",certify); completed.push_back("heterogeneous-168");
    Write(output/"completion.json",j::object{{"complete",false},{"datasets",completed}});
    for (const auto * name:synthetic_names)
    {
        auto data=SyntheticData(name); RunDataset(name,data,output/"datasets"/name,certify); completed.push_back(name);
        Write(output/"completion.json",j::object{{"complete",completed.size()==9},{"datasets",completed}});
    }
}
} // namespace second_stage_test::matched::coverage
