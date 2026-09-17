#include "support/UniqueStencilGrid.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

namespace second_stage_test::matched::unique_grid {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
constexpr double eps{std::numeric_limits<double>::epsilon()};
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("Missing unique-grid input: "+path.string());
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
j::value Number(double x) {return std::isfinite(x) ? j::value(x) : j::value(nullptr);}
Eigen::VectorXd Values(const j::value & values)
{
    Eigen::VectorXd out(static_cast<Eigen::Index>(values.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(values.at(static_cast<std::size_t>(k)));
    return out;
}
void GlobalLabels(j::value & value)
{
    if (value.is_array()) for (auto & v:value.as_array()) GlobalLabels(v);
    if (!value.is_object()) return;
    auto & o{value.as_object()};
    if (o.contains("owner")) {o.erase("owner"); o["scope"]="global";}
    if (o.contains("failure_owner"))
    {
        const bool failed{j::value_to<int>(o.at("failure_owner"))>=0};
        o.erase("failure_owner"); o["failure_block"]=failed ? j::value("global") : j::value(nullptr);
    }
    for (auto & v:o) GlobalLabels(v.value());
}
j::object Diagnostics(const Eigen::MatrixXd & x,const Eigen::VectorXd & weights)
{
    if (weights.size()!=x.rows() || !weights.allFinite() || weights.sum()<=0 || (weights.array()<0).any())
        return {{"available",false},{"reason","invalid-weights"}};
    Eigen::MatrixXd z{weights.cwiseSqrt().asDiagonal()*x};
    for (Eigen::Index k=0;k<z.cols();++k) {const double n{z.col(k).norm()}; if (n>0) z.col(k)/=n;}
    const Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::HouseholderQRPreconditioner> svd(z);
    const auto s{svd.singularValues()};
    const double threshold{eps*static_cast<double>(std::max(x.rows(),x.cols()))};
    std::vector<double> sorted(weights.data(),weights.data()+weights.size());
    std::sort(sorted.begin(),sorted.end(),std::greater<double>());
    const auto top{static_cast<std::size_t>(std::ceil(.01*static_cast<double>(sorted.size())))};
    j::array quantiles;
    for (double q:{0.0,.01,.1,.5,.9,.99,1.0})
    {
        const auto k{static_cast<std::size_t>(q*static_cast<double>(sorted.size()-1))};
        quantiles.emplace_back(j::array{q,sorted[sorted.size()-1-k]});
    }
    return {{"available",true},{"rank",static_cast<int>((s.array()>threshold*s(0)).count())},
        {"minimum_singular",Number(s.tail(1)(0))},{"condition",Number(s(0)/s.tail(1)(0))},
        {"effective_n",Number(weights.sum()*weights.sum()/weights.squaredNorm())},
        {"maximum_row_share",weights.maxCoeff()/weights.sum()},
        {"top_one_percent_rows",top},{"top_one_percent_share",std::accumulate(sorted.begin(),sorted.begin()+static_cast<std::ptrdiff_t>(top),0.0)/weights.sum()},
        {"underflow_count",static_cast<std::size_t>((weights.array()==0).count())},{"weight_quantiles",quantiles}};
}
j::object Residuals(const fs::path & output,const std::string & name,const Grid & grid,
    const std::vector<Stencil> & stencils,const Eigen::VectorXd & sample_y,const Eigen::VectorXd & y,
    const Eigen::VectorXd & prediction,const Eigen::VectorXd * weights)
{
    auto voxels{CSV(output/"residuals"/(name+".csv"),"row,prediction,residual,weight")};
    for (Eigen::Index p=0;p<y.size();++p)
    {
        voxels<<p<<','<<prediction(p)<<','<<y(p)-prediction(p)<<',';
        if (weights && weights->size()==y.size()) voxels<<(*weights)(p);
        voxels<<'\n';
    }
    const Eigen::VectorXd sampled{Project(grid,stencils,prediction)};
    auto samples{CSV(output/"residuals"/(name+"-samples.csv"),"sample,prediction,residual")};
    for (Eigen::Index p=0;p<sample_y.size();++p) samples<<p<<','<<sampled(p)<<','<<sample_y(p)-sampled(p)<<'\n';
    return {{"grid_rmse",std::sqrt((y-prediction).squaredNorm()/static_cast<double>(y.size()))},
        {"matched_sample_rmse",std::sqrt((sample_y-sampled).squaredNorm()/static_cast<double>(sample_y.size()))}};
}
} // namespace

Grid BuildGrid(const std::vector<Stencil> & stencils,const rhbm_gem::MapObject & map)
{
    Grid out; std::map<std::size_t,Voxel> nodes;
    for (const auto & stencil:stencils) for (const auto & slot:stencil.slots)
    {
        if (slot.index>=map.GetMapValueArraySize()) throw std::invalid_argument("Stencil voxel outside map.");
        auto [it,inserted]{nodes.try_emplace(slot.index,Voxel{slot.index,0,slot.position,map.GetMapValue(slot.index)})};
        if ((!inserted && it->second.position!=slot.position) || !std::isfinite(it->second.observed))
            throw std::invalid_argument("Inconsistent voxel position or nonfinite observation.");
        ++it->second.multiplicity;
    }
    std::map<std::size_t,std::size_t> rows;
    for (const auto & [index,voxel]:nodes) {rows[index]=out.voxels.size(); out.voxels.push_back(voxel);}
    for (const auto & stencil:stencils)
    {
        auto & row{out.sample_rows.emplace_back()};
        for (std::size_t k=0;k<64;++k) row[k]=rows.at(stencil.slots[k].index);
    }
    return out;
}
double Direct(const Position & p,const std::vector<Atom> & atoms,double cutoff,double * absolute)
{
    double value{},sum{};
    for (const auto & atom:atoms)
    {
        const auto b{EvaluateBasis(SquareDistance(p,atom.position),atom.width,cutoff)};
        const double contribution{atom.amplitude*b.gaussian+atom.charge*b.charge};
        value+=contribution; sum+=std::abs(atom.amplitude*b.gaussian)+std::abs(atom.charge*b.charge);
    }
    if (absolute) *absolute=sum;
    return value;
}
Eigen::MatrixXd BuildDesign(const Grid & grid,const std::vector<Atom> & atoms,double cutoff)
{
    Eigen::MatrixXd x(static_cast<Eigen::Index>(grid.voxels.size()),static_cast<Eigen::Index>(2*atoms.size()));
    for (std::size_t a=0;a<atoms.size();++a) for (std::size_t p=0;p<grid.voxels.size();++p)
    {
        const auto b{EvaluateBasis(SquareDistance(grid.voxels[p].position,atoms[a].position),atoms[a].width,cutoff)};
        x(static_cast<Eigen::Index>(p),static_cast<Eigen::Index>(2*a))=b.gaussian;
        x(static_cast<Eigen::Index>(p),static_cast<Eigen::Index>(2*a+1))=b.charge;
    }
    return x;
}
Eigen::VectorXd Project(const Grid & grid,const std::vector<Stencil> & stencils,const Eigen::VectorXd & values)
{
    if (grid.sample_rows.size()!=stencils.size() || values.size()!=static_cast<Eigen::Index>(grid.voxels.size()))
        throw std::invalid_argument("Grid projection population mismatch.");
    Eigen::VectorXd out=Eigen::VectorXd::Zero(static_cast<Eigen::Index>(stencils.size()));
    for (std::size_t p=0;p<stencils.size();++p) for (std::size_t k=0;k<64;++k)
        out(static_cast<Eigen::Index>(p))+=stencils[p].slots[k].coefficient*values(static_cast<Eigen::Index>(grid.sample_rows[p][k]));
    return out;
}
joint_ac::Blocks GlobalBlock(Eigen::Index rows,double alpha)
{
    joint_ac::Block block{0,alpha,{}}; block.rows.resize(static_cast<std::size_t>(rows));
    std::iota(block.rows.begin(),block.rows.end(),0); return {std::move(block)};
}
j::object Fit(const Eigen::MatrixXd & x,const Eigen::VectorXd & y,const Eigen::VectorXd & initial,double alpha,
    const Eigen::SparseMatrix<double> * sparse_design)
{
    const auto blocks{GlobalBlock(y.size(),alpha)};
    auto result{joint_ac::Fit(x,y,initial,blocks,iteration_budget,refinement_budget,true,sparse_design)};
    result["experiment"]="unique-stencil-grid"; result["alpha"]=alpha;
    result["svd_preconditioner"]="householder-qr";
    result["linear_solver"]=sparse_design ? "sparse-qr" : "dense-qr";
    if (sparse_design) result["sparse_row_reduction"]="householder-qr-1024";
    if (sparse_design) result["design_nonzeros"]=sparse_design->nonZeros();
    result["iteration_budget"]=iteration_budget; result["refinement_budget"]=refinement_budget;
    for (auto & branch:result.at("branches").as_array())
    {
        const auto & endpoint{branch.at("primary")};
        const auto evidence{joint_ac::Evaluate(x,y,Values(endpoint.at("beta")),Values(endpoint.at("variances")),blocks)};
        branch.as_object()["endpoint_diagnostics"]=Diagnostics(x,evidence.linear_weights);
    }
    j::value labeled{std::move(result)}; GlobalLabels(labeled); return std::move(labeled.as_object());
}

void Run(const std::string & manifest_path,const std::string & map_path,
    const std::string & index_path,const std::string & output_path)
{
    Eigen::setNbThreads(1);
    const fs::path output{output_path};
    if (fs::exists(output/"fits") || fs::exists(output/"residuals")) throw std::runtime_error("Unique-grid output already exists.");
    const auto manifest{Read(manifest_path)},index{Read(index_path)};
    const auto & settings{manifest.at("settings")};
    if (settings.at("potential_model")!="single_gaus" || index.at("schema_version")!=1 ||
        sim::FileSha256(map_path)!=j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Unsupported unique-grid manifest or map hash.");
    const auto original{rhbm_gem::ReadMap(map_path)};
    const double cutoff{j::value_to<double>(settings.at("cutoff_distance"))};
    rhbm_gem::MapObject generation(j::value_to<std::array<int,3>>(settings.at("grid_size")),
        j::value_to<Position>(settings.at("grid_spacing")),j::value_to<Position>(settings.at("origin")));
    std::vector<j::value> contexts; std::set<std::string> ids;
    for (const auto & entry:index.at("states").as_array())
    {
        const auto path{j::value_to<std::string>(entry.at("context"))};
        if (!ids.insert(j::value_to<std::string>(entry.at("id"))).second ||
            sim::FileSha256(path)!=j::value_to<std::string>(entry.at("sha256"))) throw std::runtime_error("State identity/hash mismatch.");
        contexts.push_back(Read(path)); const auto & c{contexts.back()};
        if (c.at("schema_version")!=1 || c.at("state").as_array().size()!=168 || c.at("atoms").as_array().size()!=168 ||
            c.at("phase")!=entry.at("phase") || c.at("attempt")!=entry.at("attempt")) throw std::runtime_error("Invalid checkpoint.");
    }
    if (contexts.size()!=8 || manifest.at("atoms").as_array().size()!=168) throw std::runtime_error("Require complete fold-168 population.");
    const auto selected{j::value_to<std::vector<std::string>>(index.at("selected_state_ids"))};
    if (selected.empty() || std::set<std::string>(selected.begin(),selected.end()).size()!=selected.size()) throw std::runtime_error("Invalid selected states.");
    for (const auto & id:selected) if (!ids.contains(id)) throw std::runtime_error("Unknown selected state.");
    std::vector<Atom> geometry,truth; sim::SimulationAtomPreparationResult generator_atoms;
    std::vector<Stencil> stencils; std::vector<double> observations; j::array atoms_json,samples_json;
    std::set<int> serials;
    for (std::size_t a=0;a<168;++a)
    {
        const auto & captured{contexts[0].at("atoms").at(a)}, & identity{captured.at("identity")};
        const auto serial{j::value_to<int>(identity.at("serial_id"))};
        const auto & entries{manifest.at("atoms").as_array()};
        const auto it{std::find_if(entries.begin(),entries.end(),[&](const auto & entry){return entry.at("serial_id")==identity.at("serial_id");})};
        if (it==entries.end() || !serials.insert(serial).second) throw std::runtime_error("Missing/duplicate atom identity.");
        for (const char * key:{"serial_id","chain_id","sequence_id","component_id","atom_id","alternate_indicator","position"})
            if (it->at(key)!=identity.at(key)) throw std::runtime_error("Atom identity mismatch.");
        for (const auto & c:contexts) if (c.at("atoms").at(a).at("index")!=a || c.at("atoms").at(a).at("identity")!=identity ||
            c.at("atoms").at(a).at("samples")!=captured.at("samples")) throw std::runtime_error("Changed sample population.");
        const auto position{j::value_to<Position>(identity.at("position"))};
        geometry.push_back({position,0,1,0}); atoms_json.push_back(identity);
        const auto element{static_cast<Element>(j::value_to<int>(it->at("element")))};
        truth.push_back({position,static_cast<double>(ChemicalDataHelper::GetAtomicNumber(element)),
            j::value_to<double>(settings.at("blurring_width")),j::value_to<double>(it->at("charge_used"))});
        // This independent generator reference is never passed to the fitting API.
        generator_atoms.atom_list.push_back(sim::SimulationAtom{.serial_id=serial,.element=element,.position=position,.charge_used=truth.back().charge});
        if (captured.at("samples").as_array().size()!=200) throw std::runtime_error("Changed sample count.");
        for (const auto & sample:captured.at("samples").as_array())
        {
            const auto p{j::value_to<Position>(sample.at("position"))};
            stencils.push_back(MakeStencil(generation,*original,p));
            observations.push_back(j::value_to<double>(sample.at("response")));
            samples_json.emplace_back(j::object{{"owner",a},{"sample",observations.size()-1},{"position",j::value_from(p)},
                {"response",observations.back()},{"selected",sample.at("selected")},{"boundary",stencils.back().boundary}});
        }
    }
    const auto grid{BuildGrid(stencils,*original)};
    if (grid.voxels.size()!=602995) throw std::runtime_error("Unique stencil population differs from preflight: "+std::to_string(grid.voxels.size()));
    Eigen::VectorXd y(static_cast<Eigen::Index>(grid.voxels.size()));
    for (std::size_t p=0;p<grid.voxels.size();++p) y(static_cast<Eigen::Index>(p))=grid.voxels[p].observed;
    const Eigen::Map<const Eigen::VectorXd> sample_y(observations.data(),static_cast<Eigen::Index>(observations.size()));
    const Eigen::VectorXd replay{Project(grid,stencils,y)};
    if ((replay-sample_y).cwiseAbs().maxCoeff()>512*eps*std::max(1.0,sample_y.cwiseAbs().maxCoeff())) throw std::runtime_error("Slot replay failed.");
    fs::create_directories(output/"fits"); fs::create_directories(output/"residuals");
    rhbm_gem::core::MapSimulationRequest request; request.job_count=1; request.cutoff_distance=cutoff;
    request.potential_model_choice=rhbm_gem::core::PotentialModel::SINGLE_GAUS;
    sim::PopulateMapValueArray(generation,generator_atoms,request,j::value_to<double>(settings.at("blurring_width")));
    auto voxels{CSV(output/"voxels.csv","row,index,x,y,z,observed,multiplicity,reference_double,quantization_bound")};
    double max_double{},max_quantized{},max_excess{};
    for (std::size_t p=0;p<grid.voxels.size();++p)
    {
        const auto & v{grid.voxels[p]}; double absolute{};
        const double predicted{Direct(v.position,truth,cutoff,&absolute)},reference{generation.GetMapValue(v.index)};
        const double rounded{static_cast<double>(static_cast<float>(reference))},bound{std::abs(rounded-reference)};
        const double tolerance{512*eps*std::max(1.0,absolute)};
        max_double=std::max(max_double,std::abs(predicted-reference)); max_quantized=std::max(max_quantized,std::abs(rounded-v.observed));
        max_excess=std::max(max_excess,std::abs(predicted-v.observed)-bound-tolerance);
        if (std::abs(predicted-reference)>tolerance || rounded!=v.observed || std::abs(predicted-v.observed)>bound+tolerance)
            throw std::runtime_error("Independent voxel forward check failed at "+std::to_string(v.index));
        voxels<<p<<','<<v.index<<','<<v.position[0]<<','<<v.position[1]<<','<<v.position[2]<<','<<v.observed<<','<<v.multiplicity<<','<<reference<<','<<bound<<'\n';
    }
    voxels.close();
    auto slots{CSV(output/"slots.csv","sample,slot,voxel_row,coefficient")};
    for (std::size_t p=0;p<stencils.size();++p) for (std::size_t k=0;k<64;++k)
        slots<<p<<','<<k<<','<<grid.sample_rows[p][k]<<','<<stencils[p].slots[k].coefficient<<'\n';
    slots.close();
    Write(output/"dataset.json",j::object{{"schema_version",1},{"experiment","unique-stencil-grid"},{"atoms",atoms_json},
        {"samples",samples_json},{"row_count",grid.voxels.size()},{"slot_count",64*stencils.size()},
        {"grid_size",j::value_from(original->GetGridSize())},{"generation_origin",j::value_from(generation.GetOrigin())},
        {"generation_spacing",j::value_from(generation.GetGridSpacing())},{"sampling_origin",j::value_from(original->GetOrigin())},
        {"sampling_spacing",j::value_from(original->GetGridSpacing())},{"cutoff",cutoff},{"alphas",j::value_from(alphas)},
        {"voxel_table_sha256",sim::FileSha256((output/"voxels.csv").string())},{"slot_table_sha256",sim::FileSha256((output/"slots.csv").string())}});
    Write(output/"forward-status.json",j::object{{"passed",true},{"voxel_count",grid.voxels.size()},{"sample_count",stencils.size()},
        {"maximum_double_difference",max_double},{"maximum_quantized_difference",max_quantized},{"maximum_bound_excess",max_excess},
        {"maximum_replay_difference",(replay-sample_y).cwiseAbs().maxCoeff()}});
    j::array cases;
    for (std::size_t s=0;s<contexts.size();++s)
    {
        const auto id{j::value_to<std::string>(index.at("states").at(s).at("id"))};
        if (std::find(selected.begin(),selected.end(),id)==selected.end()) continue;
        auto state{geometry}; Eigen::VectorXd initial(336);
        for (std::size_t a=0;a<168;++a)
        {
            const auto abc{j::value_to<std::array<double,3>>(contexts[s].at("state").at(a))};
            if (!std::isfinite(abc[0]) || abc[0]<0 || !std::isfinite(abc[1]) || abc[1]<=0 || !std::isfinite(abc[2])) throw std::runtime_error("Invalid checkpoint ABC.");
            state[a].amplitude=abc[0]; state[a].width=abc[1]; state[a].charge=abc[2];
            initial(static_cast<Eigen::Index>(2*a))=abc[0]; initial(static_cast<Eigen::Index>(2*a+1))=abc[2];
        }
        std::cout<<"unique-grid building state="<<id<<" rows="<<y.size()<<std::endl;
        const auto x{BuildDesign(grid,state,cutoff)}; const Eigen::VectorXd input_prediction{x*initial}; double max_design{};
        for (std::size_t p=0;p<grid.voxels.size();++p)
        {
            double absolute{}; const double direct{Direct(grid.voxels[p].position,state,cutoff,&absolute)};
            const double error{std::abs(input_prediction(static_cast<Eigen::Index>(p))-direct)};
            if (error>1024*eps*std::max(1.0,absolute)) throw std::runtime_error("Grid X beta mismatch.");
            max_design=std::max(max_design,error);
        }
        auto input{Residuals(output,id+"-input",grid,stencils,sample_y,y,input_prediction,nullptr)};
        input["abc"]=contexts[s].at("state"); input["maximum_design_difference"]=max_design;
        Write(output/(id+"-input.json"),input);
        // Exact zeros only: keep every nonzero basis coefficient, including
        // arbitrarily small values. All alpha branches share this cache.
        const Eigen::SparseMatrix<double> sparse_design{x.sparseView(0.0,0.0)};
        if (sparse_design.nonZeros()!=(x.array()!=0).count()) throw std::runtime_error("Sparse design lost nonzero coefficients.");
        std::cout<<"unique-grid sparse nonzeros="<<sparse_design.nonZeros()<<std::endl;
        for (std::size_t k=0;k<alphas.size();++k)
        {
            const std::string name{id+"-alpha-"+std::to_string(k)};
            std::cout<<"unique-grid fit="<<name<<" alpha="<<alphas[k]<<std::endl;
            auto fit{Fit(x,y,initial,alphas[k],&sparse_design)}; fit["state_id"]=id;
            const auto beta{Values(fit.at("beta"))}; const Eigen::VectorXd prediction{x*beta};
            Eigen::VectorXd weights;
            if (fit.contains("variances")) weights=joint_ac::Evaluate(x,y,beta,Values(fit.at("variances")),GlobalBlock(y.size(),alphas[k])).weights;
            fit["residuals"]=Residuals(output,name,grid,stencils,sample_y,y,prediction,&weights);
            j::array abc;
            for (std::size_t a=0;a<168;++a) abc.emplace_back(j::array{beta(static_cast<Eigen::Index>(2*a)),state[a].width,beta(static_cast<Eigen::Index>(2*a+1))});
            fit["abc"]=abc; Write(output/"fits"/(name+".json"),fit); cases.emplace_back(name);
            std::cout<<"unique-grid finished="<<name<<" reason="<<fit.at("reason")<<std::endl;
        }
    }
    Write(output/"fit-index.json",j::object{{"schema_version",1},{"experiment","unique-stencil-grid"},{"states",index.at("states")},
        {"selected_state_ids",index.at("selected_state_ids")},{"alphas",j::value_from(alphas)},{"cases",cases},
        {"iteration_budget",iteration_budget},{"refinement_budget",refinement_budget},{"joint_fits",cases.size()}});
}
} // namespace second_stage_test::matched::unique_grid
