#include "support/MatchedJointAC.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <algorithm>
#include <bit>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

namespace second_stage_test::matched::joint_ac {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if (!in) throw std::runtime_error("Missing joint AC input: "+path.string());
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},options);
}
void Write(const fs::path & path,const j::value & value)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';
}
Eigen::VectorXd Values(const j::value & v)
{
    Eigen::VectorXd out(static_cast<Eigen::Index>(v.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=v.at(static_cast<std::size_t>(k)).is_null() ?
        std::numeric_limits<double>::quiet_NaN() : j::value_to<double>(v.at(static_cast<std::size_t>(k)));
    return out;
}
bool Crosses(const Stencil & stencil,const std::vector<Atom> & atoms,double cutoff)
{
    for (const auto & atom:atoms) for (double radius:{cutoff,std::min(2.5,cutoff)})
    {
        bool in{},out{};
        for (const auto & slot:stencil.slots)
        {
            if (slot.coefficient==0) continue;
            if (SquareDistance(slot.position,atom.position)<=radius*radius) in=true; else out=true;
        }
        if (in&&out) return true;
    }
    return false;
}
} // namespace

Components BuildComponents(const std::vector<Stencil> & stencils,const std::vector<Atom> & atoms,double cutoff,
    const std::vector<std::size_t> & owners)
{
    if (owners.size()!=stencils.size()) throw std::invalid_argument("Missing sample ownership.");
    Components out; std::vector<std::size_t> parents(atoms.size()); std::iota(parents.begin(),parents.end(),0);
    auto root=[&](std::size_t k) {while (parents[k]!=k) {parents[k]=parents[parents[k]]; k=parents[k];} return k;};
    for (const auto & stencil:stencils)
    {
        auto & contributors{out.contributors.emplace_back()};
        for (std::size_t a=0;a<atoms.size();++a)
            if (std::any_of(stencil.slots.begin(),stencil.slots.end(),[&](const Slot & slot) {
                return slot.coefficient!=0 && SquareDistance(slot.position,atoms[a].position)<=cutoff*cutoff;
            })) contributors.push_back(a);
        if (contributors.empty()) throw std::runtime_error("Sample has no structural contributor.");
        for (std::size_t k=1;k<contributors.size();++k) parents[root(contributors[k])]=root(contributors[0]);
    }
    // All rows sharing a nuisance scale must belong to the same component.
    for (std::size_t p=0;p<stencils.size();++p)
    {
        if (owners[p]>=atoms.size()) throw std::invalid_argument("Invalid sample owner.");
        parents[root(out.contributors[p][0])]=root(owners[p]);
    }
    std::map<std::size_t,std::size_t> ids;
    for (std::size_t a=0;a<atoms.size();++a)
    {
        const auto r{root(a)};
        if (!ids.contains(r)) {ids[r]=out.atoms.size(); out.atoms.emplace_back(); out.rows.emplace_back();}
        out.atoms[ids.at(r)].push_back(a);
    }
    for (std::size_t p=0;p<stencils.size();++p)
    {
        const auto id{ids.at(root(out.contributors[p][0]))};
        for (const auto a:out.contributors[p]) if (ids.at(root(a))!=id)
            throw std::runtime_error("Cross-component structural contribution.");
        out.rows[id].push_back(p);
    }
    return out;
}

Eigen::MatrixXd BuildDesign(const std::vector<Stencil> & stencils,const std::vector<Atom> & atoms,
    const std::vector<std::size_t> & rows,const std::vector<std::size_t> & members,double cutoff)
{
    Eigen::MatrixXd x{Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rows.size()),static_cast<Eigen::Index>(2*members.size()))};
    for (std::size_t i=0;i<rows.size();++i) for (std::size_t a=0;a<members.size();++a)
    {
        const auto & atom{atoms.at(members[a])}; double g{},c{};
        for (const auto & slot:stencils.at(rows[i]).slots)
        {
            if (slot.coefficient==0) continue;
            const auto b{EvaluateBasis(SquareDistance(slot.position,atom.position),atom.width,cutoff)};
            g+=slot.coefficient*b.gaussian; c+=slot.coefficient*b.charge;
        }
        x(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(2*a))=g;
        x(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(2*a+1))=c;
    }
    return x;
}

void Run(const std::string & manifest_path,const std::string & map_path,
    const std::string & index_path,const std::string & output_path)
{
    Eigen::setNbThreads(1);
    const fs::path output{output_path};
    if (fs::exists(output/"fits")||fs::exists(output/"samples")) throw std::runtime_error("Joint output already exists.");
    const auto manifest{Read(manifest_path)},index{Read(index_path)};
    if (manifest.at("schema_version")!=1) throw std::runtime_error("Historical experiment only supports frozen manifest v1.");
    const auto & settings{manifest.at("settings")};
    if (settings.at("potential_model")!="single_gaus" || index.at("schema_version")!=1 ||
        sim::FileSha256(map_path)!=j::value_to<std::string>(manifest.at("output").at("map_sha256")))
        throw std::runtime_error("Unsupported manifest or map hash.");
    const auto original{rhbm_gem::ReadMap(map_path)};
    const double cutoff{j::value_to<double>(settings.at("cutoff_distance"))};
    rhbm_gem::MapObject generation(j::value_to<std::array<int,3>>(settings.at("grid_size")),
        j::value_to<Position>(settings.at("grid_spacing")),j::value_to<Position>(settings.at("origin")));
    const int jobs{std::getenv("OMP_NUM_THREADS") ? std::stoi(std::getenv("OMP_NUM_THREADS")) : 4};
    if (jobs!=1 && jobs!=4) throw std::runtime_error("Require jobs 1 or 4.");
    std::vector<j::value> contexts; std::set<std::string> ids;
    for (const auto & entry:index.at("states").as_array())
    {
        const auto path{j::value_to<std::string>(entry.at("context"))};
        if (!ids.insert(j::value_to<std::string>(entry.at("id"))).second ||
            sim::FileSha256(path)!=j::value_to<std::string>(entry.at("sha256"))) throw std::runtime_error("State hash or ID mismatch.");
        contexts.push_back(Read(path)); const auto & c{contexts.back()};
        if (c.at("schema_version")!=1 || c.at("state").as_array().size()!=168 ||
            c.at("atoms").as_array().size()!=168 || c.at("phase")!=entry.at("phase") || c.at("attempt")!=entry.at("attempt"))
            throw std::runtime_error("Invalid joint checkpoint.");
    }
    if (contexts.size()!=8) throw std::runtime_error("Require eight checkpoint states.");
    const auto selected{j::value_to<std::vector<std::string>>(index.at("selected_state_ids"))};
    if (selected.empty() || std::set<std::string>(selected.begin(),selected.end()).size()!=selected.size())
        throw std::runtime_error("Invalid selected states.");
    for (const auto & id:selected) if (!ids.contains(id)) throw std::runtime_error("Unknown selected state.");
    std::vector<double> alphas;
    std::vector<Atom> truth;
    std::vector<Stencil> stencils;
    std::vector<std::size_t> owners;
    std::vector<bool> signal,crosses;
    std::vector<double> observed,truth_prediction,truth_self;
    j::array atoms_json,rows_json;
    for (std::size_t a=0;a<168;++a)
    {
        const auto & capture{contexts[0].at("atoms").at(a)};
        const auto & identity{capture.at("identity")};
        const auto & atoms{manifest.at("atoms").as_array()};
        const auto it{std::find_if(atoms.begin(),atoms.end(),[&](const auto & atom){return atom.at("serial_id")==identity.at("serial_id");})};
        if (it==atoms.end()) throw std::runtime_error("Missing manifest identity.");
        for (const char * key:{"serial_id","chain_id","sequence_id","component_id","atom_id","alternate_indicator","position"})
            if (it->at(key)!=identity.at(key)) throw std::runtime_error("Identity mismatch.");
        for (const auto & context:contexts) if (context.at("atoms").at(a).at("index")!=a ||
            context.at("atoms").at(a).at("identity")!=identity || context.at("atoms").at(a).at("samples")!=capture.at("samples"))
            throw std::runtime_error("Changed sample population.");
        const double alpha{j::value_to<double>(capture.at("alpha"))};
        if (!std::isfinite(alpha) || alpha<0) throw std::runtime_error("Invalid checkpoint alpha.");
        for (const auto & context:contexts)
            if (std::bit_cast<std::uint64_t>(j::value_to<double>(context.at("atoms").at(a).at("alpha")))!=std::bit_cast<std::uint64_t>(alpha))
                throw std::runtime_error("Checkpoint alpha mapping differs.");
        alphas.push_back(alpha);
        // Truth is used exclusively in the independent forward/scoring paths below.
        truth.push_back({j::value_to<Position>(identity.at("position")),
            static_cast<double>(ChemicalDataHelper::GetAtomicNumber(static_cast<Element>(j::value_to<int>(it->at("element"))))),
            j::value_to<double>(settings.at("blurring_width")),j::value_to<double>(it->at("charge_used"))});
        atoms_json.push_back(identity);
        SamplingPointList points;
        for (const auto & sample:capture.at("samples").as_array())
            points.push_back({j::value_to<double>(sample.at("distance")),j::value_to<Position>(sample.at("position")),sample.at("selected").as_bool()});
        if (points.size()!=200) throw std::runtime_error("Changed sample count.");
        const auto replay{SampleExperimentPoints(*original,points)};
        for (std::size_t p=0;p<points.size();++p)
        {
            const double y{j::value_to<double>(capture.at("samples").at(p).at("response"))};
            if (!std::isfinite(y) || y!=replay[p].response) throw std::runtime_error("Map replay failed.");
            stencils.push_back(MakeStencil(generation,*original,points[p].position)); owners.push_back(a);
            signal.push_back(points[p].distance<=1); observed.push_back(y);
            rows_json.emplace_back(j::object{{"owner",a},{"sample",p},{"position",j::value_from(points[p].position)},
                {"response",y},{"signal",static_cast<bool>(signal.back())},{"selected",points[p].is_selected},{"boundary",stencils.back().boundary}});
        }
    }
    double max_forward{},max_bound_excess{};
    for (std::size_t p=0;p<stencils.size();++p)
    {
        double absolute{},bound{};
        truth_prediction.push_back(Predict(stencils[p],truth,cutoff,false,&absolute,&bound));
        const double quantized{Predict(stencils[p],truth,cutoff,true)};
        const double tau{512*std::numeric_limits<double>::epsilon()*std::max(1.0,absolute)};
        const double error{std::abs(quantized-observed[p])};
        const double excess{std::abs(truth_prediction.back()-observed[p])-bound-tau};
        if (error>tau || excess>0) throw std::runtime_error("Independent matched forward check failed.");
        max_forward=std::max(max_forward,error); max_bound_excess=std::max(max_bound_excess,excess);
        truth_self.push_back(Predict(stencils[p],{truth[owners[p]]},cutoff));
    }
    std::vector<Atom> geometry;
    for (std::size_t a=0;a<168;++a)
    {
        const auto abc{j::value_to<std::array<double,3>>(contexts[0].at("state").at(a))};
        geometry.push_back({j::value_to<Position>(atoms_json[a].at("position")),abc[0],abc[1],abc[2]});
    }
    const auto components{BuildComponents(stencils,geometry,cutoff,owners)};
    for (std::size_t p=0;p<stencils.size();++p)
    {
        crosses.push_back(Crosses(stencils[p],geometry,cutoff));
        rows_json[p].as_object()["cutoff_crossing"]=static_cast<bool>(crosses.back());
        rows_json[p].as_object()["contributors"]=j::value_from(components.contributors[p]);
    }
    fs::create_directories(output/"fits"); fs::create_directories(output/"samples");
    Write(output/"dataset.json",j::object{{"atoms",atoms_json},{"rows",rows_json},
        {"alphas",j::value_from(alphas)},{"component_atoms",j::value_from(components.atoms)},{"component_rows",j::value_from(components.rows)}});
    Write(output/"forward-status.json",j::object{{"passed",true},{"rows",stencils.size()},
        {"maximum_float_difference",max_forward},{"maximum_bound_excess",max_bound_excess}});
    std::vector<std::vector<std::string>> state_cases(contexts.size());
    std::vector<std::exception_ptr> state_errors(contexts.size());
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1) num_threads(jobs)
#endif
    for (std::size_t state_index=0;state_index<contexts.size();++state_index)
    {
        try
        {
        const auto state_id{j::value_to<std::string>(index.at("states").at(state_index).at("id"))};
        if (std::find(selected.begin(),selected.end(),state_id)==selected.end()) continue;
        auto state{geometry};
        for (std::size_t a=0;a<168;++a)
        {
            const auto abc{j::value_to<std::array<double,3>>(contexts[state_index].at("state").at(a))};
            if (!std::isfinite(abc[0]) || abc[0]<0 || !std::isfinite(abc[1]) || abc[1]<=0 || !std::isfinite(abc[2]))
                throw std::runtime_error("Invalid checkpoint ABC.");
            state[a].amplitude=abc[0]; state[a].width=abc[1]; state[a].charge=abc[2];
        }
        for (std::size_t component=0;component<components.atoms.size();++component)
        {
            const auto & members{components.atoms[component]}; const auto & rows{components.rows[component]};
            const auto x{BuildDesign(stencils,state,rows,members,cutoff)};
            Eigen::VectorXd y(static_cast<Eigen::Index>(rows.size())),input(static_cast<Eigen::Index>(2*members.size()));
            for (std::size_t p=0;p<rows.size();++p) y(static_cast<Eigen::Index>(p))=observed[rows[p]];
            for (std::size_t a=0;a<members.size();++a) {input(static_cast<Eigen::Index>(2*a))=state[members[a]].amplitude; input(static_cast<Eigen::Index>(2*a+1))=state[members[a]].charge;}
            const Eigen::VectorXd prediction{x*input};
            double max_design_error{};
            for (std::size_t p=0;p<rows.size();++p)
            {
                double absolute{}; const double direct{Predict(stencils[rows[p]],state,cutoff,false,&absolute)};
                const double delta{std::abs(direct-prediction(static_cast<Eigen::Index>(p)))};
                max_design_error=std::max(max_design_error,delta);
                if (delta>1024*std::numeric_limits<double>::epsilon()*std::max(1.0,absolute)) throw std::runtime_error("X beta prediction mismatch.");
            }
            Blocks blocks;
            for (const auto owner:members)
            {
                Block block{owner,alphas[owner],{}};
                for (std::size_t p=0;p<rows.size();++p) if (owners[rows[p]]==owner) block.rows.push_back(static_cast<Eigen::Index>(p));
                if (!block.rows.empty()) blocks.push_back(std::move(block));
            }
            const Eigen::VectorXd v0{BlockVariances(y-prediction,blocks)};
            const auto input_evidence{Evaluate(x,y,input,v0,blocks)};
            const bool common_valid{input_evidence.reason!="exact-fit-boundary" && v0.allFinite() && (v0.array()>0).all()};
            {
                const auto name{std::to_string(state_index)+"-"+std::to_string(component)};
                std::cout<<"case "<<name<<" "<<state_id<<" checkpoint atom-block alpha"<<std::endl;
                // Parallelize independent checkpoint fits in the outer loop so
                // long joint fits cannot leave the other workers idle.
                std::vector<j::object> fits(members.size()+1);
                std::vector<std::exception_ptr> errors(fits.size());
                for (std::size_t task=0;task<fits.size();++task)
                {
                    try
                    {
                        if (task==0) fits[task]=Fit(x,y,input,blocks);
                        else
                        {
                            const Eigen::Index col{static_cast<Eigen::Index>(2*(task-1))};
                            const Eigen::MatrixXd local{x.middleCols(col,2)};
                            const Eigen::VectorXd response{y-prediction+local*input.segment(col,2)};
                            fits[task]=Fit(local,response,input.segment(col,2),blocks);
                        }
                        fits[task]["state_id"]=state_id; fits[task]["component"]=component;
                        fits[task]["mode"]=task==0 ? "joint" : "frozen";
                        if (task>0) fits[task]["atom_index"]=members[task-1];
                        Write(output/"fits"/(name+"-"+std::to_string(task)+".json"),fits[task]);
                    }
                    catch (...) {errors[task]=std::current_exception();}
                }
                for (const auto & error:errors) if (error) std::rethrow_exception(error);
                Eigen::VectorXd frozen{input},joint{Values(fits[0].at("beta"))};
                for (std::size_t a=0;a<members.size();++a) frozen.segment(static_cast<Eigen::Index>(2*a),2)=Values(fits[a+1].at("beta"));
                const std::array<Eigen::VectorXd,3> betas{input,frozen,joint};
                std::array<Eigen::VectorXd,3> predictions{x*input,x*frozen,x*joint};
                j::array assembled;
                for (std::size_t mode=0;mode<3;++mode)
                {
                    j::array abc;
                    for (std::size_t a=0;a<members.size();++a) abc.emplace_back(j::array{
                        betas[mode](static_cast<Eigen::Index>(2*a)),state[members[a]].width,betas[mode](static_cast<Eigen::Index>(2*a+1))});
                    assembled.emplace_back(j::object{{"mode",std::array{"input","frozen","joint"}[mode]},
                        {"abc",abc},{"common_scale_objective",common_valid ? j::value(Objective(y-predictions[mode],v0,blocks)) : j::value(nullptr)},
                        {"common_scale_status",common_valid ? "valid" : "variance-boundary"},
                        {"residual_rmse",std::sqrt((y-predictions[mode]).squaredNorm()/static_cast<double>(y.size()))}});
                }
                std::ofstream samples(output/"samples"/(name+".csv")); samples.exceptions(std::ios::failbit|std::ios::badbit);
                samples<<"row,owner,sample,signal,cutoff_crossing,boundary,observed,truth_prediction,input_residual,frozen_residual,joint_residual,input_neighbor_error,frozen_neighbor_error,joint_neighbor_error\n"<<std::setprecision(17);
                for (std::size_t p=0;p<rows.size();++p)
                {
                    const auto row{rows[p]},owner{owners[row]};
                    const auto it{std::find(members.begin(),members.end(),owner)};
                    if (it==members.end()) throw std::runtime_error("Sample owner outside component.");
                    const auto col{static_cast<Eigen::Index>(2*static_cast<std::size_t>(it-members.begin()))};
                    const auto pi{static_cast<Eigen::Index>(p)};
                    samples<<row<<','<<owner<<','<<row%200<<','<<signal[row]<<','<<crosses[row]<<','<<stencils[row].boundary<<','<<y(pi)<<','<<truth_prediction[row];
                    for (std::size_t mode=0;mode<3;++mode) samples<<','<<y(pi)-predictions[mode](pi);
                    for (std::size_t mode=0;mode<3;++mode)
                    {
                        const double self{x(pi,col)*betas[mode](col)+x(pi,col+1)*betas[mode](col+1)};
                        samples<<','<<(predictions[mode](pi)-self)-(truth_prediction[row]-truth_self[row]);
                    }
                    samples<<'\n';
                }
                j::object info{{"id",name},{"state_id",state_id},{"component",component},{"blocks",fits[0].at("blocks")},
                    {"members",j::value_from(members)},{"row_count",rows.size()},{"common_variances",j::value_from(std::vector<double>(v0.data(),v0.data()+v0.size()))},
                    {"maximum_design_difference",max_design_error},{"assembled",assembled}};
                Write(output/("case-"+name+".json"),info); state_cases[state_index].push_back(name);
            }
        }
        }
        catch (...) {state_errors[state_index]=std::current_exception();}
    }
    j::array case_index;
    for (std::size_t s=0;s<contexts.size();++s)
    {
        if (state_errors[s]) std::rethrow_exception(state_errors[s]);
        for (const auto & name:state_cases[s]) case_index.emplace_back(name);
    }
    const auto joint_count{case_index.size()},frozen_count{selected.size()*geometry.size()};
    Write(output/"fit-index.json",j::object{{"schema_version",2},{"states",index.at("states")},{"selected_state_ids",index.at("selected_state_ids")},
        {"cases",case_index},{"joint_fits",joint_count},{"frozen_fits",frozen_count}});
}
} // namespace second_stage_test::matched::joint_ac
