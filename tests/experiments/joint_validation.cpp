#include "support/JointPartialSelection.hpp"
#include "support/JointPrecisionAudit.hpp"
#include "support/JointComponentChecks.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <fstream>
#include <iostream>

namespace {
namespace j=boost::json;
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
namespace p=second_stage_test::matched::joint_abc;
namespace r=second_stage_test::matched::runtime_json;
namespace audit=second_stage_test::matched::certification;
using V=Eigen::VectorXd;
using Path=std::filesystem::path;
j::value Read(const Path & path)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing input: "+path.string());
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},options);
}
void Write(const Path & path,const j::value & value)
{
    const auto temporary=Path(path.string()+".tmp");
    std::ofstream out(temporary); out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<j::serialize(value)<<'\n'; out.close(); std::filesystem::rename(temporary,path);
}
V Vector(const std::vector<double> & values) {return Eigen::Map<const V>(values.data(),static_cast<Eigen::Index>(values.size()));}
j::value Outcome(const c::JointFitResult & fit)
{
    j::parse_options options; options.numbers=j::number_precision::precise;
    auto result=j::parse(rhbm_gem::joint_result_io::Encode(c::CaptureJointAnalysisResult(fit)),{},options);
    if(fit.assembled_state) {
        const auto & state=result.at("assembled_state");
        if(j::value_to<std::vector<double>>(state.at("ac"))!=fit.assembled_state->ac ||
           j::value_to<std::vector<double>>(state.at("log_b"))!=fit.assembled_state->log_b)
            throw std::runtime_error("Experiment JSON changed the runtime state.");
    }
    return result;
}
j::object Census(const c::JointProblem & problem)
{
    const auto & data=c::JointProblemAccess::Get(problem);
    std::size_t memberships=0; for(const auto & s:problem.Input().support) memberships+=s.size();
    j::array sizes; for(const auto & view:data.partition.components) sizes.push_back(view.atoms.size());
    return {{"targets",problem.Input().selection_domain->target_indices.size()},
        {"contributors",problem.Input().atom_ids.size()},{"rows",problem.Input().row_ids.size()},
        {"memberships",memberships},{"component_sizes",sizes}};
}
void Weak(int start,const Path & output)
{
    auto fixture=joint_partial_test::Make("weak");
    const auto initialized=c::EstimateJointComponents(*fixture.map,*fixture.model);
    auto initial=initialized.initialization.b;
    if(start==1) initial={.5,.4}; else if(start==2) initial={.6,.5};
    else if(start!=0) throw std::invalid_argument("Weak start must be 0, 1 or 2.");
    const auto fit=start==0 ? initialized : c::FitJointComponents(*initialized.problem,initial);
    const auto & problem=*fit.problem; const auto & data=c::JointProblemAccess::Get(problem);
    const V y=data.y; const p::Domain domain(data.domain);
    j::array support; for(const auto & a:problem.Input().support) {
        j::array entries; for(const auto & s:a) entries.push_back(j::array{s.row,s.squared_distance});
        support.push_back(std::move(entries));
    }
    j::object snapshot{{"observations",j::value_from(problem.Input().observations)},
        {"support",support},{"atom_ids",j::value_from(problem.Input().atom_ids)},
        {"row_ids",j::value_from(problem.Input().row_ids)},{"scale",problem.ObservationScale()},
        {"truth_a",j::value_from(fixture.a)},{"truth_b",j::value_from(fixture.b)},
        {"truth_c",j::value_from(fixture.c)}};
    j::object out{{"start",start},{"initial_b",j::value_from(initial)},
        {"production_initial_b",j::value_from(initialized.initialization.b)},
        {"snapshot",snapshot},{"outcome",Outcome(fit)}};
    Write(output,out); // Retain the formal result if an offline audit is interrupted.
    if(!fit.assembled_state) return;
    const V eta=Vector(fit.assembled_state->log_b),beta=Vector(fit.assembled_state->ac);
    const auto e=p::AtState(domain,y,eta,beta,data.context);
    const auto assessment=n::AssessProfile(domain,y,eta,data.context,&beta);
    out["actual_state_assessment"]=r::Assessment(assessment);
    out["same_state"]=p::SameState(domain,y,eta,beta,data.partition,data.context);
    const auto dense=p::DenseDifferentiate(e,data.context.scale,&data.context);
    const auto tiled=p::MaterializeDerivative(e,data.context.scale,&data.context);
    out["dense_tiled_jacobian_relative_difference"]=(dense.jacobian-tiled.jacobian).norm()/std::max(1e-300,dense.jacobian.norm());
    n::Matrix directions(eta.size(),2);
    directions.col(0)=assessment.correction.normalized();
    directions.col(1)=assessment.weak_directions.col(0).normalized();
    out["directions"]=j::array{r::Values(directions.col(0)),r::Values(directions.col(1))};
    out["precision"]=audit::PrecisionAudit(domain,y,e,directions,&data.context);
    Write(output,out);
    j::array scans;
    const auto base=n::EvaluateProfile(domain,y,eta,false,&data.context);
    for(int d=0;d<2;++d) for(int exponent=-6;exponent<=-1;++exponent) for(double sign:{-1.,1.}) {
        const double step=sign*std::pow(10.,exponent); const V at=eta+step*directions.col(d);
        const auto trial=n::EvaluateProfile(domain,y,at,false,&data.context);
        j::object row{{"direction",d},{"step",step},{"eta",r::Values(at)},{"endpoint",r::Endpoint(trial)}};
        if(trial.valid) {
            row["trust"]=r::Trust(n::CheckTrust(domain,y,trial,data.context));
            row["objective_change"]=(trial.residual.squaredNorm()-base.residual.squaredNorm())/(2*data.context.scale*data.context.scale);
            row["prediction_change_norm"]=(trial.residual-base.residual).norm();
        }
        if(exponent==-3 || exponent==-1) row["precision_profile"]=audit::PrecisionProfileChange(domain,y,eta,at);
        scans.push_back(std::move(row)); out["scans"]=scans; Write(output,out);
    }
    auto extended=data.context; extended.profile_budget=1000; extended.update_budget=500;
    out["diagnostic_restart"]=p::FitComponent(data.partition.components.front(),y,eta.array().exp(),extended);
    out["audit_complete"]=true; Write(output,out);
}
std::unique_ptr<rhbm_gem::ModelObject> TwoAtoms(double distance)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(int i=0;i<2;++i) {
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(i+1);
        atom->SetElement(Element::CARBON); atom->SetPosition(i*distance,0,0);
        atom->SetChainID("A"); atom->SetComponentID("ALA"); atom->SetAtomID(i ? "CB" : "CA");
        atoms.push_back(std::move(atom));
    }
    auto model=std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
    model->SelectAtoms([](const auto & a){return a.GetSerialID()==1;}); return model;
}
void Statistical(const Path & input,const Path & output)
{
    const auto spec=Read(input); auto model=TwoAtoms(j::value_to<double>(spec.at("distance")));
    const auto values=j::value_to<std::vector<double>>(spec.at("values"));
    if(values.size()!=25*25*25) throw std::invalid_argument("Wrong statistical map size.");
    auto buffer=std::make_unique<double[]>(values.size()); std::copy(values.begin(),values.end(),buffer.get());
    rhbm_gem::MapObject map({25,25,25},{.5,.5,.5},{-6,-6,-6},std::move(buffer));
    const auto problem=c::BuildJointProblem(map,*model);
    const double distance=j::value_to<double>(spec.at("distance"))+j::value_to<double>(spec.at("shift"));
    const std::vector<second_stage_test::matched::Atom> truth{{{0,0,0},2,.5,.2},{{distance,0,0},2.3,.55,.15}};
    const auto clean=j::value_to<std::vector<double>>(spec.at("clean")); double error=0;
    for(std::size_t k=0;k<values.size();++k) {
        const std::array<double,3> position{.5*static_cast<double>(k%25)-6,.5*static_cast<double>((k/25)%25)-6,.5*static_cast<double>(k/625)-6};
        error=std::max(error,std::abs(clean.at(k)-second_stage_test::matched::unique_grid::Direct(position,truth,2.5)));
    }
    if(error>1e-12) throw std::runtime_error("Independent scalar generation disagreement.");
    j::object out{{"census",Census(problem)},{"generator_max_error",error},
        {"row_ids",j::value_from(problem.Input().row_ids)}};
    for(const auto & mode:spec.at("modes").as_array()) {
        const auto fit=mode=="fixed" ? c::FitJointComponents(problem,{.6,.6}) : c::EstimateJointComponents(map,*model);
        j::object record{{"outcome",Outcome(fit)},{"prediction",nullptr}};
        if(fit.prediction) record["prediction"]=j::value_from(*fit.prediction);
        out[std::string(mode.as_string())]=std::move(record);
    }
    Write(output,out);
}
void Inspect(const Path & model_path,const Path & map_path,const Path & output)
{
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path);
    model->SelectAllAtoms(); model->ApplySymmetrySelection(false); model->ApplyElementSelection(Element::HYDROGEN,true);
    auto out=Census(c::BuildJointProblem(*map,*model));
    out["model_atoms"]=model->GetAtomList().size();
    out["catalogue"]=std::count_if(model->GetAtomList().begin(),model->GetAtomList().end(),[](const auto & a){return a->GetElement()!=Element::HYDROGEN;}); out["map_size"]=j::value_from(map->GetGridSize());
    Write(output,out);
}
void Generate(const std::string & kind,int count,const Path & output)
{
    std::filesystem::create_directories(output);
    // Chain A is selected; B is the remote, unselected catalogue.
    std::ofstream cif(output/"input.cif");
    cif<<"data_control\n_em_3d_reconstruction.resolution 1.0\nloop_\n_entity.id\n_entity.type\n1 non-polymer\n#\nloop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\n";
    if(kind=="catalogue") cif<<"B 1\n";
    cif<<"#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n";
    const std::array<int,3> dims{258,258,258};
    auto values=std::make_unique<double[]>(258*258*258);
    for(int i=0;i<count;++i) {
        std::array<double,3> position{};
        if(kind=="catalogue") position=i<2 ? std::array<double,3>{8.+1.2*i,8,8} :
            std::array<double,3>{20.+3.*((i-2)%34),20.+3.*(((i-2)/34)%34),20.+3.*((i-2)/(34*34))};
        else if(kind=="single") position={20.+3.5*(i%13),20.+3.5*((i/13)%13),20.+3.5*(i/169)};
        else if(kind=="multi") {
            const int block=i/32,k=i%32;
            position={10.+20.*(block%5)+3.5*(k%4),10.+20.*((block/5)%4)+3.5*((k/4)%4),10.+20.*(block/20)+3.5*(k/16)};
        } else throw std::invalid_argument("Unknown resource geometry.");
        cif<<"HETATM "<<i+1<<" C C . UNK "<<(kind=="catalogue" && i>=2 ? "B" : "A")<<" 1 "<<i+1<<' '<<position[0]<<' '<<position[1]<<' '<<position[2]<<" 1 0 1\n";
        const std::vector<second_stage_test::matched::Atom> atom{{position,2,.5,.2}};
        std::array<int,3> lo{},hi{}; for(std::size_t k=0;k<3;++k) {lo[k]=std::max(0,static_cast<int>(std::floor((position[k]-2.5)/.5))); hi[k]=std::min(257,static_cast<int>(std::floor((position[k]+2.5)/.5)));}
        for(int z=lo[2];z<=hi[2];++z) for(int y=lo[1];y<=hi[1];++y) for(int x=lo[0];x<=hi[0];++x)
            values[static_cast<std::size_t>(x+258*(y+258*z))]+=second_stage_test::matched::unique_grid::Direct({.5*x,.5*y,.5*z},atom,2.5);
    }
    cif<<"#\n"; cif.close();
    rhbm_gem::MapObject map(dims,{.5,.5,.5},{0,0,0},std::move(values));
    rhbm_gem::WriteMap(output/"input.map",map);
}
}
int main(int argc,char ** argv)
{
    try {
        Eigen::setNbThreads(1);
        if(argc==4 && std::string(argv[1])=="weak") Weak(std::stoi(argv[2]),argv[3]);
        else if(argc==4 && std::string(argv[1])=="stat") Statistical(argv[2],argv[3]);
        else if(argc==5 && std::string(argv[1])=="inspect") Inspect(argv[2],argv[3],argv[4]);
        else if(argc==4 && (std::string(argv[1])=="catalogue" || std::string(argv[1])=="single" || std::string(argv[1])=="multi")) Generate(argv[1],std::stoi(argv[2]),argv[3]);
        else throw std::invalid_argument("Usage: joint_validation weak START OUTPUT | stat INPUT OUTPUT | inspect MODEL MAP OUTPUT | catalogue/single/multi COUNT DIRECTORY");
        return 0;
    } catch(const std::exception & e) {std::cerr<<e.what()<<'\n';return 1;}
}
