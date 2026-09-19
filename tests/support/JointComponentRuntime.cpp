#include "support/JointComponentRuntime.hpp"
#include "support/JointRuntimeJson.hpp"
#include "support/JointABCPrecision.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace core=rhbm_gem::core;
namespace sim=core::simulation;
namespace j=boost::json;
namespace fs=std::filesystem;
using Vector=Eigen::VectorXd;
void Write(const fs::path & path,const j::value & value)
{fs::create_directories(path.parent_path()); std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';}
j::array Values(const std::vector<double> & v) {j::array out; for(double x:v) out.push_back(runtime_json::Number(x)); return out;}
Vector Parse(const j::value & v) {Vector out(static_cast<Eigen::Index>(v.as_array().size())); for(Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k))); return out;}
j::object State(const core::JointState & s)
{return {{"valid",true},{"beta",Values(s.ac)},{"eta",Values(s.log_b)},{"b",Values(s.b)}};}
j::object ComponentFit(const core::JointComponentResult & c,const std::string & name,const std::string & dataset,double scale)
{
    j::object out{{"component_id",c.id},{"case",name},{"dataset",dataset},{"residual_scale",scale},
        {"execution_complete",true},{"search_stopped_without_convergence",!c.search_completed},{"search_success",c.search_completed},
        {"usable_state",c.state.has_value()},{"joint_qualified",false},{"stop_reason",c.stop_reason},
        {"profile_evaluations",c.profile_evaluations},{"search_reference_evaluations",c.reference_evaluations},
        {"accepted_updates",c.accepted_updates},{"trials",j::array{}},{"last_trusted_state",nullptr},
        {"primary",j::object{{"valid",false}}}};
    if(c.state) {out["last_trusted_state"]=State(*c.state); out["primary"]=State(*c.state);}
    return out;
}
void Snapshot(const core::JointProblem & problem,const rhbm_gem::MapObject & map,const rhbm_gem::ModelObject & model,const fs::path & root,const std::string & dataset)
{
    const auto & input=problem.Input(); j::array atoms,offsets{0};
    for(const auto & id:input.atom_ids) atoms.push_back(j::object{{"serial_id",std::stoi(id)},
        {"position",j::value_from(model.FindAtomPtr(std::stoi(id))->GetPosition())}});
    std::vector<std::vector<std::pair<std::size_t,double>>> rows(input.observations.size());
    for(std::size_t a=0;a<input.support.size();++a) for(const auto & s:input.support[a]) rows[s.row].push_back({a,s.squared_distance});
    std::ofstream voxels(root/"voxels.csv"),contributors(root/"contributors.csv");
    voxels.exceptions(std::ios::failbit|std::ios::badbit); contributors.exceptions(std::ios::failbit|std::ios::badbit);
    voxels<<"row,index,x,y,z,multiplicity,nearest_distance,reference_double,observed\n"<<std::setprecision(17);
    contributors<<"row,atom,square\n"<<std::setprecision(17); std::size_t count{};
    for(std::size_t r=0;r<rows.size();++r)
    {
        const auto index=std::stoull(input.row_ids[r]);
        const auto position=sim::GridPosition(map.GetGridIndex(index),map.GetGridSpacing(),map.GetOrigin());
        double nearest=6.25;
        for(const auto & [atom,square]:rows[r]) {contributors<<r<<','<<atom<<','<<square<<'\n'; nearest=std::min(nearest,square); ++count;}
        offsets.push_back(count);
        voxels<<r<<','<<index<<','<<position[0]<<','<<position[1]<<','<<position[2]<<','<<rows[r].size()<<','<<std::sqrt(nearest)<<','<<input.observations[r]<<','<<input.observations[r]<<'\n';
    }
    voxels.close(); contributors.close();
    Write(root/"dataset.json",j::object{{"name",dataset},{"atoms",atoms},{"geometry_source","input-MapObject"},
        {"grid_size",j::value_from(map.GetGridSize())},{"spacing",j::value_from(map.GetGridSpacing())},{"origin",j::value_from(map.GetOrigin())}});
    Write(root/"snapshot.json",j::object{{"schema_version",1},{"support_policy","sphere-fma-v1"},
        {"rows",input.observations.size()},{"atoms",input.atom_ids.size()},{"memberships",count},{"row_offsets",offsets},
        {"voxels_sha256",sim::FileSha256(root/"voxels.csv")},{"contributors_sha256",sim::FileSha256(root/"contributors.csv")}});
}
void Fresh(rhbm_gem::MapObject & map,rhbm_gem::ModelObject & model,const fs::path & root,const std::string & precision,bool physical)
{
    if(fs::exists(root)) throw std::runtime_error("Use a fresh runtime output directory.");
    fs::create_directories(root); Eigen::setNbThreads(1);
    auto first=core::EstimateJointComponents(map,model);
    j::array initialization;
    for(const auto & a:first.initialization.atoms)
        initialization.push_back(j::object{{"id",a.id},{"ols",j::value_from(a.ols)},{"mdpde",j::value_from(a.mdpde)},
            {"alpha",a.alpha},{"sample_count",a.sample_count},{"native_status",a.native_status ? j::value(*a.native_status) : j::value(nullptr)}});
    Write(root/"initialization.json",j::object{{"valid",first.initialization.valid},{"reason",first.initialization.reason},
        {"source","fresh-first-stage-mdpde"},{"uses_truth",false},{"uses_peeling",false},
        {"b0",Values(first.initialization.b)},{"atoms",initialization},{"seconds",first.costs.initialization_seconds}});
    if(!first.initialization.valid) throw std::runtime_error("Runtime first-stage initialization failed: "+first.initialization.reason);
    const auto problem=*first.problem; const auto & data=core::JointProblemAccess::Get(problem);
    const std::string dataset=physical ? "physical-two" : "fresh-map";
    Snapshot(problem,map,model,root,dataset);
    auto context=data.context; context.snapshot_hash=sim::FileSha256(root/"snapshot.json");
    auto audit_context=context;
    audit_context.audit=RegisteredAudit(static_cast<Eigen::Index>(problem.Input().atom_ids.size()),"baseline","first-stage");
    audit_context.audit.cache_precision=true;
    if(physical) {audit_context.audit.precision=true; audit_context.audit.expanded_if_unverified=true; audit_context.audit.block_precision=true;}
    Write(root/"census.json",Census(data.domain,data.partition,context));
    if(physical && data.partition.components.size()!=2) throw std::runtime_error("Physical fixture did not produce two structural components.");
    j::array complete;
    for(const std::string start:{"first-stage","narrower","wider","mixed"})
    {
        auto initial=first.initialization.b;
        for(std::size_t a=0;a<initial.size();++a)
        {
            const int serial=std::stoi(problem.Input().atom_ids[a]);
            if(start=="narrower" || (start=="mixed" && serial%2)) initial[a]*=.8;
            else if(start=="wider" || start=="mixed") initial[a]*=1.2;
        }
        const auto fit=start=="first-stage" ? first : core::FitJointComponents(problem,initial);
        const std::string name=start+"-"+precision; const auto target=root/"cases"/name;
        const Vector b=Eigen::Map<const Vector>(initial.data(),static_cast<Eigen::Index>(initial.size()));
        auto mono=Fit(data.domain,data.y,b,nullptr,"guarded",&context); mono["case"]=name; mono["dataset"]=dataset;
        j::array components;
        for(std::size_t k=0;k<fit.components.size();++k)
        {
            auto child=ComponentFit(fit.components[k],name,dataset,context.scale);
            Write(target/"components"/(std::to_string(k)+".json"),child); components.push_back(child);
            WriteLocalAudit(data.partition.components[k],data.y,child,audit_context,target/"local-components"/std::to_string(k));
        }
        auto assembled=Assemble(data.domain,data.y,data.partition,context,components); assembled["case"]=name; assembled["dataset"]=dataset;
        if(fit.prediction.has_value()!=assembled.at("prediction_available").as_bool()) throw std::runtime_error("Runtime availability differs from reference assembly.");
        double prediction_difference{};
        if(fit.prediction)
        {
            const auto state=runtime::EvaluateState(data.domain,data.y,Parse(assembled.at("assembled_state").at("eta")),Parse(assembled.at("assembled_state").at("beta")),context);
            const Vector runtime_prediction=Eigen::Map<const Vector>(fit.prediction->data(),static_cast<Eigen::Index>(fit.prediction->size()));
            prediction_difference=(runtime_prediction-(state.x*state.beta)).lpNorm<Eigen::Infinity>();
            if(prediction_difference!=0 || !fit.objective || *fit.objective!=state.certificate.objective)
                throw std::runtime_error("Runtime result changed during audit assembly.");
        }
        Write(target/"monolithic-fit.json",mono); Write(target/"assembled-fit.json",assembled);
        certification::ResetPrecisionCache(); WriteEndpointAudit(data.domain,data.y,mono,audit_context,target/"monolithic");
        certification::ResetPrecisionCache(); WriteEndpointAudit(data.domain,data.y,assembled,audit_context,target/"assembled");
        Write(target/"runtime.json",j::object{{"initial_b",Values(initial)},{"search_completed",fit.search_completed},
            {"prediction_available",fit.prediction.has_value()},{"objective",fit.objective ? runtime_json::Number(*fit.objective) : j::value(nullptr)},
            {"regular_certificate","not-run"},{"prediction_reassembly_difference",prediction_difference},
            {"initialization_seconds",fit.costs.initialization_seconds},{"search_seconds",fit.costs.search_seconds},
            {"search_reference_seconds",fit.costs.search_reference_seconds},{"assessment_seconds",fit.costs.assessment_seconds},{"assembly_seconds",fit.costs.assembly_seconds}});
        complete.emplace_back(name);
    }
    Write(root/"completion.json",j::object{{"complete",true},{"cases",complete},{"requires_regular_parity",physical},
        {"historical_files_required",false},{"geometry_source","input-MapObject"}});
}
}
void RunJointRuntime(const std::string & model_path,const std::string & map_path,const std::string & output)
{
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path); model->SelectAllAtoms();
    Fresh(*map,*model,output,"float32",false);
}
void RunPhysicalJointRuntime(const std::string & output,bool inputs_only)
{
    if(fs::exists(output)) throw std::runtime_error("Use a fresh physical runtime directory.");
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms; sim::SimulationAtomPreparationResult generator;
    for(int block=0;block<2;++block) for(int z=0;z<2;++z) for(int y=0;y<2;++y) for(int x=0;x<3;++x)
    {
        const int serial=static_cast<int>(atoms.size())+1,element=6+(serial-1)%3;
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(serial); atom->SetSequenceID(serial);
        atom->SetChainID("T"); atom->SetComponentID("ALA"); atom->SetIndicator(".");
        atom->SetAtomID(element==6 ? "C" : element==7 ? "N" : "O"); atom->SetElement(static_cast<Element>(element));
        atom->SetPosition(12*block+1.2*x,1.2*y,1.2*z);
        generator.atom_list.push_back(sim::SimulationAtom{.serial_id=serial,.element=static_cast<Element>(element),
            .position=atom->GetPosition(),.charge_used=serial%2 ? .2 : -.2}); atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    rhbm_gem::MapObject map({73,29,29},{.3,.3,.3},{-3.6,-3.6,-3.6});
    core::MapSimulationRequest request; request.job_count=1; request.cutoff_distance=2.5; request.potential_model_choice=core::PotentialModel::SINGLE_GAUS;
    if(sim::PopulateMapValueArray(map,generator,request,.5)!=1) throw std::runtime_error("Physical generation must be single threaded.");
    fs::create_directories(output);
    std::ofstream cif(fs::path(output)/"model.cif");
    cif.exceptions(std::ios::failbit|std::ios::badbit);
    cif<<"data_joint_physical\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
       <<"_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
       <<"_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
       <<"_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n"<<std::setprecision(17);
    for(const auto & atom:model.GetAtomList())
    {
        const auto p=atom->GetPosition();
        cif<<"ATOM "<<atom->GetSerialID()<<' '<<atom->GetAtomID()<<' '<<atom->GetAtomID()
           <<" . ALA T "<<atom->GetSerialID()<<' '<<p[0]<<' '<<p[1]<<' '<<p[2]<<" 1 0 1\n";
    }
    cif.close();
    rhbm_gem::WriteMap(fs::path(output)/"map.mrc",map);
    Write(fs::path(output)/"generation-geometry.json",j::object{{"grid_size",j::value_from(map.GetGridSize())},
        {"spacing",j::value_from(map.GetGridSpacing())},{"origin",j::value_from(map.GetOrigin())},
        {"support_policy","sphere-fma-v1"},{"translation_angstrom",12},{"observation_encoding","float32-MRC"}});
    if(inputs_only) return;
    Fresh(map,model,fs::path(output)/"double","double",true);
    auto values=std::make_unique<double[]>(map.GetMapValueArraySize());
    for(std::size_t k=0;k<map.GetMapValueArraySize();++k) values[k]=static_cast<double>(static_cast<float>(map.GetMapValue(k)));
    rhbm_gem::MapObject quantized(map.GetGridSize(),map.GetGridSpacing(),map.GetOrigin(),std::move(values));
    Fresh(quantized,model,fs::path(output)/"float32","float32",true);
    Write(fs::path(output)/"completion.json",j::object{{"complete",true},{"datasets",j::array{"double","float32"}},
        {"generator","production-single-gaus"},{"translation_angstrom",12},{"atoms",24}});
}
}
