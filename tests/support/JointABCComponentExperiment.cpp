#include "support/JointABCComponentExperiment.hpp"
#include "support/JointABCComponents.hpp"
#include "support/JointABCCertification.hpp"
#include "support/JointABCLocalCertification.hpp"
#include "support/JointABCPrecision.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <chrono>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
using Vector=Eigen::VectorXd;
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing component input: "+path.string());
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},options);
}
void Write(const fs::path & path,const j::value & value)
{std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';}
Vector Parse(const j::value & v)
{Vector out(static_cast<Eigen::Index>(v.as_array().size())); for(Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k))); return out;}
j::array Values(const Vector & v) {j::array out; for(double x:v) out.push_back(x); return out;}
std::vector<std::vector<double>> Table(const fs::path & path)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing snapshot table.");
    std::string line; std::getline(in,line); std::vector<std::vector<double>> rows;
    while(std::getline(in,line)) {std::istringstream stream(line); std::string cell; auto & r=rows.emplace_back(); while(std::getline(stream,cell,',')) r.push_back(std::stod(cell));}
    return rows;
}
struct Input
{
    Domain domain{0,{}};
    Vector y64,y32;
    std::string hash,name;
    std::vector<std::string> ids;
    std::vector<fs::path> sources;
};
Input Load(const fs::path & path)
{
    if(fs::exists(path/"composition.json"))
    {
        const auto recipe=Read(path/"composition.json"); Input out; out.name=j::value_to<std::string>(recipe.at("name"));
        out.hash=sim::FileSha256(path/"composition.json");
        for(const auto & entry:recipe.at("sources").as_array())
        {
            const auto source=path/j::value_to<std::string>(entry.at("path")); const auto in=Load(source);
            if(in.hash!=j::value_to<std::string>(entry.at("snapshot_sha256"))) throw std::runtime_error("Composite source hash mismatch.");
            const auto offset=out.domain.rows; out.domain.rows+=in.domain.rows;
            out.y64.conservativeResize(out.domain.rows); out.y32.conservativeResize(out.domain.rows);
            out.y64.tail(in.domain.rows)=in.y64; out.y32.tail(in.domain.rows)=in.y32;
            const auto prefix=j::value_to<std::string>(entry.at("identity_prefix"));
            for(std::size_t a=0;a<in.domain.atoms.size();++a)
            {
                out.domain.atoms.emplace_back(); for(const auto & s:in.domain.atoms[a]) out.domain.atoms.back().push_back({s.row+offset,s.square});
                out.ids.push_back(prefix+in.ids[a]);
            }
            out.sources.push_back(source);
        }
        return out;
    }
    const auto snapshot=Read(path/"snapshot.json"),dataset=Read(path/"dataset.json");
    for(const auto & [file,key]:std::vector<std::pair<std::string,std::string>>{{"voxels.csv","voxels_sha256"},{"contributors.csv","contributors_sha256"}})
        if(sim::FileSha256(path/file)!=j::value_to<std::string>(snapshot.at(key))) throw std::runtime_error("Snapshot hash mismatch.");
    Input in; in.hash=sim::FileSha256(path/"snapshot.json"); in.name=j::value_to<std::string>(dataset.at("name"));
    const auto voxels=Table(path/"voxels.csv"),contributors=Table(path/"contributors.csv");
    in.domain.rows=static_cast<Eigen::Index>(voxels.size()); in.domain.atoms.resize(j::value_to<std::size_t>(snapshot.at("atoms")));
    if(j::value_to<std::size_t>(snapshot.at("rows"))!=voxels.size() || j::value_to<std::size_t>(snapshot.at("memberships"))!=contributors.size()) throw std::runtime_error("Snapshot population mismatch.");
    std::vector<std::size_t> offsets(voxels.size()+1); std::pair<Eigen::Index,Eigen::Index> previous{-1,-1};
    for(const auto & r:contributors)
    {
        const auto row=static_cast<Eigen::Index>(r.at(0)),atom=static_cast<Eigen::Index>(r.at(1));
        if(row<0 || row>=in.domain.rows || atom<0 || static_cast<std::size_t>(atom)>=in.domain.atoms.size() ||
            r[0]!=static_cast<double>(row) || r[1]!=static_cast<double>(atom) || std::pair{row,atom}<=previous)
            throw std::runtime_error("Invalid contributor CSR.");
        previous={row,atom}; ++offsets[static_cast<std::size_t>(row)+1];
        in.domain.atoms[static_cast<std::size_t>(atom)].push_back({row,r.at(2)});
    }
    for(std::size_t k=1;k<offsets.size();++k) offsets[k]+=offsets[k-1];
    if(snapshot.at("row_offsets").as_array().size()!=offsets.size()) throw std::runtime_error("Invalid CSR offsets.");
    for(std::size_t k=0;k<offsets.size();++k) if(j::value_to<std::size_t>(snapshot.at("row_offsets").at(k))!=offsets[k]) throw std::runtime_error("Invalid CSR offsets.");
    in.y64.resize(in.domain.rows); in.y32.resize(in.domain.rows);
    for(std::size_t k=0;k<voxels.size();++k) {in.y64(static_cast<Eigen::Index>(k))=voxels[k].at(7); in.y32(static_cast<Eigen::Index>(k))=voxels[k].at(8);}
    for(const auto & a:dataset.at("atoms").as_array()) in.ids.push_back(std::to_string(j::value_to<int>(a.at("serial_id"))));
    in.sources.push_back(path); return in;
}
struct FrozenState {Vector eta,beta,initial_b;};
FrozenState Frozen(const Input & in,const std::string & name,bool endpoint)
{
    FrozenState out; const auto count=static_cast<Eigen::Index>(in.domain.atoms.size());
    out.eta.resize(count); out.beta.resize(2*count); out.initial_b.resize(count); Eigen::Index offset{};
    for(const auto & source:in.sources)
    {
        const auto fit=Read(source/"guarded"/"fits"/(name+".json")); const auto & state=fit.at(endpoint ? "primary" : "initial");
        if(j::value_to<std::string>(fit.at("observation_snapshot_sha256"))!=sim::FileSha256(source/"snapshot.json"))
            throw std::runtime_error("Frozen start belongs to a different observation snapshot.");
        const auto initial=Parse(fit.at("initial_b")); const auto n=initial.size();
        out.eta.segment(offset,n)=endpoint ? Parse(state.at("eta")) : Vector(initial.array().log());
        out.beta.segment(2*offset,2*n)=Parse(state.at("beta")); out.initial_b.segment(offset,n)=initial; offset+=n;
    }
    return out;
}
EvaluationContext Context(const Input & in,const Vector & y,const std::string & name)
{
    AuditPlan plan; Eigen::Index offset{};
    for(const auto & source:in.sources)
    {
        const auto dataset=Read(source/"dataset.json"); const auto n=static_cast<Eigen::Index>(dataset.at("atoms").as_array().size());
        const auto local=RegisteredAudit(n,j::value_to<std::string>(dataset.at("name")),name);
        plan.trial_details|=local.trial_details; plan.expanded_if_unverified|=local.expanded_if_unverified;
        plan.precision|=local.precision; plan.boundary|=local.boundary;
        for(auto a:local.boundary_atoms) plan.boundary_atoms.push_back(a+offset); offset+=n;
    }
    plan.block_precision=in.sources.size()>1; plan.cache_precision=true;
    auto c=MakeContext(y,static_cast<Eigen::Index>(in.domain.atoms.size()),in.hash,&plan); c.atom_ids=in.ids; return c;
}
EvaluationContext SearchContext(const Input & in,const Vector & y,const std::string & name,const Vector & initial_b)
{
    auto context=Context(in,y,name);
    if(in.sources.size()==1) return context;
    const auto reference=Assess(in.domain,y,Vector(initial_b.array().log()),context);
    auto & directions=context.audit.directions; directions.resize(initial_b.size(),3);
    directions.col(0)=Vector::Ones(initial_b.size()).normalized();
    for(Eigen::Index k=0;k<initial_b.size();++k) directions(k,1)=k%2 ? -1 : 1;
    directions.col(1).normalize(); directions.col(2).setZero(); directions(0,2)=1;
    if(reference.contains("width_spectrum")) directions.col(2)=Parse(reference.at("width_spectrum").at("weak_directions").at(0));
    return context;
}
EvaluationContext RestoreContext(const Input & in,const Vector & y,const std::string & name,const j::value & saved)
{
    auto context=Context(in,y,name); const auto & directions=saved.at("audit").at("directions").as_array();
    if(!directions.empty())
    {
        context.audit.directions.resize(static_cast<Eigen::Index>(in.domain.atoms.size()),static_cast<Eigen::Index>(directions.size()));
        for(std::size_t k=0;k<directions.size();++k) context.audit.directions.col(static_cast<Eigen::Index>(k))=Parse(directions[k]);
    }
    if(ContextEvidence(context)!=saved) throw std::runtime_error("Frozen context does not match input and registered policy.");
    return context;
}
void Label(j::object & fit,const Input & in,const std::string & name,const Vector & initial)
{
    fit["dataset"]=in.name; fit["case"]=name; fit["initial_b"]=Values(initial);
    fit["observation_snapshot_sha256"]=in.hash;
}
EvaluationContext EndpointContext(EvaluationContext context,const j::object & fit)
{
    const auto n=static_cast<Eigen::Index>(context.atom_ids.size());
    auto & directions=context.audit.directions; directions.resize(n,3);
    directions.col(0)=Vector::Ones(n).normalized();
    for(Eigen::Index k=0;k<n;++k) directions(k,1)=k%2 ? -1 : 1;
    directions.col(1).normalize(); directions.col(2).setZero(); directions(0,2)=1;
    if(fit.contains("width_spectrum")) directions.col(2)=Parse(fit.at("width_spectrum").at("weak_directions").at(0));
    return context;
}
void AuditFit(const Domain & domain,const Vector & y,const j::object & search_fit,const EvaluationContext & context,
    const fs::path & output,bool refresh)
{
    fs::create_directories(output); auto fit=search_fit;
    const auto start=std::chrono::steady_clock::now();
    // Search diagnostics retain their frozen initial directions. Certification
    // must freshly test the endpoint directions, including its weakest mode.
    const bool refresh_endpoint=refresh && fit.contains("primary") && fit.at("primary").at("valid").as_bool();
    if(refresh_endpoint)
    {
        const auto eta=Parse(fit.at("primary").at("eta")),beta=Parse(fit.at("primary").at("beta"));
        const bool supplied=fit.if_contains("assembled_state_preserved") && fit.at("assembled_state_preserved").as_bool();
        auto assessment=Assess(domain,y,eta,context,supplied ? &beta : nullptr);
        for(auto & field:assessment) fit[field.key()]=std::move(field.value());
        if(supplied) fit["joint_qualified"]=fit.at("joint_qualified").as_bool() &&
            fit.at("assembled_profile_agrees").as_bool() && !fit.at("search_stopped_without_convergence").as_bool();
    }
    fit["audit_assessment_refreshed"]=refresh_endpoint;
    fit["audit_assessment_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    Write(output/"fit.json",fit);
    certification::Audit(domain,y,fit,output/"audit.json",&context);
    auto audit=Read(output/"audit.json").as_object();
    const bool verified=fit.if_contains("derivative_verified") && fit.at("derivative_verified").as_bool();
    const bool requested=context.audit.expanded_if_unverified && (!verified || context.audit.precision);
    audit["scope_evidence"]=j::object{{"legacy_derivative_verified",verified},{"expanded_requested",requested},
        {"expanded_executed",audit.contains("ladders") && !audit.at("ladders").as_array().empty()},
        {"trigger",!requested ? "legacy-evidence-or-unregistered-expansion" : context.audit.precision ? "registered-high-precision-target" : "unverified-legacy-derivative"},
        {"precision_requested",context.audit.precision},{"boundary_requested",context.audit.boundary}};
    if(context.audit.precision && context.observations) audit["precision_normalization"]=certification::PrecisionNormalization(*context.observations);
    Write(output/"audit.json",audit);
    Write(output/"context.json",ContextEvidence(context));
    j::array direction_diagnostics;
    if(fit.contains("primary") && fit.at("primary").at("valid").as_bool())
    {
        const auto & state=fit.at("primary");
        const auto e=AtState(domain,y,Parse(state.at("eta")),Parse(state.at("beta")),context);
        const auto d=Differentiate(e,context.scale,&context);
        for(Eigen::Index k=0;k<context.audit.directions.cols();++k)
        {
            const double norm=context.audit.directions.col(k).norm();
            direction_diagnostics.push_back(j::object{{"direction",k},{"norm",norm},{"zero_direction",norm==0},
                {"jacobian_direction_norm",d.valid ? j::value((d.jacobian*context.audit.directions.col(k)).norm()) : j::value(nullptr)}});
        }
    }
    Write(output/"direction-diagnostics.json",direction_diagnostics);
    if(context.audit.boundary && fit.contains("primary") && fit.at("primary").at("valid").as_bool())
        Write(output/"boundary.json",certification::BoundaryAudit(domain,y,Parse(fit.at("primary").at("eta")),Parse(fit.at("primary").at("beta")),&context));
}
}
void ComponentSameState(const std::string & dataset_path,const std::string & output_path)
{
    Eigen::setNbThreads(1); const fs::path output(output_path); const auto in=Load(dataset_path);
    if(fs::exists(output)) throw std::runtime_error("Use a fresh same-state directory.");
    fs::create_directories(output); const auto partition=BuildPartition(in.domain,in.ids);
    Write(output/"census.json",Census(in.domain,partition,Context(in,in.y64,"first-stage-double")));
    bool passed=true; j::array records;
    for(const std::string precision:{"double","float32"}) for(const std::string start:{"first-stage","narrower","wider","mixed"})
    {
        const auto name=start+"-"+precision; const auto & y=precision=="double" ? in.y64 : in.y32; const auto context=Context(in,y,name);
        for(bool endpoint:{false,true})
        {
            const auto state=Frozen(in,name,endpoint); const auto id=name+(endpoint ? "-endpoint" : "-initial");
            std::cout<<in.name<<'/'<<id<<std::endl;
            auto parity=SameState(in.domain,y,state.eta,state.beta,partition,context);
            parity["dataset"]=in.name; parity["case"]=name; parity["state"]=endpoint ? "frozen-guarded-endpoint" : "frozen-initial";
            Write(output/(id+".json"),parity); passed &= parity.at("passed").as_bool();
            records.push_back(j::object{{"id",id},{"passed",parity.at("passed")},{"full_equivalence",parity.at("full_equivalence")}});
        }
    }
    Write(output/"completion.json",j::object{{"complete",true},{"passed",passed},{"states",records},{"snapshot_sha256",in.hash}});
}
void ComponentRegression(const std::string & dataset_path,const std::string & output_path)
{
    Eigen::setNbThreads(1); const fs::path source(dataset_path),output(output_path);
    if(fs::exists(output)) throw std::runtime_error("Use a fresh component regression directory.");
    const auto in=Load(source); fs::create_directories(output);
    auto context=MakeContext(in.y64,static_cast<Eigen::Index>(in.domain.atoms.size()),in.hash); context.atom_ids=in.ids;
    Write(output/"census.json",Census(in.domain,BuildPartition(in.domain,in.ids),context));
    std::map<std::string,std::pair<j::object,std::string>> audits;
    for(const std::string variant:{"legacy","guarded","guarded-log"})
    {
        const auto target=output/variant; fs::create_directories(target/"fits"); fs::create_directories(target/"audits");
        std::vector<fs::path> files; for(const auto & f:fs::directory_iterator(source/variant/"fits")) files.push_back(f.path()); std::sort(files.begin(),files.end());
        for(const auto & file:files)
        {
            const auto frozen=Read(file).as_object(); const auto name=file.stem().string(); const auto & y=name.ends_with("double") ? in.y64 : in.y32;
            if(j::value_to<std::string>(frozen.at("observation_snapshot_sha256"))!=in.hash) throw std::runtime_error("Changed fit snapshot.");
            std::cout<<in.name<<'/'<<variant<<'/'<<name<<std::endl;
            j::object resources; auto fit=Fit(in.domain,y,Parse(frozen.at("initial_b")),&resources,variant);
            for(const auto & field:frozen) if(!fit.contains(field.key())) fit[field.key()]=field.value();
            fit["resources"]=resources; Write(target/"fits"/file.filename(),fit);
            // Match the established deterministic cache; a fresh invocation recomputes it.
            const std::string key=j::serialize(j::array{name.ends_with("double"),fit.at("primary"),fit.at("joint_qualified"),fit.if_contains("width_spectrum") ? fit.at("width_spectrum") : j::value(nullptr)});
            const auto audit_path=target/"audits"/file.filename(); const auto audit_source=variant+"/audits/"+file.filename().string();
            if(audits.contains(key))
            {
                auto audit=audits.at(key).first; audit["case"]=fit.at("case"); audit["reused_identical_endpoint"]=audits.at(key).second;
                audit["seconds"]=0; if(audit.contains("precision")) audit.at("precision").as_object()["seconds"]=0; Write(audit_path,audit);
            }
            else {certification::Audit(in.domain,y,fit,audit_path); audits.emplace(key,std::make_pair(Read(audit_path).as_object(),audit_source));}
        }
    }
    Write(output/"completion.json",j::object{{"complete",true},{"cases",24},{"snapshot_sha256",in.hash}});
}
void ComponentRun(const std::string & dataset_path,const std::string & output_path)
{
    Eigen::setNbThreads(1); const auto in=Load(dataset_path); const fs::path output(output_path);
    if(fs::exists(output)) throw std::runtime_error("Use a fresh component experiment directory.");
    fs::create_directories(output); ComponentSameState(dataset_path,(output/"same-state").string());
    if(!Read(output/"same-state/completion.json").at("passed").as_bool()) throw std::runtime_error("Same-state gate failed before independent search.");
    const auto partition=BuildPartition(in.domain,in.ids); j::array cases;
    Write(output/"census.json",Census(in.domain,partition,Context(in,in.y64,"first-stage-double")));
    for(const std::string precision:{"double","float32"}) for(const std::string start:{"first-stage","narrower","wider","mixed"})
    {
        const auto name=start+"-"+precision; std::cout<<in.name<<'/'<<name<<" search"<<std::endl;
        const auto & y=precision=="double" ? in.y64 : in.y32; const auto frozen=Frozen(in,name,false);
        const auto setup=std::chrono::steady_clock::now(); const auto context=SearchContext(in,y,name,frozen.initial_b);
        const auto target=output/"cases"/name; fs::create_directories(target/"components");
        Write(target/"context.json",ContextEvidence(context));
        Write(target/"setup.json",j::object{{"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-setup).count()},
            {"direction_reference",in.sources.size()>1 ? "frozen-initial-global-profile; unit-axis fallback if unavailable" : "historical-endpoint-directions"}});
        j::object resources; auto monolithic=Fit(in.domain,y,frozen.initial_b,&resources,"guarded",&context);
        Label(monolithic,in,name,frozen.initial_b); monolithic["resources"]=resources;
        Write(target/"monolithic-fit.json",monolithic);
        auto assembled=FitComponents(in.domain,y,frozen.initial_b,partition,context);
        auto components=std::move(assembled.at("components").as_array()); j::array ids;
        for(std::size_t k=0;k<components.size();++k)
        {
            auto fit=components[k].as_object(); Label(fit,in,name,Select(frozen.initial_b,partition.components[k].atoms));
            Write(target/"components"/(std::to_string(k)+"-fit.json"),fit); ids.emplace_back(partition.components[k].id);
        }
        assembled["components"]=ids; Label(assembled,in,name,frozen.initial_b); Write(target/"assembled-fit.json",assembled);
        // Both endpoints are replayed through the same-state checker, whether
        // a search succeeds, stalls, or converges to a different stationary point.
        for(const std::string variant:{"monolithic","assembled"})
        {
            const auto & fit=variant=="monolithic" ? monolithic : assembled;
            const char * state=variant=="assembled" ? "assembled_state" : "primary";
            if(!fit.contains(state) || !fit.at(state).at("valid").as_bool())
            {Write(target/(variant+"-same-state.json"),j::object{{"available",false},{"reason","no-valid-endpoint"}}); continue;}
            Write(target/(variant+"-same-state.json"),SameState(in.domain,y,Parse(fit.at(state).at("eta")),Parse(fit.at(state).at("beta")),partition,context));
        }
        cases.emplace_back(name);
    }
    Write(output/"completion.json",j::object{{"complete",true},{"cases",cases},{"snapshot_sha256",in.hash}});
}
namespace {
void LocalAuditFit(const ComponentView & view,const Vector & parent_y,const j::object & fit,
    const EvaluationContext & parent,const fs::path & output)
{
    const auto y=Select(parent_y,view.rows);
    auto local=certification::PrepareLocalAudit(view.domain,y,fit,ComponentContext(parent,view,true));
    certification::ResetPrecisionCache();
    AuditFit(view.domain,y,local.fit,local.context,output,false);
    Write(output/"scope.json",local.scope);
}
}
void ComponentAudit(const std::string & dataset_path,const std::string & run_path,const std::string & output_path,
    const std::string & only_case,bool local_only)
{
    Eigen::setNbThreads(1); const auto in=Load(dataset_path); const fs::path run(run_path),output(output_path);
    if(fs::exists(only_case.empty() ? output : output/only_case)) throw std::runtime_error("Use a fresh component audit directory.");
    const auto partition=BuildPartition(in.domain,in.ids); const auto completed=Read(run/"completion.json");
    j::array cases;
    for(const auto & entry:completed.at("cases").as_array())
        if(only_case.empty() || entry.as_string()==only_case) cases.push_back(entry);
    if(cases.empty()) throw std::invalid_argument("Unknown frozen audit case.");
    for(const auto & entry:cases)
    {
        const auto name=j::value_to<std::string>(entry); const auto & y=name.ends_with("double") ? in.y64 : in.y32;
        const auto source=run/"cases"/name,target=output/name; const auto context=RestoreContext(in,y,name,Read(source/"context.json"));
        certification::ResetPrecisionCache();
        std::cout<<in.name<<'/'<<name<<" audit"<<std::endl;
        if(!local_only)
        {
        for(const std::string variant:{"monolithic","assembled"})
        {
            const auto fit=Read(source/(variant+"-fit.json")).as_object();
            AuditFit(in.domain,y,fit,EndpointContext(context,fit),target/variant,in.sources.size()>1);
        }
        const auto assembled=Read(source/"assembled-fit.json").as_object();
        const auto global=EndpointContext(context,assembled);
        Write(target/"direction-scope.json",j::object{{"policy","endpoint-global-directions; component restrictions without renormalization"},
            {"assembled_weak_direction_available",assembled.contains("width_spectrum")},
            {"unavailable_fallback","unit-axis; does not provide global weak-direction evidence"}});
        for(std::size_t k=0;k<partition.components.size();++k)
        {
            const auto & view=partition.components[k]; const auto child=ComponentContext(global,view,true);
            AuditFit(view.domain,Select(y,view.rows),Read(source/"components"/(std::to_string(k)+"-fit.json")).as_object(),child,target/"components"/std::to_string(k),in.sources.size()>1);
        }
        Write(target/"reference-costs.json",certification::PrecisionCacheCosts());
        }
        for(std::size_t k=0;k<partition.components.size();++k)
            LocalAuditFit(partition.components[k],y,Read(source/"components"/(std::to_string(k)+"-fit.json")).as_object(),
                context,target/"local-components"/std::to_string(k));
    }
    Write((only_case.empty() ? output : output/only_case)/"completion.json",j::object{{"complete",true},{"cases",cases},{"cache_policy","exact-input precision references within one case; cleared before every case; no cross-run cache"}});
}
void ComponentRerun(const std::string & dataset_path,const std::string & name,const std::string & id,
    const std::string & context_path,const std::string & output_path,bool local_only)
{
    Eigen::setNbThreads(1); const auto in=Load(dataset_path); const fs::path output(output_path);
    if(fs::exists(output)) throw std::runtime_error("Use a fresh isolated component directory.");
    const auto & y=name.ends_with("double") ? in.y64 : in.y32;
    const auto bundle=Read(context_path);
    const auto context=RestoreContext(in,y,name,bundle.at("search_context"));
    const auto partition=BuildPartition(in.domain,in.ids);
    const auto found=std::find_if(partition.components.begin(),partition.components.end(),[&](const auto & c){return c.id==id;});
    if(found==partition.components.end()) throw std::invalid_argument("Unknown stable component identity.");
    const auto initial=Frozen(in,name,false).initial_b; auto fit=FitComponent(*found,y,initial,context);
    Label(fit,in,name,Select(initial,found->atoms)); fs::create_directories(output); Write(output/"fit.json",fit);
    certification::ResetPrecisionCache();
    if(local_only) LocalAuditFit(*found,y,fit,context,output/"audit");
    else
    {
        const auto audit_context=RestoreContext(in,y,name,bundle.at("audit_context"));
        AuditFit(found->domain,Select(y,found->rows),fit,ComponentContext(audit_context,*found,true),output/"audit",in.sources.size()>1);
    }
    Write(output/"completion.json",j::object{{"complete",true},{"component_id",id},{"snapshot_sha256",in.hash}});
}
} // namespace second_stage_test::matched::joint_abc
