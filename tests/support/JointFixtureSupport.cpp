#include "support/JointFixtureSupport.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <fstream>
#include <sstream>
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
std::vector<std::vector<double>> Table(const fs::path & path)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing snapshot table.");
    std::string line; std::getline(in,line); std::vector<std::vector<double>> rows;
    while(std::getline(in,line)) {std::istringstream stream(line); std::string cell; auto & r=rows.emplace_back(); while(std::getline(stream,cell,',')) r.push_back(std::stod(cell));}
    return rows;
}
}
FrozenFixture LoadFixture(const fs::path & path)
{
    const auto snapshot=Read(path/"snapshot.json"),dataset=Read(path/"dataset.json");
    for(const auto & [file,key]:std::vector<std::pair<std::string,std::string>>{{"voxels.csv","voxels_sha256"},{"contributors.csv","contributors_sha256"}})
        if(sim::FileSha256(path/file)!=j::value_to<std::string>(snapshot.at(key))) throw std::runtime_error("Snapshot hash mismatch.");
    FrozenFixture in; in.hash=sim::FileSha256(path/"snapshot.json"); in.name=j::value_to<std::string>(dataset.at("name"));
    const auto voxels=Table(path/"voxels.csv"),contributors=Table(path/"contributors.csv");
    in.domain.rows=static_cast<Eigen::Index>(voxels.size()); std::vector<std::vector<Support>> support(j::value_to<std::size_t>(snapshot.at("atoms")));
    if(j::value_to<std::size_t>(snapshot.at("rows"))!=voxels.size() || j::value_to<std::size_t>(snapshot.at("memberships"))!=contributors.size()) throw std::runtime_error("Snapshot population mismatch.");
    std::vector<std::size_t> offsets(voxels.size()+1); std::pair<Eigen::Index,Eigen::Index> previous{-1,-1};
    for(const auto & r:contributors)
    {
        const auto row=static_cast<Eigen::Index>(r.at(0)),atom=static_cast<Eigen::Index>(r.at(1));
        if(row<0 || row>=in.domain.rows || atom<0 || static_cast<std::size_t>(atom)>=support.size() ||
            r[0]!=static_cast<double>(row) || r[1]!=static_cast<double>(atom) || std::pair{row,atom}<=previous)
            throw std::runtime_error("Invalid contributor CSR.");
        previous={row,atom}; ++offsets[static_cast<std::size_t>(row)+1];
        support[static_cast<std::size_t>(atom)].push_back({row,r.at(2)});
    }
    in.domain=Domain(in.domain.rows,std::move(support));
    for(std::size_t k=1;k<offsets.size();++k) offsets[k]+=offsets[k-1];
    if(snapshot.at("row_offsets").as_array().size()!=offsets.size()) throw std::runtime_error("Invalid CSR offsets.");
    for(std::size_t k=0;k<offsets.size();++k) if(j::value_to<std::size_t>(snapshot.at("row_offsets").at(k))!=offsets[k]) throw std::runtime_error("Invalid CSR offsets.");
    in.y64.resize(in.domain.rows); in.y32.resize(in.domain.rows);
    for(std::size_t k=0;k<voxels.size();++k) {in.y64(static_cast<Eigen::Index>(k))=voxels[k].at(7); in.y32(static_cast<Eigen::Index>(k))=voxels[k].at(8);}
    for(const auto & a:dataset.at("atoms").as_array()) in.ids.push_back(std::to_string(j::value_to<int>(a.at("serial_id"))));
    return in;
}
void RunFrozenFixture(const fs::path & path,const std::string & name,const fs::path & output)
{
    namespace core=rhbm_gem::core;
    Eigen::setNbThreads(1);
    const auto in=LoadFixture(path);
    const auto cases=Read(path/"cases.json");
    const auto initial=Parse(cases.at(name).at("initial_b"));
    const auto & y=name.ends_with("double") ? in.y64 : in.y32;
    const auto plan=RegisteredAudit(initial.size(),in.name,name);
    auto context=MakeContext(y,initial.size(),in.hash,&plan); context.atom_ids=in.ids;
    const auto partition=BuildPartition(in.domain,in.ids);
    if(partition.components.size()!=1) throw std::runtime_error("Frozen catalog case must be a single structural component.");
    const auto fit=FitComponent(partition.components[0],y,initial,context);
    j::object record;
    for(const char * key:{"search_success","usable_state","stop_reason","accepted_updates","profile_evaluations",
        "last_trusted_state","runtime_checks","runtime_convergence"})
        if(fit.contains(key)) record[key]=fit.at(key);
    j::array trials;
    for(const auto & t:fit.at("trials").as_array())
    {
        const auto & trial=t.as_object();
        const j::value trust=trial.if_contains("trust") ? trial.at("trust").at("passed") : j::value(nullptr);
        if(trial.at("accepted").as_bool() && (!trust.is_bool() || !trust.as_bool()))
            throw std::runtime_error("Guarded accepted an untrusted trial.");
        trials.push_back(j::object{{"accepted",trial.at("accepted")},{"valid",trial.at("valid")},{"trust_passed",trust}});
    }
    record["trials"]=trials;
    core::JointProblemInput input; input.atom_ids=in.ids; input.row_ids=context.row_ids;
    input.observations.assign(y.data(),y.data()+y.size()); input.support.resize(in.domain.atoms.size());
    for(std::size_t a=0;a<input.support.size();++a) for(const auto & point:in.domain.atoms[a])
        input.support[a].push_back({static_cast<std::size_t>(point.row),point.square});
    const core::JointProblem problem(std::move(input));
    const auto result=core::FitJointComponents(problem,{initial.data(),initial.data()+initial.size()});
    const bool available=fit.at("usable_state").as_bool();
    if(result.components.size()!=1 || result.components[0].state.has_value()!=available || result.search_completed!=fit.at("search_success").as_bool())
        throw std::runtime_error("Public runtime availability differs from frozen component search.");
    if(runtime_json::Status(result.components[0].RuntimeConvergence())!=fit.at("runtime_convergence").as_string() ||
        runtime_json::Status(result.RuntimeConvergence())!=fit.at("runtime_convergence").as_string())
        throw std::runtime_error("Public runtime convergence differs from the single-component evidence.");
    for(const auto & component:result.components) for(const auto & check:component.evidence)
        if(check.name=="two-step-derivative" && check.status!=core::JointCheckStatus::NotRun)
            throw std::runtime_error("Runtime executed an offline derivative audit.");
    if(result.regular_certificate!=core::JointCheckStatus::NotRun) throw std::runtime_error("Runtime promoted an offline certificate.");
    if(available)
    {
        const auto & state=*result.components[0].state;
        const auto & expected=fit.at("last_trusted_state");
        const auto beta=Parse(expected.at("beta")),eta=Parse(expected.at("eta"));
        for(Eigen::Index k=0;k<beta.size();++k)
            if(std::abs(state.ac[static_cast<std::size_t>(k)]-beta(k))>1e-10*(1+std::abs(beta(k))))
                throw std::runtime_error("Public coefficients differ from component search.");
        for(Eigen::Index k=0;k<eta.size();++k)
            if(std::abs(state.log_b[static_cast<std::size_t>(k)]-eta(k))>1e-10)
                throw std::runtime_error("Public widths differ from component search.");
        if(!result.prediction || !result.objective || !result.assembled_state) throw std::runtime_error("Missing assembled runtime state.");
        const Vector prediction=Eigen::Map<const Vector>(result.prediction->data(),y.size());
        const double objective=.5*(prediction-y).squaredNorm()/(context.scale*context.scale);
        for(double value:{state.objective,result.assembled_state->objective,*result.objective})
            if(std::abs(value-objective)>1e-15*(1+objective)) throw std::runtime_error("Public objective normalization differs.");
    }
    else if(result.prediction || result.objective) throw std::runtime_error("Failed component fabricated a complete prediction.");
    fs::create_directories(output.parent_path());
    Write(output,j::object{{"dataset",in.name},{"case",name},{"record",record},{"census",Census(in.domain,partition,context)},{"api_contract_passed",true}});
}
}
