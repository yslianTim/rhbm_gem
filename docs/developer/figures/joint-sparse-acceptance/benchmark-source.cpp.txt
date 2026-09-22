#include "core/detail/joint_component/Problem.hpp"
#ifndef SPARSE_BASELINE_DRIVER
#include "core/detail/joint_component/SparseFactor.hpp"
#endif
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointFixtureSupport.hpp"
#include "support/JointRuntimeJson.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <fstream>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
namespace p=second_stage_test::matched::joint_abc;
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
j::value Read(const char * path) {std::ifstream f(path); if(!f) throw std::runtime_error("Missing input"); return j::parse(std::string(std::istreambuf_iterator<char>(f),{}));}
void Write(const char * path,const j::object & v) {std::ofstream f(path);f<<j::serialize(v)<<'\n';}
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
j::array Values(const n::Vector & v) {j::array out;for(auto x:v) out.push_back(std::isfinite(x) ? j::value(x) : j::value(nullptr));return out;}
void Snapshot(const char * output,j::object & report)
{
#ifndef SPARSE_BASELINE_DRIVER
    const auto & w=n::SparseWorkForTesting();
    report["work"]=j::object{{"symbolic",w.symbolic},{"numeric",w.numeric},{"symbolic_reuses",w.symbolic_reuses},
        {"factor_reuses",w.factor_reuses},{"cancellation_reductions",w.cancellation_reductions},{"factor_nonzeros_upper_bound",w.factor_nonzeros},
        {"basis_and_csc_preparation_seconds",w.matrix_preparation_seconds},{"symbolic_seconds",w.symbolic_seconds},{"numeric_seconds",w.numeric_seconds},
        {"reference_qr_seconds",w.reference_seconds},{"reference_svd_seconds",w.reference_svd_seconds}};
#endif
    Write(output,report);
}
void Run(const n::Domain & domain,n::VectorRef y,const n::Vector & b,const n::EvaluationContext & context,const char * output)
{
    j::object report{{"rows",domain.rows},{"atoms",b.size()},{"initial_b",Values(b)},{"stage","primary"}};Snapshot(output,report);
#ifndef SPARSE_BASELINE_DRIVER
    n::SparseWorkForTesting()={}; n::LinearWorkspace workspace;
#endif
    const n::Vector eta=b.array().log();auto started=Clock::now();
#ifdef SPARSE_BASELINE_DRIVER
    auto primary=n::EvaluateProfile(domain,y,eta,false,&context);
#else
    auto primary=n::EvaluateProfile(domain,y,eta,false,&context,nullptr,&workspace);
#endif
    report["primary_seconds"]=Seconds(started);report["primary"]=second_stage_test::matched::runtime_json::Endpoint(primary);
    report["stage"]="reference";Snapshot(output,report);started=Clock::now();
    auto reference=n::EvaluateProfile(domain,y,eta,true,&context);
    report["reference_seconds"]=Seconds(started);report["reference"]=second_stage_test::matched::runtime_json::Endpoint(reference);
    report["trust"]=second_stage_test::matched::runtime_json::Trust(n::CheckTrust(domain,y,primary,context,reference));
    report["stage"]="derivative";Snapshot(output,report);started=Clock::now();
    const auto derivative=n::PrepareDerivative(primary,context.scale,&context);
    report["derivative_seconds"]=Seconds(started);report["derivative_valid"]=derivative.valid;
    report["stage"]="repeat-primary";Snapshot(output,report);started=Clock::now();
#ifdef SPARSE_BASELINE_DRIVER
    n::EvaluateProfile(domain,y,eta,false,&context);
#else
    n::EvaluateProfile(domain,y,eta,false,&context,nullptr,&workspace);

#endif
    report["repeat_primary_seconds"]=Seconds(started);rusage usage{};getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    report["peak_rss_bytes"]=usage.ru_maxrss;
#else
    report["peak_rss_bytes"]=usage.ru_maxrss*1024;
#endif
    report["stage"]="complete";Snapshot(output,report);
}
}
int main(int argc,char ** argv)
{
    try {
        Eigen::setNbThreads(1); const std::string mode=argc>1 ? argv[1] : "";
        if(mode=="fixture" && argc==5)
        {
            auto fixture=p::LoadFixture(argv[2]); const auto cases=Read((std::filesystem::path(argv[2])/"cases.json").c_str());
            const auto widths=j::value_to<std::vector<double>>(cases.at(argv[3]).at("initial_b"));
            const n::Vector & y=std::string(argv[3]).ends_with("double") ? fixture.y64 : fixture.y32;
            const auto context=n::CreateContext(y,static_cast<Eigen::Index>(widths.size()));
            Run(fixture.domain,y,Eigen::Map<const n::Vector>(widths.data(),static_cast<Eigen::Index>(widths.size())),context,argv[4]);return 0;
        }
        if((mode!="prepare" || argc!=5) && (mode!="initial" || argc!=6)) throw std::invalid_argument("prepare MODEL MAP WIDTHS | initial MODEL MAP WIDTHS OUTPUT | fixture DATASET CASE OUTPUT");
        auto model=rhbm_gem::ReadModel(argv[2]);auto map=rhbm_gem::ReadMap(argv[3]);
        model->SelectAllAtoms();model->ApplySymmetrySelection(false);model->ApplyElementSelection(Element::HYDROGEN,true);
        const auto problem=c::BuildJointProblem(*map,*model);const auto & data=c::JointProblemAccess::Get(problem);
        if(mode=="prepare")
        {
#ifndef SPARSE_BASELINE_DRIVER
            const auto init=n::InitializeContributors(*map,*model,problem);
            Write(argv[4],{{"b",j::value_from(init.b)},{"atom_ids",j::value_from(problem.Input().atom_ids)}});return 0;
#else
            throw std::invalid_argument("Prepare inputs with the candidate's unchanged production initializer");
#endif
        }
        const auto input=Read(argv[4]);const auto widths=j::value_to<std::vector<double>>(input.at("b"));
        if(j::value_to<std::vector<std::string>>(input.at("atom_ids"))!=problem.Input().atom_ids) throw std::runtime_error("Widths identity mismatch");
        const auto & component=*std::max_element(data.partition.components.begin(),data.partition.components.end(),[](const auto & a,const auto & b){return a.atoms.size()<b.atoms.size();});
        const n::Vector all=Eigen::Map<const n::Vector>(widths.data(),static_cast<Eigen::Index>(widths.size()));
        const auto y=n::SelectValues(data.y,component.rows);const auto b=n::SelectValues(all,component.atoms);
        const auto context=n::ChildContext(data.context,component,true);
        Run(component.domain,y,b,context,argv[5]);return 0;
    } catch(const std::exception & e) {std::cerr<<e.what()<<'\n';return 1;}
}
