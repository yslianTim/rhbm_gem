#include "core/detail/joint_component/Problem.hpp"
#ifndef SPARSE_BASELINE_DRIVER
#include "core/detail/joint_component/SparseFactor.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
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
#include <bit>
#include <sys/resource.h>
namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
namespace p=second_stage_test::matched::joint_abc;
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
bool audit{};
std::filesystem::path capture;
j::array svd_records;
j::value Read(const char * path,bool precise=false)
{
    std::ifstream f(path); if(!f) throw std::runtime_error("Missing input");
    j::parse_options options; if(precise) options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(f),{}),{},options);
}
void Write(const char * path,const j::object & v) {std::ofstream f(path);f<<j::serialize(v)<<'\n';}
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
j::array Values(const n::Vector & v) {j::array out;for(auto x:v) out.push_back(std::isfinite(x) ? j::value(x) : j::value(nullptr));return out;}
#ifndef SPARSE_BASELINE_DRIVER
void Mode(const std::string & name)
{
    if(name=="legacy") n::CompactSvdModeForTesting()=n::CompactSvdMode::Legacy;
    else if(name=="values") n::CompactSvdModeForTesting()=n::CompactSvdMode::ValuesOnly;
    else if(name=="auto") n::CompactSvdModeForTesting()=n::CompactSvdMode::Automatic;
    else throw std::invalid_argument("Unknown compact SVD mode");
}
void MatrixFile(const std::filesystem::path & path,const n::Matrix & a)
{
    std::ofstream f(path,std::ios::binary);
    f.write(reinterpret_cast<const char *>(a.data()),static_cast<std::streamsize>(static_cast<std::size_t>(a.size())*sizeof(double)));
    if(!f) throw std::runtime_error("Cannot write compact matrix");
}
j::object SvdRecord(const n::CompactSvdResult & r)
{
    return {{"valid",r.valid},{"rank",r.rank},{"threshold",r.threshold},{"singular_values",Values(r.singular_values)},
        {"solution",Values(r.solution)},{"used_bdc",r.used_bdc},{"jacobi_retry",r.jacobi_retry}};
}
void Capture(const n::Matrix & a,double relative,double absolute,const n::Vector * rhs,const n::CompactSvdResult & r)
{
    auto record=SvdRecord(r); const auto index=svd_records.size();
    record["role"]=rhs ? "reference" : "free-design"; record["rows"]=a.rows(); record["columns"]=a.cols();
    record["relative_threshold"]=relative; record["absolute_override"]=absolute;
    if(!capture.empty())
    {
        const auto stem=std::string(rhs ? "reference-" : "free-design-")+std::to_string(index);
        record["matrix_file"]=stem+".bin"; record["rhs"]=rhs ? j::value(Values(*rhs)) : j::value(nullptr);
        record["byte_order"]=std::endian::native==std::endian::little ? "little" : "big";
        MatrixFile(capture/(stem+".bin"),a); Write((capture/(stem+".json")).c_str(),record);
    }
    svd_records.push_back(std::move(record));
}
#endif
void Snapshot(const char * output,j::object & report)
{
#ifndef SPARSE_BASELINE_DRIVER
    const auto & w=n::SparseWorkForTesting();
    report["work"]=j::object{{"symbolic",w.symbolic},{"numeric",w.numeric},{"symbolic_reuses",w.symbolic_reuses},
        {"factor_reuses",w.factor_reuses},{"cancellation_reductions",w.cancellation_reductions},{"factor_nonzeros_upper_bound",w.factor_nonzeros},
        {"basis_and_csc_preparation_seconds",w.matrix_preparation_seconds},{"symbolic_seconds",w.symbolic_seconds},{"numeric_seconds",w.numeric_seconds},
        {"reference_qr_seconds",w.reference_seconds},{"reference_svd_seconds",w.reference_svd_seconds},
        {"derivative_preparations",w.derivative_preparations},{"derivative_compacts",w.derivative_compacts},
        {"reference_compacts",w.reference_compacts},{"free_design_svds",w.free_design_svds},{"reference_svds",w.reference_svds},
        {"reference_solves",w.reference_solves},{"bdc_svds",w.bdc_svds},{"jacobi_retries",w.jacobi_retries},
        {"derivative_compact_seconds",w.derivative_compact_seconds},{"reference_compact_seconds",w.reference_compact_seconds},
        {"free_design_svd_seconds",w.free_design_svd_seconds},{"reference_solve_seconds",w.reference_solve_seconds},
        {"cancellation_seconds",w.cancellation_seconds},{"jacobi_retry_seconds",w.jacobi_retry_seconds},
        {"derivative_inclusive_seconds",w.derivative_seconds}};
    if(audit) report["svd_records"]=svd_records;
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
    report["derivative_reason"]=derivative.reason;
#ifndef SPARSE_BASELINE_DRIVER
    if(audit && derivative.valid)
    {
        report["stage"]="derivative-audit"; Snapshot(output,report);
        const auto reduced=n::ReduceDerivative(derivative,primary.residual,true);
        report["derivative_audit"]=j::object{{"valid",reduced.valid},
            {"coefficients",Values(Eigen::Map<const n::Vector>(derivative.coefficients.data(),derivative.coefficients.size()))},
            {"correction",Values(Eigen::Map<const n::Vector>(derivative.correction.data(),derivative.correction.size()))},
            {"projected",Values(Eigen::Map<const n::Vector>(reduced.projected.data(),reduced.projected.size()))},
            {"jacobian",Values(Eigen::Map<const n::Vector>(reduced.jacobian.data(),reduced.jacobian.size()))},
            {"response",Values(reduced.response)}};
    }
#endif
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
        // Audit/capture runs are explicitly separate from timing repetitions.
#ifndef SPARSE_BASELINE_DRIVER
        for(int k=1;k<argc;++k) if(std::string(argv[k]).starts_with("--"))
        {
            const int end=argc; argc=k;
            for(;k<end;++k)
            {
                const std::string option=argv[k];
                if(option=="--svd-mode" && k+1<end) Mode(argv[++k]);
                else if(option=="--audit") audit=true;
                else if(option=="--capture" && k+1<end) {audit=true; capture=argv[++k]; std::filesystem::create_directories(capture);}
                else throw std::invalid_argument("Invalid benchmark option");
            }
        }
        if(audit) n::CompactSvdCaptureForTesting()=Capture;
#endif
        Eigen::setNbThreads(1); const std::string mode=argc>1 ? argv[1] : "";
#ifndef SPARSE_BASELINE_DRIVER
        if(mode=="replay" && argc==5)
        {
            Mode(argv[3]); const auto record=Read(argv[2],true);
            const std::string order=std::endian::native==std::endian::little ? "little" : "big";
            if(j::value_to<std::string>(record.at("byte_order"))!=order) throw std::runtime_error("Compact matrix byte order mismatch");
            n::Matrix a(j::value_to<Eigen::Index>(record.at("rows")),j::value_to<Eigen::Index>(record.at("columns")));
            const auto path=std::filesystem::path(argv[2]).parent_path()/j::value_to<std::string>(record.at("matrix_file"));
            const auto bytes=static_cast<std::size_t>(a.size())*sizeof(double);
            if(std::filesystem::file_size(path)!=bytes) throw std::runtime_error("Compact matrix size mismatch");
            std::ifstream f(path,std::ios::binary); f.read(reinterpret_cast<char *>(a.data()),static_cast<std::streamsize>(bytes));
            n::Vector rhs; const bool solve=!record.at("rhs").is_null();
            if(solve) {const auto values=j::value_to<std::vector<double>>(record.at("rhs"));rhs=Eigen::Map<const n::Vector>(values.data(),static_cast<Eigen::Index>(values.size()));}
            const auto started=Clock::now(); const auto result=n::CompactSvd(a,j::value_to<double>(record.at("relative_threshold")),
                j::value_to<double>(record.at("absolute_override")),solve ? &rhs : nullptr);
            auto report=SvdRecord(result); report["svd_wall_seconds"]=Seconds(started); Snapshot(argv[4],report); return result.valid ? 0 : 1;
        }
#endif
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
