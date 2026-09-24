#include "core/detail/joint_component/Problem.hpp"
#ifndef SPARSE_BASELINE_DRIVER
#include "core/detail/joint_component/SparseFactor.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#endif
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointFixtureSupport.hpp"
#include "support/JointRuntimeJson.hpp"
#ifndef SPARSE_BASELINE_DRIVER
#include "support/JointOperatorWorkload.hpp"
#include "core/detail/joint_component/ProfileJacobianOperator.hpp"
#include "core/detail/JointUncertainty.hpp"
#ifndef PR23_BASELINE_DRIVER
#include "core/detail/joint_component/OperatorSearch.hpp"
#endif
#endif
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
bool operator_audit{};
std::string search_kind;
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
    const auto & resource=n::ResourceWorkForTesting();
    if(resource.enabled)
    {
        j::array shapes,phases;
        for(const auto & r:resource.dense_shapes) shapes.push_back(j::object{{"phase",r.phase},{"role",r.role},
            {"rows",r.rows},{"columns",r.columns},{"observations",r.observations},{"maximum_matrix_bytes",r.maximum_bytes}});
        for(const auto & r:resource.phases) phases.push_back(j::object{{"phase",r.phase},{"calls",r.calls},{"inclusive_seconds",r.inclusive_seconds}});
        report["resources"]=j::object{{"dense_shape_probes",shapes},{"phases",phases},
            {"semantics","phase times overlap; shapes are known matrix probes, not an allocation trace"}};
    }
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

#ifndef SPARSE_BASELINE_DRIVER
j::object SearchWork()
{
    j::object out;
#ifndef PR23_BASELINE_DRIVER
    const auto & w=n::SearchWorkForTesting(); const auto & op=n::OperatorWorkForTesting();
    j::array regularizations;
    for(const auto & r:w.regularizations) regularizations.push_back(j::object{{"local_build",r.local_build},{"factor_build",r.factor_build},
        {"block",r.block},{"lambda",r.lambda},{"damping",r.damping},{"tau",r.tau},{"attempt",r.attempt}});
    out={{"linearizations",w.linearizations},{"pcg_solves",w.pcg_solves},{"pcg_iterations",w.pcg_iterations},
        {"damping_trials",w.damping_trials},{"local_builds",w.local_builds},{"factor_builds",w.factor_builds},
        {"inverse_actions",w.inverse_actions},{"topology_bytes",w.topology_bytes},{"storage_bytes",w.storage_bytes},
        {"scratch_bytes_bound",w.scratch_bytes},{"maximum_block_atoms",w.maximum_block_atoms},
        {"partition_seconds",w.partition_seconds},{"metric_seconds",w.metric_seconds},{"local_seconds",w.local_seconds},
        {"factor_seconds",w.factor_seconds},{"inverse_seconds",w.inverse_seconds},{"pcg_seconds",w.pcg_seconds},
        {"maximum_lambda",w.maximum_lambda},{"maximum_tau",w.maximum_tau},{"last_relative_residual",w.last_relative_residual},
        {"operator_normal_seconds",op.normal_seconds},{"operator_normals",op.normals},
        {"operator_design_seconds",op.design_seconds},{"operator_fixed_factor_seconds",op.factor_seconds},
        {"operator_compact_seconds",op.compact_seconds},{"operator_svd_seconds",op.svd_seconds},
        {"operator_prepare_seconds",op.preparation_seconds},{"operator_rank_seconds",op.rank_seconds},
        {"operator_apply_seconds",op.apply_seconds},{"operator_adjoint_seconds",op.adjoint_seconds},
        {"operator_applications",op.applications},{"operator_adjoints",op.adjoints},{"regularizations",regularizations}};
#endif
    return out;
}
void RunSearch(const n::Domain & domain,n::VectorRef y,const n::Vector & b,n::EvaluationContext context,const char * output)
{
#ifndef PR23_BASELINE_DRIVER
    if(search_kind!="legacy")
    {
        context.search.method=n::SearchMethod::OperatorPcg;
        if(search_kind=="identity") context.search.preconditioner=n::PreconditionerKind::Identity;
        else if(search_kind=="diagonal") context.search.preconditioner=n::PreconditionerKind::Diagonal;
        else if(search_kind=="schwarz") context.search.preconditioner=n::PreconditionerKind::Schwarz;
        else throw std::invalid_argument("Invalid search kind");
    }
#else
    if(search_kind!="legacy") throw std::invalid_argument("Baseline supports only legacy search");
#endif
    j::object report{{"stage","search"},{"search_kind",search_kind},{"atoms",b.size()},{"rows",domain.rows}};
    Snapshot(output,report);
    auto search=n::SearchProfile(domain,y,b,context);
    report["search"]=second_stage_test::matched::runtime_json::Search(search,context,domain.rows);
#ifndef PR23_BASELINE_DRIVER
    if(search_kind=="schwarz")
    {
        const auto partition=n::SearchPartition(domain,context);
        j::array memberships; for(const auto & blocks:partition->atom_blocks) if(!blocks.empty()) memberships.push_back(blocks.size());
        report["partition"]=j::object{{"blocks",partition->blocks.size()},{"atom_memberships",memberships}};
    }
#endif
    report["search_seconds"]=search.seconds; report["search_work"]=SearchWork();
    report["search_global_derivative_preparations"]=n::SparseWorkForTesting().derivative_preparations;
    report["stage"]="assessment"; Snapshot(output,report);
    const auto fit=n::AssessComponentSearch(domain,y,context,std::move(search));
    report["search_completed"]=fit.search_success;
    report["assessment_seconds"]=fit.assessment_seconds;
    report["assessment"]=second_stage_test::matched::runtime_json::Assessment(fit.assessment);
    report["returned_state"]=fit.trusted_state ? j::value(second_stage_test::matched::runtime_json::Endpoint(*fit.trusted_state)) : j::value(nullptr);
    report["returned_assessment"]=fit.trusted_assessment ? j::value(second_stage_test::matched::runtime_json::Assessment(*fit.trusted_assessment)) : j::value(nullptr);
    report["stage"]="complete"; Snapshot(output,report);
}

#endif
#if !defined(SPARSE_BASELINE_DRIVER) && !defined(PR23_BASELINE_DRIVER)
#include "support/JointFixedDiagnostic.hpp"
#ifndef PR4_BASELINE_DRIVER
#include "support/JointRankDiagnostic.hpp"
#endif
#endif
void Run(const n::Domain & domain,n::VectorRef y,const n::Vector & b,const n::EvaluationContext & context,const char * output)
{

#ifndef SPARSE_BASELINE_DRIVER
#ifndef PR23_BASELINE_DRIVER
    if(!fixed_mode.empty()) {RunFixed(domain,y,b,context,output); return;}
#endif
    if(!search_kind.empty()) {RunSearch(domain,y,b,context,output); return;}
#endif
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
    if(operator_audit && primary.valid)
    {
        report["stage"]="operator"; Snapshot(output,report);
        n::OperatorWorkForTesting()={}; const n::ProfileJacobianOperator op(primary,context);
        j::object result{{"valid",op.Valid()},{"reason",op.Reason()},{"krylov_iterations",nullptr}};
        if(!op.Valid() && !derivative.valid) result["passed"]=op.Reason()==derivative.reason;
        if(op.Valid() && derivative.valid)
        {
            const n::Vector v=n::Vector::LinSpaced(op.Columns(),-.3,.7);
            n::Vector w(op.Rows()); for(Eigen::Index r=0;r<w.size();++r) w(r)=std::sin(.17*static_cast<double>(r));
            const auto actual=op.Apply(v),adjoint=op.ApplyAdjoint(w);
            n::Vector expected(op.Rows()),expected_adjoint=n::Vector::Zero(op.Columns());
            n::Matrix projected,jacobian;
            {n::ResourcePhase phase("operator-oracle");
            for(Eigen::Index first=0;first<op.Rows();first+=1024)
            {
                const auto count=std::min<Eigen::Index>(1024,op.Rows()-first);
                derivative.Rows(first,count,projected,jacobian);
                expected.segment(first,count)=jacobian*v;
                expected_adjoint+=jacobian.transpose()*w.segment(first,count);
            }}
            const double apply_error=(actual-expected).norm()/std::max(1e-12,expected.norm());
            const double adjoint_error=(adjoint-expected_adjoint).norm()/std::max(1e-12,expected_adjoint.norm());
            const double duality=std::abs(actual.dot(w)-v.dot(adjoint))/std::max({1.,actual.norm()*w.norm(),v.norm()*adjoint.norm()});
            result["apply_relative_error"]=apply_error; result["adjoint_relative_error"]=adjoint_error;
            result["duality_scaled_error"]=duality;
            result["passed"]=apply_error<=1e-8 && adjoint_error<=1e-8 && duality<=1e-12;
            result["raw_nonzeros"]=op.RawNonZeros(); result["free_columns"]=op.FreeColumns();
        }
        const auto & ow=n::OperatorWorkForTesting();
        result["prepare_seconds"]=ow.preparation_seconds; result["rank_seconds"]=ow.rank_seconds;
        result["apply_seconds"]=ow.apply_seconds; result["adjoint_seconds"]=ow.adjoint_seconds;
        result["rank_checks"]=ow.rank_checks; report["operator"]=result;
    }
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
                else if(option=="--operator") operator_audit=true;
                else if(option=="--search" && k+1<end) search_kind=argv[++k];
#ifndef PR23_BASELINE_DRIVER
                else if(option=="--fixed" && k+1<end) fixed_mode=argv[++k];
                else if(option=="--fixed-preconditioner" && k+1<end) fixed_preconditioner=argv[++k];
                else if(option=="--state" && k+1<end) fixed_state=argv[++k];
#endif
                else if(option=="--resources") n::ResourceWorkForTesting().enabled=true;
                else if(option=="--capture" && k+1<end) {audit=true; capture=argv[++k]; std::filesystem::create_directories(capture);}
                else throw std::invalid_argument("Invalid benchmark option");
            }
        }
        if(audit) n::CompactSvdCaptureForTesting()=Capture;
#endif
        Eigen::setNbThreads(1); const std::string mode=argc>1 ? argv[1] : "";
#ifndef SPARSE_BASELINE_DRIVER
        if(mode=="synthetic" && argc==6)
        {
            const std::string topology=argv[2],phase=argv[4]; const int atoms=std::stoi(argv[3]);
            if(phase!="prepare" && phase!="fixed" && phase!="workflow" && phase!="local" && phase!="rank" && phase!="rank-oracle") throw std::invalid_argument("Invalid synthetic phase");
            if(atoms>512 && phase!="prepare" && phase!="local" && phase!="rank") throw std::invalid_argument("Large workloads are preparation/local-only");
            const auto started=Clock::now(); const c::JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
            const double construction_seconds=Seconds(started);
            const auto & data=c::JointProblemAccess::Get(problem);
            std::size_t memberships{}; for(const auto & a:problem.Input().support) memberships+=a.size();
            if(data.partition.components.size()!=1 || data.layout.full_atoms.size()!=static_cast<std::size_t>(atoms) ||
                memberships!=static_cast<std::size_t>(atoms)*515) throw std::runtime_error("Invalid synthetic support census");
            if(phase=="fixed") {Run(data.domain,data.y,n::Vector::Constant(atoms,.55),data.context,argv[5]); return 0;}
            j::object report{{"generator","frozen-lattice-v1"},{"topology",topology},{"atoms",atoms},
                {"rows",data.y.size()},{"memberships",memberships},{"components",data.partition.components.size()},
                {"spacing",3.5},{"grid",.5},{"support",2.5},{"truth_a",2.},{"truth_c",.2},{"truth_b",.5},
                {"initial_b",.55},{"noise","1e-5*sin(.13*z+.17*y+.19*x), integer half-angstrom coordinates"},
                {"construction_seconds",construction_seconds}};
            const auto hash_started=Clock::now(); report["input_sha256"]=second_stage_test::OperatorWorkloadHash(problem.Input());
            report["fingerprint_seconds"]=Seconds(hash_started);
            if(phase=="prepare" || phase=="local" || phase=="rank" || phase=="rank-oracle")
            {
                const auto basis_started=Clock::now(); n::Vector beta(2*atoms);
                for(int a=0;a<atoms;++a) {beta(2*a)=2; beta(2*a+1)=.2;}
                const auto state=n::EvaluateState(data.domain,data.y,n::Vector::Constant(atoms,std::log(.55)),beta,data.context);
                report["basis_seconds"]=Seconds(basis_started); report["raw_state_valid"]=state.valid;
                report["design_nonzeros"]=state.x.nonZeros(); report["derivative_nonzeros"]=state.derivative.nonZeros();
                report["residual_norm"]=state.residual.norm();
                report["not_run"]=j::array{"ac-solve","rank","operator","reference","search","assessment","uncertainty"};
#ifndef PR23_BASELINE_DRIVER
#ifndef PR4_BASELINE_DRIVER
                if(phase=="rank" || phase=="rank-oracle") RunRank(state,data.context,phase=="rank-oracle",report,argv[5]);
#endif
                if(phase=="local")
                {
                    report["stage"]="local"; Snapshot(argv[5],report);
                    try {
                        const auto partition=n::BuildPreconditionerPartition(data.input,data.layout);
                        {
                            const auto free=n::FreeColumnMapping(*partition,state.beta);
                            report["free_coordinates"]=free.dimension;
                            report["free_mapping_blocks"]=free.blocks.size();
                            for(const auto & block:free.blocks) for(std::size_t j=0;j<block.global.size();++j)
                                if(block.LocalIndex(block.global[j])!=static_cast<Eigen::Index>(j)) throw std::runtime_error("preconditioner-mapping-failed");
                        }
                        const auto metric=n::WidthMetric(n::WidthNorms(n::RawWidthDerivative(state),data.context.scale));
                        const n::PreconditionerContext pc{std::make_shared<const n::LinearizationIdentity>(),n::PreconditionerSpace::Width,metric,1e-3};
                        const n::SchwarzModel model(*partition,state,data.context.scale,pc); const n::SchwarzPreconditioner inverse(model,pc);
                        const n::Vector u=n::Vector::LinSpaced(atoms,-.5,.7),v=n::Vector::LinSpaced(atoms,.1,1.);
                        const auto a=inverse.ApplyInverse(u,pc),b=inverse.ApplyInverse(v,pc);
                        const double error=std::abs(u.dot(b)-v.dot(a))/std::max({1.,u.norm()*b.norm(),v.norm()*a.norm()});
                        report["local_available"]=true; report["local_passed"]=error<=1e-12 && u.dot(a)>0 && v.dot(b)>0;
                        report["symmetry_error"]=error; report["blocks"]=partition->blocks.size();
                    } catch(const std::runtime_error & e) {report["local_available"]=false; report["local_reason"]=e.what();}
                    report["search_work"]=SearchWork();
                }
#else
                if(phase=="local") throw std::invalid_argument("No local baseline implementation");
#endif

            }
            else
            {
#ifndef PR23_BASELINE_DRIVER
                n::SearchPolicy policy;
                if(!search_kind.empty() && search_kind!="legacy")
                {
                    policy.method=n::SearchMethod::OperatorPcg;
                    if(search_kind=="identity") policy.preconditioner=n::PreconditionerKind::Identity;
                    else if(search_kind=="diagonal") policy.preconditioner=n::PreconditionerKind::Diagonal;
                    else if(search_kind=="schwarz") policy.preconditioner=n::PreconditionerKind::Schwarz;
                    else throw std::invalid_argument("Invalid search kind");
                }
                const auto fit=n::FitWithSearchPolicy(problem,std::vector<double>(static_cast<std::size_t>(atoms),.55),policy);
#else
                const auto fit=c::FitJointComponents(problem,std::vector<double>(static_cast<std::size_t>(atoms),.55));
#endif
                const auto captured=c::CaptureJointAnalysisResult(fit);
                const auto uncertainty=c::detail::ComputeJointUncertainty(problem,captured);
                report["state_available"]=fit.assembled_state.has_value(); report["search_completed"]=fit.search_completed;
                report["uncertainty_outputs"]=uncertainty.size();
                report["search_seconds"]=fit.costs.search_seconds; report["assessment_seconds"]=fit.costs.assessment_seconds;
                report["assembly_seconds"]=fit.costs.assembly_seconds;
            }
            rusage usage{};getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
            report["peak_rss_bytes"]=usage.ru_maxrss;
#else
            report["peak_rss_bytes"]=usage.ru_maxrss*1024;
#endif
            report["stage"]="complete"; Snapshot(argv[5],report); return 0;
        }
#endif
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
