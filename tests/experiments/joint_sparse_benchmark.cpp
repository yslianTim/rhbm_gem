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
#include <map>
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
std::string assessment_reduction{"observation-tsqr"};
std::string projected_reduction{"observation-tiled-qr"};
bool search_only{};
std::filesystem::path capture;
j::array svd_records;
n::OperatorRankMode operator_rank_mode{n::OperatorRankMode::Auto};
n::RankBudget operator_rank_budget{};
n::SchwarzPolicy schwarz_policy{};
n::SpqrOrdering spqr_ordering{n::SpqrOrdering::Colamd};
bool spqr_ordering_requested{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
n::OperatorFactorOwnershipKindForTesting factor_ownership{n::OperatorFactorOwnershipKindForTesting::ReuseAcceptedCopyOnWrite};
#endif
j::value Read(const char * path,bool precise=false)
{
    std::ifstream f(path); if(!f) throw std::runtime_error("Missing input");
    j::parse_options options; if(precise) options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(f),{}),{},options);
}
void Write(const char * path,const j::object & v) {std::ofstream f(path);f<<j::serialize(v)<<'\n';}
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
j::array Values(const n::Vector & v) {j::array out;for(auto x:v) out.push_back(std::isfinite(x) ? j::value(x) : j::value(nullptr));return out;}
j::value OptionalSize(const std::optional<std::size_t> & value) {return value ? j::value(*value) : j::value(nullptr);}
j::object LocalWitnessRecord(const n::FreeDesignLocalWitness & witness)
{
    const auto number=[](double value)->j::value {return std::isfinite(value) ? j::value(value) : j::value(nullptr);};
    return {{"groups",witness.groups},{"covered_columns",witness.covered_columns},
        {"total_columns",witness.total_columns},{"coverage_fraction",witness.coverage_fraction},
        {"exclusive_rows",witness.exclusive_rows},{"max_group_size",witness.max_group_size},
        {"minimum_lower",witness.minimum_lower ? j::value(*witness.minimum_lower) : j::value(nullptr)},
        {"rank_threshold_upper",number(witness.threshold_upper)},
        {"exclusive_rows_disjoint",witness.exclusive_rows_disjoint},
        {"would_certify",witness.would_certify},{"reason",std::string(witness.reason)}};
}
std::size_t ParseSize(const char * raw,bool allow_zero)
{
    const std::string value=raw;
    if(value.empty() || value.front()=='-') throw std::invalid_argument("Expected a nonnegative integer policy value");
    std::size_t used{}; const auto parsed=std::stoull(value,&used);
    if(used!=value.size() || parsed>std::numeric_limits<std::size_t>::max() || (!allow_zero && parsed==0))
        throw std::invalid_argument("Invalid integer policy value");
    return static_cast<std::size_t>(parsed);
}
std::size_t ParseMiB(const char * raw,bool allow_zero=false)
{
    const auto mib=ParseSize(raw,allow_zero);
    if(mib>std::numeric_limits<std::size_t>::max()/(1024*1024)) throw std::invalid_argument("Memory limit is too large");
    return mib*1024*1024;
}
double ParseSeconds(const char * raw)
{
    const std::string value=raw; std::size_t used{}; const double seconds=std::stod(value,&used);
    if(used!=value.size() || !std::isfinite(seconds) || seconds<0) throw std::invalid_argument("Expected a finite nonnegative rank time budget");
    return seconds;
}
n::OperatorRankMode ParseOperatorRankMode(const std::string & value)
{
    if(value=="auto") return n::OperatorRankMode::Auto;
    if(value=="dense") return n::OperatorRankMode::Dense;
    if(value=="spqr-bounds") return n::OperatorRankMode::SpqrBounds;
    throw std::invalid_argument("Expected --operator-rank auto|dense|spqr-bounds");
}
n::SpqrOrdering ParseSpqrOrdering(const std::string & value)
{
    if(value=="colamd") return n::SpqrOrdering::Colamd;
    if(value=="default") return n::SpqrOrdering::Default;
    if(value=="best") return n::SpqrOrdering::Best;
    if(value=="metis") return n::SpqrOrdering::Metis;
    throw std::invalid_argument("Expected --spqr-ordering colamd|default|best|metis");
}
void ConfigureSearchPolicy(n::SearchPolicy & policy)
{
    policy.operator_rank.mode=operator_rank_mode;
    policy.operator_rank.budget=operator_rank_budget;
    policy.schwarz=schwarz_policy;
}
void ConfigureSearchPolicy(n::EvaluationContext & context) {ConfigureSearchPolicy(context.search);}
j::object PolicyRecord(const n::SearchPolicy & policy)
{
    const bool rank_active=policy.method==n::SearchMethod::OperatorPcg;
    const auto resolved=rank_active ? n::ResolveOperatorRankBackend(policy.operator_rank.mode,n::ActiveSparseBackend()) :
        std::optional<n::FreeDesignRankBackend>{};
    return {{"search_method",n::SearchMethodName(policy.method)},
        {"sparse_backend",n::SparseBackendName(n::ActiveSparseBackend())},
        {"operator_rank_active",rank_active},
        {"operator_rank_mode",n::OperatorRankModeName(policy.operator_rank.mode)},
        {"resolved_rank_backend",resolved ? j::value(n::FreeDesignRankBackendName(*resolved)) : j::value(nullptr)},
        {"operator_rank_budget_seconds",policy.operator_rank.budget.seconds},
        {"operator_rank_budget_entries",policy.operator_rank.budget.entries},
        {"operator_rank_budget_workspace_bytes",policy.operator_rank.budget.workspace_bytes},
        {"preconditioner",n::PreconditionerName(policy.preconditioner)},
        {"schwarz_core_atoms",policy.schwarz.core_atoms},{"schwarz_overlap_hops",policy.schwarz.overlap_hops},
        {"schwarz_max_block_atoms",policy.schwarz.max_block_atoms},
        {"schwarz_storage_bytes",policy.schwarz.storage_bytes},{"schwarz_scratch_bytes",policy.schwarz.scratch_bytes}};
}
j::object Summary(std::vector<std::size_t> values)
{
    if(values.empty()) return {{"min",nullptr},{"median",nullptr},{"max",nullptr}};
    std::sort(values.begin(),values.end());
    const double median=values.size()%2 ? static_cast<double>(values[values.size()/2]) :
        .5*(static_cast<double>(values[values.size()/2-1])+static_cast<double>(values[values.size()/2]));
    return {{"min",values.front()},{"median",median},{"max",values.back()}};
}
j::object PartitionRecord(const n::PreconditionerPartition & partition)
{
    std::vector<std::size_t> cores,blocks,memberships; std::map<std::size_t,std::size_t> histogram;
    for(const auto & block:partition.blocks)
    {
        cores.push_back(block.core_atoms.size()); blocks.push_back(block.core_atoms.size()+block.overlap_atoms.size());
    }
    for(auto atom:partition.layout.full_atoms)
    {
        const auto count=partition.atom_blocks[static_cast<std::size_t>(atom)].size();
        memberships.push_back(count); ++histogram[count];
    }
    j::object census;
    for(const auto & [count,atoms]:histogram) census[std::to_string(count)]=atoms;
    const auto & work=n::SearchWorkForTesting();
    return {{"blocks",partition.blocks.size()},{"core_atom_counts",Summary(std::move(cores))},
        {"block_atom_counts",Summary(std::move(blocks))},{"atom_membership",Summary(std::move(memberships))},
        {"atom_membership_histogram",census},{"overlap_hops",partition.policy.overlap_hops},
        {"core_atoms",partition.policy.core_atoms},{"max_block_atoms",partition.policy.max_block_atoms},
        {"topology_bytes",work.topology_bytes},{"storage_bytes",work.storage_bytes},
        {"scratch_bytes_bound",work.scratch_bytes}};
}
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
j::object FactorResidency()
{
    const auto & residency=n::FactorResidencyWorkForTesting(); j::array factors,constructing,peak_ids,peak_generations;
    std::size_t current_count=residency.constructing.size(),current_bytes{};
    for(const auto & factor:residency.factors)
    {
        if(factor.alive) {++current_count; current_bytes+=factor.owned_factor_bytes;}
        factors.push_back(j::object{{"factor_id",factor.factor_id},{"generation",factor.generation},
            {"kind",factor.kind},{"role",factor.role},{"created_at_stage",factor.created_at_stage},
            {"destroyed_at_stage",factor.destroyed_at_stage.empty() ? j::value(nullptr) : j::value(factor.destroyed_at_stage)},
            {"search_stage",factor.search_stage.empty() ? j::value(nullptr) : j::value(factor.search_stage)},
            {"destroyed_search_stage",factor.destroyed_search_stage.empty() ? j::value(nullptr) : j::value(factor.destroyed_search_stage)},
            {"rows",factor.rows},{"columns",factor.columns},{"nonzeros",factor.nonzeros},
            {"r_nonzeros",factor.r_nonzeros},{"h_nonzeros",factor.h_nonzeros},
            {"owned_factor_bytes_estimate",factor.owned_factor_bytes},{"created_seconds",factor.created_seconds},
            {"destroyed_seconds",factor.alive ? j::value(nullptr) : j::value(factor.destroyed_seconds)},
            {"alive_at_snapshot",factor.alive}});
    }
    for(const auto & factor:residency.constructing)
    {
        current_bytes+=factor.owned_factor_bytes_estimate;
        constructing.push_back(j::object{{"factor_id",factor.factor_id},{"generation",factor.generation},
            {"kind",factor.kind},{"role",factor.role},{"stage",factor.stage},{"rows",factor.rows},
            {"columns",factor.columns},{"nonzeros",factor.nonzeros},
            {"owned_factor_bytes_estimate",factor.owned_factor_bytes_estimate}});
    }
    for(const auto id:residency.peak_rss_factor_ids) peak_ids.push_back(id);
    for(const auto & generation:residency.peak_rss_factor_generations) peak_generations.push_back(j::value(generation));
    return {{"current_concurrent_factor_count",current_count},{"current_concurrent_owned_bytes_estimate",current_bytes},
        {"maximum_concurrent_factor_count",residency.maximum_concurrent_factor_count},
        {"maximum_concurrent_owned_bytes_estimate",residency.maximum_concurrent_owned_bytes},
        {"maximum_concurrent_factor_stage",residency.maximum_concurrent_factor_stage},
        {"maximum_concurrent_owned_bytes_stage",residency.maximum_concurrent_owned_bytes_stage},
        {"process_peak_rss_bytes_observed_at_factor_events",residency.peak_rss_bytes_observed},
        {"factor_count_at_process_peak_rss_event",residency.peak_rss_factor_count},
        {"owned_factor_bytes_at_process_peak_rss_event_estimate",residency.peak_rss_owned_bytes},
        {"stage_at_process_peak_rss_event",residency.peak_rss_stage},{"factor_ids_at_process_peak_rss_event",peak_ids},
        {"factor_generations_at_process_peak_rss_event",peak_generations},{"factors",factors},
        {"constructions_in_progress",constructing},
        {"semantics","owned bytes estimate includes CHOLMOD current allocations plus owned sparse design/R arrays; it is separate from process RSS"}};
}
#endif
void Snapshot(const char * output,j::object & report)
{
#ifndef SPARSE_BASELINE_DRIVER
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    report["factor_residency"]=FactorResidency();
#endif
    const auto & resource=n::ResourceWorkForTesting();
    if(resource.enabled)
    {
        j::array shapes,sparse_shapes,phases,completed_stages;
        j::object stage_seconds,stage_calls,stage_completed_calls,stage_dimensions,stage_nnz;
        for(const auto & r:resource.dense_shapes) shapes.push_back(j::object{{"phase",r.phase},{"role",r.role},
            {"rows",r.rows},{"columns",r.columns},{"observations",r.observations},{"maximum_matrix_bytes",r.maximum_bytes}});
        for(const auto & r:resource.sparse_shapes) sparse_shapes.push_back(j::object{{"phase",r.phase},{"role",r.role},
            {"rows",r.rows},{"columns",r.columns},{"observations",r.observations},{"maximum_nonzeros",r.maximum_nonzeros}});
        for(const auto & r:resource.phases) phases.push_back(j::object{{"phase",r.phase},{"calls",r.calls},{"inclusive_seconds",r.inclusive_seconds}});
        for(const auto & stage:resource.search_stages)
        {
            stage_seconds[stage.stage]=stage.seconds;
            stage_calls[stage.stage]=stage.calls;
            stage_completed_calls[stage.stage]=stage.completed_calls;
            stage_dimensions[stage.stage]=j::object{{"rows",stage.rows},{"columns",stage.columns}};
            stage_nnz[stage.stage]=stage.nonzeros;
        }
        completed_stages.clear();
        for(const auto & stage:resource.completed_search_stages) completed_stages.push_back(j::value(stage));
        const auto optional_stage=[](const std::string & stage)->j::value {return stage.empty() ? j::value(nullptr) : j::value(stage);};
        report["last_completed_search_stage"]=optional_stage(resource.last_completed_search_stage);
        report["active_search_stage"]=optional_stage(resource.active_search_stage);
        report["completed_search_stages"]=completed_stages;
        report["stage_seconds"]=stage_seconds;
        report["stage_calls"]=stage_calls;
        report["stage_completed_calls"]=stage_completed_calls;
        report["stage_dimensions"]=stage_dimensions;
        report["stage_nnz"]=stage_nnz;
        report["resources"]=j::object{{"dense_shape_probes",shapes},{"sparse_shape_probes",sparse_shapes},{"phases",phases},
            {"semantics","phase times overlap; shape probes describe known matrices, not every allocation"}};
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
        {"derivative_inclusive_seconds",w.derivative_seconds},
        {"spqr_symbolic_rows",w.symbolic_rows},{"spqr_symbolic_columns",w.symbolic_columns},
        {"spqr_symbolic_input_nonzeros",w.symbolic_input_nonzeros},{"spqr_numeric_rows",w.numeric_rows},
        {"spqr_numeric_columns",w.numeric_columns},{"spqr_numeric_input_nonzeros",w.numeric_input_nonzeros},
        {"spqr_fixed_factor_rows",w.fixed_factor_rows},{"spqr_fixed_factor_columns",w.fixed_factor_columns},
        {"spqr_fixed_factor_input_nonzeros",w.fixed_factor_input_nonzeros},
        {"spqr_fixed_factor_nonzeros",w.fixed_factor_nonzeros},
        {"spqr_fixed_factor_storage_bytes",w.fixed_factor_storage_bytes}};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    report["work"].as_object()["spqr_ordering"]=n::SpqrOrderingName(n::SpqrOrderingForTesting());
    report["operator_factor_representation"]=n::OperatorFactorRepresentationName(n::OperatorFactorRepresentationForTesting());
    report["operator_factor_ownership"]=n::OperatorFactorOwnershipName(factor_ownership);
    report["spqr_factorization"]=j::object{
        {"ordering",n::SpqrOrderingName(n::SpqrOrderingForTesting())},
        {"symbolic",j::object{{"calls",w.symbolic},{"seconds",w.symbolic_seconds},{"matrix_rows",w.symbolic_rows},
            {"matrix_columns",w.symbolic_columns},{"matrix_nonzeros",w.symbolic_input_nonzeros}}},
        {"numeric",j::object{{"calls",w.numeric},{"seconds",w.numeric_seconds},{"matrix_rows",w.numeric_rows},
            {"matrix_columns",w.numeric_columns},{"matrix_nonzeros",w.numeric_input_nonzeros},
            {"reported_factor_nonzeros",w.factor_nonzeros}}},
        {"fixed_factor",j::object{{"calls",w.fixed_factorizations},{"seconds",w.fixed_factor_seconds},
            {"matrix_rows",w.fixed_factor_rows},{"matrix_columns",w.fixed_factor_columns},
            {"matrix_nonzeros",w.fixed_factor_input_nonzeros},{"reported_factor_nonzeros",w.fixed_factor_nonzeros},
            {"factor_exported_storage_bytes",w.fixed_factor_storage_bytes}}},
        {"native_operator_factor",j::object{{"calls",w.native_operator_factorizations},
            {"seconds",w.native_operator_factor_seconds}}}};
#endif
    if(audit) report["svd_records"]=svd_records;
#endif
    Write(output,report);
}

struct SearchStageSnapshotContext {const char * output{}; j::object * report{};};
void SearchStageSnapshot(const n::ResourceWork & resource,void * raw_context)
{
    (void)resource;
    auto & context=*static_cast<SearchStageSnapshotContext *>(raw_context);
    (*context.report)["failure_stage"]="search";
    Snapshot(context.output,*context.report);
}
class SearchStageSnapshotRegistration
{
public:
    SearchStageSnapshotRegistration(const char * output,j::object & report):context_{output,&report}
    {n::SetResourceStageObserverForTesting(SearchStageSnapshot,&context_);}
    ~SearchStageSnapshotRegistration() {Stop();}
    SearchStageSnapshotRegistration(const SearchStageSnapshotRegistration &)=delete;
    SearchStageSnapshotRegistration & operator=(const SearchStageSnapshotRegistration &)=delete;
    void Stop()
    {
        if(active_) {n::SetResourceStageObserverForTesting(nullptr,nullptr); active_=false;}
    }
private:
    SearchStageSnapshotContext context_;
    bool active_{true};
};

#ifndef SPARSE_BASELINE_DRIVER
j::object SearchWork()
{
    j::object out;
#ifndef PR23_BASELINE_DRIVER
    const auto & w=n::SearchWorkForTesting(); const auto & op=n::OperatorWorkForTesting();
    j::array regularizations,pcg_iteration_counts;
    for(const auto count:w.pcg_iteration_counts) pcg_iteration_counts.push_back(count);
    const auto pcg_iterations=Summary(w.pcg_iteration_counts);
    for(const auto & r:w.regularizations) regularizations.push_back(j::object{{"local_build",r.local_build},{"factor_build",r.factor_build},
        {"block",r.block},{"lambda",r.lambda},{"damping",r.damping},{"tau",r.tau},{"attempt",r.attempt}});
    out={{"linearizations",w.linearizations},{"pcg_solves",w.pcg_solves},{"pcg_iterations",w.pcg_iterations},
        {"pcg_iteration_counts",pcg_iteration_counts},
        {"pcg_iterations_min",pcg_iterations.at("min")},{"pcg_iterations_median",pcg_iterations.at("median")},
        {"pcg_iterations_max",pcg_iterations.at("max")},
        {"pcg_iterations_mean",w.pcg_iteration_counts.empty() ? j::value(nullptr) :
            j::value(static_cast<double>(w.pcg_iterations)/static_cast<double>(w.pcg_iteration_counts.size()))},
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
        {"operator_rank_certificate",n::FreeDesignRankCertificateName(op.rank_certificate)},
        {"operator_rank_local_witness",LocalWitnessRecord(op.rank_local_witness)},
        {"operator_rank_checks",op.rank_checks},{"operator_rank_status",op.rank_status.empty() ? "not-run" : op.rank_status},
        {"operator_rank_reason",op.rank_reason.empty() ? "not-run" : op.rank_reason},
        {"operator_rank_rows",op.rank_rows},{"operator_rank_columns",op.rank_columns},
        {"operator_rank_entries",op.rank_entries},{"operator_rank_workspace_bytes",op.rank_workspace_bytes},
        {"operator_rank_work_stage",n::FreeDesignRankWorkStageName(op.rank_work_stage)},
        {"operator_rank_estimated_total_entries",OptionalSize(op.rank_estimated_total_entries)},
        {"operator_rank_estimated_remaining_entries",OptionalSize(op.rank_estimated_remaining_entries)},
        {"operator_rank_estimated_reconstruction_entries",OptionalSize(op.rank_estimated_reconstruction_entries)},
        {"operator_rank_design_nonzeros",op.rank_design_nonzeros},{"operator_rank_r_nonzeros",op.rank_r_nonzeros},
        {"operator_rank_reflector_nonzeros",op.rank_reflector_nonzeros},{"operator_rank_reflectors",op.rank_reflectors},
        {"operator_rank_compact_extractions",op.rank_compact_extractions},
        {"operator_rank_free_design_svds",op.rank_free_design_svds},
        {"operator_apply_seconds",op.apply_seconds},{"operator_adjoint_seconds",op.adjoint_seconds},
        {"operator_applications",op.applications},{"operator_adjoints",op.adjoints},{"regularizations",regularizations}};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    out["native_factor_accepted"]=op.native_factor_accepted;
    out["native_factor_fallbacks"]=op.native_factor_fallbacks;
    out["accepted_factor_reuse_attempts"]=op.accepted_factor_reuse_attempts;
    out["accepted_factor_reuse_accepted"]=op.accepted_factor_reuse_accepted;
    out["accepted_factor_reuse_fallbacks"]=op.accepted_factor_reuse_fallbacks;
    out["accepted_factor_reuse_fallback_reason"]=op.accepted_factor_reuse_fallback_reason;
    out["operator_factor_ownership"]=op.factor_ownership;
    out["factor_residency"]=FactorResidency();
#endif
#endif
    return out;
}
j::object AssessmentWork()
{
    const auto & work=n::AssessmentWorkForTesting();
    return {{"assessments",work.assessments},{"reference_evaluations",work.reference_evaluations},
        {"compact_attempts",work.compact_attempts},{"compact_accepted",work.compact_accepted},
        {"compact_boundary_fallbacks",work.compact_boundary_fallbacks},
        {"compact_other_fallbacks",work.compact_other_fallbacks}};
}
j::object TiledQrTelemetry(const n::TiledQrTelemetry & work)
{
    return {{"role",work.role},{"append_calls",work.append_calls},{"rows_processed",work.rows_processed},
        {"columns",work.columns},{"responses",work.responses},{"maximum_assembled_rows",work.maximum_assembled_rows},
        {"maximum_dense_design_bytes",work.maximum_dense_design_bytes},
        {"maximum_dense_response_bytes",work.maximum_dense_response_bytes},{"qr_seconds",work.qr_seconds}};
}
j::object AssessmentTelemetry()
{
    const auto & work=n::AssessmentWorkForTesting(); j::array stages,completed; j::object seconds,dimensions;
    for(const auto & stage:work.stages)
    {
        const bool is_complete=stage.calls==stage.completed_calls;
        stages.push_back(j::object{{"name",stage.name},{"calls",stage.calls},
            {"completed_calls",stage.completed_calls},{"rows",stage.rows},{"columns",stage.columns},
            {"seconds",stage.seconds},{"complete",is_complete}});
        dimensions[stage.name]=j::object{{"calls",stage.calls},{"completed_calls",stage.completed_calls},
            {"rows",stage.rows},{"columns",stage.columns}};
        if(is_complete) {completed.push_back(j::value(stage.name)); seconds[stage.name]=stage.seconds;}
    }
    const auto optional_stage=[](const std::string & name)->j::value {
        return name.empty() ? j::value(nullptr) : j::value(name);
    };
    const auto stage_seconds=[&](const std::string & name) {
        const auto found=std::find_if(work.stages.begin(),work.stages.end(),[&](const auto & stage){return stage.name==name;});
        return found==work.stages.end() ? 0. : found->seconds;
    };
    const auto & derivative=n::DerivativeWorkForTesting();
    const double rows_seconds=stage_seconds("derivative-rows");
    const double projected_seconds=stage_seconds("derivative-projected-qr");
    const double jacobian_seconds=stage_seconds("derivative-jacobian-qr");
    const double compact_jacobian_seconds=stage_seconds("derivative-compact-jacobian-qr");
    const double norms_seconds=stage_seconds("derivative-norms");
    const auto & reduction=derivative.projected_reduction;
    const double projected_reduction_seconds=reduction.kind=="observation-tiled-qr" ?
        derivative.projected_qr.qr_seconds : reduction.seconds+
        (reduction.fallbacks ? derivative.projected_qr.qr_seconds : 0.);
    j::object projected_reduction_telemetry{
        {"projected_reduction_kind",reduction.kind},
        {"projected_reduction_fallback_reason",reduction.fallback_reason},
        {"projected_reduction_attempts",reduction.attempts},
        {"projected_reduction_accepted",reduction.accepted},
        {"projected_reduction_fallbacks",reduction.fallbacks},
        {"projected_reduction_seconds",projected_reduction_seconds},
        {"projected_factor_rows",reduction.factor_rows},
        {"projected_factor_columns",reduction.factor_columns},
        {"observation_projected_rows_processed",reduction.observation_projected_rows_processed},
        {"compact_projected_rows_processed",reduction.compact_projected_rows_processed},
        {"maximum_dense_bytes",reduction.maximum_dense_bytes},
        {"observations",reduction.observations},{"free_design_columns",reduction.free_design_columns},
        {"width_columns",reduction.width_columns},{"raw_nonzeros",reduction.raw_nonzeros},
        {"raw_density",reduction.raw_density},{"q_transformed_nonzeros",reduction.q_transformed_nonzeros},
        {"q_transformed_density",reduction.q_transformed_density},{"tail_rows",reduction.tail_rows},
        {"tail_nonzeros",reduction.tail_nonzeros},{"tail_density",reduction.tail_density},
        {"q_transform_seconds",reduction.q_transform_seconds},
        {"q_transformed_storage_bytes",reduction.q_transformed_storage_bytes},
        {"tail_storage_bytes",reduction.tail_storage_bytes},
        {"sparse_rows",reduction.sparse_rows},{"sparse_columns",reduction.sparse_columns},
        {"sparse_nnz",reduction.sparse_nonzeros},
        {"symbolic_seconds",reduction.symbolic_seconds},{"numeric_seconds",reduction.numeric_seconds},
        {"ordering",reduction.ordering}};
    j::object derivative_reduction{{"active_micro_stage",optional_stage(
            work.active_stage=="derivative-rows" || work.active_stage=="derivative-projected-qr" ||
            work.active_stage=="derivative-jacobian-qr" || work.active_stage=="derivative-compact-jacobian-qr" ||
            work.active_stage=="derivative-norms" ? work.active_stage : std::string{})},
        {"reduction_inclusive_seconds",stage_seconds("derivative-reduction")},
        {"rows_seconds",rows_seconds},{"projected_qr_seconds",projected_seconds},
        {"jacobian_qr_seconds",jacobian_seconds},{"compact_jacobian_qr_seconds",compact_jacobian_seconds},
        {"projected_reduction",projected_reduction_telemetry},
        {"norms_seconds",norms_seconds},
        {"exclusive_substage_seconds",rows_seconds+projected_seconds+jacobian_seconds+compact_jacobian_seconds+norms_seconds},
        {"tile_count",derivative.tile_count},{"tiled_qr",j::object{
            {"projected",TiledQrTelemetry(derivative.projected_qr)},
            {"jacobian",TiledQrTelemetry(derivative.jacobian_qr)},
            {"compact_jacobian",TiledQrTelemetry(derivative.compact_jacobian_qr)}}}};
    return {{"assessments",work.assessments},{"reference_evaluations",work.reference_evaluations},
        {"compact_attempts",work.compact_attempts},{"compact_accepted",work.compact_accepted},
        {"compact_boundary_fallbacks",work.compact_boundary_fallbacks},
        {"compact_other_fallbacks",work.compact_other_fallbacks},
        {"last_assessment_stage",optional_stage(work.last_stage)},
        {"active_assessment_stage",optional_stage(work.active_stage)},
        {"completed_assessment_stages",completed},{"completed_stage_seconds",seconds},
        {"stage_dimensions",dimensions},{"stages",stages},
        {"derivative_reduction_micro_attribution",derivative_reduction}};
}
struct AssessmentSnapshotContext {const char * output{}; j::object * report{};};
void AssessmentSnapshot(const n::AssessmentWork &,void * raw_context)
{
    auto & context=*static_cast<AssessmentSnapshotContext *>(raw_context);
    const auto telemetry=AssessmentTelemetry();
    (*context.report)["assessment_work"]=AssessmentWork();
    (*context.report)["assessment_telemetry"]=telemetry;
    (*context.report)["last_assessment_stage"]=telemetry.at("last_assessment_stage");
    (*context.report)["active_assessment_stage"]=telemetry.at("active_assessment_stage");
    (*context.report)["assessment_stage"]=telemetry.at("active_assessment_stage");
    (*context.report)["completed_assessment_stages"]=telemetry.at("completed_assessment_stages");
    (*context.report)["completed_stage_seconds"]=telemetry.at("completed_stage_seconds");
    (*context.report)["stage_dimensions"]=telemetry.at("stage_dimensions");
    Snapshot(context.output,*context.report);
}
class AssessmentSnapshotRegistration
{
public:
    AssessmentSnapshotRegistration(const char * output,j::object & report):context_{output,&report}
    {n::SetAssessmentStageObserverForTesting(AssessmentSnapshot,&context_);}
    ~AssessmentSnapshotRegistration() {Stop();}
    AssessmentSnapshotRegistration(const AssessmentSnapshotRegistration &)=delete;
    AssessmentSnapshotRegistration & operator=(const AssessmentSnapshotRegistration &)=delete;
    void Stop()
    {
        if(active_) {n::SetAssessmentStageObserverForTesting(nullptr,nullptr); active_=false;}
    }
private:
    AssessmentSnapshotContext context_;
    bool active_{true};
};
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
    ConfigureSearchPolicy(context);
    j::object report{{"stage","search"},{"search_kind",search_kind},{"atoms",b.size()},{"rows",domain.rows},
        {"measurement_scope",search_only ? "search-only" : "joint_search_and_returned_state_assessment"},
        {"assessment_reduction",assessment_reduction},{"projected_reduction",projected_reduction}};
    report["solver_policy"]=PolicyRecord(context.search);
    auto & resource=n::ResourceWorkForTesting(); const bool resources_enabled=resource.enabled;
    resource={}; resource.enabled=resources_enabled;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    n::AssessmentWorkForTesting()={};
    n::DerivativeWorkForTesting()={};
    n::ResetFactorResidencyWorkForTesting();
#endif
    n::SparseWorkForTesting()={}; n::SearchWorkForTesting()={}; n::OperatorWorkForTesting()={};
    Snapshot(output,report);
    SearchStageSnapshotRegistration search_snapshot(output,report);
    auto search=n::SearchProfile(domain,y,b,context);
    search_snapshot.Stop(); report["failure_stage"]=nullptr;
    report["search"]=second_stage_test::matched::runtime_json::Search(search,context,domain.rows);
    auto & search_record=report["search"].as_object();
    search_record["accepted_objective"]=search.accepted_objective ? j::value(*search.accepted_objective) : j::value(nullptr);
    search_record["accepted_gradient_inf_norm"]=search.accepted_gradient_inf_norm ? j::value(*search.accepted_gradient_inf_norm) : j::value(nullptr);
    search_record["returned_search_state"]=second_stage_test::matched::runtime_json::Values(search.eta);
    if(search_only) report["search"].as_object().erase("runtime_convergence");
#ifndef PR23_BASELINE_DRIVER
    if(search_kind=="schwarz")
    {
        const auto partition=n::SearchPartition(domain,context,context.search.schwarz);
        report["partition"]=PartitionRecord(*partition);
    }
#endif
    report["search_seconds"]=search.seconds; report["search_work"]=SearchWork();
    report["search_global_derivative_preparations"]=n::SparseWorkForTesting().derivative_preparations;
#ifndef PR23_BASELINE_DRIVER
    report["free_columns"]=n::OperatorWorkForTesting().rank_columns;
    report["design_nonzeros"]=n::OperatorWorkForTesting().rank_design_nonzeros;
#endif
    if(search_only)
    {
        report["assessment_execution"]="not-run";
        report["returned_assessment"]=nullptr;
        report["assessment_work"]=AssessmentWork();
        report["scope_description"]="This profile measures nonlinear Operator-PCG search only. It does not perform returned-state assessment and does not establish runtime convergence or endpoint qualification.";
        report["stage"]="complete"; Snapshot(output,report); return;
    }
    report["assessment_execution"]="running";
    report["returned_assessment"]=nullptr;
    report["assessment_stage"]=nullptr;
    report["last_assessment_stage"]=nullptr;
    report["active_assessment_stage"]=nullptr;
    report["completed_assessment_stages"]=j::array{};
    report["completed_stage_seconds"]=j::object{};
    report["stage_dimensions"]=j::object{};
    report["assessment_telemetry"]=AssessmentTelemetry();
    report["stage"]="assessment"; Snapshot(output,report);
    AssessmentSnapshotRegistration assessment_snapshot(output,report);
    const auto fit=n::AssessComponentSearch(domain,y,context,std::move(search));
    assessment_snapshot.Stop();
    report["search_completed"]=fit.search_success;
    report["assessment_seconds"]=fit.assessment_seconds;
    report["assessment_execution"]="completed";
    report["assessment"]=second_stage_test::matched::runtime_json::Assessment(fit.assessment);
    report["returned_state"]=fit.trusted_state ? j::value(second_stage_test::matched::runtime_json::Endpoint(*fit.trusted_state)) : j::value(nullptr);
    report["returned_assessment"]=fit.trusted_assessment ? j::value(second_stage_test::matched::runtime_json::Assessment(*fit.trusted_assessment)) : j::value(nullptr);
    report["assessment_work"]=AssessmentWork();
    report["assessment_telemetry"]=AssessmentTelemetry();
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
    if(!fixed_mode.empty())
    {
        auto configured=context; ConfigureSearchPolicy(configured);
        configured.search.method=n::SearchMethod::OperatorPcg;
        if(fixed_preconditioner=="identity") configured.search.preconditioner=n::PreconditionerKind::Identity;
        else if(fixed_preconditioner=="diagonal") configured.search.preconditioner=n::PreconditionerKind::Diagonal;
        else if(fixed_preconditioner=="schwarz") configured.search.preconditioner=n::PreconditionerKind::Schwarz;
        RunFixed(domain,y,b,configured,output); return;
    }
#endif
    if(!search_kind.empty()) {RunSearch(domain,y,b,context,output); return;}
#endif
    j::object report{{"rows",domain.rows},{"atoms",b.size()},{"initial_b",Values(b)},{"stage","primary"}};Snapshot(output,report);
#ifndef SPARSE_BASELINE_DRIVER
    n::SparseWorkForTesting()={}; n::ResetFactorResidencyWorkForTesting(); n::LinearWorkspace workspace;
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
                else if(option=="--search-only") search_only=true;
                else if(option=="--assessment-reduction" && k+1<end)
                {
                    assessment_reduction=argv[++k];
                    if(assessment_reduction=="observation-tsqr")
                        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::ObservationTsqr;
                    else if(assessment_reduction=="compact-stack-qr")
                        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
                    else throw std::invalid_argument("Expected --assessment-reduction observation-tsqr|compact-stack-qr");
                }
                else if(option=="--projected-reduction" && k+1<end)
                {
                    projected_reduction=argv[++k];
                    if(projected_reduction=="observation-tiled-qr")
                        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ObservationTiledQr;
                    else if(projected_reduction=="structured-compact-qr")
                        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::StructuredCompactQr;
                    else if(projected_reduction=="projected-tail-census")
                        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ProjectedTailCensus;
                    else throw std::invalid_argument("Expected --projected-reduction observation-tiled-qr|structured-compact-qr|projected-tail-census");
                }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                else if(option=="--operator-factor-representation" && k+1<end)
                {
                    const std::string representation=argv[++k];
                    if(representation=="exported-fixed")
                        n::OperatorFactorRepresentationForTesting()=n::OperatorFactorRepresentation::ExportedFixed;
                    else if(representation=="native-qr")
                        n::OperatorFactorRepresentationForTesting()=n::OperatorFactorRepresentation::NativeQr;
                    else throw std::invalid_argument("Expected --operator-factor-representation exported-fixed|native-qr");
                }
                else if(option=="--operator-factor-ownership" && k+1<end)
                {
                    const std::string ownership=argv[++k];
                    if(ownership=="dedicated-fixed") factor_ownership=n::OperatorFactorOwnershipKindForTesting::DedicatedFixed;
                    else if(ownership=="dedicated-native") factor_ownership=n::OperatorFactorOwnershipKindForTesting::DedicatedNative;
                    else if(ownership=="reuse-accepted-copy-on-write") factor_ownership=n::OperatorFactorOwnershipKindForTesting::ReuseAcceptedCopyOnWrite;
                    else if(ownership=="reuse-accepted-handoff") factor_ownership=n::OperatorFactorOwnershipKindForTesting::ReuseAcceptedHandoff;
                    else throw std::invalid_argument("Expected --operator-factor-ownership dedicated-fixed|dedicated-native|reuse-accepted-copy-on-write|reuse-accepted-handoff");
                    n::OperatorFactorOwnershipForTesting()=factor_ownership;
                }
#endif
                else if(option=="--operator-rank" && k+1<end) operator_rank_mode=ParseOperatorRankMode(argv[++k]);
                else if(option=="--operator-rank-seconds" && k+1<end) operator_rank_budget.seconds=ParseSeconds(argv[++k]);
                else if(option=="--operator-rank-work-entries" && k+1<end) operator_rank_budget.entries=ParseSize(argv[++k],true);
                else if(option=="--operator-rank-workspace-mib" && k+1<end) operator_rank_budget.workspace_bytes=ParseMiB(argv[++k],true);
                else if(option=="--schwarz-core-atoms" && k+1<end) schwarz_policy.core_atoms=ParseSize(argv[++k],false);
                else if(option=="--schwarz-overlap-hops" && k+1<end) schwarz_policy.overlap_hops=ParseSize(argv[++k],true);
                else if(option=="--schwarz-max-block-atoms" && k+1<end) schwarz_policy.max_block_atoms=ParseSize(argv[++k],false);
                else if(option=="--schwarz-storage-mib" && k+1<end) schwarz_policy.storage_bytes=ParseMiB(argv[++k]);
                else if(option=="--schwarz-scratch-mib" && k+1<end) schwarz_policy.scratch_bytes=ParseMiB(argv[++k]);
                else if(option=="--spqr-ordering" && k+1<end)
                {spqr_ordering=ParseSpqrOrdering(argv[++k]); spqr_ordering_requested=true;}
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
        if(schwarz_policy.max_block_atoms<schwarz_policy.core_atoms)
            throw std::invalid_argument("Schwarz max block atoms must be at least core atoms");
        if(spqr_ordering_requested && !n::SpqrOrderingAvailable(spqr_ordering))
            throw std::invalid_argument(std::string("spqr-ordering-unavailable: ")+n::SpqrOrderingName(spqr_ordering)+" is not provided by this SPQR build");
        n::SpqrOrderingForTesting()=spqr_ordering;
        if(search_only && (search_kind.empty() || search_kind=="legacy"))
            throw std::invalid_argument("Search-only profile requires OperatorPcg");
        if(audit) n::CompactSvdCaptureForTesting()=Capture;
#endif
        Eigen::setNbThreads(1); const std::string mode=argc>1 ? argv[1] : "";
#ifndef SPARSE_BASELINE_DRIVER
        if(mode=="synthetic" && argc==6)
        {
            const std::string topology=argv[2],phase=argv[4]; const int atoms=std::stoi(argv[3]);
            if(phase!="prepare" && phase!="fixed" && phase!="workflow" && phase!="local" && phase!="rank" && phase!="rank-oracle") throw std::invalid_argument("Invalid synthetic phase");
            const bool operator_search=search_kind=="identity" || search_kind=="diagonal" || search_kind=="schwarz";
            if(atoms>512 && phase!="prepare" && phase!="local" && phase!="rank" && !(phase=="fixed" && operator_search))
                throw std::invalid_argument("Large solves require an explicit OperatorPcg search route");
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
                if(phase=="rank" || phase=="rank-oracle")
                {
                    auto rank_context=data.context; ConfigureSearchPolicy(rank_context);
                    report["solver_policy"]=PolicyRecord(rank_context.search);
                    RunRank(state,rank_context,phase=="rank-oracle",report,argv[5]);
                }
#endif
                if(phase=="local")
                {
                    report["stage"]="local"; Snapshot(argv[5],report);
                    try {
                        const auto partition=n::BuildPreconditionerPartition(data.input,data.layout,data.context.search.schwarz);
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
                n::SearchPolicy policy; ConfigureSearchPolicy(policy);
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
