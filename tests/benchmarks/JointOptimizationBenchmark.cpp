#include "core/detail/joint_component/ActiveSetBoundedLM.hpp"
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/JointSyntheticWorkload.hpp"

#include <boost/hash2/sha2.hpp>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <sys/resource.h>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

struct Digest
{
    boost::hash2::sha2_256 hash;
    template<class Integer> requires std::is_integral_v<Integer>
    void Add(Integer value)
    {
        std::array<unsigned char,8> bytes{};
        const auto bits=static_cast<std::uint64_t>(value);
        for(unsigned k=0;k<8;++k) bytes[k]=static_cast<unsigned char>((bits>>(8*k))&255);
        hash.update(bytes.data(),bytes.size());
    }
    void Add(double value) {Add(std::bit_cast<std::uint64_t>(value));}
    void Add(const std::string & value)
    {Add(value.size()); hash.update(value.data(),value.size());}
    template<class Derived>
    void Add(const Eigen::MatrixBase<Derived> & value)
    {
        Add(static_cast<std::uint64_t>(value.rows())); Add(static_cast<std::uint64_t>(value.cols()));
        for(Eigen::Index k=0;k<value.size();++k) Add(value.derived().coeff(k));
    }
    template<class SparseType>
    void AddSparse(const SparseType & value)
    {
        Add(static_cast<std::uint64_t>(value.rows())); Add(static_cast<std::uint64_t>(value.cols()));
        Add(static_cast<std::uint64_t>(value.nonZeros()));
        for(Eigen::Index outer=0;outer<value.outerSize();++outer)
            for(typename SparseType::InnerIterator entry(value,outer);entry;++entry)
            {Add(static_cast<std::uint64_t>(entry.row())); Add(static_cast<std::uint64_t>(entry.col())); Add(entry.value());}
    }
    std::string Finish() {return boost::hash2::to_string(hash.result());}
};

std::string BoundedDigest(const joint::ActiveSetStep & result,const joint::BoundedLmTestWork & work)
{
    Digest digest;
    digest.Add(result.valid); digest.Add(result.status); digest.Add(result.damping);
    digest.Add(result.step); digest.Add(result.gradient); digest.Add(result.multipliers);
    digest.Add(result.scaled_step_norm); digest.Add(result.kkt_inf_norm);
    digest.Add(result.iterations); digest.Add(result.releases);
    for(bool active:result.active) digest.Add(active);
    for(const auto & candidate:work.damping_trace)
    {
        digest.Add(candidate.damping); digest.Add(candidate.phase); digest.Add(candidate.status);
        digest.Add(candidate.valid); digest.Add(candidate.step); digest.Add(candidate.gradient);
        digest.Add(candidate.multipliers); digest.Add(candidate.scaled_step_norm);
        digest.Add(candidate.kkt_inf_norm); digest.Add(candidate.iterations); digest.Add(candidate.releases);
        for(bool active:candidate.active) digest.Add(active);
        digest.Add(candidate.within_trust_region);
        digest.Add(candidate.search_accepted);
        digest.Add(candidate.active_set_sequence.size());
        for(const auto & active:candidate.active_set_sequence)
        {
            digest.Add(active.size());
            for(bool value:active) digest.Add(value);
        }
        digest.Add(candidate.step_sequence.size());
        for(const auto & step:candidate.step_sequence) digest.Add(step);
        digest.Add(candidate.gradient_sequence.size());
        for(const auto & gradient:candidate.gradient_sequence) digest.Add(gradient);
        digest.Add(candidate.kkt_tolerance_sequence.size());
        for(double tolerance:candidate.kkt_tolerance_sequence) digest.Add(tolerance);
        digest.Add(candidate.free_ranks.size());
        for(Eigen::Index rank:candidate.free_ranks) digest.Add(rank);
    }
    return digest.Finish();
}

std::string DerivativeDigest(const joint::TiledDifferential & prepared,
    const joint::ReducedDifferential & reduced)
{
    Digest digest;
    digest.Add(prepared.valid); digest.Add(prepared.reference_order); digest.Add(prepared.projected_only);
    digest.Add(prepared.reason); digest.Add(prepared.scale);
    digest.AddSparse(prepared.free_design); digest.AddSparse(prepared.raw);
    digest.Add(prepared.coefficients); digest.Add(prepared.correction);
    if(prepared.free_design_rank)
    {
        digest.Add(true); digest.Add(prepared.free_design_rank->rank);
        digest.Add(prepared.free_design_rank->threshold); digest.Add(prepared.free_design_rank->singular_values);
    }
    else digest.Add(false);
    digest.Add(reduced.valid); digest.Add(reduced.reason); digest.Add(reduced.jacobian);
    digest.Add(reduced.response); digest.Add(reduced.jacobian_norms);
    return digest.Finish();
}

double Seconds(Clock::time_point started)
{return std::chrono::duration<double>(Clock::now()-started).count();}

void ReportSamples(const std::string & name,std::vector<double> samples)
{
    auto sorted=samples; std::sort(sorted.begin(),sorted.end());
    const double median=sorted[sorted.size()/2];
    std::cout<<name<<"_ms_per_call_median="<<median*1000.0
        <<" range_ms="<<sorted.front()*1000.0<<".."<<sorted.back()*1000.0<<" samples_ms=";
    for(std::size_t k=0;k<samples.size();++k)
        std::cout<<(k ? "," : "")<<samples[k]*1000.0;
    std::cout<<'\n';
}

int Run()
{
    constexpr int samples=5,repeats=5;
    Matrix j(4096,8);
    for(Eigen::Index row=0;row<j.rows();++row)
        for(Eigen::Index col=0;col<j.cols();++col)
        {
            const double phase=static_cast<double>(row+1)*static_cast<double>(col+1);
            j(row,col)=(std::sin(.00019*phase)+.4*std::cos(.00031*(phase+7.0)))/std::sqrt(4096.0);
        }
    Vector truth=Vector::Constant(j.cols(),.6),r=-(j*truth),d=Vector::Ones(j.cols());
    Vector upper=Vector::Constant(j.cols(),.15); upper(0)=.02;
    std::vector<double> bounded_times;
    joint::ActiveSetStep bounded;
    joint::BoundedLmTestWork bounded_work;
    std::string bounded_digest;
    for(int sample=0;sample<samples+1;++sample)
    {
        joint::BoundedLmTestWork aggregate;
        const auto started=Clock::now();
        for(int repeat=0;repeat<repeats;++repeat)
        {
            joint::BoundedLmWorkForTesting()={};
            bounded=joint::SolveBoundedTrustRegionStep(j,r,d,upper,.3,.5);
            if(!bounded.valid) throw std::runtime_error("bounded LM: "+bounded.status);
            const auto & work=joint::BoundedLmWorkForTesting();
            aggregate.damping_solves+=work.damping_solves;
            aggregate.bracket_solves+=work.bracket_solves;
            aggregate.bisection_solves+=work.bisection_solves;
            aggregate.active_set_iterations+=work.active_set_iterations;
            aggregate.releases+=work.releases;
            aggregate.subproblem_validation_calls+=work.subproblem_validation_calls;
            aggregate.kkt_scale_evaluations+=work.kkt_scale_evaluations;
            aggregate.jacobian_norm_evaluations+=work.jacobian_norm_evaluations;
            aggregate.diagonal_squared_norm_evaluations+=work.diagonal_squared_norm_evaluations;
            aggregate.bound_tolerance_evaluations+=work.bound_tolerance_evaluations;
            aggregate.dense_qr_factorizations+=work.dense_qr_factorizations;
            aggregate.workspace_allocation_sites+=work.workspace_allocation_sites;
            aggregate.workspace_resize_events+=work.workspace_resize_events;
            aggregate.workspace_peak_bytes=std::max(aggregate.workspace_peak_bytes,work.workspace_peak_bytes);
            aggregate.dense_qr_factorization_seconds+=work.dense_qr_factorization_seconds;
            aggregate.damping_search_seconds+=work.damping_search_seconds;
        }
        const double elapsed=Seconds(started)/repeats;
        if(sample) bounded_times.push_back(elapsed);
        bounded_work=std::move(aggregate);
    }
    joint::BoundedLmWorkForTesting()={};
    joint::BoundedLmWorkForTesting().capture_trace=true;
    bounded=joint::SolveBoundedTrustRegionStep(j,r,d,upper,.3,.5);
    bounded_digest=BoundedDigest(bounded,joint::BoundedLmWorkForTesting());
    ReportSamples("bounded_lm",std::move(bounded_times));
    std::cout<<"bounded_lm_damping_solves_per_call="<<bounded_work.damping_solves/repeats
        <<" bracket_solves_per_call="<<bounded_work.bracket_solves/repeats
        <<" bisection_solves_per_call="<<bounded_work.bisection_solves/repeats
        <<" active_set_iterations_per_call="<<bounded_work.active_set_iterations/repeats
        <<" subproblem_validations_per_call="<<bounded_work.subproblem_validation_calls/repeats
        <<" kkt_scale_evaluations_per_call="<<bounded_work.kkt_scale_evaluations/repeats
        <<" jacobian_norm_evaluations_per_call="<<bounded_work.jacobian_norm_evaluations/repeats
        <<" diagonal_squared_norm_evaluations_per_call="<<bounded_work.diagonal_squared_norm_evaluations/repeats
        <<" bound_tolerance_evaluations_per_call="<<bounded_work.bound_tolerance_evaluations/repeats
        <<" dense_qr_count_per_call="<<bounded_work.dense_qr_factorizations/repeats
        <<" dense_qr_ms_per_call="<<bounded_work.dense_qr_factorization_seconds*1000.0/repeats
        <<" workspace_allocation_sites_per_call="<<bounded_work.workspace_allocation_sites/repeats
        <<" workspace_resize_events_per_call="<<bounded_work.workspace_resize_events/repeats
        <<" workspace_peak_bytes_estimate="<<bounded_work.workspace_peak_bytes
        <<" bitwise_trace_sha256="<<bounded_digest<<'\n';

    auto input=std::make_shared<rhbm_gem::core::JointProblemInput>(
        second_stage_test::SyntheticJointWorkload("chain",30));
    joint::Domain domain(input);
    const auto context=joint::CreateContext(input);
    Vector y=Eigen::Map<const Vector>(input->observations.data(),
        static_cast<Eigen::Index>(input->observations.size()));
    Vector eta=Vector::Constant(static_cast<Eigen::Index>(input->atom_ids.size()),std::log(.55));
    const auto evaluation=joint::EvaluateProfile(domain,y,eta,&context);
    if(!evaluation.valid) throw std::runtime_error("profile preparation: "+evaluation.reason);

    std::vector<double> prepare_times,reduce_times;
    joint::TiledDifferential prepared;
    joint::ReducedDifferential reduced;
    std::string derivative_digest;
    joint::DerivativeWork derivative_work,reduction_work;
    for(int sample=0;sample<samples+1;++sample)
    {
        joint::DerivativeWork aggregate;
        auto & work=joint::DerivativeWorkForTesting();
        const auto started=Clock::now();
        for(int repeat=0;repeat<repeats;++repeat)
        {
            work={};
            prepared=joint::PrepareDerivative(evaluation,context.scale,&context);
            if(!prepared.valid) throw std::runtime_error("derivative preparation: "+prepared.reason);
            if(repeat==0) aggregate=work;
            else
            {
                aggregate.prepare_calls+=work.prepare_calls;
                aggregate.tiled_qr_invocations+=work.tiled_qr_invocations;
                aggregate.tile_count+=work.tile_count;
                aggregate.triplet_capacity_growths+=work.triplet_capacity_growths;
                aggregate.triplet_reserved_entries+=work.triplet_reserved_entries;
                aggregate.triplet_reserved_bytes+=work.triplet_reserved_bytes;
                aggregate.tiled_qr_scratch_materializations+=work.tiled_qr_scratch_materializations;
                aggregate.tiled_qr_scratch_resize_events+=work.tiled_qr_scratch_resize_events;
                aggregate.cancellation_fallbacks+=work.cancellation_fallbacks;
                aggregate.reference_order_recomputations+=work.reference_order_recomputations;
                aggregate.maximum_generated_rows=std::max(aggregate.maximum_generated_rows,work.maximum_generated_rows);
                aggregate.maximum_reduction_rows=std::max(aggregate.maximum_reduction_rows,work.maximum_reduction_rows);
                aggregate.maximum_tile_workspace_bytes=std::max(aggregate.maximum_tile_workspace_bytes,work.maximum_tile_workspace_bytes);
                aggregate.prepare_seconds+=work.prepare_seconds;
                aggregate.tiled_qr_seconds+=work.tiled_qr_seconds;
            }
        }
        const double elapsed=Seconds(started)/repeats;
        if(sample) prepare_times.push_back(elapsed);
        else derivative_digest=DerivativeDigest(prepared,reduced);
        derivative_work=std::move(aggregate);

        joint::TiledDifferential one=prepared;
        joint::DerivativeWork reduction_aggregate;
        const auto reduce_started=Clock::now();
        for(int repeat=0;repeat<repeats;++repeat)
        {
            work={};
            reduced=joint::ReduceDerivative(one,evaluation.residual);
            if(!reduced.valid) throw std::runtime_error("derivative reduction: "+reduced.reason);
            if(repeat==0) reduction_aggregate=work;
            else
            {
                reduction_aggregate.reduction_tile_scratch_materializations+=work.reduction_tile_scratch_materializations;
                reduction_aggregate.reduction_tile_scratch_resize_events+=work.reduction_tile_scratch_resize_events;
                reduction_aggregate.tiled_qr_scratch_materializations+=work.tiled_qr_scratch_materializations;
                reduction_aggregate.tiled_qr_scratch_resize_events+=work.tiled_qr_scratch_resize_events;
                reduction_aggregate.maximum_tile_workspace_bytes=std::max(
                    reduction_aggregate.maximum_tile_workspace_bytes,work.maximum_tile_workspace_bytes);
                reduction_aggregate.tiled_qr_seconds+=work.tiled_qr_seconds;
            }
        }
        const double reduce_elapsed=Seconds(reduce_started)/repeats;
        if(sample) reduce_times.push_back(reduce_elapsed);
        reduction_work=std::move(reduction_aggregate);
        if(!sample) derivative_digest=DerivativeDigest(prepared,reduced);
    }
    ReportSamples("derivative_prepare",std::move(prepare_times));
    ReportSamples("derivative_reduce",std::move(reduce_times));
    std::cout<<"derivative_reduce_tile_scratch_materializations_per_call="
        <<reduction_work.reduction_tile_scratch_materializations/repeats
        <<" derivative_reduce_tile_scratch_resize_events_per_call="
        <<reduction_work.reduction_tile_scratch_resize_events/repeats
        <<" derivative_reduce_tiled_qr_scratch_materializations_per_call="
        <<reduction_work.tiled_qr_scratch_materializations/repeats
        <<" derivative_reduce_tiled_qr_ms_per_call="<<reduction_work.tiled_qr_seconds*1000.0/repeats<<'\n';
    std::cout<<"derivative_rows="<<evaluation.residual.size()
        <<" derivative_columns="<<evaluation.x.cols()
        <<" derivative_prepare_qr_invocations="<<derivative_work.tiled_qr_invocations/repeats
        <<" derivative_prepare_tiles="<<derivative_work.tile_count/repeats
        <<" triplet_capacity_growths_per_call="<<derivative_work.triplet_capacity_growths/repeats
        <<" triplet_reserved_entries_per_call="<<derivative_work.triplet_reserved_entries/repeats
        <<" triplet_reserved_bytes_per_call="<<derivative_work.triplet_reserved_bytes/repeats
        <<" tiled_qr_scratch_materializations_per_call="<<derivative_work.tiled_qr_scratch_materializations/repeats
        <<" tiled_qr_scratch_resize_events_per_call="<<derivative_work.tiled_qr_scratch_resize_events/repeats
        <<" max_generated_rows="<<derivative_work.maximum_generated_rows
        <<" max_reduction_rows="<<derivative_work.maximum_reduction_rows
        <<" maximum_tile_workspace_bytes_estimate="<<derivative_work.maximum_tile_workspace_bytes
        <<" cancellation_fallbacks_per_call="<<derivative_work.cancellation_fallbacks/repeats
        <<" reference_order_recomputations_per_call="<<derivative_work.reference_order_recomputations/repeats
        <<" tiled_qr_ms_per_call="<<derivative_work.tiled_qr_seconds*1000.0/repeats
        <<" bitwise_trace_sha256="<<derivative_digest<<'\n';

    rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    Vector joint_eta=Vector::Constant(3,std::log(.55));
    joint::FixedNeighborPolicy policy;
    policy.core_atoms=3; policy.maximum_sweeps=3; policy.stop_after_stationarity=false;
    policy.certify_endpoint=false; policy.collect_telemetry=true; policy.eta_upper_bound=std::log(2.2);
    std::vector<double> joint_times;
    std::size_t joint_subproblem_calls{},joint_workspace_allocations{},joint_workspace_resizes{};
    joint::FixedNeighborDiagnosticResult last_joint_result;
    std::string candidate_trace_digest;
    for(int sample=0;sample<samples+1;++sample)
    {
        auto & work=joint::BoundedLmWorkForTesting(); work={};
        problem=rhbm_gem::core::JointProblem(second_stage_test::SyntheticJointWorkload("chain",3));
        Digest candidate_trace;
        if(sample==0)
            policy.candidate_observer=[&](const joint::FixedNeighborCandidateEvent & event) {
                candidate_trace.Add(event.sweep); candidate_trace.Add(event.block);
                candidate_trace.Add(event.accepted); candidate_trace.Add(event.unchanged);
                candidate_trace.Add(event.full_replay_ran); candidate_trace.Add(event.core_atoms.size());
                for(const auto atom:event.core_atoms) candidate_trace.Add(atom);
                candidate_trace.Add(event.affected_rows.size());
                for(const auto row:event.affected_rows) candidate_trace.Add(row);
                candidate_trace.Add(event.accepted_before.eta); candidate_trace.Add(event.accepted_before.beta);
                candidate_trace.Add(event.accepted_before.prediction); candidate_trace.Add(event.accepted_before.residual);
                candidate_trace.Add(event.accepted_before.objective); candidate_trace.Add(event.candidate_eta);
                candidate_trace.Add(event.candidate_beta); candidate_trace.Add(event.full_prediction);
                candidate_trace.Add(event.full_residual); candidate_trace.Add(event.full_objective);
            };
        else policy.candidate_observer={};
        const auto started=Clock::now();
        auto result=joint::RunFixedNeighborDiagnostic(problem,joint_eta,policy);
        const double elapsed=Seconds(started);
        if(sample)
        {
            joint_times.push_back(elapsed);
            joint_subproblem_calls+=work.subproblem_validation_calls;
            joint_workspace_allocations+=work.workspace_allocation_sites;
            joint_workspace_resizes+=work.workspace_resize_events;
        }
        if(sample==0) candidate_trace_digest=candidate_trace.Finish();
        if(result.reason.empty()) throw std::runtime_error("synthetic Joint diagnostic returned no status");
        if(sample==samples) last_joint_result=std::move(result);
    }
    ReportSamples("synthetic_joint_search",std::move(joint_times));
    std::cout<<"synthetic_joint_subproblems_per_run="<<joint_subproblem_calls/samples
        <<" workspace_allocations_per_run="<<joint_workspace_allocations/samples
        <<" workspace_resizes_per_run="<<joint_workspace_resizes/samples<<'\n';
    const auto & joint_data=rhbm_gem::core::JointProblemAccess::Get(problem);
    std::vector<std::uint8_t> informative(static_cast<std::size_t>(joint_data.domain.rows));
    for(const auto row:joint_data.layout.informative_rows) informative.at(row)=1;
    std::vector<std::size_t> row_contributors(informative.size());
    std::size_t support_entries{};
    for(const auto atom:joint_data.layout.full_atoms) for(const auto & support:joint_data.domain.atoms[atom])
        if(informative.at(static_cast<std::size_t>(support.row)))
        {++support_entries; ++row_contributors[static_cast<std::size_t>(support.row)];}
    std::vector<std::size_t> contributor_counts;
    for(std::size_t row=0;row<informative.size();++row)
        if(informative[row]) contributor_counts.push_back(row_contributors[row]);
    std::sort(contributor_counts.begin(),contributor_counts.end());
    std::vector<std::size_t> block_rows;
    for(const auto & block:last_joint_result.blocks) if(block.sweep==1) block_rows.push_back(block.affected_rows);
    std::sort(block_rows.begin(),block_rows.end());
    std::size_t accepted_candidates{},unchanged_candidates{};
    for(const auto & block:last_joint_result.blocks)
    {accepted_candidates+=block.accepted; unchanged_candidates+=block.status=="unchanged";}
    const auto rejected_candidates=last_joint_result.total_block_solves-accepted_candidates-unchanged_candidates;
    Digest final_state_digest;
    final_state_digest.Add(last_joint_result.state.eta); final_state_digest.Add(last_joint_result.state.beta);
    final_state_digest.Add(last_joint_result.state.prediction); final_state_digest.Add(last_joint_result.state.residual);
    final_state_digest.Add(last_joint_result.state.objective);
    const auto replay_fraction=last_joint_result.search_seconds>0 ?
        last_joint_result.work.candidate_replay_seconds/last_joint_result.search_seconds : 0.0;
    std::cout<<"synthetic_joint_census_observations="<<joint_data.domain.rows
        <<" informative_rows="<<joint_data.layout.informative_rows.size()
        <<" targets="<<joint_data.layout.full_atoms.size()
        <<" support_entries="<<support_entries
        <<" mean_contributors_per_informative_row="<<static_cast<double>(support_entries)/
            static_cast<double>(joint_data.layout.informative_rows.size())
        <<" contributors_min_median_max="<<(contributor_counts.empty()?0:contributor_counts.front())<<":"
        <<(contributor_counts.empty()?0:contributor_counts[contributor_counts.size()/2])<<":"
        <<(contributor_counts.empty()?0:contributor_counts.back())
        <<" blocks="<<last_joint_result.block_preparations
        <<" total_affected_row_entries="<<last_joint_result.memory_estimate.shape.affected_row_entries
        <<" affected_rows_min_median_max="<<(block_rows.empty()?0:block_rows.front())<<":"
        <<(block_rows.empty()?0:block_rows[block_rows.size()/2])<<":"
        <<(block_rows.empty()?0:block_rows.back())
        <<" candidates="<<last_joint_result.total_block_solves
        <<" full_candidate_replays="<<last_joint_result.work.full_candidate_replays
        <<" unchanged_replay_skips="<<last_joint_result.work.candidate_replay_skipped_unchanged
        <<" accepted_blocks="<<accepted_candidates<<" rejected_blocks="<<rejected_candidates
        <<" unchanged_blocks="<<unchanged_candidates
        <<" preparation_seconds="<<last_joint_result.work.structural_preparation_seconds
        <<" candidate_replay_seconds="<<last_joint_result.work.candidate_replay_seconds
        <<" local_search_seconds="<<last_joint_result.work.local_search_seconds
        <<" sweep_replay_seconds="<<last_joint_result.work.sweep_replay_seconds
        <<" global_state_seconds="<<last_joint_result.work.sweep_global_state_seconds
        <<" total_search_seconds="<<last_joint_result.search_seconds
        <<" candidate_replay_fraction="<<replay_fraction
        <<" amdahl_speedup_upper_bound="<<(replay_fraction<1.0 ? 1.0/(1.0-replay_fraction) : 0.0)
        <<" estimated_peak_bytes="<<last_joint_result.global_state_memory_estimate.estimated_peak_bytes
        <<" baseline_state_sha256="<<final_state_digest.Finish()
        <<" candidate_trace_sha256="<<candidate_trace_digest
        <<" final_objective="<<last_joint_result.state.objective
        <<" sweep_count="<<last_joint_result.sweep_count
        <<" global_ac_kkt="<<(last_joint_result.final_sweep ? last_joint_result.final_sweep->global_ac_kkt : 0.0)
        <<" projected_width_kkt="<<(last_joint_result.final_sweep ?
            last_joint_result.final_sweep->projected_width_gradient_inf_norm : 0.0)
        <<" stop_reason="<<last_joint_result.reason
        <<" endpoint_certified="<<last_joint_result.endpoint_certified<<'\n';
    std::cout<<"synthetic_joint_initial_eta=";
    for(Eigen::Index k=0;k<joint_eta.size();++k) std::cout<<(k?",":"")<<joint_eta(k);
    std::cout<<" final_eta=";
    for(Eigen::Index k=0;k<last_joint_result.state.eta.size();++k)
        std::cout<<(k?",":"")<<last_joint_result.state.eta(k);
    std::cout<<" final_beta=";
    for(Eigen::Index k=0;k<last_joint_result.state.beta.size();++k)
        std::cout<<(k?",":"")<<last_joint_result.state.beta(k);
    std::cout<<'\n';
    joint::BoundedLmWorkForTesting()={};
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)==0)
    {
#if defined(__APPLE__)
        const auto peak_rss_bytes=static_cast<std::uint64_t>(usage.ru_maxrss);
#else
        const auto peak_rss_bytes=static_cast<std::uint64_t>(usage.ru_maxrss)*1024;
#endif
        std::cout<<"benchmark_peak_rss_bytes="<<peak_rss_bytes<<'\n';
    }
    return 0;
}
}

int main()
{
    try {return Run();}
    catch(const std::exception & error)
    {std::cerr<<"joint optimization benchmark failed: "<<error.what()<<'\n'; return 1;}
}
