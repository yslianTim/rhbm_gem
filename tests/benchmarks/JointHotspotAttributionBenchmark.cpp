#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestNumerics.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <boost/hash2/sha2.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;

struct Sample
{
    joint::FixedNeighborWork::CandidateReplayBreakdown replay;
    joint::ProfileSearchWork profile;
    double wall_seconds{},candidate_copy_seconds{},candidate_replay_seconds{},cache_update_seconds{},
        sweep_replay_seconds{},replay_audit_seconds{};
    std::size_t peak_rss_bytes{},unchanged_fast_paths{},full_replay_fallbacks{},sweeps{},accepted_updates{},
        active_width_bounds{},candidate_events{};
    double objective{},ac_kkt{},width_kkt{};
    std::string eta_digest,beta_digest,candidate_decision_digest,stop_reason;
    bool qualified{},independent_replay_passed{};
};

double Seconds(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}

std::size_t PeakRssBytes()
{
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
}

std::string Digest(const Eigen::VectorXd & values)
{
    boost::hash2::sha2_256 hash;
    const auto add=[&hash](std::uint64_t value) {
        std::array<unsigned char,8> bytes{};
        for(unsigned k=0;k<8;++k) bytes[k]=static_cast<unsigned char>((value>>(8*k))&255);
        hash.update(bytes.data(),bytes.size());
    };
    add(static_cast<std::uint64_t>(values.size()));
    for(Eigen::Index i=0;i<values.size();++i) add(std::bit_cast<std::uint64_t>(values(i)));
    return boost::hash2::to_string(hash.result());
}

std::string Digest(const std::vector<std::uint64_t> & values)
{
    boost::hash2::sha2_256 hash;
    const auto add=[&hash](std::uint64_t value) {
        std::array<unsigned char,8> bytes{};
        for(unsigned k=0;k<8;++k) bytes[k]=static_cast<unsigned char>((value>>(8*k))&255);
        hash.update(bytes.data(),bytes.size());
    };
    add(static_cast<std::uint64_t>(values.size()));
    for(const auto value:values) add(value);
    return boost::hash2::to_string(hash.result());
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(),values.end());
    const auto middle=values.size()/2;
    return values.size()%2 ? values[middle] : .5*(values[middle-1]+values[middle]);
}

template<class Getter>
void PrintTime(const char * name,const std::vector<Sample> & samples,Getter getter)
{
    std::vector<double> values;
    for(const auto & sample:samples) values.push_back(getter(sample));
    const auto [minimum,maximum]=std::minmax_element(values.begin(),values.end());
    std::cout<<' '<<name<<"_median_s="<<Median(values)<<' '<<name<<"_range_s=["<<*minimum<<','<<*maximum<<']';
}

Sample Run(std::size_t targets,bool telemetry)
{
    joint::SetHotspotAttributionEnabled(telemetry);
    rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload(
        "chain",static_cast<int>(targets)));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const joint::Vector initial_eta=joint::Vector::Constant(
        static_cast<Eigen::Index>(targets),std::log(.55));
    joint::FixedNeighborPolicy policy;
    policy.core_atoms=12;
    policy.convergence_criteria=joint::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    policy.collect_telemetry=telemetry;
    policy.assess_local_conditioning=true;
    policy.assess_parameter_stability=true;
    std::vector<std::uint64_t> candidate_decisions;
    policy.candidate_observer=[&](const joint::FixedNeighborCandidateEvent & event) {
        candidate_decisions.push_back(static_cast<std::uint64_t>(event.sweep));
        candidate_decisions.push_back(static_cast<std::uint64_t>(event.block));
        candidate_decisions.push_back(event.accepted);
        candidate_decisions.push_back(event.unchanged);
        candidate_decisions.push_back(event.candidate_replay_ran);
    };

    const auto started=Clock::now();
    const auto result=joint::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    Sample sample;
    sample.wall_seconds=Seconds(started);
    if(!result.final_sweep || !result.certification)
        throw std::runtime_error("missing completed or certified sweep: "+result.reason);
    const auto & certification=*result.certification;
    const auto & independent=certification.assessment.independent_replay;
    if(!result.search_converged || !certification.qualified() || !independent.passed)
        throw std::runtime_error("qualification failed: "+result.reason);
    sample.replay=result.work.candidate_replay_breakdown;
    sample.sweeps=result.sweep_count;
    sample.accepted_updates=result.total_accepted_local_updates;
    sample.active_width_bounds=result.final_sweep->active_width_bound_count;
    sample.candidate_events=candidate_decisions.size()/5;
    sample.candidate_decision_digest=Digest(candidate_decisions);
    sample.profile=result.work.local_profile_work;
    sample.candidate_copy_seconds=result.work.candidate_copy_seconds;
    sample.candidate_replay_seconds=result.work.candidate_replay_seconds;
    sample.cache_update_seconds=result.work.cache_update_seconds;
    sample.sweep_replay_seconds=result.work.sweep_replay_seconds;
    sample.replay_audit_seconds=result.work.sweep_replay_audit_seconds;
    sample.unchanged_fast_paths=result.work.candidate_replay_skipped_unchanged;
    sample.full_replay_fallbacks=result.work.full_candidate_replays;
    sample.peak_rss_bytes=PeakRssBytes();
    sample.objective=independent.independent_objective/(data.context.scale*data.context.scale);
    sample.ac_kkt=result.final_sweep->global_ac_kkt;
    sample.width_kkt=result.final_sweep->projected_width_gradient_inf_norm;
    sample.eta_digest=Digest(result.state.eta);
    sample.beta_digest=Digest(result.state.beta);
    sample.stop_reason=result.reason;
    sample.qualified=certification.qualified();
    sample.independent_replay_passed=independent.passed;
    return sample;
}

void PrintHistogram(const char * label,const std::array<std::size_t,16> & histogram)
{
    std::cout<<' '<<label<<'=';
    for(std::size_t i=0;i<histogram.size();++i) std::cout<<(i ? "," : "")<<histogram[i];
}

void PrintSummary(std::size_t targets,bool telemetry,const std::vector<Sample> & samples)
{
    const auto & representative=samples[samples.size()/2];
    const auto & replay=representative.replay;
    const auto & work=representative.profile;
    const auto & linear=work.total.linear_solve;
    const auto & linear_detail=linear.breakdown;
    const auto & derivative=work.total.derivative_breakdown;
    const auto & lm=work.bounded_lm;
    const auto & lm_detail=lm.breakdown;
    const auto peak_rss=std::max_element(samples.begin(),samples.end(),[](const auto & a,const auto & b) {
        return a.peak_rss_bytes<b.peak_rss_bytes;
    })->peak_rss_bytes;
    std::cout<<std::setprecision(9)<<"targets="<<targets<<" telemetry="<<(telemetry ? "on" : "off")
        <<" warmups=1 measured_runs="<<samples.size();
    PrintTime("search_wall",samples,[](const auto & x){return x.wall_seconds;});
    PrintTime("profile_search",samples,[](const auto & x){return x.profile.total_seconds;});
    PrintTime("candidate_copy",samples,[](const auto & x){return x.candidate_copy_seconds;});
    PrintTime("candidate_replay",samples,[](const auto & x){return x.candidate_replay_seconds;});
    PrintTime("cache_update",samples,[](const auto & x){return x.cache_update_seconds;});
    PrintTime("sweep_full_replay",samples,[](const auto & x){return x.sweep_replay_seconds;});
    PrintTime("replay_audit",samples,[](const auto & x){return x.replay_audit_seconds;});
    std::cout<<" candidate_attempts="<<replay.attempts
        <<" affected_candidates="<<replay.affected_row_candidates
        <<" unchanged_fast_paths="<<representative.unchanged_fast_paths
        <<" full_replay_fallbacks="<<representative.full_replay_fallbacks
        <<" candidate_incidence_prepare_s="<<replay.incidence_preparation_seconds
        <<" candidate_incidence_ordering_s="<<replay.incidence_ordering_seconds
        <<" candidate_lookup_count_s="<<replay.incidence_lookup_count_seconds
        <<" candidate_incidence_construction_s="<<replay.incidence_construction_seconds
        <<" candidate_row_mapping_preparation_s="<<replay.row_mapping_preparation_seconds
        <<" candidate_prediction_residual_s="<<replay.prediction_residual_seconds
        <<" candidate_incremental_objective_s="<<replay.incremental_objective_seconds
        <<" candidate_scratch_peak_bytes="<<replay.peak_scratch_bytes
        <<" affected_rows_total="<<replay.affected_rows_total
        <<" affected_rows_min="<<(replay.affected_rows_min==std::numeric_limits<std::size_t>::max()?0:replay.affected_rows_min)
        <<" affected_rows_max="<<replay.affected_rows_max
        <<" contributor_entries_total="<<replay.contributor_entries_total
        <<" contributor_entries_min="<<(replay.contributor_entries_min==std::numeric_limits<std::size_t>::max()?0:replay.contributor_entries_min)
        <<" contributor_entries_max="<<replay.contributor_entries_max;
    PrintHistogram("affected_rows_log2_bins",replay.affected_row_count_bins);
    PrintHistogram("contributor_entries_log2_bins",replay.contributor_entry_count_bins);
    std::cout<<" profile_evaluations="<<work.total.evaluations
        <<" derivative_preparations="<<work.total.derivative_preparations
        <<" derivative_reductions="<<work.total.derivative_reductions
        <<" ac_solves="<<linear.solve_calls<<" ac_valid_solves="<<linear.valid_solves
        <<" ac_linear_solve_total_s="<<linear.solver_total_seconds
        <<" ac_sparse_qr_inclusive_s="<<linear.sparse_qr_seconds
        <<" ac_initial_projected_kkt_max="<<linear.initial_kkt_max
        <<" ac_final_projected_kkt_max="<<linear.final_kkt_max
        <<" ac_active_set_iterations="<<linear.active_set_iterations
        <<" ac_active_set_releases="<<linear.active_set_releases
        <<" sparse_structure_builds="<<linear.sparse_structure_build_count
        <<" sparse_qr_factorizations="<<linear_detail.sparse_qr_factorizations
        <<" free_face_changes="<<linear.free_face_changes
        <<" free_face_unchanged="<<linear.free_face_unchanged
        <<" free_face_repeat_transitions="<<linear_detail.free_face_repeats_across_builds
        <<" free_face_change_transitions="<<linear_detail.free_face_changes_across_builds
        <<" selected_pattern_repeat_transitions="<<linear_detail.selected_pattern_repeats
        <<" selected_pattern_change_transitions="<<linear_detail.selected_pattern_changes
        <<" reduced_pattern_repeat_transitions="<<linear_detail.reduced_pattern_repeats
        <<" reduced_pattern_change_transitions="<<linear_detail.reduced_pattern_changes
        <<" reduced_numeric_repeat_same_pattern="<<linear_detail.reduced_numeric_repeats_on_same_pattern
        <<" reduced_numeric_change_same_pattern="<<linear_detail.reduced_numeric_changes_on_same_pattern
        <<" selected_rows_range="<<(linear_detail.selected_matrix_builds?linear_detail.selected_rows_min:0)
        <<".."<<linear_detail.selected_rows_max
        <<" selected_columns_range="<<(linear_detail.selected_matrix_builds?linear_detail.selected_columns_min:0)
        <<".."<<linear_detail.selected_columns_max
        <<" selected_nonzeros_range="<<(linear_detail.selected_matrix_builds?linear_detail.selected_nonzeros_min:0)
        <<".."<<linear_detail.selected_nonzeros_max
        <<" reduced_rows_range="<<(linear_detail.selected_matrix_builds?linear_detail.reduced_rows_min:0)
        <<".."<<linear_detail.reduced_rows_max
        <<" reduced_nonzeros_range="<<(linear_detail.selected_matrix_builds?linear_detail.reduced_nonzeros_min:0)
        <<".."<<linear_detail.reduced_nonzeros_max
        <<" selected_nonzeros="<<linear.selected_nonzeros
        <<" reduced_nonzeros="<<linear.reduced_nonzeros
        <<" column_normalization_s="<<linear_detail.column_normalization_seconds
        <<" sparse_structure_build_s="<<linear.sparse_structure_build_seconds
        <<" row_reduction_s="<<linear_detail.row_reduction_seconds
        <<" row_major_conversion_s="<<linear.row_major_conversion_seconds
        <<" tile_mapping_build_s="<<linear.tile_mapping_build_seconds
        <<" tiled_sparse_row_qr_s="<<linear.tiled_qr_seconds
        <<" reduced_sparse_assembly_s="<<linear.reduced_sparse_assembly_seconds
        <<" sparse_qr_analysis_s="<<linear_detail.sparse_qr_analysis_seconds
        <<" sparse_qr_factorization_s="<<linear_detail.sparse_qr_factorization_seconds
        <<" sparse_qr_solve_s="<<linear_detail.sparse_qr_solve_seconds
        <<" pattern_fingerprint_s="<<linear_detail.pattern_fingerprint_seconds
        <<" reduced_pattern_fingerprint_s="<<linear_detail.reduced_pattern_fingerprint_seconds
        <<" refinement_eligible="<<linear.refinement_eligibility_count
        <<" refinement_triggered="<<linear.refinement_trigger_count
        <<" refinement_succeeded="<<linear.refinement_success_count
        <<" refinement_failed="<<linear.refinement_failure_count
        <<" refinement_iterations="<<linear.refinement_iterations
        <<" refinement_factorizations="<<linear.refinement_factorizations
        <<" refinement_s="<<linear.refinement_seconds
        <<" refinement_matrix_s="<<linear_detail.refinement_matrix_seconds
        <<" refinement_qr_s="<<linear_detail.refinement_qr_factorization_seconds
        <<" refinement_correction_solve_s="<<linear_detail.refinement_correction_solve_seconds
        <<" refinement_certification_s="<<linear_detail.refinement_certification_seconds
        <<" derivative_prepare_s="<<work.total.derivative_prepare_seconds
        <<" derivative_reduce_s="<<work.total.derivative_reduce_seconds
        <<" derivative_raw_construction_s="<<derivative.raw_construction_seconds
        <<" derivative_tiled_qr_s="<<derivative.tiled_qr_seconds
        <<" derivative_tiled_qr_calls="<<derivative.tiled_qr_calls
        <<" derivative_rank_svd_s="<<derivative.rank_svd_seconds
        <<" derivative_rank_svd_calls="<<derivative.rank_svd_calls
        <<" derivative_triangular_solves_s="<<derivative.triangular_solve_seconds
        <<" derivative_triangular_solve_calls="<<derivative.triangular_solve_calls
        <<" derivative_correction_matrix_s="<<derivative.correction_matrix_seconds
        <<" derivative_cancellation_fallbacks="<<derivative.cancellation_fallbacks
        <<" derivative_cancellation_fallback_s="<<derivative.cancellation_fallback_seconds
        <<" derivative_workspace_resize_events="<<derivative.workspace_resize_events
        <<" lm_subproblems="<<lm.subproblem_solves<<" lm_failed_subproblems="<<lm.failed_subproblems
        <<" lm_zero_damping_solves="<<lm.zero_damping_solves
        <<" lm_damping_solves="<<lm.damping_solves<<" lm_bracket_solves="<<lm.bracket_solves
        <<" lm_bisection_solves="<<lm.bisection_solves
        <<" lm_qr_factorizations="<<lm.damping_qr_factorizations
        <<" lm_active_set_iterations="<<lm.damping_active_set_iterations
        <<" lm_releases="<<lm.damping_releases
        <<" lm_active_selection_s="<<lm_detail.active_free_selection_seconds
        <<" lm_augmented_assembly_s="<<lm_detail.augmented_matrix_assembly_seconds
        <<" lm_dense_qr_s="<<lm_detail.dense_qr_factorization_seconds
        <<" lm_qr_solve_s="<<lm_detail.qr_solve_seconds
        <<" lm_kkt_validation_s="<<lm_detail.kkt_validation_seconds
        <<" lm_damping_search_s="<<lm_detail.damping_search_seconds
        <<" lm_updates_per_log2_bin=";
    for(std::size_t i=0;i<lm.damping_solves_per_update_bins.size();++i)
        std::cout<<(i ? "," : "")<<lm.damping_solves_per_update_bins[i];
    std::cout<<" lm_qr_per_update_log2_bin=";
    for(std::size_t i=0;i<lm.qr_factorizations_per_update_bins.size();++i)
        std::cout<<(i ? "," : "")<<lm.qr_factorizations_per_update_bins[i];
    std::cout<<" peak_rss_bytes="<<peak_rss
        <<" objective="<<representative.objective<<" ac_kkt="<<representative.ac_kkt
        <<" width_kkt="<<representative.width_kkt<<" stop_reason="<<representative.stop_reason
        <<" sweeps="<<representative.sweeps<<" accepted_local_updates="<<representative.accepted_updates
        <<" active_width_bounds="<<representative.active_width_bounds
        <<" candidate_events="<<representative.candidate_events
        <<" candidate_decisions_sha256="<<representative.candidate_decision_digest
        <<" eta_sha256="<<representative.eta_digest<<" beta_sha256="<<representative.beta_digest
        <<" qualified="<<representative.qualified
        <<" independent_replay="<<representative.independent_replay_passed<<'\n';
}

void PrintCancellationSummary()
{
    std::vector<std::vector<joint::Support>> support(1);
    support[0]={{0,0.0},{1,0.25},{2,0.25}};
    joint::Vector y(3);
    for(Eigen::Index row=0;row<3;++row)
    {
        const double square=row==0 ? 0.0 : 0.25;
        const auto basis=joint::EvaluateKernel(square,.55,2.5);
        y(row)=3.0*basis.gaussian+.25*basis.charge+(row==1 ? 1e-4 : row==2 ? -1e-4 : 0.0);
    }
    const joint::Domain domain(3,std::move(support));
    const auto context=second_stage_test::matched::joint_abc::MakeContext(y,1);
    const auto eta=joint::Vector::Constant(1,std::log(.55));
    const auto evaluation=joint::EvaluateProfile(domain,y,eta,&context);
    if(!evaluation.valid) throw std::runtime_error("cancellation fixture profile failed: "+evaluation.reason);
    std::vector<double> total_seconds,fallback_seconds;
    std::size_t fallback_count{};
    for(std::size_t run=0;run<6;++run)
    {
        const auto started=Clock::now();
        const auto differential=joint::PrepareDerivative(evaluation,context.scale,&context,-1,3);
        const auto elapsed=Seconds(started);
        if(!differential.valid || !differential.reference_order || differential.breakdown.cancellation_fallbacks!=1)
            throw std::runtime_error("cancellation fixture did not use one valid reference-order fallback");
        fallback_count=differential.breakdown.cancellation_fallbacks;
        if(run)
        {
            total_seconds.push_back(elapsed);
            fallback_seconds.push_back(differential.breakdown.cancellation_fallback_seconds);
        }
    }
    const auto [minimum,maximum]=std::minmax_element(fallback_seconds.begin(),fallback_seconds.end());
    std::cout<<std::setprecision(9)<<"cancellation_fixture warmups=1 measured_runs=5 fallbacks="<<fallback_count
        <<" total_median_s="<<Median(total_seconds)<<" fallback_extra_median_s="<<Median(fallback_seconds)
        <<" fallback_extra_range_s=["<<*minimum<<','<<*maximum<<"] reference_order=1\n";
}
}

int main(int argc,char ** argv)
{
    if(argc!=3)
    {
        std::cerr<<"usage: rhbm_joint_hotspot_attribution_benchmark <target-count> <on|off>\n";
        return 2;
    }
    try
    {
        const auto targets=static_cast<std::size_t>(std::stoul(argv[1]));
        if(targets==0 || targets>1000000) throw std::invalid_argument("target count is out of range");
        const std::string mode=argv[2];
        if(mode!="on" && mode!="off") throw std::invalid_argument("telemetry mode must be on or off");
        const bool telemetry=mode=="on";
        (void)Run(targets,telemetry);
        std::vector<Sample> samples;
        samples.reserve(5);
        for(std::size_t repetition=0;repetition<5;++repetition) samples.push_back(Run(targets,telemetry));
        PrintSummary(targets,telemetry,samples);
        if(telemetry) PrintCancellationSummary();
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<"joint hotspot attribution benchmark failed: "<<error.what()<<'\n';
        return 1;
    }
}
