#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;
struct Sample
{
    std::string stop_reason,mode,split_status,split_reason,null_reference_status;
    std::optional<double> max_reoptimization_prediction_drift,max_reoptimization_log_width_difference,
        max_reoptimization_scaled_ac_difference,max_seed_prediction_robustness,
        max_seed_log_width_difference,max_seed_scaled_ac_difference,max_high_budget_width_difference,
        max_high_budget_scaled_ac_difference,max_high_budget_prediction_difference;
    double search_seconds{},wall_seconds{},cpu_seconds{},parameter_stability_seconds{},
        parameter_stability_cpu_seconds{},null_search_seconds{},null_search_cpu_seconds{},
        candidate_seconds{},sweep_replay_seconds{},
        replay_audit_seconds{},reconciliation_seconds{},global_state_seconds{},
        endpoint_certification_seconds{},independent_replay_seconds{},local_profile_seconds{},
        profile_basis_seconds{},profile_linear_solve_seconds{},profile_derivative_prepare_seconds{},
        profile_derivative_reduce_seconds{};
    std::size_t affected_candidates{},full_replay_fallbacks{},sweep_replays{},
        objective_reconciliations{},sweeps{},accepted_blocks{},unchanged_blocks{},failed_blocks{},
        first_kkt_sweep{},first_kkt_eligible_sweep{},first_confirmed_sweep{},
        certification_attempts{},block_solves{},
        preflight_bytes{},peak_rss_bytes{},profile_evaluations{},profile_derivatives{},
        parameter_stability_solves{},parameter_stability_evaluations{},
        parameter_stability_blocks_selected{},parameter_stability_peak_bytes{},
        null_profile_evaluations{},null_accepted_updates{},null_replay_checks{},
        split_qualified_blocks{},split_incomplete_blocks{},split_workspace_peak_bytes{},
        bounded_subproblems{},damping_solves{},damping_bracket_solves{},damping_bisection_solves{},
        damping_qr_factorizations{},damping_active_set_iterations{},damping_releases{};
    std::uint64_t accepted_state_fingerprint{};
    double objective{},global_ac_kkt{},projected_width_kkt{},maximum_cache_replay_error{};
    bool qualified{},independent_replay_passed{},memory_preflight_passed{},full_replay_memory_fallback{};
};

std::uint64_t HashWord(std::uint64_t hash,std::uint64_t word)
{
    for(unsigned shift=0;shift<64;shift+=8)
    {
        hash^=(word>>shift)&0xffu;
        hash*=1099511628211ull;
    }
    return hash;
}

std::uint64_t AcceptedStateFingerprint(const joint::BlockCoordinateState & state)
{
    std::uint64_t hash=14695981039346656037ull;
    for(const auto * values:{&state.eta,&state.beta,&state.prediction,&state.residual})
    {
        hash=HashWord(hash,static_cast<std::uint64_t>(values->size()));
        for(Eigen::Index index=0;index<values->size();++index)
            hash=HashWord(hash,std::bit_cast<std::uint64_t>((*values)(index)));
    }
    return HashWord(hash,std::bit_cast<std::uint64_t>(state.objective));
}

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

double Median(std::vector<double> values)
{
    std::sort(values.begin(),values.end());
    const auto middle=values.size()/2;
    return values.size()%2 ? values[middle] : .5*(values[middle-1]+values[middle]);
}

void KeepMaximum(std::optional<double> & maximum,const std::optional<double> & value)
{
    if(value && (!maximum || *value>*maximum)) maximum=value;
}

template<class Getter>
void PrintTime(const char * name,const std::vector<Sample> & samples,Getter getter)
{
    std::vector<double> values;
    values.reserve(samples.size());
    for(const auto & sample:samples) values.push_back(getter(sample));
    const auto [minimum,maximum]=std::minmax_element(values.begin(),values.end());
    std::cout<<' '<<name<<"_median_s="<<Median(values)
        <<' '<<name<<"_range_s=["<<*minimum<<','<<*maximum<<']';
}

Sample Run(std::size_t targets,joint::SplitDiagnosticsMode split_mode)
{
    rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload(
        "chain",static_cast<int>(targets)));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const joint::Vector initial_eta=joint::Vector::Constant(
        static_cast<Eigen::Index>(targets),std::log(.55));
    joint::FixedNeighborPolicy policy;
    policy.core_atoms=12;
    policy.convergence_criteria=joint::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    std::size_t remaining_perturbation_blocks=joint::MaximumPerturbedBlocksPerFit;
    std::size_t remaining_general_perturbation_blocks=
        joint::MaximumPerturbedBlocksPerFit>1 ? joint::MaximumPerturbedBlocksPerFit-1 : 0;
    bool reserved_bound_active_probe=joint::MaximumPerturbedBlocksPerFit>0;
    policy.claim_parameter_stability_block=[&](int priority) {
        if(!remaining_perturbation_blocks) return false;
        if(priority==0 && reserved_bound_active_probe)
        {
            reserved_bound_active_probe=false;
            --remaining_perturbation_blocks;
            return true;
        }
        if(!remaining_general_perturbation_blocks) return false;
        --remaining_general_perturbation_blocks;
        --remaining_perturbation_blocks;
        return true;
    };
    policy.assess_local_conditioning=true;
    policy.assess_parameter_stability=true;
    policy.split_diagnostics_mode=split_mode;
    std::size_t accepted_blocks{},unchanged_blocks{},failed_blocks{};
    std::size_t first_kkt_sweep{},first_kkt_eligible_sweep{};
    const auto & criteria=policy.convergence_criteria;
    policy.production_progress=[&](const joint::FixedNeighborBlockSweep & sweep) {
        accepted_blocks+=sweep.accepted_blocks;
        unchanged_blocks+=sweep.unchanged_blocks;
        failed_blocks+=sweep.block_solves-sweep.accepted_blocks-sweep.unchanged_blocks;
        const double width_kkt=sweep.width_bound_applied ? sweep.projected_width_gradient_inf_norm :
            sweep.global_width_gradient_inf_norm;
        const bool kkt_passed=std::isfinite(sweep.global_ac_kkt) &&
            sweep.global_ac_kkt<=criteria.global_ac_kkt && std::isfinite(width_kkt) &&
            width_kkt<=criteria.global_width_gradient;
        if(kkt_passed)
        {
            if(first_kkt_sweep==0)
                first_kkt_sweep=sweep.sweep;
            const auto failed_in_sweep=sweep.block_solves-sweep.accepted_blocks-sweep.unchanged_blocks;
            if(first_kkt_eligible_sweep==0 && sweep.coordinate_confirmation_available && failed_in_sweep==0)
                first_kkt_eligible_sweep=sweep.sweep;
        }
    };

    const auto started=Clock::now();
    const auto cpu_started=std::clock();
    const auto result=joint::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    const double wall_seconds=std::chrono::duration<double>(Clock::now()-started).count();
    if(!result.final_sweep) throw std::runtime_error("Production search returned no completed sweep: "+result.reason);
    if(!result.certification) throw std::runtime_error("Production search did not certify its endpoint: "+result.reason);
    const auto & certification=*result.certification;
    const auto & independent=certification.assessment.independent_replay;
    if(!result.search_converged || !certification.qualified() || !independent.passed)
        throw std::runtime_error("Production search qualification failed: "+result.reason);

    Sample sample;
    sample.mode=split_mode==joint::SplitDiagnosticsMode::Disabled ? "B0" :
        split_mode==joint::SplitDiagnosticsMode::Minimal ? "B1" : "B2";
    sample.stop_reason=result.reason;
    sample.search_seconds=result.search_seconds;
    sample.wall_seconds=wall_seconds;
    sample.cpu_seconds=static_cast<double>(std::clock()-cpu_started)/static_cast<double>(CLOCKS_PER_SEC);
    sample.accepted_state_fingerprint=AcceptedStateFingerprint(result.state);
    sample.parameter_stability_seconds=result.work.parameter_stability_seconds;
    sample.parameter_stability_solves=result.work.perturbation_solves;
    sample.parameter_stability_evaluations=result.work.perturbation_evaluations;
    sample.parameter_stability_peak_bytes=result.work.parameter_stability_peak_memory_estimate_bytes;
    if(result.parameter_stability)
    {
        sample.parameter_stability_cpu_seconds=result.parameter_stability->cpu_seconds;
        for(const auto & block:result.parameter_stability->blocks)
            sample.parameter_stability_blocks_selected+=block.probe_selected ? 1u : 0u;
    }
    for(const auto & block:result.split_diagnostics)
    {
        const auto & split=block.diagnostics;
        sample.null_search_seconds+=split.null_search_seconds;
        sample.null_search_cpu_seconds+=split.null_search_cpu_seconds;
        sample.null_profile_evaluations+=split.null_profile_evaluations;
        sample.null_accepted_updates+=split.null_accepted_updates;
        sample.null_replay_checks+=split.null_replay_checks;
        sample.split_workspace_peak_bytes=std::max(sample.split_workspace_peak_bytes,
            split.workspace_estimate_bytes);
        KeepMaximum(sample.max_reoptimization_prediction_drift,split.reoptimization_prediction_drift);
        KeepMaximum(sample.max_reoptimization_log_width_difference,split.reoptimization_log_width_difference_inf);
        KeepMaximum(sample.max_reoptimization_scaled_ac_difference,split.reoptimization_scaled_ac_difference);
        KeepMaximum(sample.max_seed_prediction_robustness,split.seed_prediction_robustness);
        KeepMaximum(sample.max_seed_log_width_difference,split.seed_log_width_difference_inf);
        KeepMaximum(sample.max_seed_scaled_ac_difference,split.seed_scaled_ac_difference);
        KeepMaximum(sample.max_high_budget_width_difference,split.high_budget_width_difference_inf);
        KeepMaximum(sample.max_high_budget_scaled_ac_difference,split.high_budget_scaled_ac_difference);
        KeepMaximum(sample.max_high_budget_prediction_difference,split.high_budget_prediction_difference_relative);
        sample.split_qualified_blocks+=split.status=="Qualified" ? 1u : 0u;
        sample.split_incomplete_blocks+=split.status=="Incomplete" || split.status=="Unavailable" ? 1u : 0u;
        if(sample.split_status.empty())
        {
            sample.split_status=split.status;
            sample.split_reason=split.reason;
            sample.null_reference_status=split.null_reference_status;
        }
    }
    const auto & profile_work=result.work.local_profile_work;
    sample.local_profile_seconds=profile_work.total_seconds;
    sample.profile_basis_seconds=profile_work.total.basis_seconds;
    sample.profile_linear_solve_seconds=profile_work.total.linear_rhs_solve_seconds;
    sample.profile_derivative_prepare_seconds=profile_work.total.derivative_prepare_seconds;
    sample.profile_derivative_reduce_seconds=profile_work.total.derivative_reduce_seconds;
    sample.profile_evaluations=profile_work.total.evaluations;
    sample.profile_derivatives=profile_work.total.derivative_preparations;
    sample.bounded_subproblems=profile_work.bounded_lm.subproblem_solves;
    sample.damping_solves=profile_work.bounded_lm.damping_solves;
    sample.damping_bracket_solves=profile_work.bounded_lm.bracket_solves;
    sample.damping_bisection_solves=profile_work.bounded_lm.bisection_solves;
    sample.damping_qr_factorizations=profile_work.bounded_lm.damping_qr_factorizations;
    sample.damping_active_set_iterations=profile_work.bounded_lm.damping_active_set_iterations;
    sample.damping_releases=profile_work.bounded_lm.damping_releases;
    sample.candidate_seconds=result.work.candidate_replay_seconds;
    sample.sweep_replay_seconds=result.work.sweep_replay_seconds;
    sample.replay_audit_seconds=result.work.sweep_replay_audit_seconds;
    sample.reconciliation_seconds=result.work.objective_reconciliation_seconds;
    sample.global_state_seconds=result.work.sweep_global_state_seconds;
    sample.independent_replay_seconds=result.work.independent_replay_seconds;
    sample.endpoint_certification_seconds=std::max(0.0,
        result.work.endpoint_verification_seconds-sample.independent_replay_seconds);
    sample.affected_candidates=result.work.affected_candidate_replays;
    sample.full_replay_fallbacks=result.work.full_candidate_replays;
    sample.sweep_replays=result.work.sweep_replays;
    sample.objective_reconciliations=result.work.objective_reconciliations;
    sample.preflight_bytes=result.global_state_memory_estimate.estimated_peak_bytes;
    sample.peak_rss_bytes=PeakRssBytes();
    sample.accepted_blocks=accepted_blocks;
    sample.unchanged_blocks=unchanged_blocks;
    sample.failed_blocks=failed_blocks;
    sample.objective=independent.independent_objective/(data.context.scale*data.context.scale);
    sample.global_ac_kkt=result.final_sweep->global_ac_kkt;
    sample.projected_width_kkt=result.final_sweep->projected_width_gradient_inf_norm;
    sample.maximum_cache_replay_error=result.work.maximum_cache_replay_error;
    sample.sweeps=result.sweep_count;
    sample.first_kkt_sweep=first_kkt_sweep;
    sample.first_kkt_eligible_sweep=first_kkt_eligible_sweep;
    sample.first_confirmed_sweep=result.confirmed_stationarity_sweep;
    sample.certification_attempts=result.certification_attempts;
    sample.block_solves=result.total_block_solves;
    sample.qualified=certification.qualified();
    sample.independent_replay_passed=independent.passed;
    sample.memory_preflight_passed=result.global_state_memory_estimate.status==joint::MemoryStatus::Pass;
    sample.full_replay_memory_fallback=result.work.full_replay_memory_fallback;
    return sample;
}

void PrintSummary(std::size_t targets,const std::string & mode,const std::vector<Sample> & samples)
{
    std::cout<<std::setprecision(9)<<"targets="<<targets<<" mode="<<mode
        <<" warmups=1 measured_runs="<<samples.size();
    PrintTime("search",samples,[](const auto & x){return x.search_seconds;});
    PrintTime("wall",samples,[](const auto & x){return x.wall_seconds;});
    PrintTime("cpu",samples,[](const auto & x){return x.cpu_seconds;});
    PrintTime("parameter_stability_wall",samples,[](const auto & x){return x.parameter_stability_seconds;});
    PrintTime("parameter_stability_cpu",samples,[](const auto & x){return x.parameter_stability_cpu_seconds;});
    PrintTime("null_search_wall",samples,[](const auto & x){return x.null_search_seconds;});
    PrintTime("null_search_cpu",samples,[](const auto & x){return x.null_search_cpu_seconds;});
    PrintTime("local_profile_total",samples,[](const auto & x){return x.local_profile_seconds;});
    PrintTime("local_profile_basis",samples,[](const auto & x){return x.profile_basis_seconds;});
    PrintTime("local_profile_linear_solve",samples,[](const auto & x){return x.profile_linear_solve_seconds;});
    PrintTime("local_profile_derivative_prepare",samples,[](const auto & x){return x.profile_derivative_prepare_seconds;});
    PrintTime("local_profile_derivative_reduce",samples,[](const auto & x){return x.profile_derivative_reduce_seconds;});
    PrintTime("candidate_evaluation",samples,[](const auto & x){return x.candidate_seconds;});
    PrintTime("sweep_replay",samples,[](const auto & x){return x.sweep_replay_seconds;});
    PrintTime("replay_audit",samples,[](const auto & x){return x.replay_audit_seconds;});
    PrintTime("objective_reconciliation",samples,[](const auto & x){return x.reconciliation_seconds;});
    PrintTime("global_state_evaluation",samples,[](const auto & x){return x.global_state_seconds;});
    PrintTime("endpoint_certification_excluding_independent_replay",samples,
        [](const auto & x){return x.endpoint_certification_seconds;});
    PrintTime("independent_replay",samples,[](const auto & x){return x.independent_replay_seconds;});
    const auto & representative=samples[samples.size()/2];
    const auto peak_rss=std::max_element(samples.begin(),samples.end(),[](const auto & left,const auto & right) {
        return left.peak_rss_bytes<right.peak_rss_bytes;
    })->peak_rss_bytes;
    const bool state_fingerprint_equal=std::all_of(samples.begin(),samples.end(),[&](const auto & sample) {
        return sample.accepted_state_fingerprint==representative.accepted_state_fingerprint;
    });
    std::cout<<" affected_candidate_evaluations="<<representative.affected_candidates
        <<" profile_evaluations="<<representative.profile_evaluations
        <<" profile_derivatives="<<representative.profile_derivatives
        <<" bounded_subproblems="<<representative.bounded_subproblems
        <<" damping_solves="<<representative.damping_solves
        <<" damping_bracket_solves="<<representative.damping_bracket_solves
        <<" damping_bisection_solves="<<representative.damping_bisection_solves
        <<" damping_qr_factorizations="<<representative.damping_qr_factorizations
        <<" damping_active_set_iterations="<<representative.damping_active_set_iterations
        <<" damping_releases="<<representative.damping_releases
        <<" parameter_stability_solves="<<representative.parameter_stability_solves
        <<" parameter_stability_evaluations="<<representative.parameter_stability_evaluations
        <<" parameter_stability_blocks_selected="<<representative.parameter_stability_blocks_selected
        <<" parameter_stability_peak_bytes="<<representative.parameter_stability_peak_bytes
        <<" null_profile_evaluations="<<representative.null_profile_evaluations
        <<" null_accepted_updates="<<representative.null_accepted_updates
        <<" null_replay_checks="<<representative.null_replay_checks
        <<" split_qualified_blocks="<<representative.split_qualified_blocks
        <<" split_incomplete_blocks="<<representative.split_incomplete_blocks
        <<" split_workspace_peak_bytes="<<representative.split_workspace_peak_bytes
        <<" split_status="<<representative.split_status
        <<" split_reason="<<representative.split_reason
        <<" null_reference_status="<<representative.null_reference_status
        <<" max_reoptimization_prediction_drift="<<representative.max_reoptimization_prediction_drift.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_reoptimization_log_width_difference="<<representative.max_reoptimization_log_width_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_reoptimization_scaled_ac_difference="<<representative.max_reoptimization_scaled_ac_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_seed_prediction_robustness="<<representative.max_seed_prediction_robustness.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_seed_log_width_difference="<<representative.max_seed_log_width_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_seed_scaled_ac_difference="<<representative.max_seed_scaled_ac_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_high_budget_width_difference="<<representative.max_high_budget_width_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_high_budget_scaled_ac_difference="<<representative.max_high_budget_scaled_ac_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" max_high_budget_prediction_difference="<<representative.max_high_budget_prediction_difference.value_or(
            std::numeric_limits<double>::quiet_NaN())
        <<" accepted_state_fnv1a64="<<std::hex<<std::setw(16)<<std::setfill('0')
        <<representative.accepted_state_fingerprint<<std::dec<<std::setfill(' ')
        <<" accepted_state_fingerprint_equal="<<state_fingerprint_equal
        <<" full_replay_fallbacks="<<representative.full_replay_fallbacks
        <<" sweep_replays="<<representative.sweep_replays
        <<" objective_reconciliations="<<representative.objective_reconciliations
        <<" sweeps="<<representative.sweeps
        <<" accepted_blocks="<<representative.accepted_blocks
        <<" unchanged_blocks="<<representative.unchanged_blocks
        <<" failed_blocks="<<representative.failed_blocks
        <<" memory_full_replay_fallback="<<representative.full_replay_memory_fallback
        <<" memory_preflight_peak_bytes="<<representative.preflight_bytes
        <<" memory_preflight_pass="<<representative.memory_preflight_passed
        <<" peak_rss_bytes="<<peak_rss
        <<" independent_objective="<<representative.objective
        <<" global_ac_kkt="<<representative.global_ac_kkt
        <<" projected_width_kkt="<<representative.projected_width_kkt
        <<" maximum_cache_replay_error="<<representative.maximum_cache_replay_error
        <<" initial_width=0.55"
        <<" first_kkt_sweep="<<representative.first_kkt_sweep
        <<" first_kkt_eligible_sweep="<<representative.first_kkt_eligible_sweep
        <<" first_confirmed_sweep="<<representative.first_confirmed_sweep
        <<" certification_attempts="<<representative.certification_attempts
        <<" block_solves="<<representative.block_solves
        <<" stop_reason="<<representative.stop_reason
        <<" endpoint_qualified="<<representative.qualified
        <<" independent_replay_pass="<<representative.independent_replay_passed<<'\n';
}
}

int main(int argc,char ** argv)
{
    if(argc!=3)
    {
        std::cerr<<"usage: rhbm_fixed_neighbor_closed_loop_benchmark <target-count> <B0|B1|B2>\n";
        return 2;
    }
    try
    {
        const auto targets=static_cast<std::size_t>(std::stoul(argv[1]));
        if(targets==0 || targets>1000000) throw std::invalid_argument("target count is out of range");
        const std::string mode=argv[2];
        const auto split_mode=mode=="B0" ? joint::SplitDiagnosticsMode::Disabled :
            mode=="B1" ? joint::SplitDiagnosticsMode::Minimal :
            mode=="B2" ? joint::SplitDiagnosticsMode::ReuseWorkspace :
            throw std::invalid_argument("mode must be B0, B1, or B2");
        (void)Run(targets,split_mode);
        std::vector<Sample> samples;
        samples.reserve(5);
        for(std::size_t repetition=0;repetition<5;++repetition) samples.push_back(Run(targets,split_mode));
        PrintSummary(targets,mode,samples);
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<"fixed-neighbor production benchmark failed: "<<error.what()<<'\n';
        return 1;
    }
}
