#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;

struct Sample
{
    std::string stop_reason;
    double search_seconds{},wall_seconds{},candidate_seconds{},sweep_replay_seconds{},
        replay_audit_seconds{},reconciliation_seconds{},global_state_seconds{},
        endpoint_certification_seconds{},independent_replay_seconds{};
    std::size_t affected_candidates{},full_replay_fallbacks{},sweep_replays{},
        objective_reconciliations{},sweeps{},preflight_bytes{},peak_rss_bytes{};
    double objective{},global_ac_kkt{},projected_width_kkt{},maximum_cache_replay_error{},
        mean_amplitude{},mean_charge{},mean_width{};
    bool qualified{},independent_replay_passed{},memory_preflight_passed{},full_replay_memory_fallback{};
};

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

Sample Run(std::size_t targets)
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

    const auto started=Clock::now();
    const auto result=joint::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    const double wall_seconds=std::chrono::duration<double>(Clock::now()-started).count();
    if(!result.final_sweep) throw std::runtime_error("Production E1 returned no completed sweep: "+result.reason);
    if(!result.certification) throw std::runtime_error("Production E1 did not certify its endpoint: "+result.reason);
    const auto & certification=*result.certification;
    const auto & independent=certification.assessment.independent_replay;
    if(!result.search_converged || !certification.qualified() || !independent.passed)
        throw std::runtime_error("Production E1 qualification failed: "+result.reason);

    Sample sample;
    sample.stop_reason=result.reason;
    sample.search_seconds=result.search_seconds;
    sample.wall_seconds=wall_seconds;
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
    sample.objective=independent.independent_objective/(data.context.scale*data.context.scale);
    sample.global_ac_kkt=result.final_sweep->global_ac_kkt;
    sample.projected_width_kkt=result.final_sweep->projected_width_gradient_inf_norm;
    sample.maximum_cache_replay_error=result.work.maximum_cache_replay_error;
    sample.sweeps=result.sweep_count;
    const auto atom_count=result.state.eta.size();
    for(Eigen::Index atom=0;atom<atom_count;++atom)
    {
        sample.mean_amplitude+=result.state.beta(2*atom);
        sample.mean_charge+=result.state.beta(2*atom+1);
        sample.mean_width+=std::exp(result.state.eta(atom));
    }
    const double inverse_atom_count=1.0/static_cast<double>(atom_count);
    sample.mean_amplitude*=inverse_atom_count;
    sample.mean_charge*=inverse_atom_count;
    sample.mean_width*=inverse_atom_count;
    sample.qualified=certification.qualified();
    sample.independent_replay_passed=independent.passed;
    sample.memory_preflight_passed=result.global_state_memory_estimate.status==joint::MemoryStatus::Pass;
    sample.full_replay_memory_fallback=result.work.full_replay_memory_fallback;
    return sample;
}

void PrintSummary(std::size_t targets,const std::vector<Sample> & samples)
{
    std::cout<<std::setprecision(9)<<"targets="<<targets<<" mode=production-e1"
        <<" warmups=1 measured_runs="<<samples.size();
    PrintTime("search",samples,[](const auto & x){return x.search_seconds;});
    PrintTime("wall",samples,[](const auto & x){return x.wall_seconds;});
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
    std::cout<<" affected_candidate_evaluations="<<representative.affected_candidates
        <<" full_replay_fallbacks="<<representative.full_replay_fallbacks
        <<" sweep_replays="<<representative.sweep_replays
        <<" objective_reconciliations="<<representative.objective_reconciliations
        <<" sweeps="<<representative.sweeps
        <<" memory_full_replay_fallback="<<representative.full_replay_memory_fallback
        <<" memory_preflight_peak_bytes="<<representative.preflight_bytes
        <<" memory_preflight_pass="<<representative.memory_preflight_passed
        <<" peak_rss_bytes="<<peak_rss
        <<" independent_objective="<<representative.objective
        <<" global_ac_kkt="<<representative.global_ac_kkt
        <<" projected_width_kkt="<<representative.projected_width_kkt
        <<" maximum_cache_replay_error="<<representative.maximum_cache_replay_error
        <<" mean_amplitude="<<representative.mean_amplitude
        <<" mean_charge="<<representative.mean_charge
        <<" mean_width="<<representative.mean_width
        <<" stop_reason="<<representative.stop_reason
        <<" endpoint_qualified="<<representative.qualified
        <<" independent_replay_pass="<<representative.independent_replay_passed<<'\n';
}
}

int main(int argc,char ** argv)
{
    if(argc!=2)
    {
        std::cerr<<"usage: rhbm_fixed_neighbor_closed_loop_benchmark <target-count>\n";
        return 2;
    }
    try
    {
        const auto targets=static_cast<std::size_t>(std::stoul(argv[1]));
        if(targets==0 || targets>1000000) throw std::invalid_argument("target count is out of range");
        (void)Run(targets);
        std::vector<Sample> samples;
        samples.reserve(5);
        for(std::size_t repetition=0;repetition<5;++repetition) samples.push_back(Run(targets));
        PrintSummary(targets,samples);
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<"fixed-neighbor production benchmark failed: "<<error.what()<<'\n';
        return 1;
    }
}
