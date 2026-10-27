#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/JointReplayExperiment.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>
#include <utility>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;
using Vector=Eigen::VectorXd;

struct Sample
{
    double wall{},search{},structural_preparation{},candidate_replay{},local_search{},sweep_replay{},global_state{};
    std::size_t candidate_count{},accepted_count{},unchanged_count{},peak_estimate{};
};
struct Summary
{
    std::vector<double> wall,search,structural_preparation,candidate_replay,local_search,sweep_replay,global_state;
    std::size_t candidate_count{},accepted_count{},unchanged_count{},peak_estimate{},candidate_scratch{};
};

double Seconds(Clock::time_point started)
{return std::chrono::duration<double>(Clock::now()-started).count();}
bool BitwiseEqual(double left,double right)
{return std::bit_cast<std::uint64_t>(left)==std::bit_cast<std::uint64_t>(right);}
bool BitwiseEqual(const Vector & left,const Vector & right)
{
    if(left.size()!=right.size()) return false;
    for(Eigen::Index i=0;i<left.size();++i) if(!BitwiseEqual(left(i),right(i))) return false;
    return true;
}
double Median(std::vector<double> values)
{
    std::sort(values.begin(),values.end()); return values[values.size()/2];
}
std::pair<double,double> Range(std::vector<double> values)
{
    std::sort(values.begin(),values.end()); return {values.front(),values.back()};
}
void VerifyEquivalent(const joint::FixedNeighborDiagnosticResult & full,
    const joint::FixedNeighborDiagnosticResult & exact)
{
    if(full.reason!=exact.reason || full.search_converged!=exact.search_converged ||
        full.sweep_count!=exact.sweep_count || full.blocks.size()!=exact.blocks.size() ||
        full.sweeps.size()!=exact.sweeps.size() || !BitwiseEqual(full.state.eta,exact.state.eta) ||
        !BitwiseEqual(full.state.beta,exact.state.beta) || !BitwiseEqual(full.state.prediction,exact.state.prediction) ||
        !BitwiseEqual(full.state.residual,exact.state.residual) || !BitwiseEqual(full.state.objective,exact.state.objective))
        throw std::runtime_error("Closed-loop D1 changed final Full Replay search result.");
    for(std::size_t k=0;k<full.blocks.size();++k)
    {
        const auto & a=full.blocks[k]; const auto & b=exact.blocks[k];
        if(a.sweep!=b.sweep || a.block!=b.block || a.core_atoms!=b.core_atoms || a.status!=b.status ||
            a.reason!=b.reason || a.accepted!=b.accepted || a.accepted_updates!=b.accepted_updates ||
            !BitwiseEqual(a.local_objective_before,b.local_objective_before) ||
            !BitwiseEqual(a.local_objective_after,b.local_objective_after) ||
            !BitwiseEqual(a.global_replay_delta,b.global_replay_delta) ||
            !BitwiseEqual(a.local_global_delta_error,b.local_global_delta_error))
            throw std::runtime_error("Closed-loop D1 changed a block result.");
    }
    for(std::size_t k=0;k<full.sweeps.size();++k)
    {
        const auto & a=full.sweeps[k]; const auto & b=exact.sweeps[k];
        if(a.sweep!=b.sweep || a.accepted_blocks!=b.accepted_blocks || a.unchanged_blocks!=b.unchanged_blocks ||
            !BitwiseEqual(a.objective_after,b.objective_after) || !BitwiseEqual(a.global_ac_kkt,b.global_ac_kkt) ||
            !BitwiseEqual(a.projected_width_gradient_inf_norm,b.projected_width_gradient_inf_norm))
            throw std::runtime_error("Closed-loop D1 changed a sweep result.");
    }
}
Sample RunSearch(const rhbm_gem::core::JointProblem & problem,joint::VectorRef eta,
    const joint::FixedNeighborPolicy & policy,std::size_t candidate_count,
    joint::FixedNeighborDiagnosticResult & result)
{
    const auto started=Clock::now();
    result=joint::RunFixedNeighborDiagnostic(problem,eta,policy);
    const double wall=Seconds(started);
    if(result.reason.empty()) throw std::runtime_error("FixedNeighbor diagnostic returned no stop reason.");
    std::size_t accepted{},unchanged{};
    for(const auto & block:result.blocks)
    {accepted+=block.accepted; unchanged+=block.status=="unchanged";}
    return {wall,result.search_seconds,result.work.structural_preparation_seconds,
        result.work.candidate_replay_seconds,result.work.local_search_seconds,
        result.work.sweep_replay_seconds,result.work.sweep_global_state_seconds,candidate_count,
        accepted,unchanged,result.global_state_memory_estimate.estimated_peak_bytes};
}
void Run(int target_count)
{
    constexpr std::size_t measured_samples=5;
    auto problem=rhbm_gem::core::JointProblem(second_stage_test::SyntheticJointWorkload("chain",target_count));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const Vector eta=Vector::Constant(static_cast<Eigen::Index>(data.layout.full_atoms.size()),std::log(.55));
    joint::FixedNeighborPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=5;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.collect_telemetry=true; policy.eta_upper_bound=std::log(2.2);
    joint::StructuralTopology topology(data.domain,data.context.atom_ids,data.layout);
    auto partition=topology.BuildCores(policy.core_atoms); topology.PopulateMetadata(partition);
    std::vector<std::pair<joint::Indices,const joint::StructuralCore *>> cores;
    cores.reserve(partition.cores.size());
    for(const auto & core:partition.cores) cores.emplace_back(core.atoms,&core);
    std::size_t replay_calls{},maximum_candidate_scratch{};
    auto exact_policy=policy;
    exact_policy.candidate_replay_for_testing=[&](const joint::BlockCoordinateState & accepted,
        const joint::Indices & core_atoms,const joint::Indices & affected_rows,joint::VectorRef candidate_eta,
        joint::VectorRef candidate_beta) {
        const auto found=std::find_if(cores.begin(),cores.end(),[&](const auto & entry) {
            return entry.first==core_atoms;
        });
        if(found==cores.end()) throw std::logic_error("D1 benchmark could not find structural core.");
        const auto incidence=second_stage_test::PrepareReplayIncidence(data.domain,data.layout,
            core_atoms,found->second->neighbor_atoms,affected_rows);
        auto candidate=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
            accepted,candidate_eta,candidate_beta,data.context.scale,core_atoms,affected_rows,
            second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
            second_stage_test::ReplayObjectiveStrategy::FullScan,&incidence);
        ++replay_calls;
        maximum_candidate_scratch=std::max(maximum_candidate_scratch,candidate.scratch_bytes);
        return joint::FixedNeighborReplayResult{std::move(candidate.prediction),std::move(candidate.residual),
            candidate.objective,true};
    };

    std::array<Summary,2> summaries{};
    for(std::size_t sample=0;sample<=measured_samples;++sample)
    {
        replay_calls=0; maximum_candidate_scratch=0;
        joint::FixedNeighborDiagnosticResult full_result,exact_result;
        const auto full_sample=RunSearch(problem,eta,policy,0,full_result);
        const auto exact_sample=RunSearch(problem,eta,exact_policy,replay_calls,exact_result);
        VerifyEquivalent(full_result,exact_result);
        if(sample==0) continue;
        auto & exact=summaries[1];
        const std::array<Sample,2> samples{full_sample,exact_sample};
        for(std::size_t method=0;method<samples.size();++method)
        {
            const auto & current=samples[method]; auto & summary=summaries[method];
            summary.wall.push_back(current.wall); summary.search.push_back(current.search);
            summary.structural_preparation.push_back(current.structural_preparation);
            summary.candidate_replay.push_back(current.candidate_replay);
            summary.local_search.push_back(current.local_search); summary.sweep_replay.push_back(current.sweep_replay);
            summary.global_state.push_back(current.global_state); summary.candidate_count=current.candidate_count;
            summary.accepted_count=current.accepted_count; summary.unchanged_count=current.unchanged_count;
            summary.peak_estimate=current.peak_estimate;
        }
        exact.candidate_scratch=std::max(exact.candidate_scratch,maximum_candidate_scratch);
        exact.candidate_count=replay_calls;
    }
    std::size_t support_entries{};
    for(const auto atom:data.layout.full_atoms) support_entries+=data.domain.atoms[static_cast<std::size_t>(atom)].size();
    std::cout<<std::setprecision(17)<<"fixture_targets="<<target_count<<" observations="<<data.domain.rows
        <<" informative_rows="<<data.layout.informative_rows.size()<<" support_entries="<<support_entries
        <<" blocks="<<partition.cores.size()<<" samples="<<measured_samples
        <<" sweeps_per_sample="<<policy.maximum_sweeps<<" endpoint_certification=false search_only=true\n";
    constexpr std::array<const char *,2> names{"A_full_replay","D1_affected_row_exact"};
    for(std::size_t method=0;method<summaries.size();++method)
    {
        const auto & s=summaries[method]; const auto [wall_min,wall_max]=Range(s.wall);
        const auto [search_min,search_max]=Range(s.search);
        const auto [replay_min,replay_max]=Range(s.candidate_replay);
        const double candidate_share=Median(s.candidate_replay)/Median(s.search);
        const double amdahl=candidate_share<1.0 ? 1.0/(1.0-candidate_share) : 0.0;
        std::cout<<"method="<<names[method]
            <<" wall_median_seconds="<<Median(s.wall)<<" wall_range_seconds="<<wall_min<<".."<<wall_max
            <<" search_median_seconds="<<Median(s.search)<<" search_range_seconds="<<search_min<<".."<<search_max
            <<" structural_preparation_median_seconds="<<Median(s.structural_preparation)
            <<" candidate_replay_median_seconds="<<Median(s.candidate_replay)
            <<" candidate_replay_range_seconds="<<replay_min<<".."<<replay_max
            <<" candidate_replay_median_ms="<<Median(s.candidate_replay)*1000.0
            <<" local_search_median_seconds="<<Median(s.local_search)
            <<" sweep_replay_median_seconds="<<Median(s.sweep_replay)
            <<" global_state_median_seconds="<<Median(s.global_state)
            <<" candidate_replay_fraction="<<candidate_share<<" amdahl_candidate_only_upper_bound="<<amdahl
            <<" candidates="<<s.candidate_count<<" accepted="<<s.accepted_count
            <<" unchanged="<<s.unchanged_count<<" base_peak_estimate_bytes="<<s.peak_estimate
            <<" candidate_scratch_estimate_bytes="<<s.candidate_scratch
            <<" estimated_peak_with_candidate_scratch_bytes="<<s.peak_estimate+s.candidate_scratch<<'\n';
    }
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)==0)
    {
#if defined(__APPLE__)
        const auto peak=static_cast<std::uint64_t>(usage.ru_maxrss);
#else
        const auto peak=static_cast<std::uint64_t>(usage.ru_maxrss)*1024;
#endif
        std::cout<<"process_peak_rss_bytes="<<peak<<'\n';
    }
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=2) throw std::invalid_argument("usage: rhbm_fixed_neighbor_closed_loop_benchmark <target-count>");
        Run(std::stoi(argv[1]));
    }
    catch(const std::exception & error)
    {std::cerr<<"fixed-neighbor closed-loop benchmark failed: "<<error.what()<<'\n'; return 1;}
    return 0;
}
