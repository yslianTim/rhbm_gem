#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
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
    bool full_replay_memory_fallback{};
};
struct Summary
{
    std::vector<double> wall,search,structural_preparation,candidate_replay,local_search,sweep_replay,global_state;
    std::size_t candidate_count{},accepted_count{},unchanged_count{},peak_estimate{},candidate_scratch{};
    bool full_replay_memory_fallback{};
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
    const joint::FixedNeighborPolicy & policy,
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
        result.work.sweep_replay_seconds,result.work.sweep_global_state_seconds,
        result.work.full_candidate_replays+result.work.affected_candidate_replays,
        accepted,unchanged,result.global_state_memory_estimate.estimated_peak_bytes,
        result.work.full_replay_memory_fallback};
}
struct ScalarError
{
    std::size_t mismatches{};
    std::uint64_t maximum_ulp{};
    double maximum_absolute{},maximum_relative{};
};
struct ShadowMetrics
{
    std::size_t candidates{},baseline_accepted{},e1_acceptance_disagreements{},
        c_acceptance_proxy_disagreements{};
    ScalarError e1_objective,c_prediction,c_residual,c_objective;
};
std::uint64_t OrderedBits(double value)
{
    const auto bits=std::bit_cast<std::uint64_t>(value);
    return bits&(std::uint64_t{1}<<63) ? ~bits : bits|(std::uint64_t{1}<<63);
}
void RecordError(ScalarError & error,double reference,double candidate)
{
    if(!std::isfinite(reference) || !std::isfinite(candidate)) return;
    const double absolute=std::abs(reference-candidate);
    if(absolute!=0) ++error.mismatches;
    error.maximum_absolute=std::max(error.maximum_absolute,absolute);
    error.maximum_relative=std::max(error.maximum_relative,
        absolute/std::max(std::abs(reference),std::numeric_limits<double>::min()));
    if(std::isnormal(reference) && std::isnormal(candidate))
    {
        const auto a=OrderedBits(reference),b=OrderedBits(candidate);
        error.maximum_ulp=std::max(error.maximum_ulp,a>b ? a-b : b-a);
    }
}
bool ObjectiveProxyAccept(double objective_before,double objective_after,
    double local_before,double local_after,double objective_enclosure)
{
    const double replay_reference=std::max(std::abs(objective_before),std::abs(objective_after));
    const double delta=objective_after-objective_before;
    const double consistency_error=std::abs(delta+local_before-local_after);
    return joint::WithinBlockObjectiveReplay(consistency_error,replay_reference) &&
        (delta<=0 || (local_after<=local_before && delta<=objective_enclosure));
}
void RecordIncrementalObjectiveShadow(const joint::Domain & domain,
    const rhbm_gem::JointParameterLayout & layout,joint::VectorRef observations,
    double scale,const joint::FixedNeighborCandidateEvent & event,ShadowMetrics & metrics)
{
    if(!event.candidate_replay_ran) return;
    ++metrics.candidates;
    metrics.baseline_accepted+=event.accepted;
    double e1_delta_squared{};
    for(const auto row:event.affected_rows)
    {
        const double old_residual=event.accepted_before.residual(row);
        const double new_residual=event.full_residual(row);
        e1_delta_squared+=new_residual*new_residual-old_residual*old_residual;
    }
    const double e1_objective=event.accepted_before.objective+
        .5*e1_delta_squared/(scale*scale);
    RecordError(metrics.e1_objective,event.full_objective,e1_objective);
    const double replay_reference=std::max(std::abs(event.full_objective),
        std::abs(event.accepted_before.objective));
    const double enclosure=joint::BlockObjectiveReplayEnclosure(replay_reference);
    const bool e1_accept=ObjectiveProxyAccept(event.accepted_before.objective,e1_objective,
        event.local_objective_before,event.local_objective_after,enclosure) && !event.unchanged;
    metrics.e1_acceptance_disagreements+=e1_accept!=event.accepted;

    Vector c_prediction(static_cast<Eigen::Index>(event.affected_rows.size()));
    for(std::size_t k=0;k<event.affected_rows.size();++k)
        c_prediction(static_cast<Eigen::Index>(k))=event.accepted_before.prediction(event.affected_rows[k]);
    for(const auto atom_value:layout.full_atoms)
    {
        const auto atom=static_cast<Eigen::Index>(atom_value);
        if(std::find(event.core_atoms.begin(),event.core_atoms.end(),atom)==event.core_atoms.end()) continue;
        const double old_width=std::exp(event.accepted_before.eta(atom));
        const double new_width=std::exp(event.candidate_eta(atom));
        for(const auto & support:domain.atoms[atom_value])
        {
            const auto found=std::lower_bound(event.affected_rows.begin(),event.affected_rows.end(),support.row);
            if(found==event.affected_rows.end() || *found!=support.row) continue;
            const auto local=static_cast<Eigen::Index>(found-event.affected_rows.begin());
            const auto old_basis=joint::EvaluateKernel(support.square,old_width,2.5);
            const auto new_basis=joint::EvaluateKernel(support.square,new_width,2.5);
            const double old_contribution=event.accepted_before.beta(2*atom)*old_basis.gaussian+
                event.accepted_before.beta(2*atom+1)*old_basis.charge;
            const double new_contribution=event.candidate_beta(2*atom)*new_basis.gaussian+
                event.candidate_beta(2*atom+1)*new_basis.charge;
            c_prediction(local)+=new_contribution-old_contribution;
        }
    }
    double c_delta_squared{};
    for(std::size_t k=0;k<event.affected_rows.size();++k)
    {
        const auto row=event.affected_rows[k];
        const double residual=c_prediction(static_cast<Eigen::Index>(k))-observations(row);
        RecordError(metrics.c_prediction,event.full_prediction(row),c_prediction(static_cast<Eigen::Index>(k)));
        RecordError(metrics.c_residual,event.full_residual(row),residual);
        const double old_residual=event.accepted_before.residual(row);
        c_delta_squared+=residual*residual-old_residual*old_residual;
    }
    const double c_objective=event.accepted_before.objective+.5*c_delta_squared/(scale*scale);
    RecordError(metrics.c_objective,event.full_objective,c_objective);
    const bool c_accept=ObjectiveProxyAccept(event.accepted_before.objective,c_objective,
        event.local_objective_before,event.local_objective_after,
        joint::BlockObjectiveReplayEnclosure(std::max(std::abs(c_objective),
            std::abs(event.accepted_before.objective)))) && !event.unchanged;
    metrics.c_acceptance_proxy_disagreements+=c_accept!=event.accepted;
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
    auto full_policy=policy;
    full_policy.candidate_replay_for_testing=[&](const joint::BlockCoordinateState &,
        const joint::Indices &,const joint::Indices &,joint::VectorRef candidate_eta,
        joint::VectorRef candidate_beta) {
        return joint::FullReplayForTesting(data.domain,data.layout,data.y,
            candidate_eta,candidate_beta,data.context.scale);
    };

    std::array<Summary,2> summaries{};
    for(std::size_t sample=0;sample<=measured_samples;++sample)
    {
        joint::FixedNeighborDiagnosticResult full_result,exact_result;
        const auto full_sample=RunSearch(problem,eta,full_policy,full_result);
        const auto exact_sample=RunSearch(problem,eta,policy,exact_result);
        if(exact_result.work.full_replay_memory_fallback)
            throw std::runtime_error("D1 benchmark unexpectedly fell back to Full Replay.");
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
            summary.full_replay_memory_fallback=current.full_replay_memory_fallback;
        }
        exact.candidate_scratch=std::max(exact.candidate_scratch,
            exact_result.work.affected_replay_peak_scratch_bytes);
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
            <<" unchanged="<<s.unchanged_count<<" memory_preflight_peak_estimate_bytes="<<s.peak_estimate
            <<" affected_replay_peak_scratch_bytes="<<s.candidate_scratch
            <<" full_replay_memory_fallback="<<s.full_replay_memory_fallback<<'\n';
    }
    ShadowMetrics shadow;
    auto shadow_policy=policy;
    shadow_policy.candidate_observer=[&](const joint::FixedNeighborCandidateEvent & event) {
        RecordIncrementalObjectiveShadow(data.domain,data.layout,data.y,data.context.scale,event,shadow);
    };
    const auto shadow_result=joint::RunFixedNeighborDiagnostic(problem,eta,shadow_policy);
    if(shadow_result.reason.empty()) throw std::runtime_error("D1 shadow search returned no stop reason.");
    std::cout<<"shadow_candidates="<<shadow.candidates<<" baseline_accepted="<<shadow.baseline_accepted
        <<" E1_objective_mismatches="<<shadow.e1_objective.mismatches
        <<" E1_objective_max_abs="<<shadow.e1_objective.maximum_absolute
        <<" E1_objective_max_rel="<<shadow.e1_objective.maximum_relative
        <<" E1_objective_max_ulp="<<shadow.e1_objective.maximum_ulp
        <<" E1_acceptance_disagreements="<<shadow.e1_acceptance_disagreements
        <<" E2_same_as_E1=true"
        <<" C_prediction_mismatches="<<shadow.c_prediction.mismatches
        <<" C_prediction_max_abs="<<shadow.c_prediction.maximum_absolute
        <<" C_prediction_max_rel="<<shadow.c_prediction.maximum_relative
        <<" C_prediction_max_ulp="<<shadow.c_prediction.maximum_ulp
        <<" C_residual_mismatches="<<shadow.c_residual.mismatches
        <<" C_residual_max_abs="<<shadow.c_residual.maximum_absolute
        <<" C_residual_max_rel="<<shadow.c_residual.maximum_relative
        <<" C_residual_max_ulp="<<shadow.c_residual.maximum_ulp
        <<" C_objective_mismatches="<<shadow.c_objective.mismatches
        <<" C_objective_max_abs="<<shadow.c_objective.maximum_absolute
        <<" C_objective_max_rel="<<shadow.c_objective.maximum_relative
        <<" C_objective_max_ulp="<<shadow.c_objective.maximum_ulp
        <<" C_acceptance_proxy_disagreements="<<shadow.c_acceptance_proxy_disagreements
        <<" conditional_full_scan_fallback_rate=1.0"
        <<" search_only=true endpoint_certification=false independent_replay=false"
        <<" search_converged="<<shadow_result.search_converged
        <<" final_reason="<<shadow_result.reason
        <<" final_objective="<<shadow_result.state.objective
        <<" global_ac_kkt="<<(shadow_result.final_sweep ? shadow_result.final_sweep->global_ac_kkt : 0.0)
        <<" projected_width_kkt="<<(shadow_result.final_sweep ?
            shadow_result.final_sweep->projected_width_gradient_inf_norm : 0.0)
        <<" final_eta_min="<<shadow_result.state.eta.minCoeff()
        <<" final_eta_max="<<shadow_result.state.eta.maxCoeff()
        <<" final_beta_min="<<shadow_result.state.beta.minCoeff()
        <<" final_beta_max="<<shadow_result.state.beta.maxCoeff()<<'\n';
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
