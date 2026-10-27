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
#include <numeric>
#include <sstream>
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
std::string Distribution(const std::vector<double> & values)
{
    if(values.empty()) return "empty";
    auto ordered=values; std::sort(ordered.begin(),ordered.end());
    const double mean=std::accumulate(values.begin(),values.end(),0.0)/static_cast<double>(values.size());
    const auto percentile=[&](double p) {
        return ordered[static_cast<std::size_t>(std::ceil(p*static_cast<double>(ordered.size()-1)))];
    };
    std::ostringstream out; out<<std::setprecision(8)<<"min:"<<ordered.front()<<",median:"<<percentile(.5)
        <<",p95:"<<percentile(.95)<<",p99:"<<percentile(.99)<<",max:"<<ordered.back()<<",mean:"<<mean;
    return out.str();
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
struct ControlledCandidateAudit
{
    std::size_t candidates{},accepted{},rejected{},fallbacks{},estimated_reference_disagreements{},
        actual_reference_disagreements{},reference_positive_accepts{};
    ScalarError objective;
    double maximum_cache_drift{};
    std::vector<double> decision_margins;
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
void AuditControlledCandidate(const rhbm_gem::JointParameterLayout & layout,
    double scale,const joint::FixedNeighborCandidateEvent & event,ControlledCandidateAudit & audit)
{
    if(!event.candidate_replay_ran || !std::isfinite(event.incremental_objective)) return;
    ++audit.candidates; audit.accepted+=event.accepted; audit.rejected+=!event.accepted;
    audit.fallbacks+=event.objective_full_scan_fallback;
    double old_squared{},candidate_squared{};
    for(const auto row_value:layout.informative_rows)
    {
        const auto row=static_cast<Eigen::Index>(row_value);
        const double old_residual=event.accepted_before.residual(row);
        const double candidate_residual=event.full_residual(row);
        old_squared+=old_residual*old_residual;
        candidate_squared+=candidate_residual*candidate_residual;
    }
    const double old_reference=.5*old_squared/(scale*scale);
    const double candidate_reference=.5*candidate_squared/(scale*scale);
    RecordError(audit.objective,candidate_reference,event.incremental_objective);
    audit.maximum_cache_drift=std::max(audit.maximum_cache_drift,
        std::abs(event.accepted_before.objective-old_reference));
    const double estimated_enclosure=joint::BlockObjectiveReplayEnclosure(std::max(
        std::abs(event.accepted_before.objective),std::abs(event.incremental_objective)));
    const bool estimated_decision=ObjectiveProxyAccept(event.accepted_before.objective,
        event.incremental_objective,event.local_objective_before,event.local_objective_after,estimated_enclosure) &&
        !event.unchanged;
    const double reference_enclosure=joint::BlockObjectiveReplayEnclosure(std::max(
        std::abs(old_reference),std::abs(candidate_reference)));
    const bool reference_decision=ObjectiveProxyAccept(old_reference,candidate_reference,
        event.local_objective_before,event.local_objective_after,reference_enclosure) && !event.unchanged;
    audit.estimated_reference_disagreements+=estimated_decision!=reference_decision;
    audit.actual_reference_disagreements+=event.accepted!=reference_decision;
    audit.reference_positive_accepts+=reference_decision && candidate_reference>old_reference;
    const double delta=candidate_reference-old_reference;
    const double consistency=std::abs(delta+event.local_objective_before-event.local_objective_after);
    audit.decision_margins.push_back(std::min({std::abs(delta),std::abs(delta-reference_enclosure),
        std::abs(event.local_objective_after-event.local_objective_before),
        std::abs(consistency-reference_enclosure)}));
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

void RunControlled(int target_count,bool qualify,int selected_method=-1,bool performance_qualified=false)
{
    constexpr std::size_t measured_samples=5;
    const std::size_t sample_count=qualify ? 1 : measured_samples;
    auto problem=rhbm_gem::core::JointProblem(second_stage_test::SyntheticJointWorkload("chain",target_count));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const Vector initial_eta=Vector::Constant(static_cast<Eigen::Index>(data.layout.full_atoms.size()),std::log(.55));
    joint::FixedNeighborPolicy base; base.core_atoms=12; base.maximum_sweeps=qualify || performance_qualified ? 30 : 5;
    base.convergence_criteria=joint::PracticalJointConvergenceCriteria;
    base.stop_after_stationarity=qualify || performance_qualified; base.certify_endpoint=qualify || performance_qualified;
    base.assess_local_conditioning=qualify || performance_qualified; base.assess_parameter_stability=qualify || performance_qualified;
    base.capture_local_trajectory=qualify;
    base.collect_telemetry=true; base.eta_upper_bound=std::log(2.2);
    std::array<joint::FixedNeighborPolicy,3> policies{base,base,base};
    policies[1].objective_strategy_for_testing=joint::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
    policies[2].objective_strategy_for_testing=joint::FixedNeighborObjectiveStrategyForTesting::IncrementalGuarded;
    policies[2].objective_guard_band_fraction_for_testing=.1;
    policies[2].refresh_objective_after_sweep_for_testing=true;
    constexpr std::array<const char *,3> names{"D1_full_objective_scan","E1_direct","E1_guarded"};
    std::array<ControlledCandidateAudit,3> audits{};
    if(qualify)
        for(std::size_t method=1;method<policies.size();++method)
            policies[method].candidate_observer=[&,method](const joint::FixedNeighborCandidateEvent & event) {
                AuditControlledCandidate(data.layout,data.context.scale,event,audits[method]);
            };
    struct Metrics
    {
        std::vector<double> wall,search,candidate,local,sweep,global,endpoint,independent;
        std::size_t candidates{},accepted{},unchanged{},incremental{},fallbacks{},refreshes{},
            peak_estimate{},candidate_scratch{};
        double maximum_cache_replay_error{},maximum_objective_replay_error{};
    };
    std::array<Metrics,3> metrics{};
    std::array<joint::FixedNeighborDiagnosticResult,3> final_results{};
    const std::size_t total_runs=qualify ? 1 : measured_samples+1;
    for(std::size_t sample=0;sample<total_runs;++sample)
        for(std::size_t method=0;method<policies.size();++method)
        {
            if(selected_method>=0 && method!=static_cast<std::size_t>(selected_method)) continue;
            joint::FixedNeighborDiagnosticResult result;
            const auto timing=RunSearch(problem,initial_eta,policies[method],result);
            if(sample==0 && !qualify) continue;
            auto & out=metrics[method];
            out.wall.push_back(timing.wall); out.search.push_back(timing.search);
            out.candidate.push_back(timing.candidate_replay); out.local.push_back(timing.local_search);
            out.sweep.push_back(timing.sweep_replay); out.global.push_back(timing.global_state);
            out.endpoint.push_back(result.work.endpoint_verification_seconds);
            out.independent.push_back(result.work.independent_replay_seconds);
            out.candidates+=timing.candidate_count; out.accepted+=timing.accepted_count;
            out.unchanged+=timing.unchanged_count;
            out.incremental+=result.work.incremental_objective_evaluations;
            out.fallbacks+=result.work.objective_full_scan_fallbacks;
            out.refreshes+=static_cast<std::size_t>(std::count_if(result.sweeps.begin(),result.sweeps.end(),
                [](const auto & sweep){return sweep.objective_cache_refreshed;}));
            out.peak_estimate=std::max(out.peak_estimate,result.global_state_memory_estimate.estimated_peak_bytes);
            out.candidate_scratch=std::max(out.candidate_scratch,result.work.affected_replay_peak_scratch_bytes);
            out.maximum_cache_replay_error=std::max(out.maximum_cache_replay_error,result.work.maximum_cache_replay_error);
            if(result.final_sweep) out.maximum_objective_replay_error=std::max(out.maximum_objective_replay_error,
                result.final_sweep->objective_replay_error);
            final_results[method]=std::move(result);
        }
    std::size_t support_entries{};
    for(const auto atom:data.layout.full_atoms) support_entries+=data.domain.atoms[static_cast<std::size_t>(atom)].size();
    std::cout<<std::setprecision(17)<<"controlled_fixture_targets="<<target_count
        <<" observations="<<data.domain.rows<<" informative_rows="<<data.layout.informative_rows.size()
        <<" support_entries="<<support_entries<<" samples="<<sample_count
        <<" qualification_mode="<<qualify<<" performance_qualified_mode="<<performance_qualified
        <<" sweeps_requested="<<(qualify || performance_qualified ? 30 : 5)
        <<" selected_method="<<(selected_method<0 ? "all" : names[static_cast<std::size_t>(selected_method)])
        <<" initial_eta="<<std::log(.55)<<" initial_beta=0 width_upper_bound=2.2\n";
    for(std::size_t method=0;method<metrics.size();++method)
    {
        if(selected_method>=0 && method!=static_cast<std::size_t>(selected_method)) continue;
        const auto & m=metrics[method]; const auto & result=final_results[method];
        const std::size_t replay_candidates=result.work.full_candidate_replays+
            result.work.affected_candidate_replays;
        const auto accepted=static_cast<std::size_t>(std::count_if(result.blocks.begin(),result.blocks.end(),
            [](const auto & block){return block.accepted;}));
        const auto unchanged=static_cast<std::size_t>(std::count_if(result.blocks.begin(),result.blocks.end(),
            [](const auto & block){return block.status=="unchanged";}));
        std::size_t profile_trial_count{},profile_trial_accepted{};
        for(const auto & block:result.blocks) for(const auto & trial:block.profile_trials)
        {
            ++profile_trial_count;
            profile_trial_accepted+=trial.accepted;
        }
        const std::size_t rejected=replay_candidates>accepted ? replay_candidates-accepted : 0;
        double amplitude_squared{},charge_squared{},width_squared{};
        std::vector<double> amplitudes,charges,widths,amplitude_errors,charge_errors,width_errors;
        amplitudes.reserve(static_cast<std::size_t>(result.state.eta.size()));
        charges.reserve(amplitudes.capacity()); widths.reserve(amplitudes.capacity());
        amplitude_errors.reserve(amplitudes.capacity()); charge_errors.reserve(amplitudes.capacity());
        width_errors.reserve(amplitudes.capacity());
        for(Eigen::Index atom=0;atom<result.state.eta.size();++atom)
        {
            const double amplitude=result.state.beta(2*atom),charge=result.state.beta(2*atom+1),
                width=std::exp(result.state.eta(atom));
            amplitudes.push_back(amplitude); charges.push_back(charge); widths.push_back(width);
            const double error_a=std::abs(amplitude-2.0),error_c=std::abs(charge-.2),error_width=std::abs(width-.5);
            amplitude_errors.push_back(error_a); charge_errors.push_back(error_c); width_errors.push_back(error_width);
            amplitude_squared+=error_a*error_a; charge_squared+=error_c*error_c; width_squared+=error_width*error_width;
        }
        const double atom_count=static_cast<double>(result.state.eta.size());
        const double truth_rmse_a=std::sqrt(amplitude_squared/atom_count);
        const double truth_rmse_c=std::sqrt(charge_squared/atom_count);
        const double truth_rmse_width=std::sqrt(width_squared/atom_count);
        const auto refreshes=static_cast<std::size_t>(std::count_if(result.sweeps.begin(),result.sweeps.end(),
            [](const auto & sweep){return sweep.objective_cache_refreshed;}));
        std::cout<<"method="<<names[method]
            <<" search_median_seconds="<<Median(m.search)<<" search_range_seconds="<<Range(m.search).first<<".."<<Range(m.search).second
            <<" wall_median_seconds="<<Median(m.wall)
            <<" candidate_evaluation_median_seconds="<<Median(m.candidate)
            <<" local_search_median_seconds="<<Median(m.local)
            <<" sweep_replay_median_seconds="<<Median(m.sweep)
            <<" global_state_median_seconds="<<Median(m.global)
            <<" endpoint_verification_median_seconds="<<Median(m.endpoint)
            <<" independent_replay_median_seconds="<<Median(m.independent)
            <<" candidate_replays="<<replay_candidates<<" candidate_accepted="<<accepted
            <<" candidate_rejected="<<rejected<<" unchanged_blocks="<<unchanged
            <<" profile_trials="<<profile_trial_count<<" profile_trials_accepted="<<profile_trial_accepted
            <<" profile_trials_rejected="<<(profile_trial_count-profile_trial_accepted)
            <<" incremental_objective_evaluations="<<result.work.incremental_objective_evaluations
            <<" full_scan_fallbacks="<<result.work.objective_full_scan_fallbacks
            <<" full_scan_fallback_seconds="<<result.work.objective_full_scan_fallback_seconds
            <<" sweep_objective_refreshes="<<refreshes
            <<" sweeps="<<result.sweep_count<<" stop_reason="<<result.reason
            <<" search_converged="<<result.search_converged
            <<" global_ac_kkt="<<(result.final_sweep ? result.final_sweep->global_ac_kkt : 0.0)
            <<" projected_width_kkt="<<(result.final_sweep ? result.final_sweep->projected_width_gradient_inf_norm : 0.0)
            <<" cache_replay_error_max="<<m.maximum_cache_replay_error
            <<" objective_replay_error_max="<<m.maximum_objective_replay_error
            <<" endpoint_certified="<<result.endpoint_certified
            <<" independent_replay_passed="<<(result.assessment.independent_replay.available && result.assessment.independent_replay.passed)
            <<" runtime_convergence="<<static_cast<int>(result.fit.RuntimeConvergence())
            <<" final_objective_cache="<<result.state.objective
            <<" final_objective_independent="<<(result.fit.objective ? *result.fit.objective : std::numeric_limits<double>::quiet_NaN())
            <<" final_A_distribution={"<<Distribution(amplitudes)<<'}'
            <<" final_C_distribution={"<<Distribution(charges)<<'}'
            <<" final_width_distribution={"<<Distribution(widths)<<'}'
            <<" truth_abs_A_error={"<<Distribution(amplitude_errors)<<'}'
            <<" truth_abs_C_error={"<<Distribution(charge_errors)<<'}'
            <<" truth_abs_width_error={"<<Distribution(width_errors)<<'}'
            <<" truth_rmse_A="<<truth_rmse_a<<" truth_rmse_C="<<truth_rmse_c<<" truth_rmse_width="<<truth_rmse_width
            <<" memory_preflight_peak_bytes="<<m.peak_estimate
            <<" affected_replay_scratch_bytes="<<m.candidate_scratch
            <<" parameter_stability="<<(result.parameter_stability ? result.parameter_stability->trajectory_stability_status : "unavailable")
            <<" local_conditioning_blocks="<<result.local_conditioning.size()<<'\n';
        if(qualify && method>0)
        {
            const auto & audit=audits[method];
            std::cout<<"diagnostic_reference="<<names[method]<<" candidates="<<audit.candidates
                <<" accepted="<<audit.accepted<<" rejected="<<audit.rejected
                <<" full_scan_fallbacks="<<audit.fallbacks
                <<" estimated_reference_acceptance_disagreements="<<audit.estimated_reference_disagreements
                <<" actual_reference_acceptance_disagreements="<<audit.actual_reference_disagreements
                <<" reference_positive_accepts="<<audit.reference_positive_accepts
                <<" maximum_objective_absolute_error="<<audit.objective.maximum_absolute
                <<" maximum_objective_relative_error="<<audit.objective.maximum_relative
                <<" maximum_objective_ulp="<<audit.objective.maximum_ulp
                <<" maximum_accepted_cache_drift="<<audit.maximum_cache_drift
                <<" decision_margin_distribution={"<<Distribution(audit.decision_margins)<<"}\n";
        }
        if(method>0 && result.state.eta.size()==final_results[0].state.eta.size())
        {
            std::vector<double> difference_a,difference_c,difference_width;
            for(Eigen::Index atom=0;atom<result.state.eta.size();++atom)
            {
                difference_a.push_back(std::abs(result.state.beta(2*atom)-final_results[0].state.beta(2*atom)));
                difference_c.push_back(std::abs(result.state.beta(2*atom+1)-final_results[0].state.beta(2*atom+1)));
                difference_width.push_back(std::abs(std::exp(result.state.eta(atom))-
                    std::exp(final_results[0].state.eta(atom))));
            }
            std::cout<<"method="<<names[method]<<" abs_A_difference_from_D1={"<<Distribution(difference_a)
                <<"} abs_C_difference_from_D1={"<<Distribution(difference_c)
                <<"} abs_width_difference_from_D1={"<<Distribution(difference_width)<<"}\n";
        }
    }
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)==0)
    {
#if defined(__APPLE__)
        const auto peak=static_cast<std::uint64_t>(usage.ru_maxrss);
#else
        const auto peak=static_cast<std::uint64_t>(usage.ru_maxrss)*1024;
#endif
        std::cout<<"controlled_benchmark_process_peak_rss_bytes="<<peak<<'\n';
    }
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc==2) Run(std::stoi(argv[1]));
        else if(argc==3 && std::string(argv[2])=="--controlled") RunControlled(std::stoi(argv[1]),false);
        else if(argc==3 && std::string(argv[2])=="--qualify") RunControlled(std::stoi(argv[1]),true);
        else if(argc==3 && std::string(argv[2])=="--controlled-d1") RunControlled(std::stoi(argv[1]),false,0);
        else if(argc==3 && std::string(argv[2])=="--controlled-direct") RunControlled(std::stoi(argv[1]),false,1);
        else if(argc==3 && std::string(argv[2])=="--controlled-guarded") RunControlled(std::stoi(argv[1]),false,2);
        else if(argc==3 && std::string(argv[2])=="--performance-qualified-d1") RunControlled(std::stoi(argv[1]),false,0,true);
        else if(argc==3 && std::string(argv[2])=="--performance-qualified-direct") RunControlled(std::stoi(argv[1]),false,1,true);
        else if(argc==3 && std::string(argv[2])=="--performance-qualified-guarded") RunControlled(std::stoi(argv[1]),false,2,true);
        else throw std::invalid_argument("usage: rhbm_fixed_neighbor_closed_loop_benchmark <target-count> [--controlled|--qualify|--controlled-d1|--controlled-direct|--controlled-guarded|--performance-qualified-d1|--performance-qualified-direct|--performance-qualified-guarded]");
    }
    catch(const std::exception & error)
    {std::cerr<<"fixed-neighbor closed-loop benchmark failed: "<<error.what()<<'\n'; return 1;}
    return 0;
}
