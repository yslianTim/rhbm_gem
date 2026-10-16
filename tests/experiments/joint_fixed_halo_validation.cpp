#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace json=boost::json;
using Clock=std::chrono::steady_clock;
struct Truth {double a{},c{},b{};};
struct Fields
{
    std::vector<double> first,truth,truth_gaussian,truth_charge,mixed,first_amplitude_truth_width,target_truth;
    std::array<std::size_t,4> support_counts{},low_support_qualified{},low_support_failed{};
    std::array<std::size_t,4> first_stage_samples_min{},first_stage_samples_max{};
    std::array<double,4> support_error_rmse{};
};
json::value Number(double x) {return std::isfinite(x) ? json::value(x) : json::value(nullptr);}
json::value Ratio(double x,double y) {return Number(x/std::max(1.0,y));}
double Norm(const std::vector<double> & x)
{
    long double sum{}; for(const auto value:x) sum+=static_cast<long double>(value)*value;
    return std::sqrt(static_cast<double>(sum));
}
double Rmse(const std::vector<double> & x,const std::vector<double> & y)
{
    long double sum{}; for(std::size_t i=0;i<x.size();++i)
        {const auto d=x[i]-y[i]; sum+=static_cast<long double>(d)*d;}
    return x.empty() ? std::numeric_limits<double>::quiet_NaN() :
        std::sqrt(static_cast<double>(sum/static_cast<long double>(x.size())));
}
double DifferenceNorm(const std::vector<double> & x,const std::vector<double> & y)
{
    long double sum{}; for(std::size_t i=0;i<x.size();++i)
        {const auto d=x[i]-y[i]; sum+=static_cast<long double>(d)*d;}
    return std::sqrt(static_cast<double>(sum));
}
json::value ReadJson(const std::filesystem::path & path)
{
    std::ifstream input(path); if(!input) throw std::runtime_error("Cannot read "+path.string());
    json::parse_options options; options.numbers=json::number_precision::precise;
    return json::parse(std::string(std::istreambuf_iterator<char>(input),{}),{},options);
}
std::unordered_map<std::string,Truth> ReadTruth(const std::filesystem::path & path)
{
    const auto metadata=ReadJson(path); std::unordered_map<std::string,Truth> out;
    for(const auto & atom:metadata.at("atoms").as_array())
    {
        const auto id=std::to_string(json::value_to<std::uint64_t>(atom.at("serial_id")));
        const double b=json::value_to<double>(atom.at("effective_gaussian_width"));
        if(b!=json::value_to<double>(atom.at("effective_charge_width")))
            throw std::runtime_error("Simulation truth has unequal Gaussian/charge widths at atom "+id);
        out.emplace(id,Truth{static_cast<double>(json::value_to<int>(atom.at("element"))),
            json::value_to<double>(atom.at("charge_used")),b});
    }
    return out;
}
long PeakRssBytes()
{
#ifdef _WIN32
    return -1;
#else
    rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) return -1;
#ifdef __APPLE__
    return usage.ru_maxrss;
#else
    return usage.ru_maxrss*1024;
#endif
#endif
}
json::value PeakRss()
{const auto value=PeakRssBytes(); return value>=0 ? json::value(value) : json::value(nullptr);}
std::size_t RadialCount(const std::vector<core::JointSupport> & support)
{std::set<double> values; for(const auto & item:support) values.insert(item.squared_distance); return values.size();}
void AddContribution(std::vector<double> & field,const core::JointProblemInput & input,std::size_t atom,
    double a,double b,double c)
{
    for(const auto & support:input.support[atom])
    {
        const auto basis=joint::EvaluateKernel(support.squared_distance,b,2.5);
        field[support.row]+=a*basis.gaussian+c*basis.charge;
    }
}
Fields BuildFields(const core::JointProblem & problem,const core::JointInitialization & init,
    const std::unordered_map<std::string,Truth> & truth,double amplitude_scale=1.0)
{
    const auto & input=problem.Input(); const auto rows=input.observations.size();
    std::array<std::vector<double>,4> first_by_support,truth_by_support;
    for(std::size_t bucket=0;bucket<4;++bucket)
    {first_by_support[bucket].assign(rows,0.0); truth_by_support[bucket].assign(rows,0.0);}
    Fields out; out.first.assign(rows,0); out.truth.assign(rows,0); out.truth_gaussian.assign(rows,0);
    out.truth_charge.assign(rows,0); out.mixed.assign(rows,0); out.first_amplitude_truth_width.assign(rows,0);
    out.target_truth.assign(rows,0);
    out.first_stage_samples_min.fill(std::numeric_limits<std::size_t>::max());
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        const auto found=truth.find(input.atom_ids[atom]);
        if(found==truth.end()) throw std::runtime_error("Missing simulation truth for atom "+input.atom_ids[atom]);
        const auto & t=found->second;
        const double truth_a=t.a*amplitude_scale,truth_c=t.c*amplitude_scale;
        if(std::binary_search(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),atom))
        {AddContribution(out.target_truth,input,atom,truth_a,t.b,truth_c); continue;}
        const auto & record=init.atoms.at(atom); const auto support_count=input.support[atom].size();
        const auto bucket=std::min<std::size_t>(3,support_count);
        ++out.support_counts[bucket];
        out.first_stage_samples_min[bucket]=std::min(out.first_stage_samples_min[bucket],record.sample_count);
        out.first_stage_samples_max[bucket]=std::max(out.first_stage_samples_max[bucket],record.sample_count);
        if(record.first_stage_fit_qualified) ++out.low_support_qualified[bucket]; else ++out.low_support_failed[bucket];
        const double a=record.mdpde[0],b=record.mdpde[1];
        if(!std::isfinite(a) || a<0 || !std::isfinite(b) || b<=0)
            throw std::runtime_error("No finite FirstStage A/B candidate for halo "+input.atom_ids[atom]);
        AddContribution(out.first,input,atom,a,b,0.0);
        AddContribution(out.truth,input,atom,truth_a,t.b,truth_c);
        AddContribution(out.truth_gaussian,input,atom,truth_a,t.b,0.0);
        AddContribution(out.truth_charge,input,atom,0.0,t.b,truth_c);
        AddContribution(out.first_amplitude_truth_width,input,atom,a,t.b,0.0);
        AddContribution(first_by_support[bucket],input,atom,a,b,0.0);
        AddContribution(truth_by_support[bucket],input,atom,t.a,t.b,t.c);
    }
    for(std::size_t row=0;row<rows;++row) out.mixed[row]=out.first[row]+out.truth_charge[row];
    for(std::size_t bucket=0;bucket<4;++bucket)
        if(out.support_counts[bucket]) out.support_error_rmse[bucket]=Rmse(first_by_support[bucket],truth_by_support[bucket]);
    return out;
}
std::optional<double> EvidenceValue(const std::vector<rhbm_gem::JointCheck> & checks,const std::string & name)
{for(const auto & check:checks) if(check.name==name) return check.value; return std::nullopt;}
json::object RunCase(const std::string & name,const core::JointProblem & original,
    const core::JointInitialization & original_init,const std::unordered_map<std::string,Truth> & truths,
    const std::vector<double> & fixed_field,const std::vector<double> & scenario_field,bool oracle,
    const Fields & fields,double truth_amplitude_scale,
    const joint::JointConvergenceCriteria & criteria)
{
    const auto & source=original.Input(); auto input=source;
    const double original_scale=original.ObservationScale();
    std::vector<double> adjusted(source.observations.size());
    for(std::size_t row=0;row<adjusted.size();++row)
        adjusted[row]=source.observations[row]+fixed_field[row]-scenario_field[row];
    const double adjustment=original_scale/std::max(1.0,Norm(adjusted));
    for(auto & value:adjusted) value*=adjustment;
    input.observations=std::move(adjusted);
    core::JointProblem problem(std::move(input));
    if(std::abs(problem.ObservationScale()-original_scale)>1e-12*std::max(1.0,original_scale))
        throw std::runtime_error("controlled scenario changed the original observation scale");
    auto init=original_init;
    for(std::size_t atom=0;atom<init.atoms.size();++atom)
    {
        auto & record=init.atoms[atom];
        if(std::binary_search(source.selection_domain->target_indices.begin(),source.selection_domain->target_indices.end(),atom)) continue;
        const auto & t=truths.at(source.atom_ids[atom]);
        if(oracle)
        {
            record.mdpde[0]=t.a*truth_amplitude_scale; record.mdpde[1]=t.b;
            record.native_status=static_cast<int>(rhbm_gem::RHBMEstimationStatus::SUCCESS);
            record.first_stage_fit_qualified=true;
            record.first_stage_qualification=rhbm_gem::RHBMSolveQualification::NativeSuccess;
            record.first_stage_recovery.reset();
        }
        record.mdpde[0]*=adjustment;
    }
    std::map<std::size_t,std::size_t> block_count;
    const auto & targets=source.selection_domain->target_indices;
    std::unordered_map<std::size_t,std::size_t> target_slot;
    for(std::size_t i=0;i<targets.size();++i) target_slot.emplace(targets[i],i);
    struct CapturedSweep {std::vector<json::value> parameters; std::vector<bool> captured; json::object metrics;};
    std::map<std::size_t,CapturedSweep> checkpoints;
    for(const auto sweep:{79U,96U,113U,140U,200U})
    {auto & capture=checkpoints[sweep]; capture.parameters.resize(targets.size()); capture.captured.assign(targets.size(),false);}
    json::array trajectory;
    const auto observer=[&](const auto & event) {
        if(!event.fixed_neighbor) return;
        const auto index=event.component_index;
        block_count[index]=std::max(block_count[index],event.fixed_neighbor->block_solves);
        const auto & sweep=*event.fixed_neighbor;
        trajectory.emplace_back(json::object{{"sweep",sweep.sweep},{"objective",Number(sweep.objective)},
            {"target_ac_kkt",Number(sweep.global_ac_kkt)},
            {"target_width_gradient_inf_norm",Number(sweep.width_gradient_inf_norm)},
            {"complete_sweep_log_width_change",Number(event.fixed_neighbor->eta_change_inf)},
            {"maximum_accepted_log_width_step",Number(sweep.maximum_accepted_log_width_step)},
            {"maximum_target_width",Number(sweep.maximum_target_width)},
            {"accepted_blocks",sweep.accepted_blocks},{"unchanged_blocks",sweep.unchanged_blocks},
            {"failed_blocks",sweep.failed_blocks},{"accepted_local_updates",sweep.local_accepted_updates},
            {"profile_evaluations",sweep.profile_evaluations},{"recovery_attempts",sweep.recovery_attempts},
            {"recovery_successes",sweep.recovery_successes},{"sweep_seconds",Number(event.elapsed_seconds)}});
        std::cout<<"case="<<name<<" sweep="<<sweep.sweep<<" objective="<<sweep.objective
            <<" kkt="<<sweep.global_ac_kkt<<" width_gradient="<<sweep.width_gradient_inf_norm
            <<" accepted="<<sweep.accepted_blocks<<" unchanged="<<sweep.unchanged_blocks
            <<" failed="<<sweep.failed_blocks<<" profiles="<<sweep.profile_evaluations
            <<" recovery="<<sweep.recovery_attempts<<'/'<<sweep.recovery_successes
            <<" max_log_width_step="<<sweep.maximum_accepted_log_width_step
            <<" max_B="<<sweep.maximum_target_width<<" seconds="<<event.elapsed_seconds<<std::endl;
    };
    std::vector<rhbm_gem::JointFixedHaloParameter> fixed_halos;
    for(std::size_t atom=0;atom<init.atoms.size();++atom)
        if(!std::binary_search(source.selection_domain->target_indices.begin(),
            source.selection_domain->target_indices.end(),atom))
            fixed_halos.push_back({atom,source.atom_ids[atom],init.atoms[atom].mdpde[0],init.atoms[atom].mdpde[1],0.0,
                init.atoms[atom].sample_count,init.atoms[atom].native_status,
                std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),init.atoms[atom].first_stage_fit_qualified,
                init.atoms[atom].first_stage_qualification,init.atoms[atom].first_stage_recovery});
    const auto initialized_problem=problem.WithFixedHalos(std::move(fixed_halos));
    const auto & problem_data=core::JointProblemAccess::Get(initialized_problem);
    std::vector<double> target_widths; target_widths.reserve(source.selection_domain->target_indices.size());
    for(const auto target:source.selection_domain->target_indices) target_widths.push_back(init.b.at(target));
    joint::FixedNeighborSearchPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=200;
    policy.convergence_criteria=criteria;
    policy.state_observer=[&](std::size_t component_index,std::size_t sweep,const joint::BlockCoordinateState & state,
        const joint::FixedNeighborBlockSweep & record) {
        const auto found=checkpoints.find(sweep); if(found==checkpoints.end()) return;
        auto & capture=found->second;
        capture.metrics={{"objective",Number(record.objective_after)},
            {"target_ac_kkt",Number(record.global_ac_kkt)},
            {"target_width_gradient_inf_norm",Number(record.global_width_gradient_inf_norm)},
            {"complete_sweep_log_width_change",Number(record.eta_change_inf)},
            {"maximum_accepted_log_width_step",Number(record.maximum_accepted_log_width_step)},
            {"maximum_target_width",Number(record.maximum_target_width)},
            {"failed_blocks",record.block_solves-record.accepted_blocks-record.unchanged_blocks}};
        if(component_index>=problem_data.partition.components.size()) return;
        const auto & component=problem_data.partition.components[component_index];
        for(std::size_t k=0;k<component.atoms.size() && k<static_cast<std::size_t>(state.eta.size());++k)
        {
            const auto atom=static_cast<std::size_t>(component.atoms[k]);
            const auto slot=target_slot.find(atom); if(slot==target_slot.end() || 2*k+1>=static_cast<std::size_t>(state.beta.size())) continue;
            const auto index=slot->second;
            capture.parameters[index]=json::object{{"atom_id",source.atom_ids[atom]},
                {"A",state.beta(static_cast<Eigen::Index>(2*k))/adjustment/truth_amplitude_scale},
                {"B",std::exp(state.eta(static_cast<Eigen::Index>(k)))},
                {"C",state.beta(static_cast<Eigen::Index>(2*k+1))/adjustment/truth_amplitude_scale}};
            capture.captured[index]=true;
        }
    };
    const auto start=Clock::now();
    const auto fit=joint::FitFixedNeighborComponents(initialized_problem,target_widths,policy,observer);
    const double seconds=std::chrono::duration<double>(Clock::now()-start).count();
    std::vector<double> predicted_target(source.observations.size());
    std::vector<double> estimate_a(targets.size()),estimate_b(targets.size()),estimate_c(targets.size());
    std::vector<bool> available(targets.size());
    std::vector<double> candidate_target(source.observations.size());
    std::vector<double> candidate_a(targets.size()),candidate_b(targets.size()),candidate_c(targets.size());
    std::vector<bool> candidate_available(targets.size());
    for(const auto & component:fit.components)
    {
        if(component.state)
            for(std::size_t k=0;k<component.atoms.size();++k)
            {
                const auto found=target_slot.find(component.atoms[k]); if(found==target_slot.end()) continue;
                const auto slot=found->second; estimate_a[slot]=component.state->ac[2*k]/adjustment/truth_amplitude_scale;
                estimate_c[slot]=component.state->ac[2*k+1]/adjustment/truth_amplitude_scale; estimate_b[slot]=component.state->b[k];
                available[slot]=true;
                AddContribution(predicted_target,source,targets[slot],estimate_a[slot]*truth_amplitude_scale,
                    estimate_b[slot],estimate_c[slot]*truth_amplitude_scale);
            }
        if(component.candidate_diagnostic && component.candidate_diagnostic->candidate_available)
        {
            const auto & candidate=*component.candidate_diagnostic;
            for(std::size_t k=0;k<component.atoms.size() && k<candidate.candidate_eta.size() &&
                2*k+1<candidate.candidate_ac.size();++k)
            {
                const auto found=target_slot.find(component.atoms[k]); if(found==target_slot.end()) continue;
                const auto slot=found->second; candidate_a[slot]=candidate.candidate_ac[2*k]/adjustment/truth_amplitude_scale;
                candidate_c[slot]=candidate.candidate_ac[2*k+1]/adjustment/truth_amplitude_scale;
                candidate_b[slot]=std::exp(candidate.candidate_eta[k]); candidate_available[slot]=true;
                AddContribution(candidate_target,source,targets[slot],candidate_a[slot]*truth_amplitude_scale,
                    candidate_b[slot],candidate_c[slot]*truth_amplitude_scale);
            }
        }
    }
    const auto available_count=static_cast<std::size_t>(std::count(available.begin(),available.end(),true));
    const auto candidate_count=static_cast<std::size_t>(std::count(candidate_available.begin(),candidate_available.end(),true));
    std::size_t converged{},sweeps{},profiles{},accepted{},conditioning_blocks{},stability_blocks{};
    for(const auto & component:fit.components)
    {
        if(component.TargetRuntimeConvergence()==rhbm_gem::JointCheckStatus::Passed) ++converged;
        sweeps+=component.search_sweeps; profiles+=static_cast<std::size_t>(std::max(0,component.profile_evaluations));
        accepted+=static_cast<std::size_t>(std::max(0,component.accepted_updates));
        conditioning_blocks+=component.local_conditioning.size();
        if(component.parameter_stability) stability_blocks+=component.parameter_stability->blocks.size();
    }
    const bool complete=available_count==targets.size();
    double target_squared_a{},target_squared_b{},target_squared_c{};
    for(std::size_t i=0;i<targets.size();++i) if(available[i])
    {
        const auto & t=truths.at(source.atom_ids[targets[i]]);
        target_squared_a+=std::pow(estimate_a[i]-t.a,2); target_squared_b+=std::pow(estimate_b[i]-t.b,2);
        target_squared_c+=std::pow(estimate_c[i]-t.c,2);
    }
    double candidate_squared_a{},candidate_squared_b{},candidate_squared_c{};
    for(std::size_t i=0;i<targets.size();++i) if(candidate_available[i])
    {
        const auto & t=truths.at(source.atom_ids[targets[i]]);
        candidate_squared_a+=std::pow(candidate_a[i]-t.a,2); candidate_squared_b+=std::pow(candidate_b[i]-t.b,2);
        candidate_squared_c+=std::pow(candidate_c[i]-t.c,2);
    }
    double full_rmse=std::numeric_limits<double>::quiet_NaN(),objective=std::numeric_limits<double>::quiet_NaN();
    if(complete)
    {
        std::vector<double> prediction=predicted_target;
        for(std::size_t row=0;row<prediction.size();++row) prediction[row]+=scenario_field[row];
        full_rmse=Rmse(prediction,source.observations);
        const double residual=full_rmse*std::sqrt(static_cast<double>(prediction.size()));
        objective=.5*residual*residual/(original_scale*original_scale);
    }
    const bool candidate_complete=candidate_count==targets.size();
    double candidate_full_rmse=std::numeric_limits<double>::quiet_NaN(),candidate_objective=std::numeric_limits<double>::quiet_NaN();
    if(candidate_complete)
    {
        std::vector<double> prediction=candidate_target;
        for(std::size_t row=0;row<prediction.size();++row) prediction[row]+=scenario_field[row];
        candidate_full_rmse=Rmse(prediction,source.observations);
        const double residual=candidate_full_rmse*std::sqrt(static_cast<double>(prediction.size()));
        candidate_objective=.5*residual*residual/(original_scale*original_scale);
    }
    std::optional<double> candidate_kkt,candidate_width_gradient;
    for(const auto & component:fit.components) if(component.candidate_diagnostic)
    {
        if(component.candidate_diagnostic->candidate_ac_kkt)
            candidate_kkt=std::max(candidate_kkt.value_or(0.0),*component.candidate_diagnostic->candidate_ac_kkt);
        if(component.candidate_diagnostic->candidate_width_gradient_inf_norm)
            candidate_width_gradient=std::max(candidate_width_gradient.value_or(0.0),
                *component.candidate_diagnostic->candidate_width_gradient_inf_norm);
    }
    const auto target_convergence=fit.TargetRuntimeConvergence();
    json::array comp;
    for(const auto & item:fit.components) comp.emplace_back(json::object{{"id",item.id},{"target_atoms",item.atoms.size()},
        {"rows",item.rows.size()},{"search_completed",item.search_completed},{"sweeps",item.search_sweeps},
        {"first_order_stationarity_sweep",item.first_order_stationarity_sweep},
        {"confirmed_stationarity_sweep",item.confirmed_stationarity_sweep},
        {"endpoint_certification_sweep",item.endpoint_certification_sweep},
        {"profile_evaluations",item.profile_evaluations},{"accepted_updates",item.accepted_updates},
        {"target_convergence",static_cast<int>(item.TargetRuntimeConvergence())},{"stop_reason",item.stop_reason}});
    const auto kkt=EvidenceValue(fit.evidence,"target-ac-optimality");
    const auto width=EvidenceValue(fit.evidence,"target-width-stationarity");
    const auto replay=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="full-observation-replay";
    });
    json::array detailed_replay;
    const auto replay_started=Clock::now();
    for(std::size_t component_index=0;component_index<fit.components.size();++component_index)
    {
        const auto & component=fit.components[component_index];
        if(!component.candidate_diagnostic || !component.candidate_diagnostic->candidate_available) continue;
        const auto prepared=joint::MakeComponentProblemView(problem_data.layout,
            problem_data.partition.components[component_index]);
        const auto & candidate=*component.candidate_diagnostic;
        joint::Vector eta(static_cast<Eigen::Index>(candidate.candidate_eta.size()));
        joint::Vector beta(static_cast<Eigen::Index>(candidate.candidate_ac.size()));
        for(std::size_t k=0;k<candidate.candidate_eta.size();++k) eta(static_cast<Eigen::Index>(k))=candidate.candidate_eta[k];
        for(std::size_t k=0;k<candidate.candidate_ac.size();++k) beta(static_cast<Eigen::Index>(k))=candidate.candidate_ac[k];
        const auto local_y=joint::SelectValues(problem_data.y,prepared.local_to_parent_rows);
        auto context=problem_data.context; context.atom_ids=prepared.atom_ids; context.row_ids=prepared.row_ids;
        context.independent_search=true;
        const auto atom_count=static_cast<Eigen::Index>(prepared.local_layout.full_atoms.size());
        context.rank={prepared.parent_row_count,2*atom_count,atom_count};
        context.linear.rank_relative=context.rank.Relative(2*atom_count);
        const auto primary=joint::EvaluateStateStreaming(prepared.domain,local_y,eta,beta,context);
        const auto independent=joint::CheckIndependentReplay(prepared.domain,local_y,primary,context,
            nullptr,nullptr,criteria);
        const auto assessment=joint::VerifyJointEndpoint(prepared.domain,local_y,primary,independent,context,criteria);
        detailed_replay.emplace_back(json::object{{"component_id",component.id},
            {"primary_valid",independent.primary_valid},{"primary_kkt_passed",independent.primary_kkt_passed},
            {"independent_kkt_passed",independent.independent_kkt_passed},
            {"primary_ac_kkt",Number(independent.primary_ac_kkt)},
            {"independent_ac_kkt",Number(independent.independent_ac_kkt)},
            {"primary_width_gradient",Number(independent.primary_width_gradient)},
            {"independent_width_gradient",Number(independent.independent_width_gradient)},
            {"kkt_agrees",independent.kkt_agrees},{"kkt_difference",Number(independent.kkt_difference)},
            {"gradient_agrees",independent.gradient_agrees},{"gradient_difference",Number(independent.gradient_difference)},
            {"prediction_passed",independent.prediction_passed},{"prediction_difference",Number(independent.prediction_difference)},
            {"residual_passed",independent.residual_passed},{"residual_difference",Number(independent.residual_difference)},
            {"objective_passed",independent.objective_passed},{"objective_difference",Number(independent.objective_difference)},
            {"maximum_gradient_difference_atom_id",independent.maximum_gradient_difference_atom_id},
            {"maximum_gradient_difference_primary",Number(independent.maximum_gradient_difference_primary)},
            {"maximum_gradient_difference_independent",Number(independent.maximum_gradient_difference_independent)},
            {"maximum_gradient_difference",Number(independent.gradient_difference)},
            {"maximum_gradient_difference_bound",Number(independent.maximum_gradient_difference_bound)},
            {"maximum_gradient_difference_ratio",Number(independent.maximum_gradient_difference_ratio)},
            {"primary_objective",Number(independent.primary_objective)},
            {"independent_objective",Number(independent.independent_objective)},
            {"primary_residual_norm",Number(independent.primary_residual_norm)},
            {"independent_residual_norm",Number(independent.independent_residual_norm)},
            {"parameter_feasibility_passed",independent.feasible},{"replay_passed",independent.passed},
            {"endpoint_verified",assessment.endpoint_verified},{"reason",independent.reason},
            {"maximum_gradient_difference_atom_id",independent.maximum_gradient_difference_atom_id},
            {"gradient_enclosure_failure_atom_id",independent.gradient_enclosure_failure_atom_id},
            {"gradient_enclosure_failure_primary",Number(independent.gradient_enclosure_failure_primary)},
            {"gradient_enclosure_failure_independent",Number(independent.gradient_enclosure_failure_independent)},
            {"gradient_enclosure_failure_difference",Number(independent.gradient_enclosure_failure_difference)},
            {"gradient_enclosure_failure_bound",Number(independent.gradient_enclosure_failure_bound)},
            {"gradient_enclosure_failure_ratio",Number(independent.gradient_enclosure_failure_ratio)}});
    }
    const double detailed_replay_seconds=std::chrono::duration<double>(Clock::now()-replay_started).count();
    json::array checkpoint_results;
    for(auto & [sweep,capture]:checkpoints)
    {
        json::array parameters;
        for(std::size_t i=0;i<capture.parameters.size();++i)
            if(capture.captured[i]) parameters.emplace_back(std::move(capture.parameters[i]));
        capture.metrics["sweep"]=sweep;
        capture.metrics["parameter_count"]=parameters.size();
        capture.metrics["complete"]=parameters.size()==targets.size();
        capture.metrics["parameters"]=std::move(parameters);
        checkpoint_results.emplace_back(std::move(capture.metrics));
    }
    json::array candidate_parameters,candidate_target_prediction,candidate_full_prediction;
    if(candidate_complete)
    {
        for(std::size_t i=0;i<targets.size();++i)
            candidate_parameters.emplace_back(json::object{{"atom_id",source.atom_ids[targets[i]]},
                {"A",candidate_a[i]},{"B",candidate_b[i]},{"C",candidate_c[i]}});
        for(std::size_t row=0;row<candidate_target.size();++row)
        {
            candidate_target_prediction.emplace_back(candidate_target[row]);
            candidate_full_prediction.emplace_back(candidate_target[row]+scenario_field[row]);
        }
    }
    const bool verified=complete && converged==fit.components.size() &&
        target_convergence==rhbm_gem::JointCheckStatus::Passed && replay!=fit.evidence.end() &&
        replay->status==rhbm_gem::JointCheckStatus::Passed;
    json::array local_conditioning;
    for(const auto & component:fit.components) for(const auto & item:component.local_conditioning)
    {
        json::array atom_ids;
        for(const auto atom:item.core_atoms) if(atom<source.atom_ids.size()) atom_ids.emplace_back(source.atom_ids[atom]);
        local_conditioning.emplace_back(json::object{{"block_id",item.block_id},
            {"core_atom_ids",std::move(atom_ids)},{"classification",item.conditioning_classification},
            {"free_ac_status",item.free_ac_status},{"free_ac_rank",item.free_ac_rank ? json::value(*item.free_ac_rank) : json::value(nullptr)},
            {"free_ac_columns",item.free_ac_columns},{"free_ac_condition",Number(item.free_ac_condition.value_or(
                std::numeric_limits<double>::quiet_NaN()))},
            {"projected_width_status",item.projected_width_status},
            {"projected_width_rank",item.projected_width_rank ? json::value(*item.projected_width_rank) : json::value(nullptr)},
            {"projected_width_condition",Number(item.projected_width_condition.value_or(
                std::numeric_limits<double>::quiet_NaN()))},
            {"profile_status",item.profile_status},{"profile_rank",item.profile_rank ? json::value(*item.profile_rank) : json::value(nullptr)},
            {"profile_condition",Number(item.profile_condition.value_or(
                std::numeric_limits<double>::quiet_NaN()))}});
    }
    return {{"name",name},{"research_only",true},{"maximum_sweeps",policy.maximum_sweeps},
        {"criteria_id",criteria.global_ac_kkt==joint::PracticalJointConvergenceCriteria.global_ac_kkt ?
            "practical-v1" : "strict-v1"},
        {"criteria",json::object{{"global_ac_kkt",criteria.global_ac_kkt},
            {"global_width_gradient",criteria.global_width_gradient},
            {"confirmed_log_width_change",criteria.confirmed_log_width_change}}},
        {"source_first_stage_failed_halos",std::count_if(original_init.atoms.begin(),original_init.atoms.end(),
            [](const auto & record) {return !record.target && !record.first_stage_fit_qualified;})},
        {"observation_scale_preserved",std::abs(problem.ObservationScale()-original_scale)<=
            1e-12*std::max(1.0,original_scale)},{"scale_factor",adjustment},
        {"target_count",targets.size()},{"fixed_halo_count",source.atom_ids.size()-targets.size()},
        {"component_count",fit.components.size()},{"largest_component",[&] {std::size_t m{};
            for(const auto & c:fit.components) m=std::max(m,c.atoms.size()); return m;}()},
        {"target_only_blocks_per_sweep",[&] {std::size_t n{}; for(const auto & [_,v]:block_count) n+=v; return n;}()},
        {"profile_evaluations",profiles},{"accepted_updates",accepted},{"completed_sweeps",sweeps},
        {"converged_components",converged},{"available_target_estimates",available_count},
        {"target_convergence",static_cast<int>(target_convergence)},
        {"convergence_status",verified ? "converged-and-verified" : candidate_count ?
            "finite-candidate-but-unqualified" : "unavailable"},
        {"target_ac_kkt",kkt ? Number(*kkt) : json::value(nullptr)},
        {"target_width_gradient",width ? Number(*width) : json::value(nullptr)},
        {"independent_replay",replay==fit.evidence.end() ? "unavailable" :
            replay->status==rhbm_gem::JointCheckStatus::Passed ? "passed" :
            replay->status==rhbm_gem::JointCheckStatus::Failed ? "failed" : "unavailable"},
        {"target_parameter_rmse_status",verified ? "qualified" : candidate_count ?
            "finite-candidate-but-unqualified" : "unavailable"},
        {"target_rmse_A",verified ? Number(std::sqrt(target_squared_a/static_cast<double>(available_count))) : json::value(nullptr)},
        {"target_rmse_B",verified ? Number(std::sqrt(target_squared_b/static_cast<double>(available_count))) : json::value(nullptr)},
        {"target_rmse_C",verified ? Number(std::sqrt(target_squared_c/static_cast<double>(available_count))) : json::value(nullptr)},
        {"target_prediction_rmse",verified ? Number(Rmse(predicted_target,fields.target_truth)) : json::value(nullptr)},
        {"full_prediction_rmse",verified ? Number(full_rmse) : json::value(nullptr)},
        {"objective",verified ? Number(objective) : json::value(nullptr)},
        {"candidate_available_target_estimates",candidate_count},
        {"candidate_status",verified ? "qualified" : candidate_count ? "finite-candidate-unqualified" : "unavailable"},
        {"candidate_ac_kkt",candidate_kkt ? Number(*candidate_kkt) : json::value(nullptr)},
        {"candidate_width_gradient",candidate_width_gradient ? Number(*candidate_width_gradient) : json::value(nullptr)},
        {"candidate_rmse_A",candidate_count ? Number(std::sqrt(candidate_squared_a/static_cast<double>(candidate_count))) : json::value(nullptr)},
        {"candidate_rmse_B",candidate_count ? Number(std::sqrt(candidate_squared_b/static_cast<double>(candidate_count))) : json::value(nullptr)},
        {"candidate_rmse_C",candidate_count ? Number(std::sqrt(candidate_squared_c/static_cast<double>(candidate_count))) : json::value(nullptr)},
        {"candidate_target_prediction_rmse",candidate_complete ? Number(Rmse(candidate_target,fields.target_truth)) : json::value(nullptr)},
        {"candidate_full_prediction_rmse",candidate_complete ? Number(candidate_full_rmse) : json::value(nullptr)},
        {"candidate_objective",candidate_complete ? Number(candidate_objective) : json::value(nullptr)},
        {"candidate_parameters",std::move(candidate_parameters)},
        {"candidate_target_prediction",std::move(candidate_target_prediction)},
        {"candidate_full_prediction",std::move(candidate_full_prediction)},
        {"sweep_checkpoints",std::move(checkpoint_results)},
        {"detailed_independent_replay",std::move(detailed_replay)},
        {"local_conditioning",std::move(local_conditioning)},
        {"detailed_replay_seconds",detailed_replay_seconds},
        {"full_observation_replay",replay==fit.evidence.end() ? "unavailable" :
            replay->status==rhbm_gem::JointCheckStatus::Passed ? "passed" :
            replay->status==rhbm_gem::JointCheckStatus::Failed ? "failed" : "unavailable"},
        {"scenario_halo_signal_l2_error",Number(DifferenceNorm(scenario_field,fields.truth))},
        {"scenario_halo_signal_rmse",Number(Rmse(scenario_field,fields.truth))},
        {"trajectory",std::move(trajectory)},
        {"target_conditioning_blocks_assessed",conditioning_blocks},{"target_stability_blocks_assessed",stability_blocks},
        {"wall_seconds",seconds},{"process_peak_rss_bytes",PeakRss()},{"components",std::move(comp)}};
}
}
int main(int argc,char ** argv)
{
    try
    {
        if(argc<4 || argc>7) throw std::invalid_argument(
            "Usage: joint_fixed_halo_validation MODEL MAP OUTPUT.json [scenario] [strict|practical] [--cli-normalization]");
        const std::filesystem::path model_path=argv[1],map_path=argv[2],output_path=argv[3];
        std::string scenario="all",criteria_name="strict";
        bool cli_normalization=false,scenario_seen=false,criteria_seen=false;
        for(int arg=4;arg<argc;++arg)
        {
            const std::string value=argv[arg];
            if(value=="--cli-normalization") cli_normalization=true;
            else if(value=="strict" || value=="practical")
            {
                if(criteria_seen) throw std::invalid_argument("Convergence criteria specified more than once");
                criteria_name=value; criteria_seen=true;
            }
            else if(value=="all" || value=="production" || value=="gaussian-oracle" ||
                value=="charge-oracle" || value=="full-oracle")
            {
                if(scenario_seen) throw std::invalid_argument("Controlled scenario specified more than once");
                scenario=value; scenario_seen=true;
            }
            else throw std::invalid_argument("Unknown research option: "+value);
        }
        const auto selected=[&](const std::string & name) {return scenario=="all" || scenario==name;};
        if(scenario!="all" && scenario!="production" && scenario!="gaussian-oracle" &&
            scenario!="charge-oracle" && scenario!="full-oracle")
            throw std::invalid_argument("Unknown controlled scenario: "+scenario);
        auto model=rhbm_gem::ReadModel(model_path.string()); auto map=rhbm_gem::ReadMap(map_path.string());
        double normalization_divisor=1.0;
        if(cli_normalization)
            if(const auto divisor=map->MapValueArrayNormalization()) normalization_divisor=*divisor;
        const double truth_amplitude_scale=1.0/normalization_divisor;
        const auto & criteria=criteria_name=="practical" ? joint::PracticalJointConvergenceCriteria :
            joint::StrictJointConvergenceCriteria;
        model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
            return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
        });
        const auto problem=core::BuildJointProblem(*map,*model);
        const auto initialization=joint::InitializeContributors(*map,*model,problem);
        const auto truth=ReadTruth(map_path.string()+".simulation.json");
        const auto fields=BuildFields(problem,initialization,truth,truth_amplitude_scale);
        const auto & input=problem.Input(); const auto failed=std::count_if(initialization.atoms.begin(),
            initialization.atoms.end(),[](const auto & record) {return !record.target && !record.first_stage_fit_qualified;});
        const double fixed_error=Rmse(fields.first,fields.truth),truth_norm=Norm(fields.truth);
        json::array failures,case_results;
        for(std::size_t atom=0;atom<initialization.atoms.size();++atom)
        {
            const auto & record=initialization.atoms[atom];
            if(!record.target && !record.first_stage_fit_qualified)
                failures.emplace_back(json::object{{"atom_id",record.id},{"native_status",record.native_status.value_or(-1)},
                    {"reason",record.reason},{"joint_support",input.support[atom].size()},
                    {"distinct_radial_support",RadialCount(input.support[atom])},
                    {"first_stage_samples",record.sample_count},{"A",Number(record.mdpde[0])},{"B",Number(record.mdpde[1])}});
        }
        if(!initialization.valid)
            throw std::runtime_error("P3 controlled fits require a valid FirstStage fixed-halo snapshot: "+initialization.reason);
        if(selected("production"))
            case_results.emplace_back(RunCase("production-first-stage-ab-c0",problem,initialization,truth,
                fields.first,fields.first,false,fields,truth_amplitude_scale,criteria));
        if(selected("gaussian-oracle"))
            case_results.emplace_back(RunCase("gaussian-oracle-ab-c0",problem,initialization,truth,
                fields.truth_gaussian,fields.truth_gaussian,true,fields,truth_amplitude_scale,criteria));
        if(selected("charge-oracle"))
            case_results.emplace_back(RunCase("charge-oracle-first-stage-ab-truth-c",problem,initialization,truth,
                fields.first,fields.mixed,false,fields,truth_amplitude_scale,criteria));
        if(selected("full-oracle"))
            case_results.emplace_back(RunCase("full-oracle-ab-c",problem,initialization,truth,
                fields.truth_gaussian,fields.truth,true,fields,truth_amplitude_scale,criteria));
        json::array low_support;
        for(std::size_t bucket=0;bucket<4;++bucket)
        {
            const auto min=fields.first_stage_samples_min[bucket];
            low_support.emplace_back(json::object{{"joint_support",bucket==0 ? "zero" : bucket==1 ? "one" : bucket==2 ? "two" : "three-plus"},
                {"halo_count",fields.support_counts[bucket]},{"first_stage_qualified",fields.low_support_qualified[bucket]},
                {"first_stage_failed",fields.low_support_failed[bucket]},
                {"first_stage_samples_min",min==std::numeric_limits<std::size_t>::max() ? json::value(nullptr) : json::value(min)},
                {"first_stage_samples_max",fields.first_stage_samples_max[bucket]},
                {"contribution_error_rmse",fields.support_counts[bucket] ?
                    Number(fields.support_error_rmse[bucket]) : json::value(nullptr)}});
        }
        const double amplitude_error=DifferenceNorm(fields.first_amplitude_truth_width,fields.truth_gaussian);
        const double width_error=DifferenceNorm(fields.first,fields.first_amplitude_truth_width);
        const double fixed_delta_error=DifferenceNorm(fields.first,fields.truth);
        double identity_error{};
        for(std::size_t row=0;row<input.observations.size();++row)
            identity_error=std::max(identity_error,std::abs((fields.truth[row]-fields.first[row])-
                ((fields.truth_gaussian[row]-fields.first[row])+fields.truth_charge[row])));
        json::object out{{"contract","fixed-halo-target-only-model-error-study-v3"},{"scenario_filter",scenario},
            {"criteria_id",criteria_name=="practical" ? "practical-v1" : "strict-v1"},
            {"global_ac_kkt_threshold",criteria.global_ac_kkt},
            {"global_width_gradient_threshold",criteria.global_width_gradient},
            {"confirmed_log_width_change_threshold",criteria.confirmed_log_width_change},
            {"cli_map_normalization_applied",cli_normalization && normalization_divisor!=1.0},
            {"map_normalization_divisor",normalization_divisor},
            {"model",model_path.string()},
            {"map",map_path.string()},{"simulation_truth",map_path.string()+".simulation.json"},
            {"targets",input.selection_domain->target_indices.size()},
            {"contributors",input.atom_ids.size()},{"fixed_halos",input.atom_ids.size()-input.selection_domain->target_indices.size()},
            {"observations",input.observations.size()},{"observation_scale",problem.ObservationScale()},
            {"initialization_valid",initialization.valid},{"initialization_reason",initialization.reason},
            {"unqualified_first_stage_halos",failed},{"failed_halos",std::move(failures)},
            {"first_stage_halo_accuracy_status",failed ? "candidate-only" : "qualified"},
            {"first_stage_halo_A_rmse",nullptr},{"first_stage_halo_B_rmse",nullptr},
            {"first_stage_fixed_halo_prediction_rmse",Number(fixed_error)},
            {"first_stage_fixed_halo_l2_error",Number(DifferenceNorm(fields.first,fields.truth))},
            {"first_stage_truth_halo_l2_norm",Number(truth_norm)},
            {"first_stage_fixed_halo_relative_error",Ratio(DifferenceNorm(fields.first,fields.truth),truth_norm)},
            {"first_stage_gaussian_amplitude_error_l2",Number(amplitude_error)},
            {"first_stage_gaussian_amplitude_error_rmse",Number(amplitude_error/std::sqrt(static_cast<double>(input.observations.size())))},
            {"first_stage_gaussian_width_error_l2",Number(width_error)},
            {"first_stage_gaussian_width_error_rmse",Number(width_error/std::sqrt(static_cast<double>(input.observations.size())))},
            {"missing_halo_charge_l2",Number(Norm(fields.truth_charge))},
            {"missing_halo_charge_rmse",Number(Norm(fields.truth_charge)/std::sqrt(static_cast<double>(input.observations.size())))},
            {"fixed_halo_delta_h_l2",Number(fixed_delta_error)},
            {"fixed_halo_delta_h_rmse",Number(fixed_delta_error/std::sqrt(static_cast<double>(input.observations.size())))},
            {"delta_h_decomposition_max_abs_error",Number(identity_error)},
            {"first_stage_fixed_halo_max_abs_error",[&] {
                double m{};
                for(std::size_t i=0;i<fields.first.size();++i) m=std::max(m,std::abs(fields.first[i]-fields.truth[i]));
                return m;
            }()},
            {"low_support",std::move(low_support)},{"cases",std::move(case_results)},
            {"study_note","All four cases use the same observation rows, raw-data scale, target set, target widths, FixedNeighbor solver and 200-sweep budget. The Gaussian, charge and full oracle cases exist only in this research executable; oracle values never enter production."},
            {"process_peak_rss_bytes",PeakRss()}};
        double aa{},bb{}; std::size_t halos{};
        for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
            if(!std::binary_search(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),atom))
            {const auto & record=initialization.atoms[atom]; const auto & t=truth.at(input.atom_ids[atom]);
            aa+=std::pow(record.mdpde[0]/truth_amplitude_scale-t.a,2); bb+=std::pow(record.mdpde[1]-t.b,2); ++halos;}
        out["first_stage_halo_A_rmse"]=Number(std::sqrt(aa/static_cast<double>(halos)));
        out["first_stage_halo_B_rmse"]=Number(std::sqrt(bb/static_cast<double>(halos)));
        std::ofstream output(output_path); output.exceptions(std::ios::failbit|std::ios::badbit);
        output<<json::serialize(out)<<'\n';
        json::array case_summary;
        for(const auto & case_value:out.at("cases").as_array())
        {
            const auto & result=case_value.as_object();
            case_summary.emplace_back(json::object{{"name",result.at("name")},
                {"criteria_id",result.at("criteria_id")},{"completed_sweeps",result.at("completed_sweeps")},
                {"convergence_status",result.at("convergence_status")},
                {"available_target_estimates",result.at("available_target_estimates")},
                {"target_ac_kkt",result.at("target_ac_kkt")},
                {"target_width_gradient",result.at("target_width_gradient")},
                {"candidate_objective",result.at("candidate_objective")},
                {"independent_replay",result.at("independent_replay")},
                {"wall_seconds",result.at("wall_seconds")},
                {"process_peak_rss_bytes",result.at("process_peak_rss_bytes")}});
        }
        std::cout<<json::serialize(json::object{{"output",output_path.string()},{"targets",out.at("targets")},
            {"fixed_halos",out.at("fixed_halos")},{"observations",out.at("observations")},
            {"first_stage_failures",failed},{"cases",std::move(case_summary)}})<<'\n';
        return 0;
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
