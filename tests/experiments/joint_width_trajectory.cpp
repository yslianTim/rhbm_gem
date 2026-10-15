#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace json=boost::json;
using Clock=std::chrono::steady_clock;
json::value Number(double value) {return std::isfinite(value) ? json::value(value) : json::value(nullptr);}
std::optional<double> ProcessCpuSeconds()
{
#ifdef _WIN32
    return {};
#else
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) return {};
    const auto seconds=[](const timeval & value) {
        return static_cast<double>(value.tv_sec)+static_cast<double>(value.tv_usec)*1e-6;
    };
    return seconds(usage.ru_utime)+seconds(usage.ru_stime);
#endif
}
json::value ProcessPeakRssBytes()
{
#ifdef _WIN32
    return nullptr;
#else
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) return nullptr;
#ifdef __APPLE__
    return json::value(static_cast<std::uint64_t>(usage.ru_maxrss));
#else
    return json::value(static_cast<std::uint64_t>(usage.ru_maxrss)*1024);
#endif
#endif
}
const char * CheckStatusName(rhbm_gem::JointCheckStatus status)
{
    using Status=rhbm_gem::JointCheckStatus;
    switch(status)
    {
    case Status::Passed:return "Passed";
    case Status::Failed:return "Failed";
    case Status::Unavailable:return "Unavailable";
    case Status::NotRun:return "NotRun";
    }
    return "NotRun";
}
void Write(const std::filesystem::path & path,const json::value & value)
{
    std::ofstream output(path); output.exceptions(std::ios::failbit|std::ios::badbit);
    output<<json::serialize(value)<<'\n';
}
std::vector<std::string> Split(const std::string & value)
{
    std::vector<std::string> out; std::istringstream input(value); std::string part;
    while(std::getline(input,part,',')) if(!part.empty()) out.push_back(part);
    return out;
}
json::object AtomJson(const joint::FixedNeighborProfileTrial::AtomState & atom,bool target)
{
    return {{"atom_id",atom.atom_id},{"role",target ? "target" : "halo"},
        {"amplitude_active",atom.amplitude_active},
        {"first_stage_width",Number(atom.first_stage_width)},
        {"block_entry_width",Number(atom.block_entry_width)},{"lm_proposed_width",Number(atom.lm_proposed_width)},
        {"accepted_width",Number(atom.accepted_width)},{"width_before",Number(atom.width_before)},
        {"width_after",Number(atom.width_after)},{"eta_before",Number(atom.eta_before)},
        {"eta_after",Number(atom.eta_after)},{"proposed_eta_step",Number(atom.proposed_eta_step)},
        {"accepted_eta_step",Number(atom.accepted_eta_step)},{"eta_step_inf",Number(atom.eta_step_inf)},
        {"A",Number(atom.amplitude_a)},{"C",Number(atom.amplitude_c)},
        {"gaussian_column_norm",Number(atom.gaussian_column_norm)},
        {"charge_column_norm",Number(atom.charge_column_norm)},
        {"normalized_A",Number(atom.normalized_a)},{"normalized_C",Number(atom.normalized_c)},
        {"gaussian_contribution_norm",Number(atom.gaussian_contribution_norm)},
        {"charge_contribution_norm",Number(atom.charge_contribution_norm)},
        {"ac_cancellation_ratio",Number(atom.ac_cancellation_ratio)},
        {"basis_correlation",Number(atom.basis_correlation)},
        {"local_kkt",Number(atom.local_kkt)},{"width_gradient_magnitude",Number(atom.width_gradient)},
        {"local_design_rank",atom.local_design_rank<0 ? json::value(nullptr) : json::value(atom.local_design_rank)},
        {"local_design_condition",Number(atom.local_design_condition)},
        {"numerical_status",atom.numerical_status}};
}
void Run(const std::string & model_path,const std::string & map_path,const std::filesystem::path & output,
    const std::vector<std::string> & requested_atom_ids,std::size_t maximum_events,bool accepted_updates_only)
{
    if(requested_atom_ids.empty() || requested_atom_ids.size()>2048 || maximum_events>10000)
        throw std::invalid_argument("Track 1–2,048 atom IDs and at most 10,000 events.");
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a fresh trajectory output directory.");
    std::filesystem::create_directories(output);
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path);
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    const auto problem=core::BuildJointProblem(*map,*model);
    const auto initialization_started=Clock::now();
    const auto initialization=joint::InitializeContributors(*map,*model,problem);
    const double initialization_seconds=std::chrono::duration<double>(Clock::now()-initialization_started).count();
    if(!initialization.valid) throw std::runtime_error("First-stage initialization failed: "+initialization.reason);
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto & target_indices=data.input->selection_domain->target_indices;
    std::unordered_set<std::string> target_ids;
    target_ids.reserve(target_indices.size());
    for(const auto atom:target_indices) target_ids.insert(data.input->atom_ids.at(atom));
    const auto atom_ids=requested_atom_ids.size()==1 && requested_atom_ids.front()=="all" ?
        data.input->atom_ids : requested_atom_ids;
    if(atom_ids.size()>2048) throw std::invalid_argument("The FullABC component exceeds the 2,048-atom diagnostic cap.");
    joint::FixedNeighborWidthTrajectoryCapture capture;
    capture.maximum_events=maximum_events; capture.maximum_atoms=2048; capture.maximum_rows=4096;
    capture.accepted_updates_only=accepted_updates_only;
    for(const auto & id:atom_ids)
    {
        const auto found=std::find(data.input->atom_ids.begin(),data.input->atom_ids.end(),id);
        if(found==data.input->atom_ids.end()) throw std::runtime_error("Tracked atom is not a FullABC contributor: "+id);
        const auto index=static_cast<std::size_t>(found-data.input->atom_ids.begin());
        capture.atom_initial_widths.emplace_back(id,initialization.b.at(index));
    }
    std::ofstream trajectory(output/"width-trajectory.jsonl");
    trajectory.exceptions(std::ios::failbit|std::ios::badbit);
    std::ofstream sweeps(output/"sweeps.jsonl");
    sweeps.exceptions(std::ios::failbit|std::ios::badbit);
    std::size_t recovery_attempts{},recovery_successes{};
    std::size_t lm_trials_captured{},lm_rejections_captured{},accepted_width_updates_captured{};
    double maximum_accepted_width=joint::unavailable,maximum_accepted_log_width_step=joint::unavailable,
        maximum_local_ac_condition=joint::unavailable,maximum_local_ac_kkt=joint::unavailable,
        maximum_ac_cancellation_ratio=joint::unavailable,maximum_global_ac_kkt=joint::unavailable,
        maximum_global_width_gradient=joint::unavailable;
    std::unordered_map<std::string,double> maximum_width_by_atom;
    std::unordered_map<std::string,double> maximum_cancellation_by_atom,maximum_condition_by_atom;
    json::array objective_trajectory;
    const auto update_max=[](double & current,double value) {
        if(std::isfinite(value) && (!std::isfinite(current) || value>current)) current=value;
    };
    joint::FixedNeighborDiagnosticSink sink;
    sink.block=[&](const joint::FixedNeighborBlockRecord & block) {
        if(block.recovery_retry) ++recovery_attempts;
        if(block.recovery_success) ++recovery_successes;
        for(const auto & trial:block.profile_trials)
        {
            if(std::isfinite(trial.lm_predicted_decrease))
            {
                ++lm_trials_captured;
                if(!trial.lm_proposed_acceptance) ++lm_rejections_captured;
            }
            if(block.accepted && trial.accepted && trial.accepted_update.value_or(0)>0)
            {
                ++accepted_width_updates_captured;
                for(const auto & atom:trial.tracked_atoms)
                {
                    update_max(maximum_accepted_width,atom.accepted_width);
                    update_max(maximum_accepted_log_width_step,std::abs(atom.accepted_eta_step));
                    const auto found=maximum_width_by_atom.find(atom.atom_id);
                    if(found==maximum_width_by_atom.end()) maximum_width_by_atom.emplace(atom.atom_id,atom.accepted_width);
                    else update_max(found->second,atom.accepted_width);
                }
            }
            if(block.accepted && trial.accepted)
                for(const auto & atom:trial.tracked_atoms)
                {
                    update_max(maximum_local_ac_condition,atom.local_design_condition);
                    update_max(maximum_local_ac_kkt,atom.local_kkt);
                    update_max(maximum_ac_cancellation_ratio,atom.ac_cancellation_ratio);
                    update_max(maximum_condition_by_atom[atom.atom_id],atom.local_design_condition);
                    update_max(maximum_cancellation_by_atom[atom.atom_id],atom.ac_cancellation_ratio);
                }
            if(trial.tracked_atoms.empty()) continue;
            json::array atoms;
            for(const auto & atom:trial.tracked_atoms)
                atoms.push_back(AtomJson(atom,target_ids.contains(atom.atom_id)));
            trajectory<<json::serialize(json::object{{"component_index",block.component_index},
                {"component_id",block.component_id},{"sweep",block.sweep},{"block",block.block},
                {"block_accepted",block.accepted},
                {"block_status",block.status},{"block_reason",block.reason},
                {"recovery_retry",block.recovery_retry},{"recovery_attempt",block.recovery_attempt},
                {"recovery_success",block.recovery_success},
                {"trial_index",trial.trial_index},{"profile_evaluation",trial.profile_evaluation},
                {"profile_objective_before",Number(trial.local_objective_before)},
                {"profile_objective_after",Number(trial.local_objective_after)},
                {"profile_objective_reduction",Number(trial.objective_reduction)},
                {"block_global_objective_before",Number(block.objective_before)},
                {"block_global_objective_after",Number(block.objective_after)},
                {"accepted_local_objective_improvement",Number(trial.accepted_local_objective_improvement)},
                {"accepted_global_objective_improvement",Number(trial.accepted_global_objective_improvement)},
                {"lm_actual_decrease",Number(trial.lm_actual_decrease)},
                {"lm_predicted_decrease",Number(trial.lm_predicted_decrease)},
                {"lm_ratio",Number(trial.lm_ratio)},{"lm_proposed_acceptance",trial.lm_proposed_acceptance},
                {"weak_width_step_fraction",Number(trial.weak_width_step_fraction)},
                {"weak_width_guard_rejected",trial.weak_width_guard_rejected},
                {"projected_width_singular_values",[&] {
                    json::array values;
                    for(const auto value:trial.projected_width_singular_values) values.push_back(Number(value));
                    return values;
                }()},
                {"accepted",trial.accepted},{"accepted_update",trial.accepted_update ?
                    json::value(*trial.accepted_update) : json::value(nullptr)},
                {"acceptance_reason",trial.acceptance_reason},{"atoms",std::move(atoms)}})<<'\n';
        }
    };
    sink.sweep=[&](const joint::FixedNeighborBlockSweep & sweep) {
        update_max(maximum_global_ac_kkt,sweep.global_ac_kkt);
        update_max(maximum_global_width_gradient,sweep.global_width_gradient_inf_norm);
        objective_trajectory.push_back(json::object{{"sweep",sweep.sweep},
            {"objective_before",Number(sweep.objective_before)},
            {"objective_after",Number(sweep.objective_after)},
            {"global_ac_kkt",Number(sweep.global_ac_kkt)},
            {"global_width_gradient_inf_norm",Number(sweep.global_width_gradient_inf_norm)},
            {"eta_change_inf",Number(sweep.eta_change_inf)},
            {"cache_replay_error",Number(sweep.cache_replay_error)},
            {"objective_replay_error",Number(sweep.objective_replay_error)}});
        sweeps<<json::serialize(json::object{{"sweep",sweep.sweep},
            {"objective_before",Number(sweep.objective_before)},
            {"objective_after",Number(sweep.objective_after)},
            {"global_ac_kkt",Number(sweep.global_ac_kkt)},
            {"global_width_gradient_inf_norm",Number(sweep.global_width_gradient_inf_norm)},
            {"eta_change_inf",Number(sweep.eta_change_inf)},
            {"beta_scaled_change",Number(sweep.beta_scaled_change)},
            {"cache_replay_error",Number(sweep.cache_replay_error)},
            {"objective_replay_error",Number(sweep.objective_replay_error)},
            {"accepted_blocks",sweep.accepted_blocks},{"unchanged_blocks",sweep.unchanged_blocks},
            {"accepted_local_updates",sweep.accepted_local_updates}})<<'\n';
    };
    joint::FixedNeighborSearchPolicy policy;
    policy.core_atoms=12; policy.maximum_sweeps=200; policy.order=joint::FixedNeighborBlockOrder::Forward;
    policy.configure_search_for_testing=[&](joint::FixedNeighborPolicy & local) {
        local.width_trajectory=&capture; local.diagnostic_sink=&sink;
    };
    const auto cpu_before=ProcessCpuSeconds();
    const auto started=Clock::now();
    const auto fit=joint::FitFixedNeighborComponentsImpl(problem,initialization.b,policy);
    const auto cpu_after=ProcessCpuSeconds();
    trajectory.close(); sweeps.close();
    json::array initial_widths;
    for(const auto & [id,width]:capture.atom_initial_widths)
        initial_widths.push_back(json::object{{"atom_id",id},{"width",width}});
    json::array components;
    std::size_t converged_components{},available_target_estimates{};
    bool all_independent_replays_passed=!fit.components.empty(),any_independent_replay_failed{},
        any_independent_replay_unavailable{};
    for(const auto & component:fit.components)
    {
        if(component.search_completed) ++converged_components;
        if(component.state && component.state->b.size()==component.atoms.size())
            for(const auto atom:component.atoms)
                if(std::binary_search(target_indices.begin(),target_indices.end(),atom))
                    ++available_target_estimates;
        auto endpoint_status=rhbm_gem::JointCheckStatus::NotRun;
        auto replay_status=rhbm_gem::JointCheckStatus::NotRun;
        if(component.candidate_diagnostic)
        {
            endpoint_status=component.candidate_diagnostic->endpoint_verification_status;
            const auto replay=std::find_if(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
                return check.name=="independent-replay";
            });
            if(endpoint_status!=rhbm_gem::JointCheckStatus::NotRun && replay!=component.evidence.end())
                replay_status=replay->status;
        }
        if(replay_status!=rhbm_gem::JointCheckStatus::Passed) all_independent_replays_passed=false;
        if(replay_status==rhbm_gem::JointCheckStatus::Failed) any_independent_replay_failed=true;
        if(replay_status==rhbm_gem::JointCheckStatus::Unavailable) any_independent_replay_unavailable=true;
        const auto & diagnostic=component.candidate_diagnostic;
        if(diagnostic && diagnostic->candidate_ac_kkt) update_max(maximum_global_ac_kkt,*diagnostic->candidate_ac_kkt);
        if(diagnostic && diagnostic->candidate_width_gradient_inf_norm)
            update_max(maximum_global_width_gradient,*diagnostic->candidate_width_gradient_inf_norm);
        components.push_back(json::object{{"component_id",component.id},{"stop_reason",component.stop_reason},
            {"search_completed",component.search_completed},{"state_available",component.state.has_value()},
            {"search_sweeps",component.search_sweeps},{"profile_evaluations",component.profile_evaluations},
            {"accepted_updates",component.accepted_updates},
            {"candidate_objective",diagnostic && diagnostic->candidate_objective ?
                Number(*diagnostic->candidate_objective) : json::value(nullptr)},
            {"candidate_ac_kkt",diagnostic && diagnostic->candidate_ac_kkt ?
                Number(*diagnostic->candidate_ac_kkt) : json::value(nullptr)},
            {"candidate_width_gradient_inf_norm",diagnostic && diagnostic->candidate_width_gradient_inf_norm ?
                Number(*diagnostic->candidate_width_gradient_inf_norm) : json::value(nullptr)},
            {"endpoint_verification_status",CheckStatusName(endpoint_status)},
            {"independent_replay_status",CheckStatusName(replay_status)}});
    }
    const auto independent_replay_status=all_independent_replays_passed ? rhbm_gem::JointCheckStatus::Passed :
        any_independent_replay_failed ? rhbm_gem::JointCheckStatus::Failed :
        any_independent_replay_unavailable ? rhbm_gem::JointCheckStatus::Unavailable : rhbm_gem::JointCheckStatus::NotRun;
    std::vector<std::pair<std::string,double>> extreme_widths;
    for(const auto & [id,width]:maximum_width_by_atom) if(width>=10000.0) extreme_widths.emplace_back(id,width);
    std::sort(extreme_widths.begin(),extreme_widths.end(),[](const auto & left,const auto & right) {
        return left.second>right.second;
    });
    json::array extreme_atoms;
    for(const auto & [id,width]:extreme_widths)
        extreme_atoms.push_back(json::object{{"atom_id",id},{"role",target_ids.contains(id) ? "target" : "halo"},
            {"maximum_accepted_width",Number(width)}});
    std::vector<std::pair<std::string,double>> cancellation_outliers;
    for(const auto & [id,ratio]:maximum_cancellation_by_atom)
        if(std::isfinite(ratio)) cancellation_outliers.emplace_back(id,ratio);
    std::sort(cancellation_outliers.begin(),cancellation_outliers.end(),[](const auto & left,const auto & right) {
        return left.second>right.second;
    });
    json::array cancellation_atoms;
    const auto cancellation_count=std::min<std::size_t>(20,cancellation_outliers.size());
    for(std::size_t index=0;index<cancellation_count;++index)
    {
        const auto & [id,ratio]=cancellation_outliers[index];
        const auto condition=maximum_condition_by_atom.find(id);
        cancellation_atoms.push_back(json::object{{"atom_id",id},{"role",target_ids.contains(id) ? "target" : "halo"},
            {"maximum_ac_cancellation_ratio",Number(ratio)},
            {"maximum_local_ac_condition",condition==maximum_condition_by_atom.end() ?
                json::value(nullptr) : Number(condition->second)}});
    }
    const auto cpu_seconds=cpu_before && cpu_after ? json::value(std::max(0.0,*cpu_after-*cpu_before)) :
        json::value(nullptr);
    const bool profile_trial_capture_complete=!capture.accepted_updates_only &&
        capture.events_emitted<capture.maximum_events;
    const bool accepted_update_capture_complete=capture.events_emitted<capture.maximum_events;
    Write(output/"summary.json",json::object{{"source_model",model_path},{"source_map",map_path},
        {"selected_chain","A"},{"core_atoms",12},{"block_order","forward"},{"maximum_sweeps",200},
        {"maximum_events",capture.maximum_events},{"maximum_rows",capture.maximum_rows},
        {"accepted_updates_only",capture.accepted_updates_only},
        {"recovery_attempts",recovery_attempts},{"recovery_successes",recovery_successes},
        {"component_count",fit.components.size()},{"converged_component_count",converged_components},
        {"target_estimate_count",target_indices.size()},{"target_estimates_available",available_target_estimates},
        {"formal_runtime_convergence",CheckStatusName(fit.RuntimeConvergence())},
        {"independent_replay_status",CheckStatusName(independent_replay_status)},
        {"objective",fit.objective ? Number(*fit.objective) : json::value(nullptr)},
        {"maximum_accepted_width_captured_A",Number(maximum_accepted_width)},
        {"maximum_accepted_log_width_step_captured",Number(maximum_accepted_log_width_step)},
        {"maximum_local_ac_condition_captured",Number(maximum_local_ac_condition)},
        {"maximum_local_ac_kkt_captured",Number(maximum_local_ac_kkt)},
        {"maximum_ac_cancellation_ratio_captured",Number(maximum_ac_cancellation_ratio)},
        {"maximum_projected_ac_kkt",Number(maximum_global_ac_kkt)},
        {"maximum_global_width_gradient_inf_norm",Number(maximum_global_width_gradient)},
        {"lm_trials_captured",lm_trials_captured},{"lm_rejections_captured",lm_rejections_captured},
        {"accepted_width_updates_captured",accepted_width_updates_captured},
        {"profile_trial_capture_complete",profile_trial_capture_complete},
        {"accepted_update_capture_complete",accepted_update_capture_complete},
        {"extreme_width_threshold_A",10000.0},
        {"extreme_width_atoms",std::move(extreme_atoms)},
        {"top_ac_cancellation_atoms",std::move(cancellation_atoms)},
        {"objective_trajectory",std::move(objective_trajectory)},
        {"events_emitted",capture.events_emitted},{"first_stage_initial_widths",std::move(initial_widths)},
        {"initialization_seconds",initialization_seconds},{"search_completed",fit.search_completed},
        {"global_state_available",fit.assembled_state.has_value()},{"components",std::move(components)},
        {"problem_construction_preflight_peak_bytes",data.construction_memory.estimated_peak_bytes},
        {"component_preflight_peak_bytes",nullptr},{"component_preflight_reason","Unavailable through the public result"},
        {"sparse_qr_factorization_count",nullptr},
        {"sparse_qr_factorization_count_reason","Not exposed by FixedNeighbor profile telemetry"},
        {"production_high_precision_fallback_calls",0},
        {"cpu_seconds",cpu_seconds},{"process_peak_rss_bytes",ProcessPeakRssBytes()},
        {"wall_seconds",std::chrono::duration<double>(Clock::now()-started).count()}});
}
}
int main(int argc,char ** argv)
{
    try
    {
        if(argc<5 || argc>7)
            throw std::invalid_argument("Usage: joint_width_trajectory MODEL MAP OUTPUT_DIR ATOM_ID[,ATOM_ID...]|all [MAX_EVENTS] [accepted-updates-only]");
        const auto atom_ids=Split(argv[4]);
        const auto maximum_events=argc>=6 ? static_cast<std::size_t>(std::stoull(argv[5])) : 1024;
        const bool accepted_only=argc==7 && std::string(argv[6])=="accepted-updates-only";
        if(argc==7 && !accepted_only) throw std::invalid_argument("Unknown width trajectory mode.");
        Run(argv[1],argv[2],argv[3],atom_ids,maximum_events,accepted_only); return 0;
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
