#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace json=boost::json;
using Clock=std::chrono::steady_clock;
json::value Number(double value) {return std::isfinite(value) ? json::value(value) : json::value(nullptr);}
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
json::object AtomJson(const joint::FixedNeighborProfileTrial::AtomState & atom)
{
    return {{"atom_id",atom.atom_id},{"amplitude_active",atom.amplitude_active},
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
    joint::FixedNeighborDiagnosticSink sink;
    sink.block=[&](const joint::FixedNeighborBlockRecord & block) {
        for(const auto & trial:block.profile_trials)
        {
            if(trial.tracked_atoms.empty()) continue;
            json::array atoms; for(const auto & atom:trial.tracked_atoms) atoms.push_back(AtomJson(atom));
            trajectory<<json::serialize(json::object{{"component_index",block.component_index},
                {"component_id",block.component_id},{"sweep",block.sweep},{"block",block.block},
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
    const auto started=Clock::now();
    const auto fit=joint::FitFixedNeighborComponentsImpl(problem,initialization.b,policy);
    trajectory.close(); sweeps.close();
    json::array initial_widths;
    for(const auto & [id,width]:capture.atom_initial_widths)
        initial_widths.push_back(json::object{{"atom_id",id},{"width",width}});
    json::array components;
    for(const auto & component:fit.components)
        components.push_back(json::object{{"component_id",component.id},{"stop_reason",component.stop_reason},
            {"search_completed",component.search_completed},{"search_sweeps",component.search_sweeps},
            {"profile_evaluations",component.profile_evaluations},{"accepted_updates",component.accepted_updates}});
    Write(output/"summary.json",json::object{{"source_model",model_path},{"source_map",map_path},
        {"selected_chain","A"},{"core_atoms",12},{"block_order","forward"},{"maximum_sweeps",200},
        {"maximum_events",capture.maximum_events},{"maximum_rows",capture.maximum_rows},
        {"accepted_updates_only",capture.accepted_updates_only},
        {"events_emitted",capture.events_emitted},{"first_stage_initial_widths",std::move(initial_widths)},
        {"initialization_seconds",initialization_seconds},{"search_completed",fit.search_completed},
        {"global_state_available",fit.assembled_state.has_value()},{"components",std::move(components)},
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
