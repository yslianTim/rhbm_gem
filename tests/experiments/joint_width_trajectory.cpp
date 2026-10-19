#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/HaloEligibility.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/JointReproducibilityFingerprint.hpp"
#include "support/JointPrecisionAudit.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <set>
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

struct ReplaySetup
{
    core::JointProblem problem;
    std::vector<double> initial_widths;
};

ReplaySetup BuildFullABCReplaySetup(
    const rhbm_gem::ModelObject & model, const core::JointProblem & fixed_problem,
    const rhbm_gem::JointInitialization & initialization)
{
    const auto & source=fixed_problem.Input();
    if(!source.selection_domain ||
        source.selection_domain->contributor_policy!=rhbm_gem::FixedHaloContributorPolicy ||
        !source.selection_domain->fixed_halo_selection)
        throw std::invalid_argument("FullABC replay requires the fixed-halo selection domain.");
    if(initialization.b.size()!=source.atom_ids.size())
        throw std::invalid_argument("FullABC replay initialization dimensions do not match contributors.");

    std::unordered_map<std::string,double> initial_by_id;
    std::unordered_set<std::string> source_ids;
    for(std::size_t atom=0;atom<source.atom_ids.size();++atom)
    {
        initial_by_id.emplace(source.atom_ids[atom],initialization.b[atom]);
        source_ids.insert(source.atom_ids[atom]);
    }
    const auto & source_domain=*source.selection_domain;
    const auto & fixed=*source_domain.fixed_halo_selection;
    core::JointProblemInput full;
    full.row_ids=source.row_ids;
    full.observations=source.observations;
    full.selection_domain=source_domain;
    full.selection_domain->contributor_policy=std::string(rhbm_gem::FullABCHaloContributorPolicy);
    full.selection_domain->fixed_halo_selection.reset();
    full.selection_domain->halo_eligibility.emplace();
    auto & eligibility=*full.selection_domain->halo_eligibility;
    eligibility.minimum_observations=joint::FullABCHaloMinimumObservations;
    eligibility.minimum_distinct_radial_information=joint::FullABCHaloMinimumDistinctRadialInformation;
    eligibility.candidate_halos=fixed.candidate_halos;

    std::set<int> selected_serials;
    for(const auto * atom:model.GetSelectedAtoms())
        if(atom->GetElement()!=Element::HYDROGEN) selected_serials.insert(atom->GetSerialID());
    std::vector<double> initial_widths;
    for(std::size_t old=0;old<source.atom_ids.size();++old)
    {
        const bool target=std::binary_search(source_domain.target_indices.begin(),
            source_domain.target_indices.end(),old);
        if(!target)
        {
            auto support=source.support[old];
            const auto halo=joint::EvaluateFullABCHaloEligibility(
                support,joint::FullABCHaloMinimumObservations);
            if(!halo.eligible)
            {
                eligibility.excluded_halos.push_back({source.atom_ids[old],halo.support_rows,
                    halo.distinct_squared_distances,halo.reason});
                continue;
            }
        }
        const auto new_index=full.atom_ids.size();
        if(target) full.selection_domain->target_indices.push_back(new_index);
        full.atom_ids.push_back(source.atom_ids[old]);
        full.support.push_back(source.support[old]);
        initial_widths.push_back(initial_by_id.at(source.atom_ids[old]));
    }
    for(const auto & atom_ptr:model.GetAtomList())
    {
        const auto * atom=atom_ptr.get();
        if(atom->GetElement()==Element::HYDROGEN || selected_serials.contains(atom->GetSerialID())) continue;
        const auto id=std::to_string(atom->GetSerialID());
        if(!source_ids.contains(id))
            eligibility.excluded_halos.push_back({id,0,0,"insufficient-observations"});
    }
    std::sort(eligibility.excluded_halos.begin(),eligibility.excluded_halos.end(),
        [](const auto & left,const auto & right){return left.atom_id<right.atom_id;});
    eligibility.retained_halos=full.atom_ids.size()-full.selection_domain->target_indices.size();
    if(eligibility.retained_halos+eligibility.excluded_halos.size()!=eligibility.candidate_halos)
        throw std::invalid_argument("FullABC replay halo eligibility dimensions do not match candidates.");
    return {core::JointProblem(std::move(full)),std::move(initial_widths)};
}

ReplaySetup BuildReplaySetup(rhbm_gem::MapObject & map,rhbm_gem::ModelObject & model,bool pure_exclusion)
{
    if(pure_exclusion)
    {
        auto problem=core::BuildPureHaloSupportExclusionProblem(map,model);
        const auto initialization=joint::InitializeContributors(map,model,problem);
        if(!initialization.valid) throw std::runtime_error("First-stage initialization failed: "+initialization.reason);
        return {std::move(problem),initialization.b};
    }
    auto fixed_problem=core::BuildJointProblem(map,model);
    const auto initialization=joint::InitializeContributors(map,model,fixed_problem);
    if(!initialization.valid) throw std::runtime_error("First-stage initialization failed: "+initialization.reason);
    return BuildFullABCReplaySetup(model,fixed_problem,initialization);
}

json::object ProfilePairMetrics(const joint::Evaluation & evaluation,Eigen::Index atom)
{
    if(!evaluation.x.rows() || atom<0 || 2*atom+1>=evaluation.x.cols())
        return { {"available",false} };
    const auto gaussian=evaluation.x.col(2*atom),charge=evaluation.x.col(2*atom+1);
    const double gaussian_norm=gaussian.norm(),charge_norm=charge.norm();
    const double a=evaluation.beta.size()>2*atom ? evaluation.beta(2*atom) : joint::unavailable;
    const double c=evaluation.beta.size()>2*atom+1 ? evaluation.beta(2*atom+1) : joint::unavailable;
    const double normalized_a=std::isfinite(a) ? a*gaussian_norm : joint::unavailable;
    const double normalized_c=std::isfinite(c) ? c*charge_norm : joint::unavailable;
    const auto contribution=(gaussian*a+charge*c).eval();
    const double contribution_norm=contribution.norm();
    const double contribution_sum=std::abs(normalized_a)+std::abs(normalized_c);
    const double cancellation=contribution_norm>0 ? contribution_sum/contribution_norm :
        contribution_sum==0 ? 1.0 : joint::unavailable;
    const double correlation=gaussian_norm>0 && charge_norm>0 ?
        gaussian.dot(charge)/(gaussian_norm*charge_norm) : joint::unavailable;
    return { {"available",true},{"gaussian_column_norm",Number(gaussian_norm)},
        {"charge_column_norm",Number(charge_norm)},
        {"normalized_A",Number(normalized_a)}, {"normalized_C",Number(normalized_c)},
        {"gaussian_contribution_norm",Number(std::abs(normalized_a))},
        {"charge_contribution_norm",Number(std::abs(normalized_c))},
        {"ac_cancellation_ratio",Number(cancellation)}, {"basis_correlation",Number(correlation)} };
}

json::object BuildConditionalProfile(const core::joint_component::ProblemData & data,
    std::size_t component_index,const joint::BlockCoordinateState & state,
    const std::vector<double> & captured_singular_values,
    const std::map<std::size_t,joint::BlockCoordinateState> & sweep_states,
    const std::string & requested_atom_id,std::size_t sweep)
{
    const auto & component=data.partition.components.at(component_index);
    const auto prepared=joint::MakeComponentProblemView(data.layout,component);
    const auto parent_atom=std::find(data.input->atom_ids.begin(),data.input->atom_ids.end(),requested_atom_id);
    if(parent_atom==data.input->atom_ids.end()) throw std::invalid_argument("Conditional profile atom is absent.");
    const auto parent_atom_index=static_cast<Eigen::Index>(parent_atom-data.input->atom_ids.begin());
    const auto local_atom=std::find(prepared.local_to_parent_atoms.begin(),prepared.local_to_parent_atoms.end(),parent_atom_index);
    if(local_atom==prepared.local_to_parent_atoms.end()) throw std::invalid_argument("Conditional profile atom is outside its component.");
    const auto local_atom_index=static_cast<Eigen::Index>(local_atom-prepared.local_to_parent_atoms.begin());

    auto context=data.context;
    context.atom_ids=prepared.atom_ids; context.row_ids=prepared.row_ids;
    context.independent_search=true;
    const auto component_atoms=static_cast<Eigen::Index>(prepared.local_layout.full_atoms.size());
    context.rank={prepared.parent_row_count,2*component_atoms,component_atoms};
    context.linear.rank_relative=context.rank.Relative(2*component_atoms);

    joint::StructuralTopology topology(prepared.domain,context.atom_ids,prepared.local_layout);
    const auto partition=topology.BuildCores(12,true);
    std::size_t block_index{};
    const joint::StructuralCore * target_core=nullptr;
    for(std::size_t index=0;index<partition.cores.size();++index)
        if(std::find(partition.cores[index].atoms.begin(),partition.cores[index].atoms.end(),local_atom_index)!=
            partition.cores[index].atoms.end())
        {block_index=index; target_core=&partition.cores[index]; break;}
    if(!target_core) throw std::invalid_argument("Conditional profile block is absent.");
    const auto & core_atoms=target_core->atoms;
    const auto target_core_atom=static_cast<Eigen::Index>(std::find(core_atoms.begin(),core_atoms.end(),local_atom_index)-core_atoms.begin());
    const auto row_mapping=joint::CompactRowMapping::FromParentRows(target_core->affected_rows);
    const auto profile_domain=prepared.domain.Select(core_atoms,
        static_cast<Eigen::Index>(target_core->affected_rows.size()),row_mapping);
    auto profile_context=context;
    profile_context.atom_ids=context.atom_ids.Select(core_atoms);
    profile_context.row_ids=context.row_ids.Select(target_core->affected_rows);
    profile_context.rank={profile_domain.rows,2*static_cast<Eigen::Index>(core_atoms.size()),
        static_cast<Eigen::Index>(core_atoms.size())};
    profile_context.linear.rank_relative=profile_context.rank.Relative(2*static_cast<Eigen::Index>(core_atoms.size()));
    profile_context.update_budget=1;

    joint::Vector base_eta(static_cast<Eigen::Index>(core_atoms.size()));
    joint::Vector base_beta(2*static_cast<Eigen::Index>(core_atoms.size()));
    for(std::size_t k=0;k<core_atoms.size();++k)
    {
        const auto atom=core_atoms[k];
        base_eta(static_cast<Eigen::Index>(k))=state.eta(atom);
        base_beta.segment<2>(2*static_cast<Eigen::Index>(k))=state.beta.segment<2>(2*atom);
    }
    const auto predict=[&](const joint::Vector & eta,const joint::Vector & beta) {
        joint::Vector prediction=joint::Vector::Zero(profile_domain.rows);
        for(std::size_t atom=0;atom<core_atoms.size();++atom)
        {
            const auto width=std::exp(eta(static_cast<Eigen::Index>(atom)));
            for(const auto support:profile_domain.atoms[atom])
            {
                const auto basis=joint::EvaluateKernel(support.square,width,2.5);
                prediction(support.row)+=beta(2*static_cast<Eigen::Index>(atom))*basis.gaussian+
                    beta(2*static_cast<Eigen::Index>(atom)+1)*basis.charge;
            }
        }
        return prediction;
    };
    const auto base_prediction=predict(base_eta,base_beta);
    joint::Vector local_y(profile_domain.rows);
    for(std::size_t k=0;k<target_core->affected_rows.size();++k)
        local_y(static_cast<Eigen::Index>(k))=base_prediction(static_cast<Eigen::Index>(k))-
            state.residual(target_core->affected_rows[k]);
    const auto base_evaluation=joint::EvaluateState(profile_domain,local_y,base_eta,base_beta,profile_context);
    const auto base_width=std::exp(base_eta(target_core_atom));
    std::vector<double> scan_widths={0.05,0.1,0.2,0.3,0.4,base_width,0.5,0.75,1.0,2.0,4.0,8.0,
        16.0,64.0,256.0,1024.0,4096.0,8192.0,7.647783566,1365.304216,6099.810066};
    std::sort(scan_widths.begin(),scan_widths.end());
    scan_widths.erase(std::unique(scan_widths.begin(),scan_widths.end(),[](double a,double b) {
        return std::abs(a-b)<=1e-12*std::max({1.0,std::abs(a),std::abs(b)});
    }),scan_widths.end());
    json::array profile_rows;
    for(const auto width:scan_widths)
    {
        auto eta=base_eta; eta(target_core_atom)=std::log(width);
        const auto evaluation=joint::EvaluateProfile(profile_domain,local_y,eta,false,&profile_context);
        const auto spectrum=evaluation.x.rows() ? joint::ComputeSpectrum(evaluation.x,profile_context.rank,
            evaluation.x.cols(),true) : joint::Spectrum{};
        const auto pair=ProfilePairMetrics(evaluation,target_core_atom);
        const auto prediction=evaluation.residual.size()==local_y.size() ? evaluation.residual+local_y : joint::Vector{};
        const auto objective=evaluation.certificate.available ? evaluation.certificate.objective/
            (profile_context.scale*profile_context.scale) : joint::unavailable;
        const auto prediction_difference=prediction.size()==base_prediction.size() ?
            (prediction-base_prediction).norm() : joint::unavailable;
        const auto parameter_difference=evaluation.beta.size()==base_beta.size() ?
            (evaluation.beta-base_beta).lpNorm<Eigen::Infinity>() : joint::unavailable;
        profile_rows.push_back(json::object{
            {"width",Number(width)}, {"log_width",Number(std::log(width))},
            {"profile_objective",Number(objective)},
            {"objective_delta_from_sweep_state",Number(objective-base_evaluation.certificate.objective/
                (profile_context.scale*profile_context.scale))},
            {"profile_valid",evaluation.valid},{"profile_reason",evaluation.reason},
            {"local_kkt",Number(evaluation.certificate.projected_kkt)},
            {"profile_gradient_inf_norm",Number(joint::ProfileGradientInfinityNorm(evaluation).value_or(joint::unavailable))},
            {"target_width_gradient",Number(evaluation.gradient.size()>target_core_atom ?
                std::abs(evaluation.gradient(target_core_atom)) : joint::unavailable)},
            {"active_a_atom_indices",[&] {json::array values; for(const auto atom:evaluation.certificate.active_atoms) values.push_back(atom); return values;}()},
            {"local_design_rank",spectrum.available ? json::value(spectrum.rank) : json::value(nullptr)},
            {"local_design_condition",Number(spectrum.condition)},
            {"A",evaluation.beta.size()>2*target_core_atom ? Number(evaluation.beta(2*target_core_atom)) : json::value(nullptr)},
            {"C",evaluation.beta.size()>2*target_core_atom+1 ? Number(evaluation.beta(2*target_core_atom+1)) : json::value(nullptr)},
            {"prediction_difference_norm",Number(prediction_difference)},
            {"parameter_difference_inf_norm",Number(parameter_difference)}, {"pair",pair}});
    }
    json::array captured_singular_values_json;
    for(const auto value:captured_singular_values) captured_singular_values_json.push_back(Number(value));
    json::array neighbor_changes;
    const auto is_core=[&](Eigen::Index atom) {return std::find(core_atoms.begin(),core_atoms.end(),atom)!=core_atoms.end();};
    for(const auto transition:std::array<std::pair<std::size_t,std::size_t>,2>{{{6,7},{7,8}}})
    {
        const auto before=sweep_states.find(transition.first),after=sweep_states.find(transition.second);
        if(before==sweep_states.end() || after==sweep_states.end()) continue;
        double eta_change{},beta_change{}; std::size_t changed{};
        for(Eigen::Index atom=0;atom<before->second.eta.size();++atom) if(!is_core(atom))
        {
            eta_change=std::max(eta_change,std::abs(after->second.eta(atom)-before->second.eta(atom)));
            for(Eigen::Index coefficient=0;coefficient<2;++coefficient)
            {
                const auto lhs=before->second.beta(2*atom+coefficient),rhs=after->second.beta(2*atom+coefficient);
                beta_change=std::max(beta_change,std::abs(rhs-lhs)/(1.0+std::max(std::abs(lhs),std::abs(rhs))));
            }
            if(std::abs(after->second.eta(atom)-before->second.eta(atom))>1e-14 ||
                (after->second.beta.segment<2>(2*atom)-before->second.beta.segment<2>(2*atom)).lpNorm<Eigen::Infinity>()>1e-14) ++changed;
        }
        neighbor_changes.push_back(json::object{{"from_sweep",transition.first},{"to_sweep",transition.second},
            {"changed_external_atoms",changed},{"external_eta_change_inf",Number(eta_change)},
            {"external_beta_scaled_change_inf",Number(beta_change)}});
    }
    json::array high_precision;
    for(const auto width:{base_width,0.5,6099.810066})
    {
        auto eta=base_eta; eta(target_core_atom)=std::log(width);
        try
        {
            high_precision.push_back(json::object{{"width",Number(width)},
                {"reference",second_stage_test::matched::certification::PrecisionProfileComparison(
                    profile_domain,local_y,base_eta,eta,profile_context)}});
        }
        catch(const std::exception & error)
        {high_precision.push_back(json::object{{"width",Number(width)},
            {"error",error.what()}});}
    }
    json::array core_ids; for(const auto atom:core_atoms)
        core_ids.push_back(json::value(context.atom_ids[static_cast<std::size_t>(atom)]));
    return json::object{{"requested_atom_id",requested_atom_id},{"component_index",component_index},
        {"sweep",sweep},{"block",block_index+1},{"core_atom_ids",std::move(core_ids)},
        {"affected_row_count",target_core->affected_rows.size()},{"base_width",Number(base_width)},
        {"base_local_kkt",Number(base_evaluation.certificate.projected_kkt)},
        {"base_profile_objective",Number(base_evaluation.certificate.objective/(profile_context.scale*profile_context.scale))},
        {"captured_projected_width_singular_values",std::move(captured_singular_values_json)},
        {"external_neighbor_state_changes",std::move(neighbor_changes)},
        {"profile_rows",std::move(profile_rows)},{"high_precision_reference",std::move(high_precision)}};
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
    const std::vector<std::string> & requested_atom_ids,std::size_t maximum_events,
    bool accepted_updates_only,bool pure_exclusion)
{
    if(requested_atom_ids.empty() || requested_atom_ids.size()>2048 || maximum_events>10000)
        throw std::invalid_argument("Track 1–2,048 atom IDs and at most 10,000 events.");
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a fresh trajectory output directory.");
    std::filesystem::create_directories(output);
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path);
    const auto normalization=map->MapValueArrayNormalization();
    if(!normalization) throw std::runtime_error("Map normalization unavailable.");
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    model->EditAnalysis().InitializeFromSelection();
    const auto initialization_started=Clock::now();
    auto replay=BuildReplaySetup(*map,*model,pure_exclusion);
    auto problem=std::move(replay.problem);
    const double initialization_seconds=std::chrono::duration<double>(Clock::now()-initialization_started).count();
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto & target_indices=data.input->selection_domain->target_indices;
    const auto fingerprints=rhbm_gem::joint_reproducibility::WriteFingerprint(output,
        model_path,map_path,*normalization,*map,*model,problem,replay.initial_widths,
        {"diagnostic-instrumented",true});
    const auto original_observation_count=data.original_input ? data.original_input->observations.size() :
        data.input->observations.size();
    const auto fitting_observation_count=data.input->observations.size();
    const auto excluded_observation_count=original_observation_count-fitting_observation_count;
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
        capture.atom_initial_widths.emplace_back(id,replay.initial_widths.at(index));
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
    std::vector<double> sweep_eight_projected_width_singular_values;
    std::map<std::size_t,joint::BlockCoordinateState> sweep_states;
    std::optional<json::object> conditional_profile;
    std::string conditional_profile_error;
    const auto update_max=[](double & current,double value) {
        if(std::isfinite(value) && (!std::isfinite(current) || value>current)) current=value;
    };
    joint::FixedNeighborDiagnosticSink sink;
    sink.block=[&](const joint::FixedNeighborBlockRecord & block) {
        if(block.sweep==8)
            for(const auto & trial:block.profile_trials)
                if(trial.accepted && !trial.projected_width_singular_values.empty())
                    sweep_eight_projected_width_singular_values=trial.projected_width_singular_values;
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
    capture.state_observer=[&](std::size_t sweep,const joint::BlockCoordinateState & state,
        const joint::FixedNeighborBlockSweep &) {
        if(sweep==6 || sweep==7 || sweep==8) sweep_states[sweep]=state;
    };
    joint::FixedNeighborSearchPolicy policy;
    policy.core_atoms=12; policy.maximum_sweeps=200; policy.order=joint::FixedNeighborBlockOrder::Forward;
    policy.configure_search_for_testing=[&](joint::FixedNeighborPolicy & local) {
        local.width_trajectory=&capture; local.diagnostic_sink=&sink;
    };
    const auto cpu_before=ProcessCpuSeconds();
    const auto started=Clock::now();
    const auto fit=joint::FitFixedNeighborComponentsImpl(problem,replay.initial_widths,policy);
    const auto cpu_after=ProcessCpuSeconds();
    trajectory.close(); sweeps.close();
    if(const auto state=sweep_states.find(8); state!=sweep_states.end())
    {
        try
        {
            const auto profile_atom=requested_atom_ids.size()==1 && requested_atom_ids.front()!="all" ?
                requested_atom_ids.front() : std::string("1032");
            conditional_profile=BuildConditionalProfile(data,0,state->second,sweep_eight_projected_width_singular_values,sweep_states,
                profile_atom,8);
        }
        catch(const std::exception & error) {conditional_profile_error=error.what();}
    }
    if(conditional_profile)
        Write(output/"conditional-profile.json",*conditional_profile);
    else
        Write(output/"conditional-profile.json",json::object{{"status","Unavailable"},
            {"reason",conditional_profile_error.empty() ? "sweep-8-state-unavailable" : conditional_profile_error}});
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
            const auto replay_check=std::find_if(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
                return check.name=="independent-replay";
            });
            if(endpoint_status!=rhbm_gem::JointCheckStatus::NotRun && replay_check!=component.evidence.end())
                replay_status=replay_check->status;
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
        {"selected_chain","A"},{"replay_parameterization",pure_exclusion ? "PureExclusion" : "FullABC"},
        {"target_count",target_indices.size()},{"original_observation_count",original_observation_count},
        {"fitting_observation_count",fitting_observation_count},{"excluded_observation_count",excluded_observation_count},
        {"core_atoms",12},{"block_order","forward"},{"maximum_sweeps",200},
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
        {"conditional_profile_status",conditional_profile ? "Captured" : "Unavailable"},
        {"conditional_profile_reason",conditional_profile ? "sweep-8-state-profile" :
            conditional_profile_error.empty() ? "sweep-8-state-state-unavailable" : conditional_profile_error},
        {"events_emitted",capture.events_emitted},{"first_stage_initial_widths",std::move(initial_widths)},
        {"fingerprints",fingerprints},
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
            throw std::invalid_argument("Usage: joint_width_trajectory MODEL MAP OUTPUT_DIR ATOM_ID[,ATOM_ID...]|all [MAX_EVENTS] [accepted-updates-only|pure|pure-accepted-updates-only]");
        const auto atom_ids=Split(argv[4]);
        const auto maximum_events=argc>=6 ? static_cast<std::size_t>(std::stoull(argv[5])) : 1024;
        const auto mode=argc==7 ? std::string(argv[6]) : std::string{};
        const bool accepted_only=mode=="accepted-updates-only" || mode=="pure-accepted-updates-only";
        const bool pure_exclusion=mode=="pure" || mode=="pure-accepted-updates-only";
        if(argc==7 && !accepted_only && !pure_exclusion) throw std::invalid_argument("Unknown width trajectory mode.");
        Run(argv[1],argv[2],argv[3],atom_ids,maximum_events,accepted_only,pure_exclusion); return 0;
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
