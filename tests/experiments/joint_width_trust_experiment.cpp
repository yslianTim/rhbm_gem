#include "core/detail/joint_component/Problem.hpp"
#include "support/JointWidthTrustPolicy.hpp"
#include "support/JointWidthTrustTruth.hpp"

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
#include <limits>
#include <numbers>
#include <optional>
#include <random>
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
namespace trust=joint_test::width_trust;
using Clock=std::chrono::steady_clock;

using Truth=trust::Truth;

struct ExperimentInput
{
    core::JointProblem problem;
    std::vector<double> initial_widths;
    std::vector<Truth> truth;
    bool truth_available{};
    bool normalization_applied{};
    std::string normalization_operation{"not-applicable"};
    double normalization_divisor{1};
    std::string source;
    std::string simulation_mode{"not-applicable"},truth_source{"none"},truth_scale{"N/A — no truth"},
        estimate_scale{"original-domain-A/C; physical-width-B"},comparison_scale{"N/A — no truth"};
};

json::value Number(double value)
{ return std::isfinite(value) ? json::value(value) : json::value(nullptr); }

void Write(const std::filesystem::path & path,const json::value & value)
{
    std::ofstream output(path);
    output.exceptions(std::ios::failbit|std::ios::badbit);
    output<<json::serialize(value)<<'\n';
}

std::string Status(rhbm_gem::JointCheckStatus status)
{
    using S=rhbm_gem::JointCheckStatus;
    switch(status)
    {
    case S::Passed:return "Passed";
    case S::Failed:return "Failed";
    case S::Unavailable:return "Unavailable";
    case S::NotRun:return "NotRun";
    }
    return "NotRun";
}

std::uint64_t Seed(const std::string & value)
{ return static_cast<std::uint64_t>(std::stoull(value)); }

long PeakRssBytes()
{
#ifdef _WIN32
    return -1;
#else
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) return -1;
#ifdef __APPLE__
    return usage.ru_maxrss;
#else
    return usage.ru_maxrss*1024;
#endif
#endif
}

std::vector<std::array<double,3>> SyntheticAtoms()
{
    std::vector<std::array<double,3>> positions;
    for(int y=0;y<3;++y) for(int x=0;x<4;++x)
        positions.push_back({(static_cast<double>(x)-1.5)*1.5,(static_cast<double>(y)-1.0)*1.5,0.0});
    return positions;
}

bool IsSyntheticDataset(const std::string & dataset)
{
    return dataset.ends_with("-synthetic");
}

ExperimentInput BuildSynthetic(const std::string & dataset,std::uint64_t seed)
{
    const bool well=dataset=="well-synthetic",weak=dataset=="weak-synthetic",
        moderate=dataset=="moderate-synthetic",noise_high=dataset=="noise-high-synthetic",
        mixed_width=dataset=="mixed-width-synthetic",ac_correlated=dataset=="ac-correlated-synthetic";
    if(!well && !weak && !moderate && !noise_high && !mixed_width && !ac_correlated)
        throw std::invalid_argument("Unknown synthetic dataset.");
    const auto positions=SyntheticAtoms();
    const std::size_t atom_count=positions.size();
    core::JointProblemInput input;
    input.atom_ids.reserve(atom_count); input.support.resize(atom_count);
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(rhbm_gem::PureHaloSupportExclusionContributorPolicy);
    for(std::size_t atom=0;atom<atom_count;++atom)
    {
        input.atom_ids.push_back("synthetic-"+std::to_string(atom));
        input.selection_domain->target_indices.push_back(atom);
    }
    std::vector<Truth> truth(atom_count);
    for(std::size_t atom=0;atom<atom_count;++atom)
    {
        if(weak) truth[atom]={1.0+0.03*static_cast<double>(atom),8.0+0.1*static_cast<double>(atom),
            -0.5-0.01*static_cast<double>(atom)};
        else if(moderate) truth[atom]={1.4+0.04*static_cast<double>(atom),
            1.8+0.04*static_cast<double>(atom),0.12+0.005*static_cast<double>(atom)};
        else if(mixed_width) truth[atom]={1.4+0.04*static_cast<double>(atom),
            0.55+0.16*static_cast<double>(atom),0.12+0.005*static_cast<double>(atom)};
        else if(ac_correlated) truth[atom]={6.0+0.1*static_cast<double>(atom),
            3.2+0.03*static_cast<double>(atom),-0.05-0.001*static_cast<double>(atom)};
        else truth[atom]={1.4+0.04*static_cast<double>(atom),0.72+0.01*static_cast<double>(atom),
            0.12+0.005*static_cast<double>(atom)};
    }
    std::vector<std::array<double,3>> rows;
    for(int z=-5;z<=5;++z) for(int y=-8;y<=8;++y) for(int x=-8;x<=8;++x)
    {
        const std::array<double,3> position{0.5*static_cast<double>(x),0.5*static_cast<double>(y),0.5*static_cast<double>(z)};
        bool supported=false;
        for(const auto & atom:positions)
        {
            const auto dx=position[0]-atom[0],dy=position[1]-atom[1],dz=position[2]-atom[2];
            supported|=dx*dx+dy*dy+dz*dz<=6.25;
        }
        if(supported) rows.push_back(position);
    }
    const double noise_sigma=weak ? 1e-5 : noise_high ? 5e-3 : 5e-4;
    std::mt19937_64 generator(seed);
    std::normal_distribution<double> noise(0.0,noise_sigma);
    for(const auto & position:rows)
    {
        const auto row=input.observations.size();
        input.row_ids.push_back(std::to_string(row));
        double value=noise(generator);
        for(std::size_t atom=0;atom<atom_count;++atom)
        {
            const auto dx=position[0]-positions[atom][0],dy=position[1]-positions[atom][1],
                dz=position[2]-positions[atom][2];
            const double square=dx*dx+dy*dy+dz*dz;
            if(square>6.25) continue;
            input.support[atom].push_back({row,square});
            const auto basis=joint::EvaluateKernel(square,truth[atom].width,2.5);
            value+=truth[atom].amplitude*basis.gaussian+truth[atom].charge*basis.charge;
        }
        input.observations.push_back(value);
    }
    input.selection_domain->pure_halo_support_exclusion=rhbm_gem::JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    std::vector<double> initial(atom_count,weak ? 0.9 : 0.81);
    return {core::JointProblem(std::move(input)),std::move(initial),std::move(truth),true,
        false,"none",1,dataset+":"+std::to_string(seed),"synthetic-truth-controlled",
        "joint_width_trust_experiment synthetic generator",
        "original-domain-A/C; physical-width-B","original-domain-A/C; physical-width-B",
        "original-domain-A/C; physical-width-B"};
}

json::value ReadJson(const std::filesystem::path & path)
{
    std::ifstream input(path);
    if(!input) throw std::runtime_error("Cannot read JSON: "+path.string());
    json::parse_options options; options.numbers=json::number_precision::precise;
    return json::parse(std::string(std::istreambuf_iterator<char>(input),{}),{},options);
}

std::unordered_map<std::string,Truth> ReadTruth(const std::filesystem::path & path)
{
    const auto metadata=ReadJson(path);
    std::unordered_map<std::string,Truth> truth;
    for(const auto & atom:metadata.at("atoms").as_array())
    {
        const auto id=std::to_string(json::value_to<std::uint64_t>(atom.at("serial_id")));
        const auto width=json::value_to<double>(atom.at("effective_gaussian_width"));
        const auto charge_width=json::value_to<double>(atom.at("effective_charge_width"));
        if(width!=charge_width) throw std::runtime_error("Simulation truth has split widths at "+id);
        truth.emplace(id,Truth{static_cast<double>(json::value_to<int>(atom.at("element"))),width,
            json::value_to<double>(atom.at("charge_used"))});
    }
    return truth;
}

ExperimentInput BuildReal(const std::string & dataset,const std::filesystem::path & model_path,
    const std::filesystem::path & map_path)
{
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path);
    const auto normalization=map->MapValueArrayNormalization();
    if(!normalization) throw std::runtime_error("Map normalization unavailable.");
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    model->EditAnalysis().InitializeFromSelection();
    auto problem=core::BuildPureHaloSupportExclusionProblem(*map,*model);
    const auto initialization=joint::InitializeContributors(*map,*model,problem);
    if(!initialization.valid) throw std::runtime_error("First-stage initialization failed: "+initialization.reason);
    const auto & data=core::JointProblemAccess::Get(problem);
    std::vector<Truth> truth(data.input->atom_ids.size());
    bool truth_available=dataset=="6z6u-simulation";
    if(truth_available)
    {
        const auto values=ReadTruth(map_path.string()+".simulation.json");
        for(std::size_t atom=0;atom<truth.size();++atom)
        {
            const auto found=values.find(data.input->atom_ids[atom]);
            if(found==values.end()) throw std::runtime_error("Missing simulation truth for "+data.input->atom_ids[atom]);
            truth[atom]=found->second;
        }
    }
    const bool normalized=normalization.has_value();
    const std::string normalized_scale="normalized-domain-A/C; physical-width-B";
    return {std::move(problem),initialization.b,std::move(truth),truth_available,
        normalized,"MapValueArrayNormalization",*normalization,map_path.string(),
        truth_available ? "normalized-simulation-research" : "real-map-production-normalized",
        truth_available ? map_path.string()+".simulation.json" : "none",
        truth_available ? "original-domain-A/C; physical-width-B" : "N/A — no truth",
        normalized_scale,truth_available ? normalized_scale : "N/A — no truth"};
}

struct RunMetrics
{
    trust::PolicyStatistics policy_statistics;
    std::size_t policy_rejections{},accepted_width_updates{},sweeps{},profile_evaluations{};
    double maximum_width{joint::unavailable},maximum_condition{joint::unavailable},
        maximum_cancellation{joint::unavailable},minimum_projected_width_singular_value{joint::unavailable};
    double maximum_observed_width{joint::unavailable},maximum_observed_condition{joint::unavailable},
        maximum_observed_cancellation{joint::unavailable},
        minimum_observed_projected_width_singular_value{joint::unavailable},
        maximum_observed_abs_ac_basis_correlation{joint::unavailable};
    std::optional<std::size_t> first_rejection_sweep,first_rejection_block;
    json::array objective_trajectory,block_update_trace;
};

void Max(double & value,double candidate)
{
    if(std::isfinite(candidate) && (!std::isfinite(value) || candidate>value)) value=candidate;
}

void Min(double & value,double candidate)
{
    if(std::isfinite(candidate) && (!std::isfinite(value) || candidate<value)) value=candidate;
}

double Rmse(const std::vector<double> & values)
{
    if(values.empty()) return joint::unavailable;
    long double sum{};
    for(const auto value:values) sum+=static_cast<long double>(value)*value;
    return std::sqrt(static_cast<double>(sum/static_cast<long double>(values.size())));
}

json::object RunOne(const ExperimentInput & experiment,trust::Policy policy_kind)
{
    const auto & data=core::JointProblemAccess::Get(experiment.problem);
    const auto & policy_name=trust::Name(policy_kind);
    const bool synthetic=IsSyntheticDataset(experiment.source.substr(0,experiment.source.find(':')));
    std::vector<std::string> policy_target_ids;
    if(synthetic)
        for(const auto atom:data.input->selection_domain->target_indices)
            policy_target_ids.push_back(data.input->atom_ids.at(atom));
    else
        policy_target_ids.emplace_back("1032");
    const auto policy_scope=trust::LocateMutableCoreScope(experiment.problem,policy_target_ids,12);
    joint::FixedNeighborWidthTrajectoryCapture capture;
    capture.maximum_events=100000;
    capture.maximum_atoms=32;
    capture.maximum_rows=4096;
    capture.accepted_updates_only=policy_kind==trust::Policy::Baseline;
    capture.policy_evidence_only=policy_kind!=trust::Policy::Baseline;
    std::vector<std::size_t> tracked_indices;
    if(synthetic)
        tracked_indices=data.input->selection_domain->target_indices;
    else
    {
        const auto found=std::find(data.input->atom_ids.begin(),data.input->atom_ids.end(),"1032");
        if(found==data.input->atom_ids.end()) throw std::runtime_error("Serial 1032 is not a Joint contributor.");
        tracked_indices.push_back(static_cast<std::size_t>(found-data.input->atom_ids.begin()));
    }
    for(const auto atom:tracked_indices)
        capture.atom_initial_widths.emplace_back(data.input->atom_ids.at(atom),experiment.initial_widths.at(atom));

    RunMetrics metrics;
    joint::FixedNeighborDiagnosticSink sink;
    sink.block=[&](const joint::FixedNeighborBlockRecord & block) {
        metrics.block_update_trace.push_back(json::object{
            {"component_index",block.component_index},{"component_id",block.component_id},
            {"sweep",block.sweep},{"block",block.block},{"accepted",block.accepted},
            {"accepted_updates",block.accepted_updates},{"objective_before",Number(block.objective_before)},
            {"objective_after",Number(block.objective_after)},{"status",block.status}});
        for(const auto & trial:block.profile_trials)
        {
            if(trial.width_trust_policy_rejected)
            {
                ++metrics.policy_rejections;
                if(!metrics.first_rejection_sweep)
                {
                    metrics.first_rejection_sweep=block.sweep;
                    metrics.first_rejection_block=block.block;
                }
            }
            if(trial.accepted)
            {
                for(const auto & atom:trial.tracked_atoms)
                {
                    Max(metrics.maximum_observed_width,atom.accepted_width);
                    Max(metrics.maximum_observed_condition,atom.local_design_condition);
                    Max(metrics.maximum_observed_cancellation,atom.ac_cancellation_ratio);
                    Max(metrics.maximum_observed_abs_ac_basis_correlation,std::abs(atom.basis_correlation));
                }
                for(const auto value:trial.projected_width_singular_values) Min(
                    metrics.minimum_observed_projected_width_singular_value,value);
                if(block.accepted && !trial.width_trust_policy_rejected)
                {
                    if(trial.accepted_update.value_or(0)>0) ++metrics.accepted_width_updates;
                    for(const auto & atom:trial.tracked_atoms)
                    {
                        Max(metrics.maximum_width,atom.accepted_width);
                        Max(metrics.maximum_condition,atom.local_design_condition);
                        Max(metrics.maximum_cancellation,atom.ac_cancellation_ratio);
                        Max(metrics.maximum_observed_abs_ac_basis_correlation,std::abs(atom.basis_correlation));
                    }
                    for(const auto value:trial.projected_width_singular_values) Min(
                        metrics.minimum_projected_width_singular_value,value);
                }
            }
        }
    };
    sink.sweep=[&](const joint::FixedNeighborBlockSweep & sweep) {
        metrics.sweeps=sweep.sweep;
        metrics.objective_trajectory.push_back(json::object{
            {"sweep",sweep.sweep},{"objective_before",Number(sweep.objective_before)},
            {"objective_after",Number(sweep.objective_after)},
            {"global_ac_kkt",Number(sweep.global_ac_kkt)},
            {"global_width_gradient_inf_norm",Number(sweep.global_width_gradient_inf_norm)},
            {"eta_change_inf",Number(sweep.eta_change_inf)}});
    };
    joint::FixedNeighborSearchPolicy production_policy;
    production_policy.core_atoms=12; production_policy.maximum_sweeps=200;
    production_policy.order=joint::FixedNeighborBlockOrder::Forward;
    production_policy.configure_search_for_testing=[&](joint::FixedNeighborPolicy & local) {
        local.capture_local_trajectory=true;
        local.collect_telemetry=true;
        local.diagnostic_sink=&sink;
        local.width_trajectory=&capture;
        local.width_trial_trust_for_testing=trust::Make(policy_kind,policy_scope,
            local.component_index,&metrics.policy_statistics);
    };
    const auto started=Clock::now();
    const auto fit=joint::FitFixedNeighborComponentsImpl(experiment.problem,experiment.initial_widths,production_policy);
    const double wall_seconds=std::chrono::duration<double>(Clock::now()-started).count();
    if(fit.components.size()!=1) throw std::runtime_error("P2 runner expects one component.");
    if(metrics.policy_rejections!=metrics.policy_statistics.in_scope_rejected_trials)
        throw std::runtime_error("Width-trust policy rejection audit does not match block records.");
    if(metrics.policy_statistics.out_of_scope_interventions!=0)
        throw std::runtime_error("Width-trust policy intervened outside its declared scope.");
    const auto & component=fit.components.front();
    metrics.profile_evaluations=static_cast<std::size_t>(component.profile_evaluations);
    const auto & diagnostic=component.candidate_diagnostic;
    std::vector<double> parameter_error_a,parameter_error_b,parameter_error_c,prediction_error;
    bool candidate_available=diagnostic && diagnostic->candidate_available &&
        diagnostic->candidate_eta.size()==data.input->atom_ids.size() &&
        diagnostic->candidate_ac.size()==2*data.input->atom_ids.size();
    if(candidate_available)
    {
        prediction_error.assign(data.input->observations.size(),0.0);
        for(std::size_t atom=0;atom<data.input->atom_ids.size();++atom)
        {
            const auto a=diagnostic->candidate_ac.at(2*atom),c=diagnostic->candidate_ac.at(2*atom+1);
            const auto b=std::exp(diagnostic->candidate_eta.at(atom));
            for(const auto & support:data.input->support.at(atom))
            {
                const auto basis=joint::EvaluateKernel(support.squared_distance,b,2.5);
                prediction_error.at(support.row)+=a*basis.gaussian+c*basis.charge;
            }
            if(experiment.truth_available)
            {
                const auto truth=experiment.normalization_applied ?
                    trust::ToNormalizedDomain(experiment.truth.at(atom),experiment.normalization_divisor) :
                    experiment.truth.at(atom);
                parameter_error_a.push_back(a-truth.amplitude);
                parameter_error_b.push_back(b-truth.width);
                parameter_error_c.push_back(c-truth.charge);
            }
        }
        for(std::size_t row=0;row<prediction_error.size();++row)
            prediction_error.at(row)-=data.input->observations.at(row);
    }
    const auto independent_replay=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="independent-replay";
    });
    const auto runtime=fit.RuntimeConvergence();
    const std::string accuracy_metric_status=!experiment.truth_available ? "N/A — no truth" :
        runtime==rhbm_gem::JointCheckStatus::Passed ? "qualified-estimate" : "unqualified-candidate-metrics";
    const std::string prediction_metric_status=!candidate_available ? "Unavailable" :
        runtime==rhbm_gem::JointCheckStatus::Passed ? "qualified-fit" : "unqualified-candidate";
    json::array policy_scope_json;
    for(const auto & block:policy_scope)
    {
        json::array atom_ids;
        for(const auto & atom_id:block.core_atom_ids) atom_ids.emplace_back(atom_id);
        policy_scope_json.emplace_back(json::object{{"component_index",block.component_index},
            {"component_id",block.component_id},{"block_index_zero_based",block.block_index},
            {"block_number_one_based",block.block_index+1},{"mutable_core_atom_ids",std::move(atom_ids)}});
    }
    json::object unavailable_by_reason,rejections_by_reason,rejections_by_block;
    for(const auto & [reason,count]:metrics.policy_statistics.evidence_unavailable_by_reason)
        unavailable_by_reason[reason]=count;
    for(const auto & [reason,count]:metrics.policy_statistics.rejections_by_reason)
        rejections_by_reason[reason]=count;
    for(const auto & [block,count]:metrics.policy_statistics.rejections_by_block)
        rejections_by_block[block]=count;
    return {
        {"strategy",policy_name},{"dataset",experiment.source},{"source",experiment.source},
        {"simulation_mode",experiment.simulation_mode},
        {"normalization",json::object{{"applied",experiment.normalization_applied},
            {"operation",experiment.normalization_operation},
            {"divisor",Number(experiment.normalization_divisor)}}},
        {"truth_source",experiment.truth_source},{"truth_scale",experiment.truth_scale},
        {"estimate_scale",experiment.estimate_scale},{"comparison_scale",experiment.comparison_scale},
        {"seed",experiment.source.find(':')==std::string::npos ? json::value(nullptr) :
            json::value(experiment.source.substr(experiment.source.find(':')+1))},
        {"target_count",data.input->selection_domain->target_indices.size()},
        {"observation_count",data.input->observations.size()},
        {"truth_available",experiment.truth_available},
        {"policy_scope",std::move(policy_scope_json)},
        {"diagnostic_scope",synthetic ? "all-selected-targets" : "serial-1032-only"},
        {"runtime_convergence",Status(runtime)},
        {"component_qualification",Status(component.RuntimeConvergence())},
        {"global_qualification",Status(fit.RuntimeConvergence())},
        {"search_completed",fit.search_completed},
        {"component_stop_reason",component.stop_reason},
        {"independent_replay",independent_replay==fit.evidence.end() ? "NotMeasured" : Status(independent_replay->status)},
        {"candidate_available",candidate_available},
        {"accuracy_metric_status",accuracy_metric_status},
        {"prediction_metric_status",prediction_metric_status},
        {"candidate_objective",diagnostic && diagnostic->candidate_objective ? Number(*diagnostic->candidate_objective) : json::value(nullptr)},
        {"candidate_ac_kkt",diagnostic && diagnostic->candidate_ac_kkt ? Number(*diagnostic->candidate_ac_kkt) : json::value(nullptr)},
        {"candidate_width_gradient_inf_norm",diagnostic && diagnostic->candidate_width_gradient_inf_norm ? Number(*diagnostic->candidate_width_gradient_inf_norm) : json::value(nullptr)},
        {"A_rmse",experiment.truth_available ? Number(Rmse(parameter_error_a)) : json::value(nullptr)},
        {"B_rmse_angstrom",experiment.truth_available ? Number(Rmse(parameter_error_b)) : json::value(nullptr)},
        {"C_rmse",experiment.truth_available ? Number(Rmse(parameter_error_c)) : json::value(nullptr)},
        {"prediction_rmse",candidate_available ? Number(Rmse(prediction_error)) : json::value(nullptr)},
        {"maximum_accepted_width_A",Number(metrics.maximum_width)},
        {"maximum_local_ac_condition",Number(metrics.maximum_condition)},
        {"maximum_ac_cancellation_ratio",Number(metrics.maximum_cancellation)},
        {"minimum_accepted_projected_width_singular_value",Number(metrics.minimum_projected_width_singular_value)},
        {"maximum_observed_width_A",Number(metrics.maximum_observed_width)},
        {"maximum_observed_local_ac_condition",Number(metrics.maximum_observed_condition)},
        {"maximum_observed_ac_cancellation_ratio",Number(metrics.maximum_observed_cancellation)},
        {"maximum_observed_abs_ac_basis_correlation",Number(metrics.maximum_observed_abs_ac_basis_correlation)},
        {"minimum_observed_projected_width_singular_value",Number(
            metrics.minimum_observed_projected_width_singular_value)},
        {"extreme_width_incident",metrics.maximum_observed_width>=10000.0},
        {"policy_rejected_trials",metrics.policy_rejections},
        {"policy_statistics",json::object{
            {"in_scope_evaluated_trials",metrics.policy_statistics.in_scope_evaluated_trials},
            {"in_scope_rejected_trials",metrics.policy_statistics.in_scope_rejected_trials},
            {"out_of_scope_evaluated_trials",metrics.policy_statistics.out_of_scope_evaluated_trials},
            {"out_of_scope_interventions",metrics.policy_statistics.out_of_scope_interventions},
            {"evidence_unavailable_trials",metrics.policy_statistics.evidence_unavailable_trials},
            {"evidence_unavailable_by_reason",std::move(unavailable_by_reason)},
            {"rejections_by_reason",std::move(rejections_by_reason)},
            {"rejections_by_block",std::move(rejections_by_block)},
            {"first_rejection_sweep",metrics.policy_statistics.first_rejection_sweep ?
                json::value(*metrics.policy_statistics.first_rejection_sweep) : json::value(nullptr)},
            {"first_rejection_block",metrics.policy_statistics.first_rejection_block ?
                json::value(*metrics.policy_statistics.first_rejection_block) : json::value(nullptr)},
            {"first_rejection_component",metrics.policy_statistics.first_rejection_component ?
                json::value(*metrics.policy_statistics.first_rejection_component) : json::value(nullptr)}}},
        {"diagnostic_capture_limits",json::object{
            {"maximum_events",capture.maximum_events},{"maximum_atoms",capture.maximum_atoms},
            {"maximum_rows",capture.maximum_rows},{"events_emitted",capture.events_emitted},
            {"events_skipped_due_to_event_limit",capture.events_skipped_due_to_event_limit},
            {"events_skipped_due_to_atom_limit",capture.events_skipped_due_to_atom_limit},
            {"events_unavailable_due_to_row_limit",capture.events_unavailable_due_to_row_limit},
            {"events_unavailable_due_to_spectrum",capture.events_unavailable_due_to_spectrum}}},
        {"first_policy_rejection_sweep",metrics.first_rejection_sweep ? json::value(*metrics.first_rejection_sweep) : json::value(nullptr)},
        {"first_policy_rejection_block",metrics.first_rejection_block ? json::value(*metrics.first_rejection_block) : json::value(nullptr)},
        {"accepted_width_updates",metrics.accepted_width_updates},
        {"search_sweeps",metrics.sweeps},
        {"profile_evaluations",metrics.profile_evaluations},
        {"objective_trajectory",std::move(metrics.objective_trajectory)},
        {"block_update_trace",std::move(metrics.block_update_trace)},
        {"wall_seconds",wall_seconds},
        {"process_peak_rss_bytes",PeakRssBytes()}
    };
}

void Run(const std::string & dataset,const std::filesystem::path & output,const std::string & policy_name,
    std::uint64_t seed,const std::filesystem::path & model_path,const std::filesystem::path & map_path)
{
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a fresh P2 output directory.");
    std::filesystem::create_directories(output);
    const auto policy=policy_name=="A" ? trust::Policy::Baseline :
        policy_name=="B1" ? trust::Policy::ConditioningAware :
        policy_name=="B2" ? trust::Policy::ConditioningAndCancellationAware :
        throw std::invalid_argument("Policy must be A, B1, or B2.");
    ExperimentInput experiment=IsSyntheticDataset(dataset) ?
        BuildSynthetic(dataset,seed) : BuildReal(dataset,model_path,map_path);
    Write(output/"result.json",RunOne(experiment,policy));
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=5 && argc!=7)
            throw std::invalid_argument("Usage: joint_width_trust_experiment DATASET OUTPUT A|B1|B2 SEED [MODEL MAP]");
        const auto dataset=std::string(argv[1]);
        const auto seed=Seed(argv[4]);
        Run(dataset,argv[2],argv[3],seed,argc==7 ? std::filesystem::path(argv[5]) : std::filesystem::path{},
            argc==7 ? std::filesystem::path(argv[6]) : std::filesystem::path{});
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
