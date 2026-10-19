#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "support/JointReproducibilityFingerprint.hpp"

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace core=rhbm_gem::core;
namespace json=boost::json;

void Write(const std::filesystem::path & path,const json::value & value)
{
    std::ofstream output(path);
    output.exceptions(std::ios::failbit|std::ios::badbit);
    output<<json::serialize(value)<<'\n';
}

const char * Status(rhbm_gem::JointCheckStatus status);

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

json::object ParameterMetrics(const core::JointFitResult & fit,const core::JointProblem & problem,
    const std::filesystem::path & map_path,double normalization_divisor)
{
    json::object out;
    if(!fit.assembled_state || !fit.assembled_state->b.size()) return out;
    const auto & input=core::JointProblemAccess::Get(problem).input;
    const auto calibration=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="global-width-calibration";
    });
    if(calibration==fit.evidence.end() || !calibration->threshold) return out;
    const double bmax=*calibration->threshold,eta_max=std::log(bmax);
    double min_width=std::numeric_limits<double>::infinity(),max_width=0.0,maximum_violation{};
    std::size_t active{};
    double serial_1032_width{};
    for(std::size_t atom=0;atom<fit.assembled_state->b.size();++atom)
    {
        const double width=fit.assembled_state->b[atom];
        const double eta=fit.assembled_state->log_b[atom];
        min_width=std::min(min_width,width); max_width=std::max(max_width,width);
        maximum_violation=std::max(maximum_violation,std::max(0.0,eta-eta_max));
        active+=eta>=eta_max-64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(eta_max));
        if(input->atom_ids[atom]=="1032") serial_1032_width=width;
    }
    out["global_bmax_angstrom"]=bmax;
    out["minimum_width_angstrom"]=min_width;
    out["maximum_width_angstrom"]=max_width;
    out["serial_1032_width_angstrom"]=serial_1032_width>0.0 ? json::value(serial_1032_width) : json::value(nullptr);
    out["active_bound_atoms"]=active;
    out["maximum_bound_violation_log_width"]=maximum_violation;
    out["truth_used_for_calibration"]=false;

    const auto truth_path=std::filesystem::path(map_path.string()+".simulation.json");
    if(!std::filesystem::exists(truth_path)) return out;
    std::ifstream truth_file(truth_path);
    json::parse_options parse_options; parse_options.numbers=json::number_precision::precise;
    const auto truth_root=json::parse(std::string(std::istreambuf_iterator<char>(truth_file),{}),{},parse_options);
    struct Truth {double amplitude{},width{},charge{};};
    std::unordered_map<std::string,Truth> truth_by_id;
    for(const auto & atom:truth_root.at("atoms").as_array())
        truth_by_id.emplace(std::to_string(json::value_to<std::uint64_t>(atom.at("serial_id"))),
            Truth{static_cast<double>(json::value_to<int>(atom.at("element"))),
                json::value_to<double>(atom.at("effective_gaussian_width")),
                json::value_to<double>(atom.at("charge_used"))});
    long double sum_a{},sum_b{},sum_c{};
    std::size_t covered{};
    for(std::size_t atom=0;atom<fit.assembled_state->b.size();++atom)
    {
        const auto & truth=truth_by_id.at(input->atom_ids[atom]);
        const double da=fit.assembled_state->ac[2*atom]*normalization_divisor-truth.amplitude;
        const double db=fit.assembled_state->b[atom]-truth.width;
        const double dc=fit.assembled_state->ac[2*atom+1]*normalization_divisor-truth.charge;
        sum_a+=static_cast<long double>(da)*da;
        sum_b+=static_cast<long double>(db)*db;
        sum_c+=static_cast<long double>(dc)*dc;
        covered+=truth.width<=bmax;
    }
    const double count=static_cast<double>(fit.assembled_state->b.size());
    out["truth_coverage"]=static_cast<double>(covered)/count;
    out["truth_outside_bound_atoms"]=fit.assembled_state->b.size()-covered;
    out["A_rmse"]=std::sqrt(static_cast<double>(sum_a/count));
    out["B_rmse_angstrom"]=std::sqrt(static_cast<double>(sum_b/count));
    out["C_rmse"]=std::sqrt(static_cast<double>(sum_c/count));
    return out;
}

json::object TargetQualityMetrics(const core::JointFitResult & fit)
{
    std::map<std::string,std::size_t> quality_counts,quality_reason_counts,stability_counts,conditioning_counts;
    std::size_t available{};
    if(fit.target_evidence)
        for(const auto & point:fit.target_evidence->point_atoms)
        {
            ++quality_counts[point.point_quality];
            ++quality_reason_counts[point.reason];
            ++stability_counts[point.parameter_stability];
            ++conditioning_counts[point.local_conditioning];
            available+=point.point_quality=="Available";
        }
    const auto Counts=[](const std::map<std::string,std::size_t> & counts) {
        json::object result;
        for(const auto & [name,count]:counts) result[name]=count;
        return result;
    };
    json::array component_stability;
    for(const auto & component:fit.components)
    {
        if(!component.parameter_stability)
        {
            component_stability.push_back(json::object{{"component_id",component.id},
                {"perturbation_status","NotAssessed"},{"trajectory_status","NotAssessed"}});
            continue;
        }
        std::map<std::string,std::size_t> block_counts,block_reason_counts;
        json::array block_audit;
        std::size_t selected_probe_count{},successfully_assessed_count{};
        for(const auto & block:component.parameter_stability->blocks)
        {
            ++block_counts[block.status];
            ++block_reason_counts[block.reason];
            selected_probe_count+=block.probe_selected ? 1u : 0u;
            successfully_assessed_count+=(block.status=="Stable" || block.status=="Sensitive" ||
                block.status=="Unstable") ? 1u : 0u;
            const auto number=[](const std::optional<double> & value)->json::value {
                return value ? json::value(*value) : json::value(nullptr);
            };
            const auto size=[](const std::optional<std::size_t> & value)->json::value {
                return value ? json::value(*value) : json::value(nullptr);
            };
            const auto conditioning=std::find_if(component.local_conditioning.begin(),component.local_conditioning.end(),
                [&](const auto & value){return value.block_id==block.block_id;});
            json::array probe_audit;
            for(const auto & probe:block.probes)
                probe_audit.emplace_back(json::object{{"direction",probe.direction},{"status",probe.status},
                    {"minus_seed_status",probe.minus_seed_status},{"minus_seed_reason",probe.minus_seed_reason},
                    {"minus_requested_displacement_inf",number(probe.minus_requested_displacement_inf)},
                    {"minus_effective_displacement_inf",number(probe.minus_effective_displacement_inf)},
                    {"plus_seed_status",probe.plus_seed_status},{"plus_seed_reason",probe.plus_seed_reason},
                    {"plus_requested_displacement_inf",number(probe.plus_requested_displacement_inf)},
                    {"plus_effective_displacement_inf",number(probe.plus_effective_displacement_inf)}});
            block_audit.emplace_back(json::object{{"block_id",block.block_id},{"status",block.status},
                {"reason",block.reason},{"probe_selected",block.probe_selected},
                {"baseline_local_ac_kkt",number(block.baseline_local_ac_kkt)},
                {"baseline_projected_width_kkt",number(block.baseline_projected_width_kkt)},
                {"baseline_width_bound_slack",number(block.baseline_width_bound_slack)},
                {"baseline_local_objective",number(block.baseline_local_objective)},
                {"baseline_local_residual_rmse_normalized",number(block.baseline_local_residual_rmse_normalized)},
                {"baseline_active_ac_constraints",size(block.baseline_active_ac_constraints)},
                {"baseline_active_width_bound_targets",size(block.baseline_active_width_bound_targets)},
                {"core_atoms",json::value_from(block.core_atoms)},
                {"local_conditioning_classification",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : json::value(conditioning->conditioning_classification)},
                {"local_free_ac_rank",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : size(conditioning->free_ac_rank)},
                {"local_free_ac_columns",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : json::value(conditioning->free_ac_columns)},
                {"local_free_ac_condition",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : number(conditioning->free_ac_condition)},
                {"local_projected_width_condition",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : number(conditioning->projected_width_condition)},
                {"local_profile_condition",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : number(conditioning->profile_condition)},
                {"local_active_amplitude_count",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : json::value(conditioning->active_amplitude_count)},
                {"local_active_width_bound_targets",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : size(conditioning->active_width_bound_targets)},
                {"local_ac_projected_kkt",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : number(conditioning->local_ac_projected_kkt)},
                {"local_projected_width_kkt",conditioning==component.local_conditioning.end() ?
                    json::value(nullptr) : number(conditioning->local_projected_width_kkt)},
                {"baseline_independent_replay_status",block.baseline_independent_replay_status},
                {"baseline_independent_replay_reason",block.baseline_independent_replay_reason},
                {"baseline_numerical_validity_status",block.baseline_numerical_validity_status},
                {"baseline_conditional_optimality_status",block.baseline_conditional_optimality_status},
                {"baseline_eligibility_status",block.baseline_eligibility_status},
                {"baseline_eligibility_reason",block.baseline_eligibility_reason},
                {"conditional_reprofile_status",block.conditional_reprofile_status},
                {"conditional_reprofile_local_ac_kkt",number(block.conditional_reprofile_local_ac_kkt)},
                {"conditional_reprofile_projected_width_kkt",number(block.conditional_reprofile_projected_width_kkt)},
                {"conditional_reprofile_objective_difference",number(block.conditional_reprofile_objective_difference)},
                {"conditional_reprofile_scaled_ac_parameter_difference",
                    number(block.conditional_reprofile_scaled_ac_parameter_difference)},
                {"conditional_reprofile_prediction_change_relative",
                    number(block.conditional_reprofile_prediction_change_relative)},
                {"conditional_reprofile_independent_replay_status",
                    block.conditional_reprofile_independent_replay_status},
                {"conditional_reprofile_independent_replay_reason",
                    block.conditional_reprofile_independent_replay_reason},
                {"probes",std::move(probe_audit)},
                {"cpu_seconds",block.cpu_seconds},{"perturbation_solves",block.perturbation_solves},
                {"perturbation_evaluations",block.perturbation_evaluations}});
        }
        component_stability.push_back(json::object{{"component_id",component.id},
            {"perturbation_status",component.parameter_stability->perturbation_status},
            {"trajectory_status",component.parameter_stability->trajectory_stability_status},
            {"blocks",component.parameter_stability->blocks.size()},
            {"block_status_counts",Counts(block_counts)},
            {"block_reason_counts",Counts(block_reason_counts)},
            {"selected_probe_blocks",selected_probe_count},
            {"successfully_assessed_blocks",successfully_assessed_count},
            {"block_audit",std::move(block_audit)},
            {"perturbation_solves",component.parameter_stability->perturbation_solves},
            {"perturbation_evaluations",component.parameter_stability->perturbation_evaluations},
            {"cpu_seconds",component.parameter_stability->cpu_seconds},
            {"wall_seconds",component.parameter_stability->seconds},
            {"workspace_peak_memory_estimate_bytes",
                component.parameter_stability->workspace_peak_memory_estimate_bytes ?
                    json::value(*component.parameter_stability->workspace_peak_memory_estimate_bytes) :
                    json::value(nullptr)}});
    }
    return {{"target_runtime_convergence",Status(fit.TargetRuntimeConvergence())},
        {"target_count",fit.target_evidence ? fit.target_evidence->point_atoms.size() : 0},
        {"point_quality_available",available},{"point_quality_counts",Counts(quality_counts)},
        {"point_quality_reason_counts",Counts(quality_reason_counts)},
        {"parameter_stability_counts",Counts(stability_counts)},
        {"local_conditioning_counts",Counts(conditioning_counts)},
        {"components",std::move(component_stability)}};
}

const char * Status(rhbm_gem::JointCheckStatus status)
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

void Run(const std::filesystem::path & model_path,const std::filesystem::path & map_path,
    const std::filesystem::path & output)
{
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a fresh reproducibility output directory.");
    std::filesystem::create_directories(output);
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path);
    const auto normalization=map->MapValueArrayNormalization();
    if(!normalization) throw std::runtime_error("Map normalization unavailable.");
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    model->EditAnalysis().InitializeFromSelection();
    const auto started=std::chrono::steady_clock::now();
    const auto problem=core::BuildPureHaloSupportExclusionProblem(*map,*model);
    auto initialization=core::joint_component::InitializeContributors(*map,*model,problem);
    auto fit=core::detail::FitJointComponentsInitializedBoundedWithProgress(problem,initialization,{});
    if(!fit.problem) throw std::runtime_error("Production fit did not retain its problem snapshot.");
    const auto & data=core::JointProblemAccess::Get(*fit.problem);
    const auto fingerprints=rhbm_gem::joint_reproducibility::WriteFingerprint(output,
        model_path,map_path,*normalization,*map,*model,*fit.problem,fit.initialization.b,
        {"production-bounded-candidate-no-instrumentation",false});
    json::array components;
    for(const auto & component:fit.components)
    {
        const auto bounded_lm=std::find_if(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
            return check.name=="active-set-lm-subproblems";
        });
        components.push_back(json::object{
            {"component_id",component.id},
            {"search_sweeps",component.search_sweeps},
            {"profile_evaluations",component.profile_evaluations},
            {"accepted_updates",component.accepted_updates},
            {"search_completed",component.search_completed},
            {"stop_reason",component.stop_reason},
            {"runtime_convergence",Status(component.RuntimeConvergence())},
            {"candidate_objective",component.candidate_diagnostic && component.candidate_diagnostic->candidate_objective ?
                json::value(*component.candidate_diagnostic->candidate_objective) : json::value(nullptr)},
            {"candidate_ac_kkt",component.candidate_diagnostic && component.candidate_diagnostic->candidate_ac_kkt ?
                json::value(*component.candidate_diagnostic->candidate_ac_kkt) : json::value(nullptr)},
            {"candidate_width_gradient_inf_norm",component.candidate_diagnostic && component.candidate_diagnostic->candidate_width_gradient_inf_norm ?
                json::value(*component.candidate_diagnostic->candidate_width_gradient_inf_norm) : json::value(nullptr)},
            {"active_set_lm_status",bounded_lm==component.evidence.end() ? "NotRun" : Status(bounded_lm->status)},
            {"active_set_lm_kkt",bounded_lm!=component.evidence.end() && bounded_lm->value ?
                json::value(*bounded_lm->value) : json::value(nullptr)},
            {"active_set_lm_work",bounded_lm==component.evidence.end() ? "not-recorded" : bounded_lm->reason}});
    }
    json::array evidence;
    for(const auto & check:fit.evidence)
        evidence.push_back(json::object{{"name",check.name},{"status",Status(check.status)},
            {"scope",static_cast<int>(check.scope)},
            {"value",check.value && std::isfinite(*check.value) ? json::value(*check.value) : json::value(nullptr)},
            {"threshold",check.threshold && std::isfinite(*check.threshold) ?
                json::value(*check.threshold) : json::value(nullptr)},
            {"reason",check.reason}});
    double prediction_rmse=std::numeric_limits<double>::quiet_NaN();
    if(fit.prediction && fit.prediction->size()==data.input->observations.size())
    {
        long double squared_error{};
        for(std::size_t row=0;row<fit.prediction->size();++row)
        {
            const double error=fit.prediction->at(row)-data.input->observations[row];
            squared_error+=static_cast<long double>(error)*error;
        }
        prediction_rmse=std::sqrt(static_cast<double>(squared_error/
            static_cast<long double>(fit.prediction->size())));
    }
    Write(output/"summary.json",json::object{
        {"execution_role","production-bounded-candidate-no-instrumentation"},
        {"instrumentation_enabled",false},
        {"target_count",data.input->selection_domain->target_indices.size()},
        {"observation_count",data.input->observations.size()},
        {"formal_runtime_convergence",Status(fit.RuntimeConvergence())},
        {"target_quality",TargetQualityMetrics(fit)},
        {"search_completed",fit.search_completed},
        {"objective",fit.objective ? json::value(*fit.objective) : json::value(nullptr)},
        {"prediction_rmse",std::isfinite(prediction_rmse) ? json::value(prediction_rmse) : json::value(nullptr)},
        {"initialization_valid",fit.initialization.valid},
        {"initialization_reason",fit.initialization.reason},
        {"components",std::move(components)},
        {"evidence",std::move(evidence)},
        {"parameter_metrics",ParameterMetrics(fit,*fit.problem,map_path,*normalization)},
        {"fingerprints",fingerprints},
        {"peak_rss_bytes",PeakRssBytes()},
        {"wall_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()}});
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=4) throw std::invalid_argument("Usage: joint_width_repro_no_instrumentation MODEL MAP OUTPUT_DIR");
        Run(argv[1],argv[2],argv[3]);
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
