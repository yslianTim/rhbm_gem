#include "support/JointFixedNeighborDiagnostic.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif
#if defined(__APPLE__)
#include <mach/mach.h>
#include <mach/task_info.h>
#endif
#if defined(__linux__)
#include <unistd.h>
#endif

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
using Input=rhbm_gem::core::JointProblemInput;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
std::optional<std::size_t> PeakRssBytes()
{
 #if defined(__APPLE__) || defined(__linux__)
    rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) return {};
#if defined(__APPLE__)
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
 #else
    return {};
 #endif
}
std::optional<std::size_t> CurrentRssBytes()
{
#if defined(__APPLE__)
    mach_task_basic_info_data_t info{}; mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
    if(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)!=KERN_SUCCESS) return {};
    return static_cast<std::size_t>(info.resident_size);
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm"); std::size_t total_pages{},resident_pages{};
    const auto page_bytes=sysconf(_SC_PAGESIZE);
    if(!(statm>>total_pages>>resident_pages) || page_bytes<=0 ||
        resident_pages>std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(page_bytes)) return {};
    return resident_pages*static_cast<std::size_t>(page_bytes);
#else
    return {};
#endif
}
const rhbm_gem::JointCheck * FindCheck(const rhbm_gem::core::JointFitResult & fit,const std::string & name)
{
    const auto found=std::find_if(fit.evidence.begin(),fit.evidence.end(),[&](const auto & check) {
        return check.name==name && check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    return found==fit.evidence.end() ? nullptr : &*found;
}
std::string CheckName(rhbm_gem::JointCheckStatus status)
{
    using S=rhbm_gem::JointCheckStatus;
    switch(status) {case S::Passed:return "Passed"; case S::Failed:return "Failed";
        case S::Unavailable:return "Unavailable"; case S::NotRun:return "NotRun";}
    return "NotRun";
}
j::value CheckValue(const rhbm_gem::core::JointFitResult & fit,const std::string & name)
{
    const auto * check=FindCheck(fit,name);
    return check ? j::value{{"status",CheckName(check->status)},
        {"value",check->value ? j::value(*check->value) : j::value(nullptr)},
        {"threshold",check->threshold ? j::value(*check->threshold) : j::value(nullptr)}} : j::value(nullptr);
}
double WidthGradient(const rhbm_gem::core::JointFitResult & fit)
{
    if(!fit.assembled_state || fit.assembled_state->width_gradient.empty()) return std::numeric_limits<double>::infinity();
    double out{}; for(double value:fit.assembled_state->width_gradient) out=std::max(out,std::abs(value)); return out;
}
j::value Number(double value)
{return std::isfinite(value) ? j::value(value) : j::value(nullptr);}
j::value BytesValue(std::optional<std::size_t> value)
{return value ? j::value(*value) : j::value(nullptr);}
j::object MemoryEstimateJson(const n::MemoryEstimate & estimate)
{
    return {{"stage",n::MemoryStageName(estimate.stage)},
        {"status",n::MemoryStatusName(estimate.status)},
        {"reason",estimate.reason},
        {"resident_bytes",estimate.resident_bytes},
        {"temporary_bytes",estimate.temporary_bytes},
        {"estimated_peak_bytes",estimate.estimated_peak_bytes},
        {"budget_bytes",estimate.budget_bytes}};
}
j::object ResourceCensusJson(const n::FixedNeighborDiagnosticResult & result)
{
    const auto & shape=result.memory_estimate.shape;
    return {{"component_atoms",shape.component_atoms},{"component_rows",shape.component_observations},
        {"support_entries",shape.component_support_entries},{"component_blocks",shape.component_blocks},
        {"total_affected_rows",shape.affected_row_entries},{"maximum_block_rows",shape.maximum_block_rows},
        {"total_block_support_entries",shape.total_block_support_entries},
        {"maximum_block_support_entries",shape.maximum_block_support_entries},
        {"total_neighbor_atoms",shape.component_neighbor_entries},
        {"maximum_neighbor_atoms",shape.maximum_neighbor_atoms},
        {"prepared_mapping_entries",result.prepared_mapping_entries},
        {"prepared_mapping_bytes",result.prepared_mapping_bytes},
        {"estimated_preparation_peak_bytes",result.memory_estimate.estimated_peak_bytes},
        {"estimated_global_state_peak_bytes",result.global_state_memory_estimate.estimated_peak_bytes},
        {"preparation_preflight",MemoryEstimateJson(result.memory_estimate)},
        {"global_state_preflight",MemoryEstimateJson(result.global_state_memory_estimate)}};
}
j::object RssCheckpointsJson(std::optional<std::size_t> before,std::optional<std::size_t> after_preparation,
    std::optional<std::size_t> after_first_sweep,std::optional<std::size_t> after_certification)
{
    return {{"measurement",CurrentRssBytes() ? "available" : "unavailable"},
        {"rss_before_preparation",BytesValue(before)},
        {"rss_after_preparation",BytesValue(after_preparation)},
        {"rss_after_first_sweep",BytesValue(after_first_sweep)},
        {"rss_after_certification",BytesValue(after_certification)}};
}
std::size_t DerivativeTileRows()
{
    return static_cast<std::size_t>(n::derivative_tile_rows);
}
j::array ProfileTrialsJson(const std::vector<n::FixedNeighborProfileTrial> & trials)
{
    j::array out;
    for(const auto & trial:trials)
        out.push_back({{"trial_index",trial.trial_index},{"profile_evaluation",trial.profile_evaluation},
            {"accepted",trial.accepted},{"accepted_update",trial.accepted_update ?
                j::value(*trial.accepted_update) : j::value(nullptr)},
            {"local_objective_before",Number(trial.local_objective_before)},
            {"local_objective_after",Number(trial.local_objective_after)},
            {"objective_reduction",Number(trial.objective_reduction)},
            {"eta_change_inf",Number(trial.eta_change_inf)},
            {"gradient_inf_norm",Number(trial.gradient_inf_norm)},
            {"profile_seconds",Number(trial.profile_seconds)},
            {"factor_seconds",Number(trial.factor_seconds)},
            {"cumulative_factor_seconds",Number(trial.cumulative_factor_seconds)}});
    return out;
}
j::value SpectrumJson(const std::optional<n::Spectrum> & spectrum)
{
    if(!spectrum) return nullptr;
    return j::object{{"rank",spectrum->rank},{"minimum_singular_value",Number(spectrum->minimum)},
        {"rank_threshold",Number(spectrum->threshold)},{"condition_estimate",Number(spectrum->condition)}};
}
j::object AssessmentJson(const n::Assessment & assessment,const n::TrustEvidence & trust)
{
    const double correction=assessment.correction.size()==0 ? std::numeric_limits<double>::infinity() :
        assessment.correction.lpNorm<Eigen::Infinity>();
    const auto & replay=assessment.independent_replay;
    return {{"endpoint_verified",assessment.endpoint_verified},{"numerically_qualified",assessment.qualified},
        {"primary_ac_kkt_passed",replay.primary_kkt_passed},
        {"independent_ac_kkt_passed",replay.independent_kkt_passed},
        {"ac_feasible",replay.feasible},{"independent_replay_passed",replay.passed},
        {"primary_ac_kkt",Number(replay.primary_ac_kkt)},
        {"independent_ac_kkt",Number(replay.independent_ac_kkt)},
        {"independent_kkt_difference",Number(replay.kkt_difference)},
        {"primary_width_gradient_inf_norm",Number(replay.primary_width_gradient)},
        {"independent_width_gradient_inf_norm",Number(replay.independent_width_gradient)},
        {"width_gradient_difference",Number(replay.gradient_difference)},
        {"prediction_difference",Number(replay.prediction_difference)},
        {"objective_difference",Number(replay.objective_difference)},
        {"local_correction_passed",assessment.local},{"numerical_identifiability_passed",assessment.identified},
        {"failure",assessment.failure},{"verification_reason",replay.reason},
        {"primary_endpoint_seconds",Number(assessment.primary_endpoint_seconds)},
        {"independent_replay_seconds",Number(assessment.independent_replay_seconds)},
        {"rank_assessment_seconds",Number(assessment.rank_assessment_seconds)},
        {"endpoint_verification_memory",MemoryEstimateJson(assessment.verification_memory_estimate)},
        {"rank_assessment_memory",MemoryEstimateJson(assessment.rank_memory_estimate)},
        {"primary_gradient_inf_norm",assessment.primary.gradient.size() ?
            Number(assessment.primary.gradient.lpNorm<Eigen::Infinity>()) : j::value(nullptr)},
        {"local_correction_inf_norm",Number(correction)},
        {"projected_width",SpectrumJson(assessment.widths)},
        {"corrected_jacobian",SpectrumJson(assessment.jacobian)},
        {"normalized_width",SpectrumJson(assessment.normalized_widths)},
        {"endpoint_verification",j::object{{"passed",trust.passed},{"reason",trust.reason}}},
        {"reference_solve_calls",0}};
}
j::object SweepJson(const n::FixedNeighborBlockSweep & sweep)
{
    return {{"sweep",sweep.sweep},{"objective",sweep.objective_after},
        {"objective_before",sweep.objective_before},{"objective_after",sweep.objective_after},
        {"global_a_feasibility",sweep.global_a_feasibility},{"global_ac_kkt",sweep.global_ac_kkt},
        {"global_width_gradient_inf_norm",sweep.global_width_gradient_inf_norm},
        {"eta_change_inf",sweep.eta_change_inf},{"beta_scaled_change",sweep.beta_scaled_change},
        {"coordinate_confirmation_available",sweep.coordinate_confirmation_available},
        {"accepted_blocks",sweep.accepted_blocks},{"unchanged_blocks",sweep.unchanged_blocks},
        {"cache_replay_error",sweep.cache_replay_error},{"objective_replay_error",sweep.objective_replay_error},
        {"block_solves",sweep.block_solves},{"profile_evaluations",sweep.profile_evaluations},
        {"maximum_block_rows",sweep.maximum_block_rows},{"maximum_block_columns",sweep.maximum_block_columns},
        {"wall_seconds",sweep.wall_seconds}};
}
j::object ProfileWorkJson(const n::ProfileSearchWork &);
j::object WorkJson(const n::FixedNeighborWork & work)
{
    return {{"structural_preparation_seconds",Number(work.structural_preparation_seconds)},
        {"endpoint_verification_seconds",Number(work.endpoint_verification_seconds)},
        {"primary_endpoint_seconds",Number(work.primary_endpoint_seconds)},
        {"independent_replay_seconds",Number(work.independent_replay_seconds)},
        {"rank_assessment_seconds",Number(work.rank_assessment_seconds)},
        {"old_core_seconds",Number(work.old_core_seconds)},
        {"effective_response_seconds",Number(work.effective_response_seconds)},
        {"local_profile_seconds",Number(work.local_search_seconds)},
        {"local_search_seconds",Number(work.local_search_seconds)},
        {"local_state_seconds",Number(work.local_state_seconds)},
        {"candidate_copy_seconds",Number(work.candidate_copy_seconds)},
        {"candidate_replay_seconds",Number(work.candidate_replay_seconds)},
        {"incremental_candidate_replay_seconds",Number(work.incremental_candidate_replay_seconds)},
        {"incremental_candidate_replays",work.incremental_candidate_replays},
        {"incremental_acceptance_mismatches",work.incremental_acceptance_mismatches},
        {"maximum_incremental_objective_difference",Number(work.maximum_incremental_objective_difference)},
        {"maximum_incremental_relative_objective_difference",Number(work.maximum_incremental_relative_objective_difference)},
        {"maximum_incremental_residual_difference",Number(work.maximum_incremental_residual_difference)},
        {"maximum_incremental_delta_condition",Number(work.maximum_incremental_delta_condition)},
        {"last_incremental_candidate_objective",Number(work.last_incremental_candidate_objective)},
        {"last_full_candidate_objective",Number(work.last_full_candidate_objective)},
        {"maximum_cache_replay_error",Number(work.maximum_cache_replay_error)},
        {"cache_update_seconds",Number(work.cache_update_seconds)},
        {"sweep_replay_seconds",Number(work.sweep_replay_seconds)},
        {"sweep_global_state_seconds",Number(work.sweep_global_state_seconds)},
        {"maximum_block_seconds",Number(work.maximum_block_seconds)},
        {"block_visits",work.block_visits},
        {"block_solves",work.block_solves},
        {"full_candidate_replays",work.full_candidate_replays},
        {"candidate_state_full_copies",work.candidate_state_full_copies},
        {"affected_row_updates",work.affected_row_updates},
        {"old_core_basis_builds",work.old_core_basis_builds},
        {"local_problem_count",work.local_problem_count},
        {"minimum_local_atoms",work.minimum_local_atoms},
        {"mean_local_atoms",Number(work.local_problem_count ? work.total_local_atoms/
            static_cast<double>(work.local_problem_count) : 0.)},
        {"maximum_local_atoms",work.maximum_local_atoms},
        {"local_profile_work",ProfileWorkJson(work.local_profile_work)}};
}
j::object IncrementalReplayJson(const n::FixedNeighborWork & work,bool enabled)
{
    return {{"enabled",enabled},{"candidate_comparisons",work.incremental_candidate_replays},
        {"incremental_candidate_replays",work.incremental_candidate_replays},
        {"full_replay_oracle_replays",work.full_candidate_replays},
        {"acceptance_mismatches",work.incremental_acceptance_mismatches},
        {"acceptance_agreements",work.incremental_candidate_replays-work.incremental_acceptance_mismatches},
        {"incremental_seconds",Number(work.incremental_candidate_replay_seconds)},
        {"full_candidate_replay_seconds",Number(work.candidate_replay_seconds)},
        {"maximum_absolute_objective_difference",Number(work.maximum_incremental_objective_difference)},
        {"maximum_relative_objective_difference",Number(work.maximum_incremental_relative_objective_difference)},
        {"maximum_residual_difference",Number(work.maximum_incremental_residual_difference)},
        {"maximum_delta_condition",Number(work.maximum_incremental_delta_condition)},
        {"last_incremental_objective",Number(work.last_incremental_candidate_objective)},
        {"last_full_replay_objective",Number(work.last_full_candidate_objective)},
        {"maximum_cache_replay_error",Number(work.maximum_cache_replay_error)}};
}
j::object ProfileRoleWorkJson(const n::ProfileRoleWork & work)
{
    return {{"evaluations",work.evaluations},{"evaluation_seconds",Number(work.evaluation_seconds)},
        {"profile_basis_seconds",Number(work.basis_seconds)},
        {"linear_matrix_preparation_seconds",Number(work.linear_matrix_preparation_seconds)},
        {"linear_rhs_solve_seconds",Number(work.linear_rhs_solve_seconds)},
        {"linear_certificate_seconds",Number(work.linear_certificate_seconds)},
        {"derivative_preparations",work.derivative_preparations},
        {"derivative_prepare_seconds",Number(work.derivative_prepare_seconds)},
        {"derivative_raw_assembly_seconds",Number(work.derivative_raw_assembly_seconds)},
        {"derivative_free_design_assembly_seconds",Number(work.derivative_free_design_assembly_seconds)},
        {"derivative_factor_match_seconds",Number(work.derivative_factor_match_seconds)},
        {"derivative_factor_build_seconds",Number(work.derivative_factor_build_seconds)},
        {"derivative_factor_compact_seconds",Number(work.derivative_factor_compact_seconds)},
        {"derivative_rank_seconds",Number(work.derivative_rank_seconds)},
        {"derivative_least_squares_seconds",Number(work.derivative_least_squares_seconds)},
        {"derivative_normal_solve_seconds",Number(work.derivative_normal_solve_seconds)},
        {"derivative_cancellation_check_seconds",Number(work.derivative_cancellation_check_seconds)},
        {"derivative_cancellation_fallback_seconds",Number(work.derivative_cancellation_fallback_seconds)},
        {"derivative_reductions",work.derivative_reductions},
        {"derivative_reduce_seconds",Number(work.derivative_reduce_seconds)},
        {"derivative_rows_seconds",Number(work.derivative_rows_seconds)},
        {"derivative_jacobian_qr_seconds",Number(work.derivative_jacobian_qr_seconds)},
        {"derivative_norms_seconds",Number(work.derivative_norms_seconds)},
        {"derivative_outer_overhead_seconds",Number(work.derivative_outer_overhead_seconds)},
        {"tiled_qr_assembly_copy_seconds",Number(work.tiled_qr_assembly_copy_seconds)},
        {"tiled_qr_householder_seconds",Number(work.tiled_qr_householder_seconds)},
        {"tiled_qr_rhs_transform_seconds",Number(work.tiled_qr_rhs_transform_seconds)},
        {"replay_checks",work.replay_checks},{"replay_trust_seconds",Number(work.replay_trust_seconds)}};
}
j::object ProfileWorkJson(const n::ProfileSearchWork & work)
{
    return {{"total_seconds",Number(work.total_seconds)},
        {"lm_overhead_seconds",Number(work.lm_overhead_seconds)},
        {"total",ProfileRoleWorkJson(work.total)},
        {"initial_profile",ProfileRoleWorkJson(work.initial_profile)},
        {"trial_profile",ProfileRoleWorkJson(work.trial_profile)},
        {"accepted_endpoint",ProfileRoleWorkJson(work.accepted_endpoint)},
        {"reference_evaluation",ProfileRoleWorkJson(work.reference_evaluation)}};
}
void Write(const std::filesystem::path &,const j::value &);
j::object Run(const std::string & topology,int atoms,const std::filesystem::path & output_path,
    bool scaling_only=false,std::size_t maximum_sweeps=30,bool detailed=false,
    bool compare_incremental_replay=false)
{
    auto input=std::make_shared<Input>(second_stage_test::SyntheticJointWorkload(topology,atoms));
    const rhbm_gem::core::JointProblem problem(*input);
    n::Vector initial_eta= n::Vector::Constant(atoms,std::log(.55)); n::FixedNeighborPolicy neighbor_policy;
    neighbor_policy.core_atoms=12;
    neighbor_policy.maximum_sweeps=maximum_sweeps;
    neighbor_policy.capture_local_trajectory=detailed;
    neighbor_policy.collect_diagnostics=detailed;
    neighbor_policy.retain_diagnostic_records=detailed;
    neighbor_policy.compare_incremental_replay=compare_incremental_replay;
    neighbor_policy.certify_endpoint=!scaling_only && maximum_sweeps>0;
    neighbor_policy.collect_telemetry=true;
    const auto rss_before=CurrentRssBytes();
    std::optional<std::size_t> rss_after_preparation,rss_after_first_sweep,rss_after_certification;
    neighbor_policy.preparation_observer=[&] {rss_after_preparation=CurrentRssBytes();};
    j::array progress_sweeps;
    const auto progress_path=std::filesystem::path(output_path.string()+".progress.json");
    std::size_t first_sweep_block_solves{},accepted_blocks{},maximum_local_rows{},maximum_local_columns{};
    neighbor_policy.sweep_observer=[&](const auto & sweep) {
        if(sweep.sweep==1) rss_after_first_sweep=CurrentRssBytes();
        if(sweep.sweep==1) first_sweep_block_solves=sweep.block_solves;
        accepted_blocks+=sweep.accepted_blocks;
        maximum_local_rows=std::max(maximum_local_rows,sweep.maximum_block_rows);
        maximum_local_columns=std::max(maximum_local_columns,sweep.maximum_block_columns);
        progress_sweeps.push_back(SweepJson(sweep));
        j::object progress{{"topology",topology},{"atoms",atoms},
            {"rows",problem.Input().observations.size()},{"parameter_count",3*atoms},
            {"sweeps_completed",progress_sweeps.size()},{"sweep_telemetry",progress_sweeps}};
        Write(progress_path,progress);
        std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor sweep "<<progress_sweeps.size()
            <<" KKT="<<sweep.global_ac_kkt<<" width-grad="<<sweep.global_width_gradient_inf_norm
            <<" seconds="<<sweep.wall_seconds<<'\n';
    };
    std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor started\n";
    auto started=Clock::now(); const auto neighbor=n::RunFixedNeighborDiagnostic(problem,initial_eta,neighbor_policy);
    if(neighbor.certification_attempted) rss_after_certification=CurrentRssBytes();
    const double neighbor_seconds=Seconds(started);
    const double neighbor_search_seconds=neighbor.search_seconds;
    std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor finished in "<<neighbor_seconds<<" s; converged="
        <<neighbor.search_converged<<" sweeps="<<neighbor.sweep_count<<" KKT="
        <<(neighbor.final_sweep ? neighbor.final_sweep->global_ac_kkt : 0.0);
    if(!scaling_only) std::cerr<<" endpoint="<<neighbor.endpoint_certified
        <<" runtime="<<CheckName(neighbor.fit.RuntimeConvergence());
    std::cerr<<'\n';
    j::array sweeps=progress_sweeps;
    j::array blocks;
    for(const auto & block:neighbor.blocks)
    {
        j::object block_json{{"sweep",block.sweep},{"block",block.block},{"atoms",block.core_atoms.size()},
            {"affected_rows",block.affected_rows},{"objective_before",block.objective_before},
            {"objective_after",block.objective_after},{"local_objective_before",block.local_objective_before},
            {"local_objective_after",block.local_objective_after},{"global_replay_delta",block.global_replay_delta},
            {"local_global_delta_error",block.local_global_delta_error},
            {"objective_replay_enclosure",block.objective_replay_enclosure},
            {"profile_evaluations",block.profile_evaluations},{"accepted_updates",block.accepted_updates},
            {"profile_factor_seconds",Number(block.profile_factor_seconds)},
            {"profile_trials",ProfileTrialsJson(block.profile_trials)},
            {"search_seconds",block.search_seconds},{"block_seconds",block.block_seconds},{"accepted",block.accepted},
            {"status",block.status},{"reason",block.reason},{"local_search_stop_reason",block.local_search_stop_reason},
            {"local_final_gradient_inf_norm",Number(block.local_final_gradient_inf_norm)},
            {"local_final_ac_kkt",Number(block.local_final_ac_kkt)}};
        blocks.push_back(std::move(block_json));
    }
    if(scaling_only)
    {
        const double local_factor_seconds=neighbor.work.local_search_seconds;
        const double profile_factor_seconds=std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),0.0,
            [](double total,const auto & block){return total+block.profile_factor_seconds;});
        const auto accepted_local_updates=neighbor.total_accepted_local_updates;
        const auto profile_evaluations=neighbor.total_profile_evaluations;
        const auto blocks_per_sweep=first_sweep_block_solves;
        const auto & last=neighbor.final_sweep ? *neighbor.final_sweep : n::FixedNeighborBlockSweep{};
        j::object scaling_json{{"search_converged",neighbor.search_converged},
                {"measurement_scope",maximum_sweeps==0 ? "preparation-only" :
                    maximum_sweeps==1 ? "one-full-sweep" : "bounded-multi-sweep"},
                {"search_reason",neighbor.reason},{"blocks_per_sweep",blocks_per_sweep},
                {"first_order_stationarity_sweep",neighbor.first_order_stationarity_sweep},
                {"confirmed_stationarity_sweep",neighbor.confirmed_stationarity_sweep},
                {"confirmation_extra_sweeps",neighbor.first_order_stationarity_sweep && neighbor.confirmed_stationarity_sweep ?
                    j::value(neighbor.confirmed_stationarity_sweep-neighbor.first_order_stationarity_sweep) : j::value(nullptr)},
                {"sweeps",neighbor.sweep_count},{"objective",neighbor.state.objective},
                {"final_global_ac_kkt",neighbor.final_sweep ? Number(last.global_ac_kkt) : j::value(nullptr)},
                {"final_raw_width_gradient_inf_norm",neighbor.final_sweep ? Number(last.global_width_gradient_inf_norm) : j::value(nullptr)},
                {"block_solves",neighbor.total_block_solves},
                {"profile_evaluations",profile_evaluations},
                {"accepted_blocks",accepted_blocks},
                {"accepted_updates",accepted_local_updates},
                {"accepted_local_updates",accepted_local_updates},
                {"prepared_block_count",neighbor.block_preparations},
                {"block_preparations",neighbor.block_preparations},
                {"domain_preparations",neighbor.domain_preparations},
                {"mapping_preparations",neighbor.mapping_preparations},
                {"local_factor_seconds",local_factor_seconds},
                {"profile_factor_seconds",detailed ? j::value(profile_factor_seconds) : j::value(nullptr)},
                {"resource_census",ResourceCensusJson(neighbor)},
                {"incremental_replay_validation",IncrementalReplayJson(neighbor.work,neighbor_policy.compare_incremental_replay)},
                {"rss_checkpoints",RssCheckpointsJson(rss_before,rss_after_preparation,rss_after_first_sweep,rss_after_certification)},
                {"fixed_neighbor_work",WorkJson(neighbor.work)},
                {"search_seconds",neighbor_search_seconds},
                {"candidate_replay_fraction",Number(neighbor_search_seconds>0 ?
                    neighbor.work.candidate_replay_seconds/neighbor_search_seconds : 0.)},
                {"local_profile_fraction",Number(neighbor_search_seconds>0 ?
                    neighbor.work.local_search_seconds/neighbor_search_seconds : 0.)},
                {"local_atoms",j::object{{"minimum",neighbor.work.minimum_local_atoms},
                    {"mean",Number(neighbor.work.local_problem_count ? neighbor.work.total_local_atoms/
                        static_cast<double>(neighbor.work.local_problem_count) : 0.)},
                    {"maximum",neighbor.work.maximum_local_atoms}}},
                {"derivative_tile_rows",DerivativeTileRows()},
                {"total_elapsed_seconds",neighbor_seconds},
                {"maximum_local_rows",maximum_local_rows},
                {"maximum_local_columns",maximum_local_columns},
                {"sweep_telemetry",sweeps},{"block_telemetry",blocks}};
        scaling_json["local_profile_work"]=ProfileWorkJson(neighbor.work.local_profile_work);
        const auto peak_rss=PeakRssBytes();
        return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
            {"parameter_count",3*atoms},{"core_atoms",12},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},
            {"block_order","Forward"},
            {"measurement_scope",maximum_sweeps==0 ? "preparation-only" :
                maximum_sweeps==1 ? "one-full-sweep" : "bounded-multi-sweep"},
            {"fixed_neighbor",scaling_json},
            {"peak_rss_bytes",BytesValue(peak_rss)},
            {"peak_rss_mb",peak_rss ? j::value(static_cast<double>(*peak_rss)/(1024.0*1024.0)) : j::value(nullptr)}};
    }
    auto endpoint_assessment=AssessmentJson(neighbor.assessment,neighbor.endpoint_trust);
    j::object neighbor_json{{"method","FixedNeighbor"},
        {"local_atoms",j::object{{"minimum",neighbor.work.minimum_local_atoms},
            {"mean",Number(neighbor.work.local_problem_count ? neighbor.work.total_local_atoms/
                static_cast<double>(neighbor.work.local_problem_count) : 0.)},
            {"maximum",neighbor.work.maximum_local_atoms}}},
        {"derivative_tile_rows",DerivativeTileRows()},
        {"search_converged",neighbor.search_converged},
        {"block_order","Forward"},
        {"search_reason",neighbor.reason},{"sweeps",neighbor.sweep_count},
        {"measurement_scope",maximum_sweeps==0 ? "preparation-only" :
            maximum_sweeps==1 ? "one-full-sweep" : "bounded-multi-sweep"},
        {"first_order_stationarity_sweep",neighbor.first_order_stationarity_sweep},
        {"confirmed_stationarity_sweep",neighbor.confirmed_stationarity_sweep},
        {"endpoint_certified",neighbor.endpoint_certified},
        {"endpoint_trust_reason",neighbor.endpoint_trust.reason},
        {"objective",neighbor.fit.objective ? j::value(*neighbor.fit.objective) : j::value(nullptr)},
        {"global_kkt",CheckValue(neighbor.fit,"conditional-ac-optimality")},
        {"independent_kkt",CheckValue(neighbor.fit,"independent-ac-optimality")},
        {"independent_kkt_consistency",CheckValue(neighbor.fit,"independent-kkt-consistency")},
        {"independent_replay",CheckValue(neighbor.fit,"independent-replay")},
        {"width_gradient_consistency",CheckValue(neighbor.fit,"width-gradient-consistency")},
        {"ac_feasibility",CheckValue(neighbor.fit,"ac-feasibility")},
        {"width_gradient_inf_norm",Number(WidthGradient(neighbor.fit))},
        {"assessment_local",CheckValue(neighbor.fit,"local-correction")},
        {"assessment_identified",CheckValue(neighbor.fit,"numerical-identifiability")},
        {"endpoint_assessment",std::move(endpoint_assessment)},
        {"runtime_convergence",CheckName(neighbor.fit.RuntimeConvergence())},
        {"sweep_telemetry",sweeps},{"block_telemetry",blocks},
        {"search_seconds",neighbor_search_seconds},
        {"assessment_seconds",Number(neighbor.work.endpoint_verification_seconds)},
        {"reference_solve_seconds",0.0},
        {"reference_solve_calls",neighbor.work.local_profile_work.reference_evaluation.evaluations},
        {"endpoint_verification_seconds",Number(neighbor.work.endpoint_verification_seconds)},
        {"independent_replay_seconds",Number(neighbor.work.independent_replay_seconds)},
        {"rank_assessment_seconds",Number(neighbor.work.rank_assessment_seconds)},
        {"estimated_peak_bytes",std::max(neighbor.assessment.verification_memory_estimate.estimated_peak_bytes,
            neighbor.assessment.rank_memory_estimate.estimated_peak_bytes)},
        {"observed_peak_rss_bytes",BytesValue(PeakRssBytes())},
        {"total_search_seconds",neighbor_search_seconds},
        {"total_assessment_seconds",Number(neighbor.work.endpoint_verification_seconds)},
        {"prepared_block_count",neighbor.block_preparations},
        {"block_preparations",neighbor.block_preparations},
        {"domain_preparations",neighbor.domain_preparations},
        {"mapping_preparations",neighbor.mapping_preparations},
        {"fixed_neighbor_work",WorkJson(neighbor.work)},
        {"profile_factor_seconds",detailed ? j::value(std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),0.0,
            [](double total,const auto & block){return total+block.profile_factor_seconds;})) : j::value(nullptr)},
        {"accepted_updates",neighbor.total_accepted_local_updates},
        {"profile_evaluations",neighbor.total_profile_evaluations},
        {"candidate_replay_fraction",Number(neighbor_search_seconds>0 ?
            neighbor.work.candidate_replay_seconds/neighbor_search_seconds : 0.)},
        {"local_profile_fraction",Number(neighbor_search_seconds>0 ?
            neighbor.work.local_search_seconds/neighbor_search_seconds : 0.)},
        {"total_elapsed_seconds",neighbor_seconds}};
    neighbor_json["local_profile_work"]=ProfileWorkJson(neighbor.work.local_profile_work);
    neighbor_json["resource_census"]=ResourceCensusJson(neighbor);
    neighbor_json["incremental_replay_validation"]=IncrementalReplayJson(neighbor.work,neighbor_policy.compare_incremental_replay);
    neighbor_json["rss_checkpoints"]=RssCheckpointsJson(rss_before,rss_after_preparation,rss_after_first_sweep,rss_after_certification);
    const auto peak_rss=PeakRssBytes();
    return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
        {"parameter_count",3*atoms},{"core_atoms",12},
        {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
        {"fixed_neighbor",neighbor_json},{"peak_rss_bytes",BytesValue(peak_rss)},
        {"peak_rss_mb",peak_rss ? j::value(static_cast<double>(*peak_rss)/(1024.0*1024.0)) : j::value(nullptr)}};
}
void Write(const std::filesystem::path & path,const j::value & value)
{std::ofstream output(path); if(!output) throw std::runtime_error("Could not open experiment output."); output<<j::serialize(value)<<'\n';}
}
int main(int argc,char ** argv)
{
    try {
        if(argc<5 || (std::string(argv[1])!="--case" && std::string(argv[1])!="--scaling-only"))
            throw std::invalid_argument("Usage: joint_fixed_neighbor_experiment --case|--scaling-only OUTPUT_FILE TOPOLOGY ATOMS [--maximum-sweeps N] [--detailed] [--compare-incremental-replay]");
        Eigen::setNbThreads(1);
        const std::filesystem::path output_path(argv[2]);
        if(output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());
        const std::string mode(argv[1]);
        std::size_t maximum_sweeps=30; bool detailed=false,compare_incremental_replay=false;
        for(int i=5;i<argc;++i)
        {
            const std::string option(argv[i]);
            if(option=="--detailed") detailed=true;
            else if(option=="--compare-incremental-replay") compare_incremental_replay=true;
            else if(option=="--maximum-sweeps" && i+1<argc)
            {
                const auto value=std::stoull(argv[++i]);
                if(value>200) throw std::invalid_argument("Benchmark sweep bound must not exceed 200.");
                maximum_sweeps=static_cast<std::size_t>(value);
            }
            else throw std::invalid_argument("Unknown FixedNeighbor experiment option: "+option);
        }
        Write(output_path,Run(argv[3],std::stoi(argv[4]),output_path,mode=="--scaling-only",maximum_sweeps,
            detailed,compare_incremental_replay));
        std::cout<<argv[3]<<'-'<<argv[4]<<" fixed-neighbor experiment complete\n";
        return 0;
    } catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
