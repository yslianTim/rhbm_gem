#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/SparseFactor.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointOperatorWorkload.hpp"
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
#include <sys/resource.h>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
using Input=rhbm_gem::core::JointProblemInput;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
double PeakRssMb()
{
    rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("Could not read peak RSS.");
#if defined(__APPLE__)
    return static_cast<double>(usage.ru_maxrss)/(1024.0*1024.0);
#else
    return static_cast<double>(usage.ru_maxrss)/1024.0;
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
std::size_t DerivativeTileRows()
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    return static_cast<std::size_t>(n::DerivativeTileRowsForTesting());
#else
    return static_cast<std::size_t>(n::derivative_tile_rows);
#endif
}
j::array NumberArray(const n::Vector & values)
{j::array out; for(Eigen::Index k=0;k<values.size();++k) out.push_back(Number(values(k))); return out;}
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
    return {{"inner",assessment.inner},{"gradient",assessment.gradient},{"local",assessment.local},
        {"identified",assessment.identified},{"failure",assessment.failure},
        {"coefficient_difference",Number(assessment.coefficient_difference)},
        {"primary_gradient_inf_norm",assessment.primary.gradient.size() ?
            Number(assessment.primary.gradient.lpNorm<Eigen::Infinity>()) : j::value(nullptr)},
        {"reference_gradient_inf_norm",assessment.reference.gradient.size() ?
            Number(assessment.reference.gradient.lpNorm<Eigen::Infinity>()) : j::value(nullptr)},
        {"local_correction_inf_norm",Number(correction)},
        {"projected_width",SpectrumJson(assessment.widths)},
        {"corrected_jacobian",SpectrumJson(assessment.jacobian)},
        {"normalized_width",SpectrumJson(assessment.normalized_widths)},
        {"endpoint_trust",j::object{{"passed",trust.passed},{"reason",trust.reason}}}};
}
j::object CorrectionJson(const n::Assessment & assessment,const n::Vector & eta)
{
    const auto & correction=assessment.correction;
    j::object out{{"inf_norm",correction.size() ? Number(correction.lpNorm<Eigen::Infinity>()) : j::value(nullptr)},
        {"vector",NumberArray(correction)}};
    if(correction.size()==0)
    {
        out["max_coordinate"]=nullptr;
        return out;
    }
    Eigen::Index index{}; correction.cwiseAbs().maxCoeff(&index);
    const double log_width=eta.size()==correction.size() ? eta(index) : std::numeric_limits<double>::quiet_NaN();
    out["max_coordinate"]={{"atom",index},{"value",Number(correction(index))},
        {"eta",Number(log_width)},{"width",Number(std::exp(log_width))}};
    return out;
}
j::array AcScalingWeights(const rhbm_gem::core::JointProblem & problem,const n::Vector & eta)
{
    const auto & input=problem.Input(); const auto & layout=problem.ParameterLayout();
    const double scale=problem.ObservationScale(); j::array weights;
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
    {
        const auto atom=layout.full_atoms[k]; const double width=std::exp(eta(static_cast<Eigen::Index>(atom)));
        double gaussian{},charge{};
        for(const auto & support:input.support.at(atom))
        {
            if(!std::binary_search(layout.informative_rows.begin(),layout.informative_rows.end(),support.row)) continue;
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            gaussian+=basis.gaussian*basis.gaussian; charge+=basis.charge*basis.charge;
        }
        weights.push_back(Number(std::sqrt(gaussian)/scale));
        weights.push_back(Number(std::sqrt(charge)/scale));
    }
    return weights;
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
    return {{"old_core_seconds",Number(work.old_core_seconds)},
        {"effective_response_seconds",Number(work.effective_response_seconds)},
        {"local_search_seconds",Number(work.local_search_seconds)},
        {"local_state_seconds",Number(work.local_state_seconds)},
        {"candidate_copy_seconds",Number(work.candidate_copy_seconds)},
        {"candidate_replay_seconds",Number(work.candidate_replay_seconds)},
        {"cache_update_seconds",Number(work.cache_update_seconds)},
        {"sweep_replay_seconds",Number(work.sweep_replay_seconds)},
        {"sweep_global_state_seconds",Number(work.sweep_global_state_seconds)},
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
j::object ProfileRoleWorkJson(const n::ProfileRoleWork & work)
{
    return {{"evaluations",work.evaluations},{"evaluation_seconds",Number(work.evaluation_seconds)},
        {"profile_basis_seconds",Number(work.basis_seconds)},
        {"linear_matrix_preparation_seconds",Number(work.linear_matrix_preparation_seconds)},
        {"linear_symbolic_seconds",Number(work.linear_symbolic_seconds)},
        {"linear_numeric_seconds",Number(work.linear_numeric_seconds)},
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
void AddSparseAttributionJson(j::object & output,const n::SparseWork & work)
{
    output["numeric_factor_requests"]=work.numeric_factor_requests;
    output["numeric_factor_exact_reuse_opportunities"]=work.numeric_factor_exact_reuse_opportunities;
    output["numeric_factor_pattern_only_matches"]=work.numeric_factor_pattern_only_matches;
    output["numeric_factor_value_mismatches"]=work.numeric_factor_value_mismatches;
    output["numeric_factor_column_mismatches"]=work.numeric_factor_column_mismatches;
    output["numeric_factor_policy_mismatches"]=work.numeric_factor_policy_mismatches;
    output["numeric_factor_pattern_mismatches"]=work.numeric_factor_pattern_mismatches;
    output["initial_profile_exact_reuse_opportunities"]=work.initial_profile_exact_reuse_opportunities;
    output["trial_profile_exact_reuse_opportunities"]=work.trial_profile_exact_reuse_opportunities;
    output["reference_exact_reuse_opportunities"]=work.reference_exact_reuse_opportunities;
    output["accepted_endpoint_exact_reuse_opportunities"]=work.accepted_endpoint_exact_reuse_opportunities;
}
void Write(const std::filesystem::path &,const j::value &);
j::object Run(const std::string & topology,int atoms,const std::filesystem::path & output_path,
    bool scaling_only=false,bool reverse_order=false,bool record_final_state=false,
    std::size_t core_atoms=12)
{
    auto input=std::make_shared<Input>(second_stage_test::OperatorWorkload(topology,atoms));
    const rhbm_gem::core::JointProblem problem(*input);
    n::Vector initial_eta= n::Vector::Constant(atoms,std::log(.55)); n::FixedNeighborPolicy neighbor_policy;
    neighbor_policy.core_atoms=core_atoms;
    neighbor_policy.order=reverse_order ? n::FixedNeighborBlockOrder::Reverse : n::FixedNeighborBlockOrder::Forward;
    neighbor_policy.capture_local_trajectory=true;
    neighbor_policy.assess_final_endpoint=!scaling_only;
    neighbor_policy.collect_telemetry=true;
    j::array progress_sweeps;
    const auto progress_path=std::filesystem::path(output_path.string()+".progress.json");
    neighbor_policy.sweep_observer=[&](const auto & sweep) {
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
    n::SparseWorkForTesting()={};
    auto started=Clock::now(); const auto neighbor=n::SearchFixedNeighbor(problem,initial_eta,neighbor_policy);
    const double neighbor_seconds=Seconds(started);
    const auto sparse_work=n::SparseWorkForTesting();
    const double sweep_seconds=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),0.0,
        [](double total,const auto & sweep){return total+sweep.wall_seconds;});
    const double neighbor_search_seconds=sweep_seconds;
    std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor finished in "<<neighbor_seconds<<" s; converged="
        <<neighbor.search_converged<<" sweeps="<<neighbor.sweeps.size()<<" KKT="
        <<(neighbor.sweeps.empty() ? 0.0 : neighbor.sweeps.back().global_ac_kkt);
    if(!scaling_only) std::cerr<<" endpoint="<<neighbor.endpoint_certified
        <<" runtime="<<CheckName(neighbor.fit.RuntimeConvergence());
    std::cerr<<'\n';
    j::array sweeps;
    for(const auto & sweep:neighbor.sweeps)
        sweeps.push_back(SweepJson(sweep));
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
            {"search_seconds",block.search_seconds},{"accepted",block.accepted},
            {"status",block.status},{"reason",block.reason},{"local_search_stop_reason",block.local_search_stop_reason},
            {"local_final_gradient_inf_norm",Number(block.local_final_gradient_inf_norm)},
            {"local_final_ac_kkt",Number(block.local_final_ac_kkt)}};
        blocks.push_back(std::move(block_json));
    }
    if(scaling_only)
    {
        const auto sum_sweeps=[&](auto member) {
            return std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
                [&](std::size_t total,const auto & sweep){return total+sweep.*member;});
        };
        const double local_factor_seconds=std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),0.0,
            [](double total,const auto & block){return total+block.search_seconds;});
        const double profile_factor_seconds=std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),0.0,
            [](double total,const auto & block){return total+block.profile_factor_seconds;});
        const auto accepted_local_updates=std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),std::size_t{},
            [](std::size_t total,const auto & block){return total+static_cast<std::size_t>(block.accepted_updates);});
        const auto blocks_per_sweep=static_cast<std::size_t>(std::count_if(neighbor.blocks.begin(),neighbor.blocks.end(),
            [](const auto & block){return block.sweep==1;}));
        const auto & last=neighbor.sweeps.empty() ? n::FixedNeighborBlockSweep{} : neighbor.sweeps.back();
        j::object scaling_json{{"search_converged",neighbor.search_converged},
                {"search_reason",neighbor.reason},{"blocks_per_sweep",blocks_per_sweep},
                {"first_order_stationarity_sweep",neighbor.first_order_stationarity_sweep},
                {"confirmed_stationarity_sweep",neighbor.confirmed_stationarity_sweep},
                {"confirmation_extra_sweeps",neighbor.first_order_stationarity_sweep && neighbor.confirmed_stationarity_sweep ?
                    j::value(neighbor.confirmed_stationarity_sweep-neighbor.first_order_stationarity_sweep) : j::value(nullptr)},
                {"sweeps",neighbor.sweeps.size()},{"objective",neighbor.state.objective},
                {"final_global_ac_kkt",neighbor.sweeps.empty() ? j::value(nullptr) : Number(last.global_ac_kkt)},
                {"final_raw_width_gradient_inf_norm",neighbor.sweeps.empty() ? j::value(nullptr) : Number(last.global_width_gradient_inf_norm)},
                {"block_solves",sum_sweeps(&n::FixedNeighborBlockSweep::block_solves)},
                {"profile_evaluations",sum_sweeps(&n::FixedNeighborBlockSweep::profile_evaluations)},
                {"accepted_blocks",sum_sweeps(&n::FixedNeighborBlockSweep::accepted_blocks)},
                {"accepted_local_updates",accepted_local_updates},
                {"prepared_block_count",neighbor.block_preparations},
                {"block_preparations",neighbor.block_preparations},
                {"domain_preparations",neighbor.domain_preparations},
                {"mapping_preparations",neighbor.mapping_preparations},
                {"symbolic_factorizations",sparse_work.symbolic},
                {"symbolic_reuses",sparse_work.symbolic_reuses},
                {"numeric_factorizations",sparse_work.numeric},
                {"fresh_workspace_symbolic_factorizations",sparse_work.numeric},
                {"symbolic_seconds",sparse_work.symbolic_seconds},
                {"numeric_seconds",sparse_work.numeric_seconds},
                {"matrix_preparation_seconds",sparse_work.matrix_preparation_seconds},
                {"factor_storage_bytes",sparse_work.factor_storage_bytes},
                {"local_factor_seconds",local_factor_seconds},{"profile_factor_seconds",profile_factor_seconds},
                {"fixed_neighbor_work",WorkJson(neighbor.work)},
                {"search_seconds",neighbor_search_seconds},
                {"outer_core_atoms",neighbor_policy.core_atoms},
                {"local_atoms",j::object{{"minimum",neighbor.work.minimum_local_atoms},
                    {"mean",Number(neighbor.work.local_problem_count ? neighbor.work.total_local_atoms/
                        static_cast<double>(neighbor.work.local_problem_count) : 0.)},
                    {"maximum",neighbor.work.maximum_local_atoms}}},
                {"derivative_tile_rows",DerivativeTileRows()},
                {"total_elapsed_seconds",neighbor_seconds},
                {"maximum_local_rows",std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
                    [](std::size_t maximum,const auto & sweep){return std::max(maximum,sweep.maximum_block_rows);})},
                {"maximum_local_columns",std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
                    [](std::size_t maximum,const auto & sweep){return std::max(maximum,sweep.maximum_block_columns);})},
                {"sweep_telemetry",sweeps},{"block_telemetry",blocks}};
        scaling_json["local_profile_work"]=ProfileWorkJson(neighbor.work.local_profile_work);
        AddSparseAttributionJson(scaling_json,sparse_work);
        if(record_final_state)
        {
            scaling_json["final_eta"]=NumberArray(neighbor.state.eta);
            scaling_json["final_beta"]=NumberArray(neighbor.state.beta);
            scaling_json["final_ac_scaling_weights"]=AcScalingWeights(problem,neighbor.state.eta);
        }
        return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
            {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
            {"outer_core_atoms",neighbor_policy.core_atoms},
            {"workspace_mode","persistent"},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},
            {"block_order",reverse_order ? "reverse" : "forward"},
            {"measurement_scope","fixed-neighbor-search-only"},
            {"fixed_neighbor",scaling_json},
            {"peak_rss_mb",PeakRssMb()}};
    }
    auto endpoint_assessment=AssessmentJson(neighbor.assessment,neighbor.endpoint_trust);
    if(record_final_state)
        endpoint_assessment["correction"]=CorrectionJson(neighbor.assessment,neighbor.state.eta);
    j::object neighbor_json{{"method","FixedNeighbor"},
        {"outer_core_atoms",neighbor_policy.core_atoms},
        {"local_atoms",j::object{{"minimum",neighbor.work.minimum_local_atoms},
            {"mean",Number(neighbor.work.local_problem_count ? neighbor.work.total_local_atoms/
                static_cast<double>(neighbor.work.local_problem_count) : 0.)},
            {"maximum",neighbor.work.maximum_local_atoms}}},
        {"derivative_tile_rows",DerivativeTileRows()},
        {"search_converged",neighbor.search_converged},
        {"workspace_mode","persistent"},
        {"block_order",reverse_order ? "reverse" : "forward"},
        {"search_reason",neighbor.reason},{"sweeps",neighbor.sweeps.size()},
        {"first_order_stationarity_sweep",neighbor.first_order_stationarity_sweep},
        {"confirmed_stationarity_sweep",neighbor.confirmed_stationarity_sweep},
        {"endpoint_certified",neighbor.endpoint_certified},
        {"endpoint_trust_reason",neighbor.endpoint_trust.reason},
        {"objective",neighbor.fit.objective ? j::value(*neighbor.fit.objective) : j::value(nullptr)},
        {"global_kkt",CheckValue(neighbor.fit,"kkt")},{"width_gradient_inf_norm",WidthGradient(neighbor.fit)},
        {"assessment_inner",CheckValue(neighbor.fit,"inner")},{"assessment_gradient",CheckValue(neighbor.fit,"width-stationarity")},
        {"assessment_local",CheckValue(neighbor.fit,"local-correction")},{"assessment_identified",CheckValue(neighbor.fit,"numerical-identifiability")},
        {"endpoint_assessment",std::move(endpoint_assessment)},
        {"runtime_convergence",CheckName(neighbor.fit.RuntimeConvergence())},
        {"sweep_telemetry",sweeps},{"block_telemetry",blocks},
        {"search_seconds",neighbor_search_seconds},
        {"assessment_seconds",std::max(0.0,neighbor_seconds-neighbor_search_seconds)},
        {"prepared_block_count",neighbor.block_preparations},
        {"block_preparations",neighbor.block_preparations},
        {"domain_preparations",neighbor.domain_preparations},
        {"mapping_preparations",neighbor.mapping_preparations},
        {"symbolic_factorizations",sparse_work.symbolic},
        {"symbolic_reuses",sparse_work.symbolic_reuses},
        {"numeric_factorizations",sparse_work.numeric},
        {"fresh_workspace_symbolic_factorizations",sparse_work.numeric},
        {"symbolic_seconds",sparse_work.symbolic_seconds},
        {"numeric_seconds",sparse_work.numeric_seconds},
        {"matrix_preparation_seconds",sparse_work.matrix_preparation_seconds},
        {"factor_storage_bytes",sparse_work.factor_storage_bytes},
        {"fixed_neighbor_work",WorkJson(neighbor.work)},
        {"profile_factor_seconds",std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),0.0,
            [](double total,const auto & block){return total+block.profile_factor_seconds;})},
        {"total_elapsed_seconds",neighbor_seconds}};
    neighbor_json["local_profile_work"]=ProfileWorkJson(neighbor.work.local_profile_work);
    AddSparseAttributionJson(neighbor_json,sparse_work);
    if(record_final_state)
    {
        neighbor_json["final_eta"]=NumberArray(neighbor.state.eta);
        neighbor_json["final_beta"]=NumberArray(neighbor.state.beta);
    }
    if(record_final_state)
        neighbor_json["final_ac_scaling_weights"]=AcScalingWeights(problem,neighbor.state.eta);
    return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
        {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
        {"outer_core_atoms",neighbor_policy.core_atoms},
        {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
        {"fixed_neighbor",neighbor_json},{"peak_rss_mb",PeakRssMb()}};
}
void Write(const std::filesystem::path & path,const j::value & value)
{std::ofstream output(path); if(!output) throw std::runtime_error("Could not open experiment output."); output<<j::serialize(value)<<'\n';}
}
int main(int argc,char ** argv)
{
    try {
        if((argc<5 || argc>6) || (std::string(argv[1])!="--case" && std::string(argv[1])!="--neighbor-only" &&
            std::string(argv[1])!="--scaling-only" &&
            std::string(argv[1])!="--neighbor-forward" && std::string(argv[1])!="--neighbor-reverse" &&
            std::string(argv[1])!="--scaling-forward" && std::string(argv[1])!="--scaling-reverse"))
            throw std::invalid_argument("Usage: joint_fixed_neighbor_experiment --case|--neighbor-only|--scaling-only|--neighbor-forward|--neighbor-reverse|--scaling-forward|--scaling-reverse OUTPUT_FILE TOPOLOGY ATOMS [OUTER_CORE_ATOMS]");
        Eigen::setNbThreads(1);
        const std::filesystem::path output_path(argv[2]);
        if(output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());
        const std::string mode(argv[1]);
        if(argc>=6 && mode!="--neighbor-only" && mode!="--scaling-only" &&
            mode!="--neighbor-forward" && mode!="--neighbor-reverse" && mode!="--scaling-forward" &&
            mode!="--scaling-reverse")
            throw std::invalid_argument("OUTER_CORE_ATOMS is supported only with a FixedNeighbor search mode.");
        const std::size_t core_atoms=argc>=6 ? static_cast<std::size_t>(std::stoul(argv[5])) : 12;
        if(core_atoms==0) throw std::invalid_argument("CORE_ATOMS must be positive.");
        Write(output_path,Run(argv[3],std::stoi(argv[4]),output_path,
            mode=="--scaling-only" || mode=="--scaling-forward" || mode=="--scaling-reverse",
            mode=="--neighbor-reverse" || mode=="--scaling-reverse",
            mode=="--neighbor-forward" || mode=="--neighbor-reverse" ||
                mode=="--scaling-forward" || mode=="--scaling-reverse",
            core_atoms));
        std::cout<<argv[3]<<'-'<<argv[4]<<" fixed-neighbor experiment complete\n";
        return 0;
    } catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
