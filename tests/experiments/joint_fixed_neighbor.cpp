#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
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
const char * LocalWorkName(n::FixedNeighborLocalWork work)
{
    switch(work)
    {
    case n::FixedNeighborLocalWork::Full: return "FullLocalSearch";
    case n::FixedNeighborLocalWork::OneAcceptedUpdate: return "OneAcceptedLocalUpdate";
    case n::FixedNeighborLocalWork::TwoAcceptedUpdates: return "TwoAcceptedLocalUpdates";
    }
    return "FullLocalSearch";
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
double CoefficientDifference(const n::Vector & lhs,const n::Vector & rhs)
{
    if(lhs.size()!=rhs.size() || lhs.size()==0) return std::numeric_limits<double>::infinity();
    return ((lhs-rhs).array().abs()/(1+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
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
        {"local_correction_inf_norm",Number(correction)},
        {"projected_width",SpectrumJson(assessment.widths)},
        {"corrected_jacobian",SpectrumJson(assessment.jacobian)},
        {"normalized_width",SpectrumJson(assessment.normalized_widths)},
        {"endpoint_trust",j::object{{"passed",trust.passed},{"reason",trust.reason}}}};
}
j::object EndpointStateJson(const std::string & name,const n::Evaluation & endpoint,
    const n::Assessment & assessment,const n::TrustEvidence & trust,double scale)
{
    return {{"state",name},{"valid",endpoint.valid},{"evaluation_reason",endpoint.reason},
        {"objective",endpoint.certificate.evaluated ? Number(endpoint.certificate.objective/(scale*scale)) : j::value(nullptr)},
        {"a_feasible",endpoint.certificate.feasible},
        {"ac_kkt",Number(endpoint.certificate.projected_kkt)},
        {"ac_kkt_passed",endpoint.certificate.available && endpoint.certificate.feasible &&
            endpoint.certificate.kkt_passed && endpoint.certificate.projected_kkt<=1e-10},
        {"width_gradient_inf_norm",endpoint.gradient.size() ? Number(endpoint.gradient.lpNorm<Eigen::Infinity>()) : j::value(nullptr)},
        {"assessment",AssessmentJson(assessment,trust)}};
}
j::object DecomposeEndpoint(const n::Domain & domain,const n::Vector & y,const n::EvaluationContext & context,
    const n::BlockCoordinateState & state,std::size_t sweep,const n::FixedNeighborBlockSweep & sweep_record,
    const std::vector<n::FixedNeighborBlockRecord> & block_records)
{
    const n::Vector eta=state.eta,beta=state.beta;
    const auto raw=n::EvaluateState(domain,y,eta,beta,context);
    const auto primary=n::EvaluateProfile(domain,y,eta,false,&context);
    const auto reference=n::EvaluateProfile(domain,y,eta,true,&context);
    const auto raw_assessment=n::AssessEvaluated(domain,y,raw,reference,context,true);
    const auto primary_assessment=n::AssessEvaluated(domain,y,primary,reference,context,false);
    const auto reference_assessment=n::AssessEvaluated(domain,y,reference,primary,context,false);
    const auto raw_trust=n::CheckTrust(domain,y,raw,context,reference);
    const auto primary_trust=n::CheckTrust(domain,y,primary,context,reference);
    const auto reference_trust=n::CheckTrust(domain,y,reference,context,primary);
    const bool eta_unchanged=(state.eta.array()==eta.array()).all();
    j::array local_blocks;
    for(const auto & block:block_records) if(block.sweep==sweep)
        local_blocks.push_back({{"block",block.block},{"core_atoms",block.core_atoms.size()},
            {"local_search_stop_reason",block.local_search_stop_reason},
            {"profile_evaluations",block.profile_evaluations},{"accepted_updates",block.accepted_updates},
            {"local_final_gradient_inf_norm",Number(block.local_final_gradient_inf_norm)},
            {"local_final_ac_kkt",Number(block.local_final_ac_kkt)}});
    return {{"sweep",sweep},{"eta_unchanged",eta_unchanged},
        {"eta",NumberArray(eta)},{"beta",NumberArray(beta)},
        {"search_state",j::object{{"objective",state.objective},{"global_ac_kkt",sweep_record.global_ac_kkt},
            {"raw_width_gradient_inf_norm",sweep_record.global_width_gradient_inf_norm}}},
        {"coefficient_difference_block_primary",Number(CoefficientDifference(beta,primary.beta))},
        {"coefficient_difference_primary_reference",Number(CoefficientDifference(primary.beta,reference.beta))},
        {"coefficient_difference_block_reference",Number(CoefficientDifference(beta,reference.beta))},
        {"raw",EndpointStateJson("raw",raw,raw_assessment,raw_trust,context.scale)},
        {"same_eta_primary",EndpointStateJson("same-eta-primary",primary,primary_assessment,primary_trust,context.scale)},
        {"same_eta_reference",EndpointStateJson("same-eta-reference",reference,reference_assessment,reference_trust,context.scale)},
        {"local_blocks",local_blocks}};
}
j::object AssessSweepEndpoint(const rhbm_gem::JointParameterLayout & layout,
    const n::Domain & domain,const n::Vector & y,const n::EvaluationContext & context,
    const n::BlockCoordinateState & state,std::size_t sweep)
{
    const n::Vector eta=state.eta,beta=state.beta;
    const auto endpoint=n::EvaluateState(domain,y,eta,beta,context);
    const auto reference=n::EvaluateProfile(domain,y,eta,true,&context);
    const auto assessment=n::AssessEvaluated(domain,y,endpoint,reference,context,true);
    const auto trust=n::CheckTrust(domain,y,endpoint,context,reference);
    const auto control=n::EvaluateProfile(domain,y,eta,false,&context);
    bool agrees=control.valid;
    double coefficient_difference=std::numeric_limits<double>::infinity();
    n::Vector expected=n::Vector::Zero(y.size());
    if(agrees)
    {
        expected=control.x*beta;
        coefficient_difference=CoefficientDifference(beta,control.beta);
        agrees=coefficient_difference<=1e-10;
    }
    n::Vector actual(y.size());
    for(std::size_t k=0;k<layout.informative_rows.size();++k)
        actual(static_cast<Eigen::Index>(k))=state.prediction(static_cast<Eigen::Index>(layout.informative_rows[k]));
    const double prediction_difference=((actual-expected).array().abs() /
        (1.0+actual.array().abs().max(expected.array().abs()))).maxCoeff();
    const double profile_objective=.5*(expected-y).squaredNorm()/(context.scale*context.scale);
    const bool reconstruction=actual.allFinite() && expected.allFinite() && prediction_difference<=1e-10 &&
        std::abs(state.objective-profile_objective)<=1e-12;
    agrees &= reconstruction;

    using Scope=rhbm_gem::JointEvidenceScope;
    using Status=rhbm_gem::JointCheckStatus;
    auto global_evidence=n::AssessmentEvidence(assessment,Scope::AssembledGlobal);
    global_evidence.push_back({"assembled-profile",agrees ? Status::Passed : Status::Failed,
        Scope::AssembledGlobal,std::max(coefficient_difference,prediction_difference),1e-10,
        reconstruction ? "" : "full-domain-reconstruction-failed"});
    const auto global_status=n::ConvergenceStatus(global_evidence,Scope::AssembledGlobal,true);
    const auto component_status=n::ConvergenceStatus(
        n::AssessmentEvidence(assessment,Scope::ComponentLocal),Scope::ComponentLocal);
    const auto runtime_status=n::MergeConvergenceStatus(global_status,component_status);
    const bool endpoint_certified=n::IsCertifiedLocalEndpoint(assessment,trust);
    return {{"sweep",sweep},{"objective",state.objective},
        {"global_ac_kkt",Number(endpoint.certificate.projected_kkt)},
        {"width_gradient_inf_norm",Number(endpoint.gradient.lpNorm<Eigen::Infinity>())},
        {"inner",assessment.inner},{"gradient",assessment.gradient},{"local",assessment.local},
        {"identified",assessment.identified},{"endpoint_trust",trust.passed},
        {"endpoint_trust_reason",trust.reason},{"endpoint_certified",endpoint_certified},
        {"runtime_convergence",CheckName(runtime_status)},
        {"assessment",AssessmentJson(assessment,trust)}};
}
j::object PackFit(const std::string & method,const rhbm_gem::core::JointFitResult & fit,double seconds)
{
    return {{"method",method},{"search_completed",fit.search_completed},
        {"search_seconds",fit.costs.search_seconds},{"total_elapsed_seconds",seconds},
        {"objective",fit.objective ? j::value(*fit.objective) : j::value(nullptr)},
        {"global_kkt",CheckValue(fit,"kkt")},{"width_gradient_inf_norm",WidthGradient(fit)},
        {"assessment_inner",CheckValue(fit,"inner")},{"assessment_gradient",CheckValue(fit,"width-stationarity")},
        {"assessment_local",CheckValue(fit,"local-correction")},{"assessment_identified",CheckValue(fit,"numerical-identifiability")},
        {"runtime_convergence",CheckName(fit.RuntimeConvergence())}};
}
struct MethodState {std::string name; const rhbm_gem::core::JointFitResult * fit{};};
double ScaledAcDifference(const rhbm_gem::core::JointProblem & problem,const MethodState & lhs,const MethodState & rhs)
{
    if(!lhs.fit || !rhs.fit || !lhs.fit->assembled_state || !rhs.fit->assembled_state) return std::numeric_limits<double>::infinity();
    const auto & input=problem.Input(); const auto & layout=problem.ParameterLayout();
    const auto & a=lhs.fit->assembled_state->ac; const auto & b=rhs.fit->assembled_state->ac;
    if(a.size()!=b.size()) return std::numeric_limits<double>::infinity();
    const double scale=problem.ObservationScale(); double difference{};
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
    {
        const auto atom=layout.full_atoms[k]; const double width=rhs.fit->assembled_state->b.at(k);
        double gaussian{},charge{};
        for(const auto & support:input.support.at(atom))
        {
            if(!std::binary_search(layout.informative_rows.begin(),layout.informative_rows.end(),support.row)) continue;
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            gaussian+=basis.gaussian*basis.gaussian; charge+=basis.charge*basis.charge;
        }
        const double scales[]{std::sqrt(gaussian)/scale,std::sqrt(charge)/scale};
        for(std::size_t kind=0;kind<2;++kind)
            difference=std::max(difference,std::abs(a[2*k+kind]-b[2*k+kind])*scales[kind]);
    }
    return difference;
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
        {"local_assessments",sweep.local_assessments},{"certified_local_candidates",sweep.certified_local_candidates},
        {"local_assessment_seconds",sweep.local_assessment_seconds},
        {"maximum_local_assessment_rows",sweep.maximum_local_assessment_rows},
        {"maximum_local_assessment_columns",sweep.maximum_local_assessment_columns},
        {"wall_seconds",sweep.wall_seconds}};
}
void Write(const std::filesystem::path &,const j::value &);
j::object Run(const std::string & topology,int atoms,const std::filesystem::path & output_path,bool compare_global,
    bool decompose=false,bool certify_local=false,bool qualification=false,bool scaling_only=false,
    bool reverse_order=false,bool record_final_state=false,bool attribution=false,
    n::FixedNeighborLocalWork local_work=n::FixedNeighborLocalWork::Full)
{
    auto input=std::make_shared<Input>(second_stage_test::OperatorWorkload(topology,atoms));
    const rhbm_gem::core::JointProblem problem(*input);
    const std::vector<double> initial_b(static_cast<std::size_t>(atoms),.55);
    n::Vector initial_eta= n::Vector::Constant(atoms,std::log(.55)); n::FixedNeighborPolicy neighbor_policy;
    neighbor_policy.core_atoms=128;
    neighbor_policy.order=reverse_order ? n::FixedNeighborBlockOrder::Reverse : n::FixedNeighborBlockOrder::Forward;
    neighbor_policy.local_work=local_work;
    neighbor_policy.certify_local_candidates=certify_local;
    neighbor_policy.capture_local_trajectory=attribution || local_work!=n::FixedNeighborLocalWork::Full;
    neighbor_policy.stop_after_no_certified_update=certify_local;
    neighbor_policy.assess_final_endpoint=!scaling_only;
    if(qualification)
    {
        neighbor_policy.stop_after_stationarity=false;
        neighbor_policy.maximum_sweeps=atoms==256 ? (topology=="cube" ? 15 : 10) : 5;
    }
    if(decompose && topology=="cube" && atoms==256)
    {neighbor_policy.maximum_sweeps=15; neighbor_policy.stop_after_stationarity=false;}
    j::array progress_sweeps;
    j::array endpoint_snapshots;
    j::array local_block_snapshots;
    j::array endpoint_assessments;
    const auto progress_path=std::filesystem::path(output_path.string()+".progress.json");
    neighbor_policy.sweep_observer=[&](const auto & sweep) {
        progress_sweeps.push_back(SweepJson(sweep));
        j::object progress{{"topology",topology},{"atoms",atoms},
            {"rows",problem.Input().observations.size()},{"parameter_count",3*atoms},
            {"sweeps_completed",progress_sweeps.size()},{"sweep_telemetry",progress_sweeps}};
        if(decompose) progress["endpoint_decomposition"]=endpoint_snapshots;
        if(certify_local || attribution) progress["local_block_telemetry"]=local_block_snapshots;
        if(qualification) progress["endpoint_assessment_by_sweep"]=endpoint_assessments;
        Write(progress_path,progress);
        std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor sweep "<<progress_sweeps.size()
            <<" KKT="<<sweep.global_ac_kkt<<" width-grad="<<sweep.global_width_gradient_inf_norm
            <<" seconds="<<sweep.wall_seconds<<'\n';
    };
    if(decompose || certify_local || qualification || attribution)
        neighbor_policy.state_observer=[&](std::size_t sweep,const n::BlockCoordinateState & state,
            const n::FixedNeighborBlockSweep & sweep_record,const std::vector<n::FixedNeighborBlockRecord> & blocks) {
            if(qualification)
            {
                const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
                const auto domain=n::ProfileDomain(data.domain,data.layout);
                const auto context=n::ProfileContext(data.context,data.layout,data.domain.rows);
                endpoint_assessments.push_back(AssessSweepEndpoint(data.layout,domain,data.y,
                    context,state,sweep));
                Write(progress_path,j::object{{"topology",topology},{"atoms",atoms},
                    {"sweeps_completed",progress_sweeps.size()},{"sweep_telemetry",progress_sweeps},
                    {"endpoint_assessment_by_sweep",endpoint_assessments}});
            }
            if(certify_local || attribution)
            {
                for(const auto & block:blocks) if(block.sweep==sweep)
                {
                    j::object snapshot{{"sweep",block.sweep},{"block",block.block},
                        {"affected_rows",block.affected_rows},{"local_assessment_rows",block.local_assessment_rows},
                        {"local_assessment_columns",block.local_assessment_columns},
                        {"local_assessment_attempted",block.local_assessment_attempted},
                        {"local_assessment_passed",block.local_assessment_passed},
                        {"local_inner_passed",block.local_inner_passed},{"local_gradient_passed",block.local_gradient_passed},
                        {"local_correction_passed",block.local_correction_passed},{"local_identified",block.local_identified},
                        {"local_trust_passed",block.local_trust_passed},
                        {"local_profile_gradient_inf_norm",Number(block.local_profile_gradient_inf_norm)},
                        {"local_reference_gradient_inf_norm",Number(block.local_reference_gradient_inf_norm)},
                        {"local_correction_inf_norm",Number(block.local_correction_inf_norm)},
                        {"local_projected_width_minimum",Number(block.local_projected_width_minimum)},
                        {"local_projected_width_rank",block.local_projected_width_rank},
                        {"local_corrected_jacobian_minimum",Number(block.local_corrected_jacobian_minimum)},
                        {"local_corrected_jacobian_rank",block.local_corrected_jacobian_rank},
                        {"local_normalized_width_minimum",Number(block.local_normalized_width_minimum)},
                        {"local_normalized_width_rank",block.local_normalized_width_rank},
                        {"local_assessment_failure",block.local_assessment_failure},
                        {"local_trust_reason",block.local_trust_reason},
                        {"local_assessment_seconds",block.local_assessment_seconds},
                        {"profile_factor_seconds",Number(block.profile_factor_seconds)},
                        {"profile_trials",ProfileTrialsJson(block.profile_trials)},
                        {"accepted",block.accepted},{"status",block.status},{"reason",block.reason}};
                    local_block_snapshots.push_back(std::move(snapshot));
                }
                Write(progress_path,j::object{{"topology",topology},{"atoms",atoms},
                    {"sweeps_completed",progress_sweeps.size()},{"sweep_telemetry",progress_sweeps},
                    {"local_block_telemetry",local_block_snapshots}});
            }
            if(!decompose) return;
            const bool cube_snapshot=topology=="cube" && atoms==256 && sweep>=8 && sweep<=15;
            const bool chain_control=topology=="chain" && atoms==256 &&
                sweep_record.global_ac_kkt<=1e-10 && sweep_record.global_width_gradient_inf_norm<=1e-12;
            if(!cube_snapshot && !chain_control) return;
            const auto diagnostic_started=Clock::now();
            const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
            const auto domain=n::ProfileDomain(data.domain,data.layout);
            const auto context=n::ProfileContext(data.context,data.layout,data.domain.rows);
            const n::Vector y=data.y;
            auto snapshot=DecomposeEndpoint(domain,y,context,state,sweep,sweep_record,blocks);
            snapshot["diagnostic_seconds"]=Seconds(diagnostic_started);
            endpoint_snapshots.push_back(std::move(snapshot));
            Write(progress_path,j::object{{"topology",topology},{"atoms",atoms},
                {"sweeps_completed",progress_sweeps.size()},{"sweep_telemetry",progress_sweeps},
                {"endpoint_decomposition",endpoint_snapshots}});
        };
    std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor started\n";
    auto started=Clock::now(); const auto neighbor=n::SearchFixedNeighbor(problem,initial_eta,neighbor_policy);
    const double neighbor_seconds=Seconds(started);
    const double sweep_seconds=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),0.0,
        [](double total,const auto & sweep){return total+sweep.wall_seconds;});
    const double local_assessment_seconds=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),0.0,
        [](double total,const auto & sweep){return total+sweep.local_assessment_seconds;});
    const double neighbor_search_seconds=sweep_seconds-local_assessment_seconds;
    const auto local_assessments=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
        [](std::size_t total,const auto & sweep){return total+sweep.local_assessments;});
    const auto certified_local_candidates=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
        [](std::size_t total,const auto & sweep){return total+sweep.certified_local_candidates;});
    const auto maximum_local_assessment_rows=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
        [](std::size_t maximum,const auto & sweep){return std::max(maximum,sweep.maximum_local_assessment_rows);});
    const auto maximum_local_assessment_columns=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
        [](std::size_t maximum,const auto & sweep){return std::max(maximum,sweep.maximum_local_assessment_columns);});
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
            {"local_assessment_attempted",block.local_assessment_attempted},
            {"local_assessment_passed",block.local_assessment_passed},
            {"local_inner_passed",block.local_inner_passed},{"local_gradient_passed",block.local_gradient_passed},
            {"local_correction_passed",block.local_correction_passed},{"local_identified",block.local_identified},
            {"local_trust_passed",block.local_trust_passed},{"local_assessment_failure",block.local_assessment_failure},
            {"local_profile_gradient_inf_norm",Number(block.local_profile_gradient_inf_norm)},
            {"local_reference_gradient_inf_norm",Number(block.local_reference_gradient_inf_norm)},
            {"local_correction_inf_norm",Number(block.local_correction_inf_norm)},
            {"local_projected_width_minimum",Number(block.local_projected_width_minimum)},
            {"local_projected_width_rank",block.local_projected_width_rank},
            {"local_corrected_jacobian_minimum",Number(block.local_corrected_jacobian_minimum)},
            {"local_corrected_jacobian_rank",block.local_corrected_jacobian_rank},
            {"local_normalized_width_minimum",Number(block.local_normalized_width_minimum)},
            {"local_normalized_width_rank",block.local_normalized_width_rank},
            {"local_trust_reason",block.local_trust_reason},
            {"local_assessment_rows",block.local_assessment_rows},{"local_assessment_columns",block.local_assessment_columns},
            {"local_assessment_seconds",block.local_assessment_seconds},
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
                {"local_factor_seconds",local_factor_seconds},{"profile_factor_seconds",profile_factor_seconds},
                {"search_seconds",neighbor_search_seconds},
                {"local_trajectory_telemetry",attribution},
                {"local_work_policy",LocalWorkName(local_work)},
                {"total_elapsed_seconds",neighbor_seconds},
                {"maximum_local_rows",std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
                    [](std::size_t maximum,const auto & sweep){return std::max(maximum,sweep.maximum_block_rows);})},
                {"maximum_local_columns",std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),std::size_t{},
                    [](std::size_t maximum,const auto & sweep){return std::max(maximum,sweep.maximum_block_columns);})},
                {"sweep_telemetry",sweeps},{"block_telemetry",blocks}};
        if(record_final_state)
        {
            scaling_json["final_eta"]=NumberArray(neighbor.state.eta);
            scaling_json["final_beta"]=NumberArray(neighbor.state.beta);
            scaling_json["final_ac_scaling_weights"]=AcScalingWeights(problem,neighbor.state.eta);
        }
        return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
            {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},
            {"block_order",reverse_order ? "reverse" : "forward"},
            {"local_work_policy",LocalWorkName(local_work)},
            {"measurement_scope","fixed-neighbor-search-only"},
            {"fixed_neighbor",scaling_json},
            {"peak_rss_mb",PeakRssMb()}};
    }
    j::object neighbor_json{{"method","FixedNeighbor"},{"search_converged",neighbor.search_converged},
        {"block_order",reverse_order ? "reverse" : "forward"},
        {"local_work_policy",LocalWorkName(local_work)},
        {"search_reason",neighbor.reason},{"sweeps",neighbor.sweeps.size()},
        {"first_order_stationarity_sweep",neighbor.first_order_stationarity_sweep},
        {"confirmed_stationarity_sweep",neighbor.confirmed_stationarity_sweep},
        {"endpoint_certified",neighbor.endpoint_certified},
        {"endpoint_trust_reason",neighbor.endpoint_trust.reason},
        {"objective",neighbor.fit.objective ? j::value(*neighbor.fit.objective) : j::value(nullptr)},
        {"global_kkt",CheckValue(neighbor.fit,"kkt")},{"width_gradient_inf_norm",WidthGradient(neighbor.fit)},
        {"assessment_inner",CheckValue(neighbor.fit,"inner")},{"assessment_gradient",CheckValue(neighbor.fit,"width-stationarity")},
        {"assessment_local",CheckValue(neighbor.fit,"local-correction")},{"assessment_identified",CheckValue(neighbor.fit,"numerical-identifiability")},
        {"endpoint_assessment",AssessmentJson(neighbor.assessment,neighbor.endpoint_trust)},
        {"runtime_convergence",CheckName(neighbor.fit.RuntimeConvergence())},
        {"sweep_telemetry",sweeps},{"block_telemetry",blocks},
        {"search_seconds",neighbor_search_seconds},
        {"assessment_seconds",std::max(0.0,neighbor_seconds-neighbor_search_seconds)},
        {"local_assessment_count",local_assessments},
        {"certified_local_candidates",certified_local_candidates},
        {"local_assessment_seconds",local_assessment_seconds},
        {"profile_factor_seconds",std::accumulate(neighbor.blocks.begin(),neighbor.blocks.end(),0.0,
            [](double total,const auto & block){return total+block.profile_factor_seconds;})},
        {"maximum_local_assessment_rows",maximum_local_assessment_rows},
        {"maximum_local_assessment_columns",maximum_local_assessment_columns},
        {"total_elapsed_seconds",neighbor_seconds}};
    if(certify_local || record_final_state)
    {
        neighbor_json["final_eta"]=NumberArray(neighbor.state.eta);
        neighbor_json["final_beta"]=NumberArray(neighbor.state.beta);
    }
    if(record_final_state)
        neighbor_json["final_ac_scaling_weights"]=AcScalingWeights(problem,neighbor.state.eta);
    if(decompose)
        return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
            {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
            {"fixed_neighbor",neighbor_json},{"endpoint_decomposition",endpoint_snapshots},
            {"peak_rss_mb",PeakRssMb()}};
    if(qualification)
        return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
            {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
            {"fixed_neighbor",neighbor_json},{"endpoint_assessment_by_sweep",endpoint_assessments},
            {"peak_rss_mb",PeakRssMb()}};
    if(!compare_global)
        return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
            {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
            {"fixed_neighbor",neighbor_json},{"peak_rss_mb",PeakRssMb()}};
    n::SearchPolicy legacy_policy;
    std::cerr<<topology<<'-'<<atoms<<" global LegacyCompact started\n";
    started=Clock::now(); const auto legacy=n::FitWithSearchPolicy(problem,initial_b,legacy_policy);
    const double legacy_seconds=Seconds(started);
    std::cerr<<topology<<'-'<<atoms<<" global LegacyCompact finished in "<<legacy_seconds<<" s\n";
    n::SearchPolicy operator_policy; operator_policy.method=n::SearchMethod::OperatorPcg;
    std::cerr<<topology<<'-'<<atoms<<" global OperatorPcg started\n";
    started=Clock::now(); const auto operator_fit=n::FitWithSearchPolicy(problem,initial_b,operator_policy);
    const double operator_seconds=Seconds(started);
    std::cerr<<topology<<'-'<<atoms<<" global OperatorPcg finished in "<<operator_seconds<<" s\n";
    const MethodState reference{"LegacyCompact",&legacy};
    j::array comparisons;
    for(const MethodState & candidate:std::array<MethodState,2>{{{"OperatorPcg",&operator_fit},{"FixedNeighbor",&neighbor.fit}}})
    {
        const double objective_difference=legacy.objective && candidate.fit->objective ?
            std::abs(*legacy.objective-*candidate.fit->objective) : std::numeric_limits<double>::infinity();
        const double eta_difference=legacy.assembled_state && candidate.fit->assembled_state ?
            (n::Vector::Map(legacy.assembled_state->log_b.data(),static_cast<Eigen::Index>(legacy.assembled_state->log_b.size()))-
             n::Vector::Map(candidate.fit->assembled_state->log_b.data(),static_cast<Eigen::Index>(candidate.fit->assembled_state->log_b.size()))).lpNorm<Eigen::Infinity>() :
            std::numeric_limits<double>::infinity();
        comparisons.push_back({{"against","LegacyCompact"},{"method",candidate.name},
            {"objective_difference",objective_difference},{"eta_inf_difference",eta_difference},
            {"scaled_ac_inf_difference",ScaledAcDifference(problem,reference,candidate)}});
    }
    return {{"topology",topology},{"atoms",atoms},{"rows",problem.Input().observations.size()},
        {"parameter_count",3*atoms},{"core_atoms",neighbor_policy.core_atoms},
        {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
        {"global_legacy_compact",PackFit("LegacyCompact",legacy,legacy_seconds)},
        {"global_operator_pcg",PackFit("OperatorPcg",operator_fit,operator_seconds)},
        {"fixed_neighbor",neighbor_json},{"comparisons",comparisons},{"peak_rss_mb",PeakRssMb()}};
}
void Write(const std::filesystem::path & path,const j::value & value)
{std::ofstream output(path); if(!output) throw std::runtime_error("Could not open experiment output."); output<<j::serialize(value)<<'\n';}
}
int main(int argc,char ** argv)
{
    try {
        if(argc!=5 || (std::string(argv[1])!="--case" && std::string(argv[1])!="--neighbor-only" &&
            std::string(argv[1])!="--decompose" && std::string(argv[1])!="--certified-local" &&
            std::string(argv[1])!="--qualification" && std::string(argv[1])!="--scaling-only" &&
            std::string(argv[1])!="--attribution" &&
            std::string(argv[1])!="--full-attribution" &&
            std::string(argv[1])!="--inexact-one" && std::string(argv[1])!="--inexact-two" &&
            std::string(argv[1])!="--neighbor-forward" && std::string(argv[1])!="--neighbor-reverse" &&
            std::string(argv[1])!="--scaling-forward" && std::string(argv[1])!="--scaling-reverse"))
            throw std::invalid_argument("Usage: joint_fixed_neighbor_experiment --case|--neighbor-only|--decompose|--certified-local|--qualification|--scaling-only|--attribution|--full-attribution|--inexact-one|--inexact-two|--neighbor-forward|--neighbor-reverse|--scaling-forward|--scaling-reverse OUTPUT_FILE TOPOLOGY ATOMS");
        Eigen::setNbThreads(1);
        const std::filesystem::path output_path(argv[2]);
        if(output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());
        const std::string mode(argv[1]);
        Write(output_path,Run(argv[3],std::stoi(argv[4]),output_path,mode=="--case",mode=="--decompose",
            mode=="--certified-local",mode=="--qualification",
            mode=="--scaling-only" || mode=="--attribution" || mode=="--scaling-forward" || mode=="--scaling-reverse",
            mode=="--neighbor-reverse" || mode=="--scaling-reverse",
            mode=="--neighbor-forward" || mode=="--neighbor-reverse" ||
                mode=="--scaling-forward" || mode=="--scaling-reverse",
            mode=="--attribution" || mode=="--full-attribution",
            mode=="--inexact-one" ? n::FixedNeighborLocalWork::OneAcceptedUpdate :
                mode=="--inexact-two" ? n::FixedNeighborLocalWork::TwoAcceptedUpdates :
                n::FixedNeighborLocalWork::Full));
        std::cout<<argv[3]<<'-'<<argv[4]<<" fixed-neighbor experiment complete\n";
        return 0;
    } catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
