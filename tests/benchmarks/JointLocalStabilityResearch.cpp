#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <boost/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace j=rhbm_gem::core::joint_component;
namespace json=boost::json;
using Clock=std::chrono::steady_clock;

json::value Number(double value)
{return std::isfinite(value) ? json::value(value) : json::value(nullptr);}
json::array Values(const j::Vector & values)
{
    json::array out; out.reserve(static_cast<std::size_t>(values.size()));
    for(Eigen::Index i=0;i<values.size();++i) out.push_back(Number(values(i)));
    return out;
}
template<class T> json::array Values(const std::vector<T> & values)
{
    json::array out; out.reserve(values.size());
    for(const auto & value:values) out.push_back(json::value(value));
    return out;
}
json::array Values(const j::Indices & values)
{
    json::array out; out.reserve(values.size());
    for(const auto value:values) out.push_back(static_cast<std::uint64_t>(value));
    return out;
}
json::object Certificate(const j::Certificate & value)
{
    json::array active;
    for(const auto atom:value.active_atoms) active.push_back(static_cast<std::uint64_t>(atom));
    return json::object{{"evaluated",value.evaluated},{"available",value.available},
        {"feasible",value.feasible},{"kkt_passed_strict_ac_1e-10",value.kkt_passed},
        {"linear_solver_valid",value.linear_solver_valid},{"linear_solver_reason",value.linear_solver_reason},
        {"primary_ac_kkt",Number(value.projected_kkt)},{"objective_raw",Number(value.objective)},
        {"residual_rmse",Number(value.residual_rmse)},{"active_ac_indices",active},
        {"free_rank",value.free_rank ? json::value(*value.free_rank) : json::value(nullptr)},
        {"width_bound_applied",value.width_bound_applied},{"width_bound_feasible",value.width_bound_feasible},
        {"primary_width_projected_kkt",Number(value.width_bound_projected_kkt)},
        {"active_width_bound_targets",value.active_width_bound_count},
        {"maximum_width_bound_violation",Number(value.maximum_width_bound_violation)},
        {"boundary_gradient_max",Number(value.boundary_gradient_max)},
        {"boundary_multiplier_max",Number(value.boundary_multiplier_max)},
        {"width_bound_complementarity",Number(value.width_bound_complementarity)}};
}
json::value Replay(const std::optional<j::IndependentReplayEvidence> & value)
{
    if(!value) return nullptr;
    const auto & e=*value;
    return json::object{{"available",e.available},{"valid",e.valid},{"feasible",e.feasible},
        {"passed_at_run_criteria",e.passed},{"reason",e.reason},
        {"primary_kkt_passed_at_run_criteria",e.primary_kkt_passed},
        {"independent_kkt_passed_at_run_criteria",e.independent_kkt_passed},
        {"prediction_passed",e.prediction_passed},{"residual_passed",e.residual_passed},
        {"objective_passed",e.objective_passed},{"kkt_agrees",e.kkt_agrees},
        {"gradient_agrees",e.gradient_agrees},{"primary_ac_kkt",Number(e.primary_ac_kkt)},
        {"independent_ac_kkt",Number(e.independent_ac_kkt)},
        {"primary_width_kkt",Number(e.primary_width_projected_kkt)},
        {"independent_width_kkt",Number(e.independent_width_projected_kkt)},
        {"prediction_difference",Number(e.prediction_difference)},
        {"residual_difference",Number(e.residual_difference)},
        {"objective_difference",Number(e.objective_difference)},
        {"independent_objective",Number(e.independent_objective)},
        {"gradient_difference",Number(e.gradient_difference)},
        {"maximum_gradient_difference_primary",Number(e.maximum_gradient_difference_primary)},
        {"maximum_gradient_difference_independent",Number(e.maximum_gradient_difference_independent)},
        {"maximum_gradient_difference_bound",Number(e.maximum_gradient_difference_bound)},
        {"maximum_gradient_difference_ratio",Number(e.maximum_gradient_difference_ratio)},
        {"width_bound_applied",e.width_bound_applied},{"width_bound_feasible",e.width_bound_feasible},
        {"maximum_width_bound_violation",Number(e.maximum_width_bound_violation)},
        {"active_width_bounds",e.active_width_bounds}};
}
json::array QualificationMatrix(const j::LocalStabilityTraceForTesting & trace)
{
    constexpr std::array<j::JointConvergenceCriteria,4> policies{{
        {1e-10,1e-12},{1e-10,1e-9},{1e-8,1e-9},{1e-7,1e-9}}};
    constexpr std::array<const char *,4> names{"Q0","Q1","Q2","Q3"};
    json::array out;
    for(std::size_t i=0;i<policies.size();++i)
    {
        const auto & q=policies[i];
        const auto & r=trace.replay;
        const bool replay_consistency=r && r->available && r->valid && r->feasible &&
            r->prediction_passed && r->residual_passed && r->objective_passed && r->gradient_agrees;
        const bool ac=trace.certificate.available && trace.certificate.feasible &&
            trace.certificate.projected_kkt<=q.global_ac_kkt && r && r->primary_valid && r->feasible &&
            r->independent_ac_kkt<=q.global_ac_kkt;
        const bool width=!r || !r->width_bound_applied || (r->width_bound_feasible &&
            r->primary_width_projected_kkt<=q.global_width_gradient &&
            r->independent_width_projected_kkt<=q.global_width_gradient);
        out.push_back(json::object{{"policy",names[i]},
            {"a_c_threshold",q.global_ac_kkt},{"width_threshold",q.global_width_gradient},
            {"primary_and_independent_optimality_pass",ac && width},
            {"replay_consistency_pass",replay_consistency},
            {"qualified",ac && width && replay_consistency}});
    }
    return out;
}
json::object TraceObject(const j::LocalStabilityTraceForTesting & t)
{
    json::object out{{"type","local-trace"},{"stage",t.stage},{"block_id",t.block_id},
        {"core_target_ids",Values(t.core_atom_ids)},
        {"direction",t.direction},{"sign",t.sign},{"direction_vector",Values(t.direction_vector)},
        {"baseline_eta",Values(t.baseline_eta)},{"baseline_beta",Values(t.baseline_beta)},
        {"seed_eta",Values(t.seed_eta)},{"seed_beta",Values(t.seed_beta)},
        {"requested_displacement",Values(t.requested_displacement)},
        {"effective_displacement",Values(t.effective_displacement)},
        {"endpoint_eta",Values(t.endpoint_eta)},{"endpoint_beta",Values(t.endpoint_beta)},
        {"baseline_certificate",Certificate(t.baseline_certificate)},
        {"baseline_replay",Replay(t.baseline_replay)},{"baseline_width_kkt",t.baseline_width_kkt ? Number(*t.baseline_width_kkt) : json::value(nullptr)},
        {"certificate",Certificate(t.certificate)},{"independent_replay",Replay(t.replay)},
        {"width_kkt",t.width_kkt ? Number(*t.width_kkt) : json::value(nullptr)},
        {"parameter_change",t.parameter_change ? Number(*t.parameter_change) : json::value(nullptr)},
        {"width_parameter_change",t.width_parameter_change ? Number(*t.width_parameter_change) : json::value(nullptr)},
        {"scaled_ac_parameter_change",t.scaled_ac_parameter_change ? Number(*t.scaled_ac_parameter_change) : json::value(nullptr)},
        {"prediction_change",t.prediction_change ? Number(*t.prediction_change) : json::value(nullptr)},
        {"objective_change",t.objective_change ? Number(*t.objective_change) : json::value(nullptr)},
        {"objective_difference",t.objective_difference ? Number(*t.objective_difference) : json::value(nullptr)},
        {"objective_trajectory",Values(t.objective_trajectory)},
        {"solver_stop_reason",t.stop_reason},{"solver_updates",t.accepted_updates},
        {"profile_evaluations",t.evaluations},{"budget_exhausted",t.budget_exhausted},
        {"feasible",t.feasible},{"width_bound_feasible",t.width_bound_feasible},
        {"qualification_status",t.qualification_status},{"wall_seconds",t.seconds},
        {"cpu_seconds",t.cpu_seconds}};
    out["qualification_matrix"]=QualificationMatrix(t);
    return out;
}
void Write(std::ofstream & out,const json::object & object)
{out<<json::serialize(object)<<'\n';}

struct Snapshot
{
    j::Domain domain;
    j::Vector y,eta,beta;
    j::EvaluationContext context;
    std::string block_id;
    std::vector<std::string> core_ids;
    j::Indices core_indices;
};
struct Scenario
{
    std::string name,topology{"chain"};
    int targets{};
    j::JointConvergenceCriteria search{j::PracticalJointConvergenceCriteria};
    j::JointConvergenceCriteria eligibility{j::PracticalJointConvergenceCriteria};
    bool custom_signal{};
    double signal_scale{1},noise_amplitude{},true_width{.5},true_a{2},true_c{.2};
};
Scenario SelectScenario(const std::string & name)
{
    if(name=="t0_3") return {name,"chain",3};
    if(name=="t0_30") return {name,"chain",30};
    if(name=="t4_solver_30") return {name,"chain",30,{1e-7,1e-5},j::PracticalJointConvergenceCriteria};
    if(name=="t4_full_30") return {name,"chain",30,{1e-7,1e-5},{1e-7,1e-5}};
    if(name=="t0_300") return {name,"chain",300};
    if(name=="t0_1000") return {name,"chain",1000};
    if(name=="t0_cube8") return {name,"cube",8};
    if(name=="t0_exact_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true};
    if(name=="t0_scale_small_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,.001};
    if(name=="t0_scale_large_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,1000};
    if(name=="t0_noise_low_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,1,1e-6};
    if(name=="t0_noise_high_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,1,1e-5};
    if(name=="noise_only_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,0,1e-5,.5,0,0};
    if(name=="t0_width08_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,1,0,.8};
    if(name=="t0_amp_small_3") return {name,"chain",3,j::PracticalJointConvergenceCriteria,
        j::PracticalJointConvergenceCriteria,true,1,0,.5,.2,.02};
    throw std::invalid_argument("unknown scenario");
}

double Seconds(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}
double CpuSeconds(std::clock_t start)
{return static_cast<double>(std::clock()-start)/static_cast<double>(CLOCKS_PER_SEC);}
std::size_t PeakRssBytes()
{
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
}
j::Vector Prediction(const j::Evaluation & e)
{return (e.x*e.beta).eval();}
j::Vector ColumnNorms(const j::Evaluation & e)
{
    j::Vector out(e.x.cols());
    for(Eigen::Index i=0;i<out.size();++i) out(i)=e.x.col(i).norm();
    return out;
}
json::object PredictionDifferenceStats(const j::Vector & lhs,const j::Vector & rhs,double denominator)
{
    const j::Vector difference=(lhs-rhs).cwiseAbs().eval();
    std::vector<double> ordered(difference.data(),difference.data()+difference.size());
    std::sort(ordered.begin(),ordered.end());
    const auto percentile=[&](double p) {
        if(ordered.empty()) return 0.0;
        const auto index=static_cast<std::size_t>(std::ceil(p*static_cast<double>(ordered.size()))-1.0);
        return ordered[std::min(index,ordered.size()-1)];
    };
    return {{"D_prediction",Number((lhs-rhs).norm()/std::max(1.0,denominator))},
        {"max_abs_prediction_difference",ordered.empty() ? json::value(nullptr) : Number(ordered.back())},
        {"rms_prediction_difference",Number(difference.norm()/std::sqrt(static_cast<double>(difference.size())))},
        {"p95_abs_prediction_difference",Number(percentile(.95))},
        {"p99_abs_prediction_difference",Number(percentile(.99))}};
}
json::object StateObject(const std::string & type,const Snapshot & s,const j::Evaluation & baseline,
    const j::Evaluation & endpoint,const std::optional<j::SearchResult> & search,
    const j::Evaluation * reference=nullptr,const j::JointConvergenceCriteria & eligibility=
        j::PracticalJointConvergenceCriteria)
{
    const auto baseline_prediction=Prediction(baseline),prediction=Prediction(endpoint);
    const double normalized_objective=endpoint.certificate.objective/(s.context.scale*s.context.scale);
    const double baseline_objective=baseline.certificate.objective/(s.context.scale*s.context.scale);
    const auto baseline_column_norms=ColumnNorms(baseline);
    const double width_change=(endpoint.eta-baseline.eta).lpNorm<Eigen::Infinity>();
    const double ac_change=(((endpoint.beta-baseline.beta).array()*baseline_column_norms.array()).abs()/
        s.context.scale).maxCoeff();
    auto prediction_stats=PredictionDifferenceStats(prediction,baseline_prediction,baseline_prediction.norm());
    const auto replay=j::CheckIndependentReplay(s.domain,s.y,endpoint,s.context,nullptr,nullptr,
        j::StrictJointConvergenceCriteria);
    const auto qualification=j::CheckIndependentReplay(s.domain,s.y,endpoint,s.context,nullptr,nullptr,eligibility);
    json::object out{{"type",type},{"block_id",s.block_id},{"core_target_ids",Values(s.core_ids)},
        {"baseline_eta",Values(baseline.eta)},{"baseline_beta",Values(baseline.beta)},
        {"starting_eta",Values(search ? search->initial.eta : endpoint.eta)},
        {"starting_profile_beta",Values(search ? search->initial.beta : endpoint.beta)},
        {"endpoint_eta",Values(endpoint.eta)},{"endpoint_beta",Values(endpoint.beta)},
        {"objective",Number(normalized_objective)},
        {"objective_change_vs_baseline",Number(std::abs(normalized_objective-baseline_objective))},
        {"D_B_vs_baseline",Number(width_change)},{"D_AC_vs_baseline",Number(ac_change)},
        {"prediction_change_vs_baseline",prediction_stats.at("D_prediction")},
        {"prediction_difference_stats_vs_baseline",prediction_stats},
        {"primary_ac_kkt",Number(endpoint.certificate.projected_kkt)},
        {"primary_width_kkt",Number(j::ProfileGradientInfinityNorm(endpoint,s.context).value_or(j::unavailable))},
        {"feasible",endpoint.certificate.feasible},{"active_ac_constraints",endpoint.certificate.active_atoms.size()},
        {"active_width_bound_targets",endpoint.certificate.active_width_bound_count},
        {"strict_independent_replay",Replay(replay)},{"eligibility_replay",Replay(qualification)}};
    if(search)
    {
        json::array trajectory;
        for(const auto & trial:search->trials)
            if(trial.endpoint.certificate.available)
                trajectory.push_back(Number(trial.endpoint.certificate.objective/(s.context.scale*s.context.scale)));
        out["solver_stop_reason"]=search->stop_reason;
        out["solver_updates"]=search->accepted;
        out["profile_evaluations"]=search->evaluations;
        out["budget_exhausted"]=search->stop_reason=="profile-budget" ||
            search->stop_reason=="accepted-update-budget" || search->stop_reason=="cpu-time-budget";
        out["objective_trajectory"]=trajectory;
    }
    if(reference && reference->valid && endpoint.valid)
    {
        const auto norms=ColumnNorms(*reference);
        const double d_ac=(((endpoint.beta-reference->beta).array()*norms.array()).abs()/s.context.scale).maxCoeff();
        const double d_b=(endpoint.eta-reference->eta).lpNorm<Eigen::Infinity>();
        const auto ref_prediction=Prediction(*reference);
        out["prediction_difference_stats_vs_reference"]=
            PredictionDifferenceStats(prediction,ref_prediction,ref_prediction.norm());
        out["D_B_solver_vs_reference"]=Number(d_b);
        out["D_AC_solver_vs_reference"]=Number(d_ac);
        out["D_prediction_solver_vs_reference"]=Number((prediction-ref_prediction).norm()/
            std::max(1.0,ref_prediction.norm()));
    }
    return out;
}

j::Evaluation RunSearch(Snapshot & s,const j::Vector & start_eta,double gtol,
    const std::string & type,const j::Evaluation & baseline,std::ofstream & out,
    const j::JointConvergenceCriteria & eligibility,j::ProfileSearchTolerancesForTesting controls={},
    const j::Evaluation * reference=nullptr,const std::string & direction={},
    const std::string & sign={},double delta=0.)
{
    const auto wall=Clock::now(); const auto cpu=std::clock();
    auto context=s.context; context.profile_budget=11; context.update_budget=6;
    context.profile_cpu_budget_seconds=0;
    controls.gtol=gtol;
    const auto search=j::SearchProfileForTesting(s.domain,s.y,start_eta.array().exp(),context,controls);
    j::Evaluation endpoint;
    if(search.eta.size()==start_eta.size())
        endpoint=j::EvaluateProfile(s.domain,s.y,search.eta,&context,j::ProfileEvaluationRole::AcceptedEndpoint);
    auto row=StateObject(type,s,baseline,endpoint,search.eta.size()==start_eta.size() ?
        std::optional<j::SearchResult>(search) : std::nullopt,reference,eligibility);
    row["lm_ftol"]=controls.ftol; row["lm_xtol"]=controls.xtol; row["lm_gtol"]=controls.gtol;
    if(type.rfind("N2-",0)==0)
    {
        row["initial_beta_passed_to_search"]=false;
        row["search_parameterization"]="log-width eta; beta is conditionally reprofiled at every evaluation";
    }
    row["direction"]=direction; row["sign"]=sign; row["delta"]=delta;
    row["requested_displacement"]=Values(start_eta-s.eta);
    row["wall_seconds"]=Seconds(wall); row["cpu_seconds"]=CpuSeconds(cpu);
    Write(out,row);
    return endpoint;
}

void WriteSnapshot(std::ofstream & out,const Snapshot & s,const std::string & scenario)
{
    json::array support;
    for(std::size_t atom=0;atom<s.domain.atoms.size();++atom)
    {
        json::array entries;
        for(const auto & item:s.domain.atoms[atom])
            entries.push_back(json::object{{"row",item.row},{"squared_distance",item.square}});
        support.push_back(entries);
    }
    Write(out,json::object{{"type","baseline-snapshot"},{"scenario",scenario},
        {"block_id",s.block_id},{"core_target_ids",Values(s.core_ids)},
        {"core_parent_indices",Values(s.core_indices)},{"rows",s.domain.rows},
        {"observation_scale",s.context.scale},{"eta_upper_bound",s.context.eta_upper_bound ?
            Number(*s.context.eta_upper_bound) : json::value(nullptr)},
        {"profile_budget",s.context.profile_budget},{"update_budget",s.context.update_budget},
        {"effective_response",Values(s.y)},{"eta",Values(s.eta)},{"beta",Values(s.beta)},
        {"atom_support",support}});
}

void RunControls(Snapshot & s,const rhbm_gem::JointBlockParameterStability * evidence,
    const j::JointConvergenceCriteria & eligibility,std::ofstream & out)
{
    const auto baseline=j::EvaluateState(s.domain,s.y,s.eta,s.beta,s.context,
        j::ProfileEvaluationRole::AcceptedEndpoint);
    const auto base_prediction=Prediction(baseline);
    const auto n0=j::EvaluateState(s.domain,s.y,s.eta,s.beta,s.context,
        j::ProfileEvaluationRole::AcceptedEndpoint);
    Write(out,StateObject("N0-identity",s,baseline,n0,std::nullopt,nullptr,eligibility));

    const auto n1=j::EvaluateProfile(s.domain,s.y,s.eta,&s.context,
        j::ProfileEvaluationRole::AcceptedEndpoint);
    Write(out,StateObject("N1-fixed-width-reprofile",s,baseline,n1,std::nullopt,nullptr,eligibility));

    auto strict_context=s.context; strict_context.profile_budget=200; strict_context.update_budget=100;
    const auto ref_search=j::SearchProfile(s.domain,s.y,s.eta.array().exp(),strict_context);
    const auto reference=j::EvaluateProfile(s.domain,s.y,ref_search.eta,&strict_context,
        j::ProfileEvaluationRole::AcceptedEndpoint);
    auto strict_row=StateObject("strict-local-reference",s,baseline,reference,ref_search,nullptr,eligibility);
    strict_row["solver_stop_reason"]=ref_search.stop_reason;
    strict_row["wall_seconds"]=ref_search.seconds; strict_row["accepted_updates"]=ref_search.accepted;
    Write(out,strict_row);
    const auto null_s0=RunSearch(s,s.eta,1e-12,"N2-profiled-null-S0",baseline,out,eligibility,{},&reference);
    const auto n1_prediction=Prediction(n1),n2_prediction=Prediction(null_s0);
    const double n1_objective=n1.certificate.objective/(s.context.scale*s.context.scale);
    const double n2_objective=null_s0.certificate.objective/(s.context.scale*s.context.scale);
    Write(out,json::object{{"type","N2-versus-N1"},{"block_id",s.block_id},
        {"D_prediction",Number((n2_prediction-n1_prediction).norm()/std::max(1.0,n1_prediction.norm()))},
        {"objective_difference",Number(std::abs(n2_objective-n1_objective))},
        {"D_B",Number((null_s0.eta-n1.eta).lpNorm<Eigen::Infinity>())},
        {"D_AC",Number((((null_s0.beta-n1.beta).array()*ColumnNorms(n1).array()).abs()/
            s.context.scale).maxCoeff())}});

    constexpr std::array<double,4> gtols{1e-12,1e-9,1e-8,1e-7};
    constexpr std::array<const char *,4> labels{"S0","S1","S2","S3"};
    for(std::size_t i=0;i<gtols.size();++i)
        RunSearch(s,s.eta,gtols[i],std::string("N2-null-")+labels[i],baseline,out,eligibility,{},&reference);
    (void)null_s0;
    for(const double offset:{-.25,.25})
    {
        auto start=s.eta.array()+offset;
        auto multi_context=s.context; multi_context.profile_budget=200; multi_context.update_budget=100;
        const auto multi_search=j::SearchProfile(s.domain,s.y,start.exp().eval(),multi_context);
        const auto multi_endpoint=j::EvaluateProfile(s.domain,s.y,multi_search.eta,&multi_context,
            j::ProfileEvaluationRole::AcceptedEndpoint);
        auto multi_row=StateObject(offset<0 ? "multi-start-low" : "multi-start-high",s,baseline,
            multi_endpoint,multi_search,&reference,eligibility);
        multi_row["initial_eta"]=Values(start.eval());
        multi_row["solver_stop_reason"]=multi_search.stop_reason;
        multi_row["accepted_updates"]=multi_search.accepted;
        Write(out,multi_row);
    }

    if(!evidence) return;
    struct Pair
    {
        std::optional<j::Evaluation> minus,plus;
        double minus_effective{},plus_effective{};
    };
    for(const auto & probe:evidence->probes)
    {
        if(probe.unit_direction.empty()) continue;
        const j::Vector direction=Eigen::Map<const j::Vector>(probe.unit_direction.data(),
            static_cast<Eigen::Index>(probe.unit_direction.size()));
        constexpr std::array<double,3> displacements{1e-4,1e-3,1e-2};
        for(const auto delta:displacements)
        {
            Pair strict_pair;
            for(const double sign:{-1.,1.})
            {
                j::Vector seed=s.eta+sign*delta*direction;
                if(s.context.eta_upper_bound) seed=seed.cwiseMin(
                    j::Vector::Constant(seed.size(),*s.context.eta_upper_bound));
                const std::string tag=sign<0 ? "minus" : "plus";
                std::optional<j::Evaluation> endpoint;
                if(delta==1e-3)
                {
                    auto reference_context=s.context; reference_context.profile_budget=200;
                    reference_context.update_budget=100;
                    const auto reference_search=j::SearchProfile(s.domain,s.y,seed.array().exp(),reference_context);
                    const auto perturb_reference=j::EvaluateProfile(s.domain,s.y,reference_search.eta,
                        &reference_context,j::ProfileEvaluationRole::AcceptedEndpoint);
                    auto reference_row=StateObject("strict-perturbation-reference-"+tag,s,baseline,
                        perturb_reference,reference_search,nullptr,eligibility);
                    reference_row["direction"]=probe.direction;
                    reference_row["delta"]=delta;
                    Write(out,reference_row);
                    endpoint=RunSearch(s,seed,1e-12,"perturb-"+tag+"-S0",baseline,out,
                        eligibility,{},&perturb_reference,probe.direction,tag,delta);
                    for(std::size_t i=1;i<gtols.size();++i)
                        RunSearch(s,seed,gtols[i],"perturb-"+tag+"-"+labels[i],baseline,out,
                            eligibility,{},&perturb_reference,probe.direction,tag,delta);
                }
                else endpoint=RunSearch(s,seed,1e-12,"perturb-"+tag+"-S0-delta-"+
                    std::to_string(delta),baseline,out,eligibility,{},nullptr,probe.direction,tag,delta);
                if(endpoint && endpoint->valid)
                {
                    if(sign<0)
                    {
                        strict_pair.minus=*endpoint;
                        strict_pair.minus_effective=(seed-s.eta).lpNorm<Eigen::Infinity>();
                    }
                    else
                    {
                        strict_pair.plus=*endpoint;
                        strict_pair.plus_effective=(seed-s.eta).lpNorm<Eigen::Infinity>();
                    }
                }
            }
            if(strict_pair.minus && strict_pair.plus)
            {
                const auto ym=Prediction(*strict_pair.minus),yp=Prediction(*strict_pair.plus);
                const double norm=std::max(1.0,base_prediction.norm());
                const double ratio=std::max(strict_pair.minus_effective,strict_pair.plus_effective)/
                    std::max(std::numeric_limits<double>::min(),
                        std::min(strict_pair.minus_effective,strict_pair.plus_effective));
                Write(out,json::object{{"type","paired-response"},{"block_id",s.block_id},
                    {"direction",probe.direction},{"delta",delta},
                    {"minus_effective_seed_displacement_inf",strict_pair.minus_effective},
                    {"plus_effective_seed_displacement_inf",strict_pair.plus_effective},
                    {"effective_displacements_comparable",ratio<=1.05},
                    {"D_common",Number(((yp+ym)/2.-base_prediction).norm()/norm)},
                    {"D_odd",Number(((yp-ym)/2.).norm()/norm)}});
            }
        }
    }
}

void Run(const Scenario & scenario,const std::string & output_path,const std::string & head)
{
    std::ofstream out(output_path,std::ios::trunc);
    if(!out) throw std::runtime_error("cannot open evidence output");
    out<<std::setprecision(17);
    auto input=second_stage_test::SyntheticJointWorkload(scenario.topology,scenario.targets);
    if(scenario.custom_signal)
    {
        std::fill(input.observations.begin(),input.observations.end(),0.0);
        for(const auto & atom:input.support)
            for(const auto & support:atom)
            {
                const auto basis=j::EvaluateKernel(support.squared_distance,scenario.true_width,2.5);
                input.observations[support.row]+=scenario.signal_scale*
                    (scenario.true_a*basis.gaussian+scenario.true_c*basis.charge);
            }
        for(std::size_t row=0;row<input.observations.size();++row)
            input.observations[row]+=scenario.noise_amplitude*
                std::sin(.371*static_cast<double>(row)+.23);
    }
    const auto fixture_hash=second_stage_test::SyntheticJointWorkloadHash(input);
    const rhbm_gem::core::JointProblem problem(input);
    const j::Vector initial_eta=j::Vector::Constant(scenario.targets,std::log(.55));
    j::FixedNeighborPolicy policy;
    policy.core_atoms=12; policy.maximum_sweeps=j::ProductionFixedNeighborMaximumSweeps;
    policy.convergence_criteria=scenario.search; policy.eta_upper_bound=std::log(2.2);
    policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
    policy.local_stability_eligibility_for_testing=scenario.eligibility;
    policy.claim_parameter_stability_block=[](int) {return true;};
    std::vector<Snapshot> snapshots;
    policy.local_stability_baseline_observer=[&](const j::Domain & domain,j::VectorRef y,
        const j::Vector & eta,const j::Vector & beta,const j::EvaluationContext & context,
        const std::string & block_id,const j::Indices & core_indices) {
        Snapshot snapshot{domain,j::Vector(y),eta,beta,context,block_id,context.atom_ids,core_indices};
        WriteSnapshot(out,snapshot,scenario.name);
        snapshots.push_back(std::move(snapshot));
    };
    policy.local_stability_trace_observer=[&](const j::LocalStabilityTraceForTesting & trace) {
        Write(out,TraceObject(trace));
    };
    const auto wall=Clock::now(); const auto cpu=std::clock();
    const auto result=j::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    auto no_stability_policy=policy;
    no_stability_policy.assess_parameter_stability=false;
    no_stability_policy.local_stability_baseline_observer={};
    no_stability_policy.local_stability_trace_observer={};
    no_stability_policy.claim_parameter_stability_block={};
    const auto no_stability=j::RunFixedNeighborDiagnostic(problem,initial_eta,no_stability_policy);
    const bool eta_equal=(result.state.eta.array()==no_stability.state.eta.array()).all();
    const bool beta_equal=(result.state.beta.array()==no_stability.state.beta.array()).all();
    const bool prediction_equal=(result.state.prediction.array()==no_stability.state.prediction.array()).all();
    const bool residual_equal=(result.state.residual.array()==no_stability.state.residual.array()).all();
    json::object summary{{"type","global-summary"},{"scenario",scenario.name},
        {"fixture_hash",fixture_hash},{"head",head},{"topology",scenario.topology},
        {"targets",scenario.targets},{"custom_signal",scenario.custom_signal},
        {"signal_scale",scenario.signal_scale},{"noise_amplitude",scenario.noise_amplitude},
        {"true_width",scenario.true_width},{"true_a",scenario.true_a*scenario.signal_scale},
        {"true_c",scenario.true_c*scenario.signal_scale},
        {"search_ac_kkt_threshold",scenario.search.global_ac_kkt},
        {"search_width_kkt_threshold",scenario.search.global_width_gradient},
        {"local_eligibility_ac_threshold",scenario.eligibility.global_ac_kkt},
        {"local_eligibility_width_threshold",scenario.eligibility.global_width_gradient},
        {"sweeps",result.sweep_count},{"accepted_blocks",result.final_sweep ? result.final_sweep->accepted_blocks : 0},
        {"global_ac_kkt",result.final_sweep ? Number(result.final_sweep->global_ac_kkt) : json::value(nullptr)},
        {"global_width_kkt",result.final_sweep ? Number(result.final_sweep->projected_width_gradient_inf_norm) : json::value(nullptr)},
        {"objective",Number(result.state.objective)},{"search_converged",result.search_converged},
        {"endpoint_certified",result.endpoint_certified},{"independent_replay_passed",result.assessment.independent_replay.passed},
        {"stop_reason",result.reason},{"selected_snapshots",snapshots.size()},
        {"wall_seconds",Seconds(wall)},{"cpu_seconds",CpuSeconds(cpu)},
        {"process_peak_rss_bytes",PeakRssBytes()},{"build_configuration","Release; -O3 -DNDEBUG"},
        {"thread_policy","OMP/BLAS/MKL/VECLIB=1"},
        {"state_without_stability_exactly_equal",eta_equal && beta_equal && prediction_equal && residual_equal &&
            result.state.objective==no_stability.state.objective},
        {"state_without_stability_eta_equal",eta_equal},{"state_without_stability_beta_equal",beta_equal},
        {"state_without_stability_prediction_equal",prediction_equal},
        {"state_without_stability_residual_equal",residual_equal},
        {"sweeps_without_stability",no_stability.sweep_count},
        {"objective_without_stability",Number(no_stability.state.objective)},
        {"stability_seconds",result.work.parameter_stability_seconds},
        {"perturbation_solves",result.work.perturbation_solves},
        {"perturbation_evaluations",result.work.perturbation_evaluations},
        {"stability_workspace_peak_bytes",result.work.parameter_stability_peak_memory_estimate_bytes},
        {"local_conditioning_blocks",result.local_conditioning.size()}};
    if(result.endpoint_certified && result.assessment.independent_replay.passed)
    {
        summary["accepted_eta"]=Values(result.state.eta);
        summary["accepted_beta"]=Values(result.state.beta);
        summary["endpoint_certification_reason"]=result.endpoint_trust.reason;
        summary["independent_prediction_difference"]=Number(result.assessment.independent_replay.prediction_difference);
        summary["independent_residual_difference"]=Number(result.assessment.independent_replay.residual_difference);
        summary["independent_objective_difference"]=Number(result.assessment.independent_replay.objective_difference);
        summary["global_endpoint_objective"]=Number(result.assessment.independent_replay.independent_objective);
        j::Vector truth_prediction= j::Vector::Zero(result.state.prediction.size());
        for(const auto & atom:input.support)
            for(const auto & support:atom)
            {
                const auto basis=j::EvaluateKernel(support.squared_distance,scenario.true_width,2.5);
                truth_prediction(static_cast<Eigen::Index>(support.row))+=scenario.signal_scale*
                    (scenario.true_a*basis.gaussian+scenario.true_c*basis.charge);
            }
        const j::Vector prediction_error=result.state.prediction-truth_prediction;
        summary["max_width_truth_error"]=Number((result.state.eta.array()-std::log(scenario.true_width)).abs().maxCoeff());
        double max_a_error{},max_c_error{};
        for(Eigen::Index atom=0;atom<scenario.targets;++atom)
        {
            max_a_error=std::max(max_a_error,
                std::abs(result.state.beta(2*atom)-scenario.true_a*scenario.signal_scale));
            max_c_error=std::max(max_c_error,
                std::abs(result.state.beta(2*atom+1)-scenario.true_c*scenario.signal_scale));
        }
        summary["max_a_truth_error"]=Number(max_a_error);
        summary["max_c_truth_error"]=Number(max_c_error);
        summary["prediction_truth_rmse"]=Number(prediction_error.norm()/
            std::sqrt(static_cast<double>(prediction_error.size())));
    }
    json::array blocks;
    for(const auto & block:result.local_conditioning)
        blocks.push_back(json::object{{"block_id",block.block_id},{"classification",block.conditioning_classification},
            {"free_ac_rank",block.free_ac_rank ? json::value(*block.free_ac_rank) : json::value(nullptr)},
            {"free_ac_columns",block.free_ac_columns},{"projected_width_rank",block.projected_width_rank ?
                json::value(*block.projected_width_rank) : json::value(nullptr)},
            {"projected_width_condition",block.projected_width_condition ? Number(*block.projected_width_condition) : json::value(nullptr)},
            {"active_width_bound_targets",block.active_width_bound_targets ?
                json::value(*block.active_width_bound_targets) : json::value(nullptr)}});
    summary["local_conditioning"]=blocks;
    Write(out,summary);
    if(result.parameter_stability)
        for(auto & snapshot:snapshots)
        {
            const auto found=std::find_if(result.parameter_stability->blocks.begin(),
                result.parameter_stability->blocks.end(),[&](const auto & block) {
                    return block.block_id==snapshot.block_id;
                });
            const auto * evidence=found==result.parameter_stability->blocks.end() ? nullptr : &*found;
            if(evidence)
                Write(out,json::object{{"type","block-classification"},{"block_id",snapshot.block_id},
                    {"status",evidence->status},{"reason",evidence->reason},
                    {"prediction_status",evidence->prediction_status},
                    {"baseline_eligibility_status",evidence->baseline_eligibility_status},
                    {"baseline_eligibility_reason",evidence->baseline_eligibility_reason},
                    {"baseline_ac_kkt",evidence->baseline_local_ac_kkt ? Number(*evidence->baseline_local_ac_kkt) : json::value(nullptr)},
                    {"baseline_width_kkt",evidence->baseline_projected_width_kkt ? Number(*evidence->baseline_projected_width_kkt) : json::value(nullptr)},
                    {"probes",evidence->probes.size()},{"perturbation_seeds",evidence->perturbation_seeds},
                    {"perturbation_solves",evidence->perturbation_solves},
                    {"perturbation_evaluations",evidence->perturbation_evaluations},
                    {"parameter_change",evidence->log_width_parameter_change_inf ? Number(*evidence->log_width_parameter_change_inf) : json::value(nullptr)},
                    {"scaled_ac_change",evidence->scaled_ac_parameter_change_inf ? Number(*evidence->scaled_ac_parameter_change_inf) : json::value(nullptr)},
                    {"prediction_change",evidence->local_prediction_change_relative ? Number(*evidence->local_prediction_change_relative) : json::value(nullptr)}});
            RunControls(snapshot,evidence,scenario.eligibility,out);
        }
    out.flush();
    if(!out) throw std::runtime_error("failed to write evidence output");
}
}

int main(int argc,char ** argv)
{
    if(argc!=4)
    {
        std::cerr<<"usage: rhbm_local_stability_research <scenario> <jsonl-path> <head>\n";
        return 2;
    }
    try
    {
        Run(SelectScenario(argv[1]),argv[2],argv[3]);
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<"local stability research failed: "<<error.what()<<'\n';
        return 1;
    }
}
