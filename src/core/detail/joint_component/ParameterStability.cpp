#include "ParameterStability.hpp"
#include "ResourceWork.hpp"
#include "TiledQR.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}
double CpuSeconds(std::clock_t start)
{return static_cast<double>(std::clock()-start)/static_cast<double>(CLOCKS_PER_SEC);}

void Unavailable(JointBlockParameterStability & evidence,const std::string & reason)
{
    evidence.status="Unavailable"; evidence.reason=reason;
    evidence.prediction_status="Unavailable"; evidence.objective_consistency_status="Unavailable";
}
}

LocalParameterStabilityAssessment AssessLocalParameterStability(
    const Domain & domain,VectorRef y,const Vector & eta,const Vector & beta,const EvaluationContext & context,
    std::string block_id,std::size_t resident_bytes)
{
    const auto started=Clock::now(); const auto cpu_started=std::clock();
    LocalParameterStabilityAssessment out;
    auto & evidence=out.evidence; evidence.block_id=std::move(block_id);
    const auto finish=[&]() {
        out.seconds=Seconds(started); out.cpu_seconds=CpuSeconds(cpu_started);
        evidence.cpu_seconds=out.cpu_seconds; return out;
    };
    if(domain.atoms.size()>std::numeric_limits<std::size_t>::max()/2)
    {Unavailable(evidence,"size-overflow"); return finish();}
    MemoryShape shape; shape.component_atoms=domain.atoms.size();
    shape.component_observations=static_cast<std::size_t>(std::max<Eigen::Index>(domain.rows,0));
    shape.component_free_columns=2*domain.atoms.size(); shape.component_width_columns=domain.atoms.size();
    shape.maximum_block_rows=shape.component_observations; shape.maximum_neighbor_atoms=domain.atoms.size();
    shape.maximum_tile_rows=static_cast<std::size_t>(derivative_tile_rows);
    shape.component_stability_evaluations=2*MaximumLocalPerturbationProfileEvaluations;
    shape.resident_bytes=resident_bytes;
    for(const auto & support:domain.atoms)
        if(!CheckedAddBytes(shape.component_support_entries,support.size(),shape.component_support_entries))
        {shape.component_support_entries=std::numeric_limits<std::size_t>::max(); break;}
    shape.maximum_block_support_entries=shape.component_support_entries;
    out.memory_estimate=PlanMemory(shape,MemoryStage::ParameterStability);
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {
        Unavailable(evidence,std::string("resource-limit:parameter-stability:")+out.memory_estimate.reason);
        return finish();
    }
    if(domain.rows<=0 || y.size()!=domain.rows || eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) ||
        beta.size()!=2*eta.size())
    {Unavailable(evidence,"invalid-local-dimensions"); return finish();}
    if(!eta.allFinite() || !beta.allFinite() || !y.allFinite() || !(context.scale>0) || !std::isfinite(context.scale))
    {Unavailable(evidence,"nonfinite-local-input"); return finish();}

    const auto baseline=EvaluateState(domain,y,eta,beta,context,ProfileEvaluationRole::AcceptedEndpoint);
    if(!baseline.valid || !baseline.certificate.available || !baseline.certificate.feasible ||
        !baseline.certificate.kkt_passed)
    {Unavailable(evidence,"baseline-local-kkt-failed"); return finish();}
    const auto baseline_trust=CheckTrust(domain,y,baseline,context);
    if(!baseline_trust.passed || !baseline_trust.independent || !baseline_trust.independent->objective_passed)
    {Unavailable(evidence,"baseline-independent-replay-failed"); return finish();}
    const auto baseline_gradient=ProfileGradientInfinityNorm(baseline);
    if(!baseline_gradient || !(*baseline_gradient<=1e-12))
    {Unavailable(evidence,"baseline-width-stationarity-failed"); return finish();}
    if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
    {Unavailable(evidence,"cpu-time-budget-exhausted"); return finish();}
    const Vector baseline_prediction=(baseline.x*baseline.beta).eval();
    const double baseline_objective=baseline.certificate.objective/(context.scale*context.scale);
    Vector column_norms(baseline.x.cols());
    for(Eigen::Index k=0;k<column_norms.size();++k) column_norms(k)=baseline.x.col(k).norm();

    double maximum_width_change{},maximum_ac_change{},maximum_prediction_change{},maximum_objective_difference{};
    for(const double sign:{-1.0,1.0})
    {
        const Vector seed=eta.array()+sign*LocalWidthPerturbation;
        const Vector widths=seed.array().exp();
        if(!seed.allFinite() || !widths.allFinite() || (widths.array()<=0).any())
        {Unavailable(evidence,"unrepresentable-perturbation-width"); return finish();}
        const double remaining_cpu=MaximumLocalPerturbationCpuSeconds-CpuSeconds(cpu_started);
        if(!(remaining_cpu>0))
        {Unavailable(evidence,"cpu-time-budget-exhausted"); return finish();}
        auto bounded_context=context;
        bounded_context.profile_budget=static_cast<int>(MaximumLocalPerturbationProfileEvaluations);
        bounded_context.update_budget=6;
        bounded_context.profile_cpu_budget_seconds=remaining_cpu;
        const auto search=SearchProfile(domain,y,widths,bounded_context);
        ++evidence.perturbation_solves;
        evidence.perturbation_evaluations+=static_cast<std::size_t>(std::max(search.evaluations,0));
        if(search.stop_reason=="cpu-time-budget" ||
            CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
        {Unavailable(evidence,"cpu-time-budget-exhausted"); return finish();}
        if(search.eta.size()!=eta.size() || !search.eta.allFinite())
        {Unavailable(evidence,"perturbation-search-invalid"); return finish();}
        const auto candidate=EvaluateProfile(domain,y,search.eta,false,&bounded_context,nullptr,
            ProfileEvaluationRole::AcceptedEndpoint);
        ++evidence.perturbation_evaluations;
        if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
        {Unavailable(evidence,"cpu-time-budget-exhausted"); return finish();}
        if(!candidate.valid || !candidate.certificate.available || !candidate.certificate.feasible ||
            !candidate.certificate.kkt_passed)
        {Unavailable(evidence,"perturbation-local-kkt-failed"); return finish();}
        const auto trust=CheckTrust(domain,y,candidate,bounded_context);
        if(!trust.passed || !trust.independent || !trust.independent->objective_passed)
        {
            evidence.objective_consistency_status="Failed";
            Unavailable(evidence,"perturbation-independent-replay-failed");
            evidence.objective_consistency_status="Failed";
            return finish();
        }
        const auto gradient=ProfileGradientInfinityNorm(candidate);
        if(!gradient || !(*gradient<=1e-12))
        {Unavailable(evidence,"perturbation-width-stationarity-failed"); return finish();}
        if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
        {Unavailable(evidence,"cpu-time-budget-exhausted"); return finish();}
        const double width_change=(candidate.eta-eta).lpNorm<Eigen::Infinity>();
        const double ac_change=(((candidate.beta-beta).array()*column_norms.array()).abs()/context.scale).maxCoeff();
        const Vector candidate_prediction=(candidate.x*candidate.beta).eval();
        const double prediction_change=(candidate_prediction-baseline_prediction).norm()/
            std::max(1.0,baseline_prediction.norm());
        const double objective=candidate.certificate.objective/(context.scale*context.scale);
        maximum_width_change=std::max(maximum_width_change,width_change);
        maximum_ac_change=std::max(maximum_ac_change,ac_change);
        maximum_prediction_change=std::max(maximum_prediction_change,prediction_change);
        maximum_objective_difference=std::max(maximum_objective_difference,std::abs(objective-baseline_objective));
    }
    evidence.log_width_parameter_change_inf=maximum_width_change;
    evidence.scaled_ac_parameter_change_inf=maximum_ac_change;
    evidence.local_prediction_change_relative=maximum_prediction_change;
    evidence.objective_difference=maximum_objective_difference;
    evidence.prediction_status=maximum_prediction_change<=1e-8 ? "Stable" : "Changed";
    evidence.objective_consistency_status="Passed";
    if(maximum_prediction_change>1e-8)
        evidence.status="Unstable";
    else if(maximum_width_change<=1e-8 && maximum_ac_change<=1e-8)
        evidence.status="Stable";
    else evidence.status="Sensitive";
    evidence.reason="bounded-local-perturbation-assessed";
    return finish();
}
}
