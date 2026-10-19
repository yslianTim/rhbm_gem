#include "ParameterStability.hpp"
#include "ResourceWork.hpp"
#include "TiledDerivative.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <Eigen/SVD>
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

struct ProbeDirection
{
    std::string direction,source;
    Vector unit;
};

void Canonicalize(Vector & direction)
{
    Eigen::Index pivot{};
    direction.cwiseAbs().maxCoeff(&pivot);
    if(direction(pivot)<0) direction=-direction;
}

void AddDirection(std::vector<ProbeDirection> & directions,std::string name,std::string source,Vector direction)
{
    const double norm=direction.norm();
    if(!(norm>0) || !std::isfinite(norm) || !direction.allFinite()) return;
    direction/=norm;
    Canonicalize(direction);
    for(auto & existing:directions)
        if(std::abs(existing.unit.dot(direction))>=1.-1e-12)
        {
            existing.source+="; equivalent-to-"+name+" ("+source+")";
            return;
        }
    directions.push_back({std::move(name),std::move(source),std::move(direction)});
}

struct DirectionSet
{
    std::vector<ProbeDirection> available;
    bool weak_direction_available{};
    std::string weak_direction_reason{"projected-width-derivative-unavailable"};
};

DirectionSet BuildDirections(const Evaluation & baseline,const Vector & eta,const EvaluationContext & context)
{
    DirectionSet out;
    const auto p=eta.size();
    if(p<=0) return out;

    Vector weak;
    Vector column_norms;
    const auto prepared=PrepareProjectedWidthDerivative(baseline,context.scale,&context,-1,
        derivative_tile_rows,true);
    if(prepared.valid)
    {
        const auto reduced=ReduceProjectedWidthDerivative(prepared);
        if(reduced.valid && reduced.projected.cols()==p && reduced.projected.rows()>0)
        {
            column_norms=reduced.projected_norms;
            Matrix normalized=reduced.projected;
            for(Eigen::Index column=0;column<p;++column)
            {
                const double norm=normalized.col(column).norm();
                if(norm>0) normalized.col(column)/=norm;
            }
            if(normalized.allFinite() && normalized.norm()>0)
            {
                Eigen::JacobiSVD<Matrix> svd(normalized,Eigen::ComputeFullV);
                if(svd.info()==Eigen::Success && svd.matrixV().rows()==p && svd.matrixV().cols()==p &&
                    svd.matrixV().allFinite())
                {
                    weak=svd.matrixV().col(p-1);
                    out.weak_direction_available=true;
                    out.weak_direction_reason="smallest-right-singular-vector-of-column-normalized-projected-width-jacobian";
                }
            }
            if(column_norms.size()!=p || !column_norms.allFinite()) column_norms=Vector::Zero(p);
        }
    }
    if(weak.size()==p)
    {
        Eigen::Index coordinate{};
        weak.cwiseAbs().maxCoeff(&coordinate);
        Vector axis=Vector::Zero(p); axis(coordinate)=1.;
        AddDirection(out.available,"weak-singular",out.weak_direction_reason,std::move(weak));
        AddDirection(out.available,"coordinate","largest-weak-vector-component",std::move(axis));
    }
    else
    {
        Eigen::Index coordinate{};
        if(column_norms.size()==p && column_norms.allFinite()) column_norms.minCoeff(&coordinate);
        Vector axis=Vector::Zero(p); axis(coordinate)=1.;
        AddDirection(out.available,"coordinate","smallest-projected-width-column-norm",std::move(axis));
    }
    AddDirection(out.available,"common-mode","normalized-all-width-direction",Vector::Ones(p));
    return out;
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
    shape.component_stability_evaluations=LocalPerturbationSeedsPerDirection*
        MaximumLocalPerturbationDirectionsPerBlock*MaximumLocalPerturbationProfileEvaluations;
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
    const auto baseline_gradient=ProfileGradientInfinityNorm(baseline,context);
    if(!baseline_gradient || !(*baseline_gradient<=1e-12))
    {Unavailable(evidence,"baseline-width-stationarity-failed"); return finish();}
    if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
    {Unavailable(evidence,"cpu-time-budget-exhausted"); return finish();}
    const Vector baseline_prediction=(baseline.x*baseline.beta).eval();
    const double baseline_objective=baseline.certificate.objective/(context.scale*context.scale);
    Vector column_norms(baseline.x.cols());
    for(Eigen::Index k=0;k<column_norms.size();++k) column_norms(k)=baseline.x.col(k).norm();

    const auto directions=BuildDirections(baseline,eta,context);
    if(!directions.weak_direction_available && eta.size()>1)
        evidence.probes.push_back({"weak-singular",directions.weak_direction_reason,"Unavailable",{},0,{},{},{}});
    double maximum_width_change{},maximum_ac_change{},maximum_prediction_change{},maximum_objective_difference{};
    bool all_directions_assessed=directions.weak_direction_available || eta.size()==1;
    bool any_search_failure{};
    for(const auto & direction:directions.available)
    {
        JointStabilityProbe probe;
        probe.direction=direction.direction; probe.source=direction.source;
        probe.unit_direction.assign(direction.unit.data(),direction.unit.data()+direction.unit.size());
        double direction_parameter_change{},direction_prediction_change{},direction_objective_change{};
        bool direction_assessed=true;
        std::string failure_reason;
        for(const double sign:{-1.0,1.0})
        {
            if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
            {direction_assessed=false; failure_reason="cpu-time-budget-exhausted"; break;}
            Vector seed=eta+sign*LocalWidthPerturbation*direction.unit;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(context.eta_upper_bound_for_research)
                seed=seed.cwiseMin(Vector::Constant(seed.size(),*context.eta_upper_bound_for_research));
#endif
            const Vector widths=seed.array().exp();
            if(!seed.allFinite() || !widths.allFinite() || (widths.array()<=0).any())
            {direction_assessed=false; failure_reason="unrepresentable-perturbation-width"; break;}
            const double remaining_cpu=MaximumLocalPerturbationCpuSeconds-CpuSeconds(cpu_started);
            if(!(remaining_cpu>0))
            {direction_assessed=false; failure_reason="cpu-time-budget-exhausted"; break;}
            auto bounded_context=context;
            bounded_context.profile_budget=static_cast<int>(MaximumLocalPerturbationProfileEvaluations-1);
            bounded_context.update_budget=6;
            bounded_context.profile_cpu_budget_seconds=remaining_cpu;
            const auto search=SearchProfile(domain,y,widths,bounded_context);
            ++evidence.perturbation_solves;
            const auto search_evaluations=static_cast<std::size_t>(std::max(search.evaluations,0));
            evidence.perturbation_evaluations+=search_evaluations;
            probe.evaluations+=search_evaluations;
            if(search.stop_reason=="cpu-time-budget" ||
                CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
            {direction_assessed=false; failure_reason="cpu-time-budget-exhausted"; break;}
            if(search.eta.size()!=eta.size() || !search.eta.allFinite())
            {direction_assessed=false; failure_reason="perturbation-search-invalid"; break;}
            const auto candidate=EvaluateProfile(domain,y,search.eta,false,&bounded_context,nullptr,
                ProfileEvaluationRole::AcceptedEndpoint);
            ++evidence.perturbation_evaluations; ++probe.evaluations;
            if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
            {direction_assessed=false; failure_reason="cpu-time-budget-exhausted"; break;}
            if(!candidate.valid || !candidate.certificate.available || !candidate.certificate.feasible ||
                !candidate.certificate.kkt_passed)
            {direction_assessed=false; failure_reason="perturbation-local-kkt-failed"; break;}
            const auto trust=CheckTrust(domain,y,candidate,bounded_context);
            if(!trust.passed || !trust.independent || !trust.independent->objective_passed)
            {
                direction_assessed=false; failure_reason="perturbation-independent-replay-failed";
                evidence.objective_consistency_status="Failed";
                break;
            }
            const auto gradient=ProfileGradientInfinityNorm(candidate,bounded_context);
            if(!gradient || !(*gradient<=1e-12))
            {direction_assessed=false; failure_reason="perturbation-width-stationarity-failed"; break;}
            if(CpuSeconds(cpu_started)>=MaximumLocalPerturbationCpuSeconds)
            {direction_assessed=false; failure_reason="cpu-time-budget-exhausted"; break;}
            const double width_change=(candidate.eta-eta).lpNorm<Eigen::Infinity>();
            const double ac_change=(((candidate.beta-beta).array()*column_norms.array()).abs()/context.scale).maxCoeff();
            const Vector candidate_prediction=(candidate.x*candidate.beta).eval();
            const double prediction_change=(candidate_prediction-baseline_prediction).norm()/
                std::max(1.0,baseline_prediction.norm());
            const double objective=candidate.certificate.objective/(context.scale*context.scale);
            direction_parameter_change=std::max(direction_parameter_change,std::max(width_change,ac_change));
            direction_prediction_change=std::max(direction_prediction_change,prediction_change);
            direction_objective_change=std::max(direction_objective_change,std::abs(objective-baseline_objective));
            maximum_width_change=std::max(maximum_width_change,width_change);
            maximum_ac_change=std::max(maximum_ac_change,ac_change);
            maximum_prediction_change=std::max(maximum_prediction_change,prediction_change);
            maximum_objective_difference=std::max(maximum_objective_difference,std::abs(objective-baseline_objective));
        }
        if(direction_assessed)
        {
            probe.status="Assessed";
            probe.parameter_change=direction_parameter_change;
            probe.prediction_change=direction_prediction_change;
            probe.objective_change=direction_objective_change;
        }
        else
        {
            any_search_failure=true; all_directions_assessed=false;
            probe.status=failure_reason=="cpu-time-budget-exhausted" ? "BudgetExhausted" : "SolverFailed";
            probe.source+="; "+failure_reason;
            if(evidence.reason.empty()) evidence.reason=failure_reason;
        }
        evidence.probes.push_back(std::move(probe));
    }
    evidence.perturbation_directions=static_cast<std::size_t>(std::count_if(evidence.probes.begin(),
        evidence.probes.end(),[](const JointStabilityProbe & probe){return !probe.unit_direction.empty();}));
    evidence.perturbation_seeds=evidence.perturbation_solves;
    if(!directions.weak_direction_available && evidence.reason.empty())
        evidence.reason=directions.weak_direction_reason;
    const bool objective_failed=evidence.objective_consistency_status=="Failed";
    evidence.objective_consistency_status=objective_failed ? "Failed" : "Passed";
    if(evidence.perturbation_solves==0)
    {Unavailable(evidence,evidence.reason.empty() ? "no-local-directions-available" : evidence.reason); return finish();}
    evidence.log_width_parameter_change_inf=maximum_width_change;
    evidence.scaled_ac_parameter_change_inf=maximum_ac_change;
    evidence.local_prediction_change_relative=maximum_prediction_change;
    evidence.objective_difference=maximum_objective_difference;
    evidence.prediction_status=maximum_prediction_change<=LocalPredictionChangeWarning ? "Stable" : "Changed";
    if(maximum_prediction_change>LocalPredictionChangeWarning)
        evidence.status="Unstable";
    else if(any_search_failure || !all_directions_assessed)
        evidence.status="Unavailable";
    else if(maximum_width_change<=LocalParameterChangeWarning && maximum_ac_change<=LocalParameterChangeWarning)
        evidence.status="Stable";
    else evidence.status="Sensitive";
    if(evidence.reason.empty()) evidence.reason="bounded-local-directional-perturbation-assessed";
    return finish();
}
}
