#include "Numerics.hpp"
#include "LinearSolve.hpp"
#include "CompactSvd.hpp"
#include "InstrumentedLM.hpp"
#include "TiledDerivative.hpp"
#include "ActiveSetBoundedLM.hpp"
#include <algorithm>
#include <cmath>

namespace rhbm_gem::core::joint_component {
WeakWidthStepAssessment AssessWeakWidthStep(const Matrix & jacobian,VectorRef step,const RankPolicy & policy)
{
    WeakWidthStepAssessment out;
    if(jacobian.rows()!=jacobian.cols() || jacobian.cols()!=step.size() ||
        jacobian.rows()==0 || !jacobian.allFinite() || !step.allFinite()) return out;
    const auto svd=EvaluateRank(jacobian,{policy,jacobian.cols(),-1,RankBoundary::StrictGreater},
        CompactSvdVectors::Right);
    if(!svd.valid || svd.right_vectors.cols()!=step.size() || svd.singular_values.size()!=step.size()) return out;
    const auto weak=step.size()-svd.rank;
    const double norm=step.norm();
    if(!(norm>0) || !std::isfinite(norm)) return out;
    const double weak_norm=weak ? (svd.right_vectors.rightCols(weak).transpose()*step).norm() : 0.0;
    out.weak_step_fraction=weak_norm/norm; out.available=std::isfinite(out.weak_step_fraction);
    return out;
}
bool RejectWeakWidthTrial(double step_inf,double weak_step_fraction,double objective_before,double objective_after)
{
    constexpr double large_log_width_step=2.302585092994046; // one decade in width
    constexpr double weak_direction_fraction=.95;
    if(!std::isfinite(step_inf) || step_inf<=large_log_width_step ||
        !std::isfinite(weak_step_fraction) || weak_step_fraction<weak_direction_fraction ||
        !std::isfinite(objective_before) || !std::isfinite(objective_after)) return false;
    const double improvement=objective_before-objective_after;
    const double replay_enclosure=1e-12+2e-12*std::abs(objective_before);
    return improvement>=0 && improvement<=replay_enclosure;
}
namespace {
struct Profile
{
    const Domain & domain;
    VectorRef y;
    const EvaluationContext & context;
    Evaluation cached;
    std::vector<joint_component::Trial> trace;
    int evaluations{},derivatives{};
    std::string failure;
    std::optional<double> accepted_objective;
    std::optional<ActiveSetStep> last_bounded_subproblem;
    BoundedLmWork bounded_lm_work;
    LinearSolveWorkspace linear_workspace;
    void RecordBoundedSubproblem(const ActiveSetStep & step)
    {
        ++bounded_lm_work.subproblem_solves;
        bounded_lm_work.damping_solves+=step.damping_solves;
        bounded_lm_work.bracket_solves+=step.bracket_solves;
        bounded_lm_work.damping_qr_factorizations+=step.damping_qr_factorizations;
        if(!step.valid) {++bounded_lm_work.failed_subproblems; return;}
        bounded_lm_work.active_set_iterations+=step.iterations;
        bounded_lm_work.releases+=step.releases;
        for(std::size_t i=0;i<step.active.size();++i)
            if(step.active[i]) ++bounded_lm_work.active_bound_occurrences;
        if(std::isfinite(step.kkt_inf_norm))
            bounded_lm_work.maximum_subproblem_kkt=std::isfinite(bounded_lm_work.maximum_subproblem_kkt) ?
                std::max(bounded_lm_work.maximum_subproblem_kkt,step.kkt_inf_norm) : step.kkt_inf_norm;
    }
    ActiveSetStep solve_bounded_step(const Matrix & jacobian,VectorRef response,VectorRef eta,
        VectorRef diagonal,double radius,double damping)
    {
        ActiveSetStep failed_step;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(context.force_width_subproblem_failure_for_testing)
        {failed_step.status="injected-width-subproblem-failure"; last_bounded_subproblem=failed_step;
            RecordBoundedSubproblem(failed_step); return failed_step;}
#endif
        if(!context.eta_upper_bound || !std::isfinite(*context.eta_upper_bound))
        {failed_step.status="invalid-width-upper-bound"; last_bounded_subproblem=failed_step;
            RecordBoundedSubproblem(failed_step); return failed_step;}
        for(Eigen::Index i=0;i<eta.size();++i)
            if(eta(i)>*context.eta_upper_bound)
            {failed_step.status="infeasible-current-width"; last_bounded_subproblem=failed_step;
                RecordBoundedSubproblem(failed_step); return failed_step;}
        const Vector upper=Vector::Constant(eta.size(),*context.eta_upper_bound)-eta;
        auto result=SolveBoundedTrustRegionStep(jacobian,response,diagonal,upper,radius,damping);
        last_bounded_subproblem=result;
        RecordBoundedSubproblem(result);
        return result;
    }
    bool WidthStationary() const
    {
        return cached.valid && cached.certificate.available && cached.certificate.feasible &&
            cached.certificate.kkt_passed && cached.certificate.width_bound_applied &&
            cached.certificate.width_bound_feasible &&
            cached.certificate.width_bound_projected_kkt<=PracticalJointConvergenceCriteria.global_width_gradient;
    }
    bool retry() const {return evaluations<context.profile_budget && failure!="unrepresentable-step";}
    bool Trial(const Vector & accepted,const Vector & step,const Vector & diagonal,
        double radius,double damping,double actual,double predicted,double ratio,bool proposed)
    {
        if(trace.empty()) return false;
        auto & row=trace.back();
        LmTrial lm;
        lm.accepted_eta=accepted; lm.step=step; lm.diagonal=diagonal;
        lm.radius=radius; lm.damping=damping;
        lm.actual_decrease=cached.valid ? actual : unavailable;
        lm.predicted_decrease=predicted; lm.ratio=cached.valid ? ratio : unavailable;
        lm.proposed_acceptance=proposed;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(last_bounded_subproblem)
        {
            lm.bounded_subproblem=true;
            lm.subproblem_status=last_bounded_subproblem->status;
            lm.subproblem_kkt_inf_norm=last_bounded_subproblem->kkt_inf_norm;
            lm.active_set_iterations=last_bounded_subproblem->iterations;
            lm.active_set_releases=last_bounded_subproblem->releases;
            for(std::size_t i=0;i<last_bounded_subproblem->active.size();++i)
            {
                if(!last_bounded_subproblem->active[i]) continue;
                lm.active_bound_indices.push_back(i);
                lm.active_bound_multipliers.push_back(last_bounded_subproblem->multipliers(
                    static_cast<Eigen::Index>(i)));
            }
        }
#endif
        last_bounded_subproblem.reset();
        row.lm=std::move(lm);
        bool trusted=cached.valid;
        if(proposed)
        {
            auto evidence=CheckReplay(domain,y,cached,context);
            trusted=evidence.passed; row.trust=std::move(evidence);
        }
        if(proposed && trusted && step.lpNorm<Eigen::Infinity>()>2.302585092994046)
        {
            const auto prepared=PrepareDerivative(cached,context.scale,&context);
            const auto differential=ReduceDerivative(prepared,cached.residual);
            if(differential.valid)
            {
                const auto assessment=AssessWeakWidthStep(differential.jacobian,step,context.rank);
                if(assessment.available)
                {
                    row.weak_width_step_fraction=assessment.weak_step_fraction;
                    const double candidate_objective=cached.certificate.objective/(context.scale*context.scale);
                    const double old_objective=accepted_objective.value_or(candidate_objective);
                    if(RejectWeakWidthTrial(step.lpNorm<Eigen::Infinity>(),assessment.weak_step_fraction,
                        old_objective,candidate_objective))
                    {
                        row.rejection_reason="weak-width-identifiability-limited";
                        trusted=false;
                    }
                }
            }
        }
        if(!proposed || !trusted)
            UpdateAcceptedProfileObjective(accepted_objective,cached,context.scale,false);
        if(!trusted) failure=row.rejection_reason.empty() ? "untrusted-trial" : row.rejection_reason;
        return trusted;
    }
    int values() const {return static_cast<int>(domain.rows);}
    bool Get(const Vector & eta)
    {
        if(cached.valid && cached.eta.size()==eta.size() && (cached.eta.array()==eta.array()).all()) return true;
        if(evaluations>=context.profile_budget) {failure="profile-budget"; return false;}
        failure.clear();
        cached=EvaluateProfile(domain,y,eta,&context,&linear_workspace); ++evaluations;
        joint_component::Trial row; row.endpoint=cached; row.evaluation=evaluations;
        trace.push_back(std::move(row));
        if(!cached.valid) failure="inner-"+cached.reason;
        return cached.valid;
    }
    int operator()(const Vector & eta,Vector & residual)
    {
        if(!context.eta_upper_bound || !std::isfinite(*context.eta_upper_bound))
        {failure="invalid-width-upper-bound"; return -1;}
        if((eta.array()>*context.eta_upper_bound).any())
        {failure="bounded-step-infeasible"; return -1;}
        if(!Get(eta)) return -1;
        residual=cached.residual/context.scale; return 0;
    }
    int linearize(const Vector & eta,const Vector &,Matrix & factor,Vector & response,Vector & norms)
    {
        if(!Get(eta)) return -1;
        const auto prepared=PrepareDerivative(cached,context.scale,&context);
        auto differential=ReduceDerivative(prepared,cached.residual); ++derivatives;
        if(!differential.valid)
        {failure=differential.reason; return -1;}
        factor=std::move(differential.jacobian); response=std::move(differential.response);
        norms=std::move(differential.jacobian_norms); return 0;
    }
    void Accept(const Vector & eta,int update)
    {
        UpdateAcceptedProfileObjective(accepted_objective,cached,context.scale,true);
        for(auto it=trace.rbegin();it!=trace.rend();++it)
            if((it->endpoint.eta.array()==eta.array()).all())
            {it->accepted=true; it->accepted_update=update; break;}
    }
};
}
SearchResult SearchProfile(const Domain & domain,VectorRef y,const Vector & initial_b,
    const EvaluationContext & context)
{
    const auto active_eta_max=context.eta_upper_bound;
    Profile profile{domain,y,context,{}, {},0,0,{}, {},{},{},{}};
    Vector eta=initial_b.array().log(); int accepted{};
    auto search=[&](auto & lm) {
        lm.parameters.factor=.1; lm.parameters.ftol=1e-14; lm.parameters.xtol=1e-12;
        lm.parameters.gtol=1e-12; lm.parameters.maxfev=context.profile_budget;
        auto status=lm.minimizeInit(eta);
        if(profile.cached.valid)
        {
            const bool trusted=profile.Trial(eta,Vector::Zero(eta.size()),Vector::Ones(eta.size()),0,0,0,0,0,true);
            if(!trusted) return Eigen::LevenbergMarquardtSpace::UserAsked;
        }
        if(profile.cached.valid) profile.Accept(eta,0);
        while(status==Eigen::LevenbergMarquardtSpace::NotStarted || status==Eigen::LevenbergMarquardtSpace::Running)
        {
            if(accepted>=context.update_budget) {profile.failure="accepted-update-budget"; break;}
            const auto before=lm.iterations(); status=lm.minimizeOneStep(eta);
            if(lm.iterations()>before) {++accepted; profile.Accept(eta,accepted);}
        }
        return status;
    };
    Eigen::LevenbergMarquardtSpace::Status status;
    if(!active_eta_max || !std::isfinite(*active_eta_max))
    {
        profile.failure="invalid-width-upper-bound";
        status=Eigen::LevenbergMarquardtSpace::UserAsked;
    }
    else if((eta.array()>*active_eta_max).any())
    {
        profile.failure="bounded-initial-width-infeasible";
        status=Eigen::LevenbergMarquardtSpace::UserAsked;
    }
    else
    {
        InstrumentedLM<Profile> lm(profile);
        status=search(lm);
    }
    SearchResult out; out.initial=profile.trace.empty() ? static_cast<Endpoint>(profile.cached) : profile.trace.front().endpoint;
    out.initial_accepted=!profile.trace.empty() && profile.trace.front().accepted;
    out.trials=std::move(profile.trace); out.eta=eta; out.lm_status=static_cast<int>(status);
    out.stop_reason=profile.failure.empty() ? "native-lm-stop" : profile.failure;
    out.width_failed=profile.bounded_lm_work.failed_subproblems>0 ||
        out.stop_reason=="bounded-tiny-step-nonstationary" ||
        (out.stop_reason=="unrepresentable-step" && !profile.WidthStationary());
    out.evaluations=profile.evaluations; out.derivatives=profile.derivatives; out.accepted=accepted;
    out.stopped=status==Eigen::LevenbergMarquardtSpace::UserAsked || status==Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation || !profile.failure.empty();
    out.bounded_lm_work=profile.bounded_lm_work;
    return out;
}
}
