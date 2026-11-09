#include "Numerics.hpp"
#include "ConstrainedAcLeastSquares.hpp"
#include "CompactSvd.hpp"
#include "BoundedLmDriver.hpp"
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
    constexpr double large_log_tau_step=2.302585092994046; // one decade in width
    constexpr double weak_direction_fraction=.95;
    if(!std::isfinite(step_inf) || step_inf<=large_log_tau_step ||
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
    Evaluation cached; // Latest evaluation, including unaccepted candidates.
    std::vector<Trial> trace;
    int evaluations{},derivatives{};
    std::string failure;
    std::optional<double> accepted_objective;
    bool accepted_width_stationary{};
    BoundedLmWork bounded_lm_work;
    ConstrainedAcWorkspace ac_workspace;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::optional<ActiveSetStep> last_bounded_subproblem{};
#endif
    void RecordBoundedSubproblem(const ActiveSetStep & step)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        last_bounded_subproblem=step;
#endif
        ++bounded_lm_work.subproblem_solves;
        bounded_lm_work.damping_solves+=step.damping_solves;
        bounded_lm_work.bracket_solves+=step.bracket_solves;
        bounded_lm_work.damping_qr_factorizations+=step.damping_qr_factorizations;
        if(!step.valid) {++bounded_lm_work.failed_subproblems; return;}
        bounded_lm_work.active_set_iterations+=step.iterations;
        bounded_lm_work.working_set_releases+=step.working_set_releases;
        for(std::size_t i=0;i<step.working_bound_mask.size();++i)
            if(step.working_bound_mask[i]) ++bounded_lm_work.active_bound_occurrences;
        if(std::isfinite(step.kkt_inf_norm))
            bounded_lm_work.maximum_subproblem_kkt=std::isfinite(bounded_lm_work.maximum_subproblem_kkt) ?
                std::max(bounded_lm_work.maximum_subproblem_kkt,step.kkt_inf_norm) : step.kkt_inf_norm;
    }
    ActiveSetStep SolveBoundedStep(const Matrix & jacobian,VectorRef response,VectorRef eta,
        VectorRef diagonal,double radius,double damping)
    {
        ActiveSetStep failed_step;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(context.force_width_subproblem_failure_for_testing)
        {failed_step.status="injected-width-subproblem-failure";
            RecordBoundedSubproblem(failed_step); return failed_step;}
#endif
        if(!context.eta_upper_bound || !std::isfinite(*context.eta_upper_bound))
        {failed_step.status="invalid-width-upper-bound";
            RecordBoundedSubproblem(failed_step); return failed_step;}
        for(Eigen::Index i=0;i<eta.size();++i)
            if(eta(i)>*context.eta_upper_bound)
            {failed_step.status="infeasible-current-width";
                RecordBoundedSubproblem(failed_step); return failed_step;}
        const Vector upper=Vector::Constant(eta.size(),*context.eta_upper_bound)-eta;
        auto result=SolveBoundedTrustRegionStep(jacobian,response,diagonal,upper,radius,damping);
        RecordBoundedSubproblem(result);
        return result;
    }
    bool WidthStationary() const
    {
        return cached.valid && cached.certificate.available && cached.certificate.feasible &&
            cached.certificate.kkt_passed && cached.certificate.width_bound_applied &&
            cached.certificate.width_bound_feasible &&
            cached.certificate.width_bound_projected_kkt<=PracticalJointConvergenceCriteria.global_eta_gradient;
    }
    // A retry may have replaced the cached candidate while eta is still accepted.
    bool IsStationary(const Vector &) const {return accepted_width_stationary;}
    bool CanRetry() const {return evaluations<context.profile_budget && failure!="unrepresentable-step";}
    bool AssessTrial([[maybe_unused]] const Vector & accepted,const Vector & step,
        [[maybe_unused]] const Vector & diagonal,[[maybe_unused]] double radius,
        [[maybe_unused]] double damping,double actual,double predicted,double ratio,bool proposed)
    {
        // A rejected evaluation can exit before creating a row.
        if(trace.empty() || trace.back().accepted) return false;
        auto & row=trace.back();
        LmTrial lm;
        lm.step=step;
        lm.actual_decrease=cached.valid ? actual : unavailable;
        lm.predicted_decrease=predicted; lm.ratio=cached.valid ? ratio : unavailable;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        lm.accepted_eta=accepted; lm.diagonal=diagonal;
        lm.radius=radius; lm.damping=damping;
        lm.proposed_acceptance=proposed;
        if(last_bounded_subproblem)
        {
            lm.bounded_subproblem=true;
            lm.subproblem_status=last_bounded_subproblem->status;
            lm.subproblem_kkt_inf_norm=last_bounded_subproblem->kkt_inf_norm;
            lm.active_set_iterations=last_bounded_subproblem->iterations;
            lm.working_set_releases=last_bounded_subproblem->working_set_releases;
            for(std::size_t i=0;i<last_bounded_subproblem->working_bound_mask.size();++i)
            {
                if(!last_bounded_subproblem->working_bound_mask[i]) continue;
                lm.active_bound_indices.push_back(i);
                lm.active_bound_multipliers.push_back(last_bounded_subproblem->multipliers(
                    static_cast<Eigen::Index>(i)));
            }
        }
        last_bounded_subproblem.reset();
#endif
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
    Eigen::Index ResidualCount() const {return domain.rows;}
    bool Get(const Vector & eta)
    {
        if(cached.valid && cached.eta.size()==eta.size() && (cached.eta.array()==eta.array()).all()) return true;
        if(evaluations>=context.profile_budget) {failure="profile-budget"; return false;}
        failure.clear();
        cached=EvaluateProfile(domain,y,eta,&context,&ac_workspace); ++evaluations;
        Trial row; row.endpoint=cached; row.evaluation=evaluations;
        trace.push_back(std::move(row));
        if(!cached.valid) failure="inner-"+cached.reason;
        return cached.valid;
    }
    bool Evaluate(const Vector & eta,double & residual_norm)
    {
        if(!context.eta_upper_bound || !std::isfinite(*context.eta_upper_bound))
        {failure="invalid-width-upper-bound"; return false;}
        if((eta.array()>*context.eta_upper_bound).any())
        {failure="bounded-step-infeasible"; return false;}
        if(!Get(eta)) return false;
        // Normalize each coefficient before taking the norm, as in the full residual model.
        residual_norm=(cached.residual/context.scale).stableNorm(); return true;
    }
    bool Linearize(const Vector & eta,Matrix & factor,Vector & response,Vector & norms)
    {
        if(!Get(eta)) return false;
        const auto prepared=PrepareDerivative(cached,context.scale,&context);
        auto differential=ReduceDerivative(prepared,cached.residual); ++derivatives;
        if(!differential.valid)
        {failure=differential.reason; return false;}
        factor=std::move(differential.jacobian); response=std::move(differential.response);
        norms=std::move(differential.jacobian_norms); return true;
    }
    void Accept(const Vector & eta,int update)
    {
        UpdateAcceptedProfileObjective(accepted_objective,cached,context.scale,true);
        accepted_width_stationary=WidthStationary();
        for(auto it=trace.rbegin();it!=trace.rend();++it)
            if((it->endpoint.eta.array()==eta.array()).all())
            {it->accepted=true; it->accepted_update=update; break;}
    }
};
}
SearchResult SearchProfile(const Domain & domain,VectorRef y,const Vector & initial_tau,
    const EvaluationContext & context)
{
    const auto eta_upper_bound=context.eta_upper_bound;
    Profile profile{domain,y,context,{}, {},0,0,{}, {},false,{},{}};
    Vector eta=initial_tau.array().log(); int accepted{};
    auto search=[&] {
        BoundedLmDriver lm(profile,{.initial_radius_factor=.1,.evaluation_budget=context.profile_budget,
            .ftol=1e-14,.xtol=1e-12,.gtol=1e-12});
        auto status=lm.Initialize(eta);
        if(profile.cached.valid)
        {
            const bool trusted=profile.AssessTrial(eta,Vector::Zero(eta.size()),Vector::Ones(eta.size()),0,0,0,0,0,true);
            if(!trusted) return BoundedLmStatus::ProblemStopped;
        }
        if(profile.cached.valid) profile.Accept(eta,0);
        while(status==BoundedLmStatus::Ready || status==BoundedLmStatus::Running)
        {
            if(accepted>=context.update_budget) {profile.failure="accepted-update-budget"; break;}
            status=lm.Advance(eta);
            if(lm.AcceptedIterations()>accepted)
            {
                accepted=static_cast<int>(lm.AcceptedIterations());
                profile.Accept(eta,accepted);
            }
        }
        return status;
    };
    BoundedLmStatus status;
    if(!eta_upper_bound || !std::isfinite(*eta_upper_bound))
    {
        profile.failure="invalid-width-upper-bound";
        status=BoundedLmStatus::ProblemStopped;
    }
    else if((eta.array()>*eta_upper_bound).any())
    {
        profile.failure="bounded-initial-width-infeasible";
        status=BoundedLmStatus::ProblemStopped;
    }
    else status=search();
    SearchResult out; out.initial=profile.trace.empty() ? static_cast<Endpoint>(profile.cached) : profile.trace.front().endpoint;
    out.initial_accepted=!profile.trace.empty() && profile.trace.front().accepted;
    out.trials=std::move(profile.trace); out.eta=eta;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    out.lm_status=status;
#endif
    out.stop_reason=profile.failure.empty() ? "native-lm-stop" : profile.failure;
    out.width_failed=profile.bounded_lm_work.failed_subproblems>0 ||
        out.stop_reason=="bounded-tiny-step-nonstationary" ||
        (out.stop_reason=="unrepresentable-step" && !profile.IsStationary(eta));
    out.evaluations=profile.evaluations; out.derivatives=profile.derivatives; out.accepted=accepted;
    out.stopped=status==BoundedLmStatus::ProblemStopped || status==BoundedLmStatus::EvaluationBudget || !profile.failure.empty();
    out.bounded_lm_work=profile.bounded_lm_work;
    return out;
}
}
