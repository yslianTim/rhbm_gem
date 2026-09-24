#include "OperatorSearch.hpp"
#include <chrono>

namespace rhbm_gem::core::joint_component {
SearchWork & SearchWorkForTesting() {static thread_local SearchWork work; return work;}
WidthStepResult SolvePcg(const VectorAction & action,const VectorAction & inverse,VectorRef rhs,VectorRef metric,int limit)
{
    ResourcePhase phase("pcg"); auto & work=SearchWorkForTesting(); WorkTimer timer(work.pcg_seconds); ++work.pcg_solves;
    WidthStepResult out; out.step=Vector::Zero(rhs.size()); out.reason="pcg-invalid-input";
    if(rhs.size()==0 || metric.size()!=rhs.size() || !rhs.allFinite() || !metric.allFinite() || (metric.array()<=0).any()) return out;
    const auto norm=[&](VectorRef r){return (r.array()/metric.array()).matrix().stableNorm();};
    const double initial=norm(rhs); if(initial==0) {out.valid=true; out.reason="pcg-zero-rhs"; out.relative_residual=0; return out;}
    const int maximum=static_cast<int>(std::min<Eigen::Index>(4*rhs.size(),1000));
    limit=limit<0 ? maximum : std::min(limit,maximum);
    try {
        Vector r=rhs,z=inverse(r),direction=z;
        double rz=r.dot(z);
        for(int k=0;k<limit;++k)
        {
            if(!z.allFinite() || !std::isfinite(rz)) {out.reason="pcg-nonfinite"; return out;}
            if(!(rz>0)) {out.reason="pcg-nonpositive-preconditioner"; return out;}
            const Vector hd=action(direction); const double curvature=direction.dot(hd);
            if(!hd.allFinite() || !std::isfinite(curvature)) {out.reason="pcg-nonfinite"; return out;}
            if(!(curvature>0)) {out.reason="pcg-nonpositive-curvature"; return out;}
            const double alpha=rz/curvature; out.step+=alpha*direction; r-=alpha*hd;
            ++out.iterations; ++work.pcg_iterations;
            const bool refresh=out.iterations%32==0 || norm(r)<=1e-10*initial || out.iterations==limit;
            if(refresh) r=rhs-action(out.step);
            out.relative_residual=norm(r)/initial; work.last_relative_residual=out.relative_residual;
            if(!out.step.allFinite() || !r.allFinite() || !std::isfinite(out.relative_residual)) {out.reason="pcg-nonfinite"; return out;}
            if(refresh && out.relative_residual<=1e-10) {out.valid=true; out.reason="pcg-converged"; return out;}
            z=inverse(r); const double next=r.dot(z);
            // True-residual replacement restarts the recurrence; never combine
            // a replaced residual with a direction conjugate to the old one.
            if(refresh) direction=z; else direction=z+(next/rz)*direction;
            rz=next;
        }
        out.reason="pcg-iteration-budget";
    } catch(const std::logic_error &) {out.reason="pcg-stale-or-invalid-context";}
      catch(const std::runtime_error &) {out.reason="pcg-action-failed";}
    return out;
}
WidthStepResult WidthStepSolver(const ProfileJacobianOperator & op,VectorRef gradient,const PreconditionerContext & context,
    const VectorAction & inverse,int limit)
{
    if(!context.Valid() || context.linearization!=op.Identity() || context.space!=PreconditionerSpace::Width || context.metric.size()!=op.Columns() || !(context.damping>0))
    {WidthStepResult out; out.reason="pcg-stale-or-invalid-context"; return out;}
    const Vector diagonal=context.damping*context.metric.array().square().matrix();
    const auto action=[&](VectorRef v)->Vector {return op.ApplyNormal(v)+(diagonal.array()*v.array()).matrix();};
    auto result=SolvePcg(action,inverse,-gradient,context.metric,limit);
    if(result.valid)
    {
        const Vector jstep=op.Apply(result.step);
        result.predicted=-gradient.dot(result.step)-.5*jstep.squaredNorm();
        if(!std::isfinite(result.predicted)) {result.valid=false; result.reason="pcg-nonfinite";}
    }
    return result;
}
SearchResult SearchOperatorProfile(const Domain & domain,VectorRef y,const Vector & initial_b,const EvaluationContext & context)
{
    ResourcePhase phase("search"); const auto start=std::chrono::steady_clock::now();
    SearchResult out; out.eta=initial_b.array().log(); out.stopped=true; out.lm_status=9;
    auto finish=[&](const std::string & reason,bool stopped=true,int status=9) {
        out.stop_reason=reason; out.stopped=stopped; out.lm_status=status;
        out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); return out;
    };
    LinearWorkspace trial_workspace;
    auto evaluate=[&](const Vector & eta) {
        const auto t=std::chrono::steady_clock::now();
        auto e=EvaluateProfile(domain,y,eta,false,&context,nullptr,&trial_workspace); ++out.evaluations;
        Trial trial; trial.endpoint=e; trial.evaluation=out.evaluations;
        trial.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count(); out.trials.push_back(std::move(trial)); return e;
    };
    if(context.profile_budget<=0) return finish("profile-budget");
    auto accepted=evaluate(out.eta); out.initial=accepted;
    if(!accepted.valid) return finish("inner-"+accepted.reason);
    out.trials.back().trust=CheckReplay(domain,y,accepted,context);
    if(!out.trials.back().trust->passed) return finish("untrusted-trial");
    out.initial_accepted=true; out.trials.back().accepted=true; out.trials.back().accepted_update=0;
    std::shared_ptr<const PreconditionerPartition> partition;
    Vector metric; double radius{},mu=1e-3;
    try {
        if(context.search.preconditioner==PreconditionerKind::Schwarz) partition=SearchPartition(domain,context);
        for(;;)
        {
            if(out.accepted>=context.update_budget) return finish("accepted-update-budget");
            const ProfileJacobianOperator op(accepted,context); ++out.derivatives; ++SearchWorkForTesting().linearizations;
            if(!op.Valid()) return finish(op.Reason());
            Vector norms;
            {
                ResourcePhase phase_metric("width-metric"); WorkTimer timer(SearchWorkForTesting().metric_seconds);
                norms=WidthNorms(RawWidthDerivative(accepted),context.scale);
                const Vector current=WidthMetric(norms); metric=metric.size() ? metric.cwiseMax(current).eval() : current;
            }
            if(out.accepted==0) {radius=.1*metric.cwiseProduct(out.eta).stableNorm(); if(radius==0) radius=.1;}
            const Vector gradient=op.ApplyAdjoint(accepted.residual/context.scale);
            if(gradient.lpNorm<Eigen::Infinity>()<=1e-12) return finish("operator-gradient-stop",false,4);
            PreconditionerContext pc{op.Identity(),PreconditionerSpace::Width,metric,mu};
            std::unique_ptr<SchwarzModel> local;
            if(partition) local=std::make_unique<SchwarzModel>(*partition,accepted,context.scale,pc);
            bool advanced=false,rejected=false; double lower=0,upper=0;
            for(int attempt=0;attempt<std::min(20,context.search.damping_trials);++attempt)
            {
                ++SearchWorkForTesting().damping_trials; mu=std::max(1e-12,mu); pc.damping=mu;
                std::unique_ptr<SchwarzPreconditioner> schwarz;
                if(local) schwarz=std::make_unique<SchwarzPreconditioner>(*local,pc);
                const Vector diagonal=norms.array().square()+mu*metric.array().square();
                const auto inverse=[&](VectorRef r)->Vector {
                    if(schwarz) return schwarz->ApplyInverse(r,pc);
                    if(context.search.preconditioner==PreconditionerKind::Diagonal) return (r.array()/diagonal.array()).matrix();
                    return r;
                };
                const auto step=WidthStepSolver(op,gradient,pc,inverse,context.search.pcg_iterations);
                if(!step.valid) return finish(step.reason);
                const double length=metric.cwiseProduct(step.step).stableNorm();
                if(length>radius)
                {lower=mu; mu=upper>0 ? std::sqrt(lower*upper) : mu*4; continue;}
                if(length<.9*radius && mu>1e-12 && !rejected)
                {upper=mu; mu=lower>0 ? std::sqrt(lower*upper) : std::max(1e-12,mu/4); continue;}
                if(!(step.predicted>0)) return finish("nonpositive-predicted-reduction");
                const Vector candidate_eta=out.eta+step.step;
                if((candidate_eta.array()==out.eta.array()).all()) return finish("unrepresentable-step");
                if(!candidate_eta.allFinite()) return finish("unrepresentable-step");
                if(out.evaluations>=context.profile_budget) return finish("profile-budget");
                auto candidate=evaluate(candidate_eta);
                const double objective=.5*(accepted.residual/context.scale).squaredNorm();
                const double actual=candidate.valid ? objective-.5*(candidate.residual/context.scale).squaredNorm() : unavailable;
                const double ratio=actual/step.predicted;
                const bool proposed=candidate.valid && std::isfinite(ratio) && ratio>=1e-4;
                auto & trial=out.trials.back(); trial.lm=LmTrial{out.eta,step.step,metric,radius,mu,actual,step.predicted,ratio,proposed};
                if(proposed || context.audit.trial_details) trial.trust=CheckReplay(domain,y,candidate,context);
                const bool trusted=candidate.valid && (!trial.trust || trial.trust->passed);
                if(!trusted)
                {rejected=true; radius*=.25; mu*=4; lower=upper=0; continue;}
                if(ratio<=.25) {radius*=.25; mu*=4;}
                else if(ratio>=.75) {radius=std::max(radius,2*length); mu=std::max(1e-12,mu*.5);}
                if(proposed)
                {
                    accepted=std::move(candidate); out.eta=candidate_eta; ++out.accepted;
                    trial.accepted=true; trial.accepted_update=out.accepted; advanced=true;
                    const bool small_reduction=std::abs(actual)<=1e-14*objective && step.predicted<=1e-14*objective;
                    const bool small_step=radius<=1e-12*metric.cwiseProduct(out.eta).stableNorm();
                    if(small_reduction || small_step) return finish(small_step ? "operator-step-stop" : "operator-reduction-stop",false,small_step ? 2 : 1);
                    break;
                }
                rejected=true; lower=upper=0;
                if(radius<=1e-12*metric.cwiseProduct(out.eta).stableNorm()) return finish("no-trustworthy-descent-step");
            }
            if(!advanced) return finish("damping-trial-budget");
        }
    } catch(const std::runtime_error & e) {return finish(e.what());}
}
}
