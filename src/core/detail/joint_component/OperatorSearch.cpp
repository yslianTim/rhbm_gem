#include "OperatorSearch.hpp"
#include <chrono>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <atomic>
#include <fstream>
#include <thread>
#include <sys/resource.h>
#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/task_info.h>
#elif defined(__linux__)
#include <unistd.h>
#endif
#endif

namespace rhbm_gem::core::joint_component {
SearchWork & SearchWorkForTesting() {static thread_local SearchWork work; return work;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
namespace {
std::size_t CurrentResidentBytesForTrial()
{
#ifdef __APPLE__
    mach_task_basic_info_data_t info{}; mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
    if(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)==KERN_SUCCESS)
        return static_cast<std::size_t>(info.resident_size);
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm"); std::size_t total{},resident{};
    if(statm>>total>>resident) return resident*static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
#endif
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)==0)
    {
#ifdef __APPLE__
        return static_cast<std::size_t>(usage.ru_maxrss);
#else
        return static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
    }
    return 0;
}
class TrialResidentSampler
{
    std::atomic<bool> stopped_{};
    std::atomic<std::size_t> peak_{};
    std::thread thread_;
    void Observe()
    {
        const auto bytes=CurrentResidentBytesForTrial(); auto previous=peak_.load(std::memory_order_relaxed);
        while(bytes>previous && !peak_.compare_exchange_weak(previous,bytes,std::memory_order_relaxed)) {}
    }
public:
    explicit TrialResidentSampler(bool enabled)
    {
        if(!enabled) return;
        Observe();
        thread_=std::thread([this] {
            while(!stopped_.load(std::memory_order_relaxed))
            {Observe(); std::this_thread::sleep_for(std::chrono::milliseconds(25));}
            Observe();
        });
    }
    std::size_t Stop()
    {
        if(thread_.joinable()) {stopped_.store(true,std::memory_order_relaxed); thread_.join();}
        return peak_.load(std::memory_order_relaxed);
    }
    ~TrialResidentSampler() {Stop();}
};
}
#endif
WidthStepResult SolvePcg(const VectorAction & action,const VectorAction & inverse,VectorRef rhs,VectorRef metric,int limit)
{
    ResourcePhase phase("pcg",true,rhs.size(),rhs.size()); auto & work=SearchWorkForTesting(); WorkTimer timer(work.pcg_seconds); ++work.pcg_solves;
    WidthStepResult out; out.step=Vector::Zero(rhs.size()); out.reason="pcg-invalid-input";
    struct IterationRecorder
    {
        SearchWork & work;
        const WidthStepResult & result;
        ~IterationRecorder() {work.pcg_iteration_counts.push_back(static_cast<std::size_t>(result.iterations));}
    } recorder{work,out};
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
SearchResult SearchOperatorProfile(const Domain & domain,VectorRef y,const Vector & initial_b,const EvaluationContext & context,
    const JointProgressObserver & observer,const JointProgressComponent * progress_component,ProfileSearchWork * telemetry)
{
    ResourcePhase phase("search",true,domain.rows,initial_b.size()); const auto start=std::chrono::steady_clock::now();
    SearchResult out; out.eta=initial_b.array().log(); out.stopped=true; out.lm_status=9;
    ProfileSearchWork local_work;
    ProfileEvaluationRole last_role{ProfileEvaluationRole::Unspecified};
    auto report=[&] {
        if(!progress_component) return;
        NotifyJointProgress(observer,JointProgressPhase::SearchProgress,*progress_component,
            out.evaluations,context.profile_budget,out.accepted,context.update_budget,
            std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(),{},false,
            out.accepted_objective,out.accepted_gradient_inf_norm);
    };
    auto finish=[&](const std::string & reason,bool stopped=true,int status=9) {
        out.stop_reason=reason; out.stopped=stopped; out.lm_status=status;
        out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(telemetry)
        {
            local_work.total_seconds=out.seconds;
            local_work.lm_overhead_seconds=std::max(0.0,out.seconds-local_work.AttributedSeconds());
            out.profile_work=std::move(local_work);
        }
        return out;
    };
    LinearWorkspace trial_workspace;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool evict_before_trial=OperatorFactorOwnershipForTesting()==OperatorFactorOwnershipKindForTesting::ReuseAcceptedEvictBeforeTrial;
    if(OperatorFactorOwnershipForTesting()==OperatorFactorOwnershipKindForTesting::ReuseAcceptedCopyOnWrite)
        trial_workspace.EnableCopyOnWrite();
#else
    trial_workspace.EnableCopyOnWrite();
#endif
    auto check_replay=[&](const Evaluation & evaluation) {
        const auto replay_started=std::chrono::steady_clock::now();
        auto evidence=CheckReplay(domain,y,evaluation,context);
        if(telemetry)
        {
            ProfileRoleWork work; work.replay_checks=1;
            work.replay_trust_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-replay_started).count();
            local_work.Add(last_role,work);
        }
        return evidence;
    };
    auto evaluate=[&](const Vector & eta
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ,SearchTrialDiagnostic * diagnostic
#endif
        ) {
        ResourcePhase evaluation(out.evaluations==0 ? "profile-evaluation" : "trial-evaluation",true,domain.rows,eta.size());
        const auto t=std::chrono::steady_clock::now();
        auto & sparse_work=SparseWorkForTesting();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto symbolic_before=sparse_work.symbolic_seconds,numeric_before=sparse_work.numeric_seconds;
        const auto factor_id_before=FactorResidencyWorkForTesting().next_factor_id;
        TrialResidentSampler rss_sampler(diagnostic!=nullptr);
#endif
        last_role=out.evaluations==0 ? ProfileEvaluationRole::InitialProfile : ProfileEvaluationRole::TrialProfile;
        Evaluation e;
        if(telemetry)
            e=EvaluateProfile(domain,y,eta,false,&context,nullptr,&trial_workspace,nullptr,last_role,&local_work);
        else
            e=EvaluateProfile(domain,y,eta,false,&context,nullptr,&trial_workspace);
        ++out.evaluations;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(diagnostic)
        {
            diagnostic->profile_evaluation=static_cast<std::size_t>(out.evaluations);
            diagnostic->factor_seconds=(sparse_work.symbolic_seconds-symbolic_before)+
                (sparse_work.numeric_seconds-numeric_before);
            diagnostic->factor_constructions=FactorResidencyWorkForTesting().next_factor_id-factor_id_before;
            diagnostic->peak_rss_bytes=rss_sampler.Stop();
            diagnostic->candidate_evaluated=true; diagnostic->candidate_valid=e.valid;
            if(!e.valid) diagnostic->rejection_reason=e.reason;
        }
#endif
        Trial trial; trial.endpoint=e; trial.evaluation=out.evaluations;
        trial.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count(); out.trials.push_back(std::move(trial));
        report(); return e;
    };
    if(context.profile_budget<=0) return finish("profile-budget");
    auto accepted=evaluate(out.eta
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ,nullptr
#endif
        ); out.initial=accepted;
    if(!accepted.valid) return finish("inner-"+accepted.reason);
    out.trials.back().trust=check_replay(accepted);
    if(!out.trials.back().trust->passed) return finish("untrusted-trial");
    out.initial_accepted=true; out.trials.back().accepted=true; out.trials.back().accepted_update=0;
    UpdateAcceptedProfileObjective(out.accepted_objective,accepted,context.scale,true);
    out.accepted_gradient_inf_norm=ProfileGradientInfinityNorm(accepted);
    report();
    std::shared_ptr<const PreconditionerPartition> partition;
    Vector metric; double radius{},mu=1e-3;
    try {
        if(context.search.preconditioner==PreconditionerKind::Schwarz)
            partition=SearchPartition(domain,context,context.search.schwarz);
        for(;;)
        {
            if(out.accepted>=context.update_budget) return finish("accepted-update-budget");
            const auto rank_backend=ResolveOperatorRankBackend(context.search.operator_rank.mode,ActiveSparseBackend());
            if(!rank_backend)
            {
                OperatorWorkForTesting().rank_status="unavailable";
                OperatorWorkForTesting().rank_reason="rank-backend-unavailable";
                return finish("rank-backend-unavailable");
            }
            ProfileJacobianOperator op(accepted,context,-1,*rank_backend); ++out.derivatives; ++SearchWorkForTesting().linearizations;
            if(!op.Valid()) return finish(op.Reason());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(OperatorFactorOwnershipForTesting()==OperatorFactorOwnershipKindForTesting::ReuseAcceptedHandoff)
                trial_workspace.HandoffForTesting();
#endif
            Vector norms;
            {
                ResourcePhase phase_metric("width-metric",true,accepted.residual.size(),accepted.eta.size(),
                    static_cast<std::size_t>(accepted.derivative.nonZeros()));
                WorkTimer timer(SearchWorkForTesting().metric_seconds);
                norms=WidthNorms(RawWidthDerivative(accepted),context.scale);
                const Vector current=WidthMetric(norms); metric=metric.size() ? metric.cwiseMax(current).eval() : current;
            }
            if(out.accepted==0) {radius=.1*metric.cwiseProduct(out.eta).stableNorm(); if(radius==0) radius=.1;}
            const Vector gradient=op.ApplyAdjoint(accepted.residual/context.scale);
            if(gradient.lpNorm<Eigen::Infinity>()<=1e-12) return finish("operator-gradient-stop",false,4);
            PreconditionerContext pc{op.Identity(),PreconditionerSpace::Width,metric,mu};
            std::unique_ptr<SchwarzModel> local;
            if(partition) local=std::make_unique<SchwarzModel>(*partition,accepted,context.scale,pc);
            const auto rebuild_accepted_factor=[&](Evaluation & candidate) {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(evict_before_trial)
                {
                    candidate.factor.reset();
                    trial_workspace.HandoffForTesting();
                    op.RebuildAcceptedFactorForTesting();
                    accepted.factor=op.FactorForTesting();
                }
#else
                (void)candidate;
#endif
            };
            bool advanced=false,rejected=false; double lower=0,upper=0;
            for(int attempt=0;attempt<std::min(20,context.search.damping_trials);++attempt)
            {
                ++SearchWorkForTesting().damping_trials; mu=std::max(1e-12,mu); pc.damping=mu;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                SearchTrialDiagnostic * diagnostic=nullptr;
                if(SearchWorkForTesting().capture_trial_telemetry)
                {
                    SearchTrialDiagnostic record; record.trial_index=SearchWorkForTesting().trial_diagnostics.size()+1;
                    record.accepted_update=out.accepted; record.damping_attempt=attempt+1; record.mu=mu; record.radius=radius;
                    SearchWorkForTesting().trial_diagnostics.push_back(std::move(record));
                    diagnostic=&SearchWorkForTesting().trial_diagnostics.back();
                }
#endif
                std::unique_ptr<SchwarzPreconditioner> schwarz;
                if(local) schwarz=std::make_unique<SchwarzPreconditioner>(*local,pc);
                const Vector diagonal=norms.array().square()+mu*metric.array().square();
                const auto inverse=[&](VectorRef r)->Vector {
                    if(schwarz) return schwarz->ApplyInverse(r,pc);
                    if(context.search.preconditioner==PreconditionerKind::Diagonal) return (r.array()/diagonal.array()).matrix();
                    return r;
                };
                const auto step=WidthStepSolver(op,gradient,pc,inverse,context.search.pcg_iterations);
                if(!step.valid)
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason=step.reason;
#endif
                    return finish(step.reason);
                }
                const double length=metric.cwiseProduct(step.step).stableNorm();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(diagnostic) diagnostic->step_length=length;
#endif
                if(length>radius)
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason="step-exceeds-radius";
#endif
                    lower=mu; mu=upper>0 ? std::sqrt(lower*upper) : mu*4; continue;
                }
                if(length<.9*radius && mu>1e-12 && !rejected)
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason="step-below-radius-target";
#endif
                    upper=mu; mu=lower>0 ? std::sqrt(lower*upper) : std::max(1e-12,mu/4); continue;
                }
                if(!(step.predicted>0))
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason="nonpositive-predicted-reduction";
#endif
                    return finish("nonpositive-predicted-reduction");
                }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(diagnostic) diagnostic->predicted_reduction=step.predicted;
#endif
                const Vector candidate_eta=out.eta+step.step;
                if((candidate_eta.array()==out.eta.array()).all())
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason="unrepresentable-step";
#endif
                    return finish("unrepresentable-step");
                }
                if(!candidate_eta.allFinite())
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason="unrepresentable-step";
#endif
                    return finish("unrepresentable-step");
                }
                if(out.evaluations>=context.profile_budget)
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->rejection_reason="profile-budget";
#endif
                    return finish("profile-budget");
                }
                const double objective=.5*(accepted.residual/context.scale).squaredNorm();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(evict_before_trial)
                {
                    op.EvictAcceptedFactorForTesting();
                    accepted.factor.reset();
                    trial_workspace.HandoffForTesting();
                }
#endif
                auto candidate=evaluate(candidate_eta
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    ,diagnostic
#endif
                    );
                const double actual=candidate.valid ? objective-.5*(candidate.residual/context.scale).squaredNorm() : unavailable;
                const double ratio=actual/step.predicted;
                const bool proposed=candidate.valid && std::isfinite(ratio) && ratio>=1e-4;
                auto & trial=out.trials.back(); trial.lm=LmTrial{out.eta,step.step,metric,radius,mu,actual,step.predicted,ratio,proposed};
                if(proposed || context.audit.trial_details) trial.trust=check_replay(candidate);
                const bool trusted=candidate.valid && (!trial.trust || trial.trust->passed);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(diagnostic)
                {
                    diagnostic->actual_reduction=actual; diagnostic->predicted_reduction=step.predicted;
                    diagnostic->ratio=ratio; diagnostic->trust_evaluated=trial.trust.has_value(); diagnostic->trusted=trusted;
                    if(!candidate.valid) diagnostic->rejection_reason=candidate.reason;
                    else if(!trusted) diagnostic->rejection_reason=trial.trust ? trial.trust->reason : "untrusted-candidate";
                    else if(!proposed) diagnostic->rejection_reason="trust-ratio-below-threshold";
                }
#endif
                if(!trusted)
                {
                    UpdateAcceptedProfileObjective(out.accepted_objective,candidate,context.scale,false);
                    rejected=true; radius*=.25; mu*=4; lower=upper=0;
                    rebuild_accepted_factor(candidate); continue;
                }
                if(ratio<=.25) {radius*=.25; mu*=4;}
                else if(ratio>=.75) {radius=std::max(radius,2*length); mu=std::max(1e-12,mu*.5);}
                if(proposed)
                {
                    accepted=std::move(candidate); out.eta=candidate_eta; ++out.accepted;
                    UpdateAcceptedProfileObjective(out.accepted_objective,accepted,context.scale,true);
                    out.accepted_gradient_inf_norm=ProfileGradientInfinityNorm(accepted);
                    trial.accepted=true; trial.accepted_update=out.accepted; advanced=true;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(diagnostic) diagnostic->accepted=true;
#endif
                    report();
                    const bool small_reduction=std::abs(actual)<=1e-14*objective && step.predicted<=1e-14*objective;
                    const bool small_step=radius<=1e-12*metric.cwiseProduct(out.eta).stableNorm();
                    if(small_reduction || small_step) return finish(small_step ? "operator-step-stop" : "operator-reduction-stop",false,small_step ? 2 : 1);
                    break;
                }
                UpdateAcceptedProfileObjective(out.accepted_objective,candidate,context.scale,false);
                rejected=true; lower=upper=0;
                if(radius<=1e-12*metric.cwiseProduct(out.eta).stableNorm()) return finish("no-trustworthy-descent-step");
                rebuild_accepted_factor(candidate);
            }
            if(!advanced) return finish("damping-trial-budget");
        }
    } catch(const std::runtime_error & e) {return finish(e.what());}
}
}
