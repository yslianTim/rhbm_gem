#include "Numerics.hpp"
#include "OperatorSearch.hpp"
#include "SparseFactor.hpp"
#include "InstrumentedLM.hpp"
#include "TiledDerivative.hpp"
#include <algorithm>
#include <chrono>

namespace rhbm_gem::core::joint_component {
namespace {
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
struct Profile
{
    const Domain & domain;
    VectorRef y;
    double scale;
    const EvaluationContext & context;
    Evaluation cached;
    std::vector<joint_component::Trial> trace;
    int evaluations{},derivatives{};
    std::string failure;
    LinearWorkspace workspace;
    LinearWorkspace * supplied_workspace{};
    const void * workspace_identity{};
    ProfileSearchWork * telemetry{};
    ProfileEvaluationRole last_role{ProfileEvaluationRole::Unspecified};
    const JointProgressObserver & observer;
    const JointProgressComponent * progress_component{};
    std::chrono::steady_clock::time_point search_start;
    int accepted_updates{};
    std::optional<double> accepted_objective;
    std::optional<double> accepted_gradient_inf_norm;
    void Report() const
    {
        if(!progress_component) return;
        NotifyJointProgress(observer,JointProgressPhase::SearchProgress,*progress_component,
            evaluations,context.profile_budget,accepted_updates,context.update_budget,Seconds(search_start),{},false,
            accepted_objective,accepted_gradient_inf_norm);
    }
    bool retry() const {return evaluations<context.profile_budget && failure!="unrepresentable-step";}
    bool Trial(const Vector & accepted,const Vector & step,const Vector & diagonal,
        double radius,double damping,double actual,double predicted,double ratio,bool proposed)
    {
        if(trace.empty()) return false;
        auto & row=trace.back();
        row.lm=LmTrial{accepted,step,diagonal,radius,damping,cached.valid ? actual : unavailable,predicted,
            cached.valid ? ratio : unavailable,proposed};
        bool trusted=cached.valid;
        if(proposed || context.audit.trial_details)
        {
            const auto replay_started=std::chrono::steady_clock::now();
            auto evidence=CheckReplay(domain,y,cached,context);
            if(telemetry)
            {
                ProfileRoleWork work; work.replay_checks=1;
                work.replay_trust_seconds=Seconds(replay_started);
                telemetry->Add(last_role,work);
            }
            trusted=evidence.passed; row.trust=std::move(evidence);
        }
        if(!proposed || !trusted)
            UpdateAcceptedProfileObjective(accepted_objective,cached,context.scale,false);
        if(!trusted) failure="untrusted-trial";
        return trusted;
    }
    int values() const {return static_cast<int>(domain.rows);}
    bool Get(const Vector & eta)
    {
        if(cached.valid && cached.eta.size()==eta.size() && (cached.eta.array()==eta.array()).all()) return true;
        if(evaluations>=context.profile_budget) {failure="profile-budget"; return false;}
        failure.clear();
        ResourcePhase evaluation(evaluations==0 ? "profile-evaluation" : "trial-evaluation",true,domain.rows,eta.size());
        const auto start=std::chrono::steady_clock::now();
        const auto & sparse_work=SparseWorkForTesting();
        const double factor_seconds_before=sparse_work.symbolic_seconds+sparse_work.numeric_seconds+
            sparse_work.fixed_factor_seconds;
        auto * active_workspace=supplied_workspace ? supplied_workspace : &workspace;
        const void * identity=supplied_workspace ? workspace_identity : nullptr;
        last_role=evaluations==0 ? ProfileEvaluationRole::InitialProfile : ProfileEvaluationRole::TrialProfile;
        cached=EvaluateProfile(domain,y,eta,false,&context,nullptr,active_workspace,identity,last_role,telemetry); ++evaluations;
        joint_component::Trial row; row.endpoint=cached; row.evaluation=evaluations; row.seconds=Seconds(start);
        row.factor_seconds=SparseWorkForTesting().symbolic_seconds+SparseWorkForTesting().numeric_seconds+
            SparseWorkForTesting().fixed_factor_seconds-factor_seconds_before;
        trace.push_back(std::move(row));
        Report();
        if(!cached.valid) failure="inner-"+cached.reason;
        return cached.valid;
    }
    int operator()(const Vector & eta,Vector & residual)
    {if(!Get(eta)) return -1; residual=cached.residual/scale; return 0;}
    int linearize(const Vector & eta,const Vector &,Matrix & factor,Vector & response,Vector & norms)
    {
        if(!Get(eta)) return -1;
        const auto prepare_started=std::chrono::steady_clock::now();
        const auto prepared=PrepareDerivative(cached,scale,&context);
        const auto prepare_seconds=Seconds(prepare_started);
        const auto reduce_started=std::chrono::steady_clock::now();
        auto differential=ReduceDerivative(prepared,cached.residual,false); ++derivatives;
        const auto reduce_seconds=Seconds(reduce_started);
        if(telemetry)
        {
            ProfileRoleWork work; work.derivative_preparations=1; work.derivative_prepare_seconds=prepare_seconds;
            telemetry->Add(last_role,work);
            work={}; work.derivative_reductions=1; work.derivative_reduce_seconds=reduce_seconds;
            telemetry->Add(last_role,work);
        }
        if(!differential.valid) {failure=differential.reason; return -1;}
        factor=std::move(differential.jacobian); response=std::move(differential.response);
        norms=std::move(differential.jacobian_norms); return 0;
    }
    void Accept(const Vector & eta,int update)
    {
        const bool changed=update>accepted_updates;
        accepted_updates=update;
        UpdateAcceptedProfileObjective(accepted_objective,cached,context.scale,true);
        accepted_gradient_inf_norm=ProfileGradientInfinityNorm(cached);
        for(auto it=trace.rbegin();it!=trace.rend();++it)
            if((it->endpoint.eta.array()==eta.array()).all())
            {it->accepted=true; it->accepted_update=update; break;}
        if(changed || update==0) Report();
    }
};
}
SearchResult SearchProfile(const Domain & domain,VectorRef y,const Vector & initial_b,
    const EvaluationContext & context,const JointProgressObserver & observer,
    const JointProgressComponent * progress_component,LinearWorkspace * workspace,const void * workspace_identity,
    ProfileSearchWork * telemetry)
{
    if(context.search.method==SearchMethod::OperatorPcg)
        return SearchOperatorProfile(domain,y,initial_b,context,observer,progress_component);
    ResourcePhase phase("search",true,domain.rows,initial_b.size());
    const auto start=std::chrono::steady_clock::now();
    ProfileSearchWork local_work;
    Profile profile{domain,y,context.scale,context,{}, {},0,0,{}, {},workspace,workspace_identity,
        telemetry ? &local_work : nullptr,ProfileEvaluationRole::Unspecified,observer,progress_component,start,0,{}, {}};
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
            const auto before=lm.iter; status=lm.minimizeOneStep(eta);
            if(lm.iter>before) {++accepted; profile.Accept(eta,accepted);}
        }
        return status;
    };
    InstrumentedLM<Profile> lm(profile); const auto status=search(lm);
    SearchResult out; out.initial=profile.trace.empty() ? static_cast<Endpoint>(profile.cached) : profile.trace.front().endpoint;
    out.initial_accepted=!profile.trace.empty() && profile.trace.front().accepted;
    out.trials=std::move(profile.trace); out.eta=eta; out.lm_status=static_cast<int>(status);
    out.accepted_objective=profile.accepted_objective;
    out.accepted_gradient_inf_norm=profile.accepted_gradient_inf_norm;
    out.stop_reason=profile.failure.empty() ? "native-lm-stop" : profile.failure;
    out.evaluations=profile.evaluations; out.derivatives=profile.derivatives; out.accepted=accepted;
    out.stopped=status==Eigen::LevenbergMarquardtSpace::UserAsked || status==Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation || !profile.failure.empty();
    out.seconds=Seconds(start);
    if(telemetry)
    {
        local_work.total_seconds=out.seconds;
        local_work.lm_overhead_seconds=std::max(0.0,out.seconds-local_work.AttributedSeconds());
        out.profile_work=std::move(local_work);
    }
    return out;
}
}
