#include "Numerics.hpp"
#include "SparseFactor.hpp"
#include "InstrumentedLM.hpp"
#include "TiledDerivative.hpp"
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
    int references{};
    double reference_seconds{};
    LinearWorkspace workspace;
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
            const auto start=std::chrono::steady_clock::now();
            auto evidence=CheckTrust(domain,y,cached,context); ++references;
            reference_seconds+=Seconds(start); trusted=evidence.passed; row.trust=std::move(evidence);
        }
        if(!trusted) failure="untrusted-trial";
        return trusted;
    }
    int values() const {return static_cast<int>(domain.rows);}
    bool Get(const Vector & eta)
    {
        if(cached.valid && cached.eta.size()==eta.size() && (cached.eta.array()==eta.array()).all()) return true;
        if(evaluations>=context.profile_budget) {failure="profile-budget"; return false;}
        failure.clear();
        const auto start=std::chrono::steady_clock::now();
        cached=EvaluateProfile(domain,y,eta,false,&context,nullptr,&workspace); ++evaluations;
        joint_component::Trial row; row.endpoint=cached; row.evaluation=evaluations; row.seconds=Seconds(start); trace.push_back(std::move(row));
        if(!cached.valid) failure="inner-"+cached.reason;
        return cached.valid;
    }
    int operator()(const Vector & eta,Vector & residual)
    {if(!Get(eta)) return -1; residual=cached.residual/scale; return 0;}
    int linearize(const Vector & eta,const Vector &,Matrix & factor,Vector & response,Vector & norms)
    {
        if(!Get(eta)) return -1;
        const auto prepared=PrepareDerivative(cached,scale,&context);
        auto differential=ReduceDerivative(prepared,cached.residual,false); ++derivatives;
        if(!differential.valid) {failure=differential.reason; return -1;}
        factor=std::move(differential.jacobian); response=std::move(differential.response);
        norms=std::move(differential.jacobian_norms); return 0;
    }
    void Accept(const Vector & eta,int update)
    {
        for(auto it=trace.rbegin();it!=trace.rend();++it)
            if((it->endpoint.eta.array()==eta.array()).all())
            {it->accepted=true; it->accepted_update=update; break;}
    }
};
}
SearchResult SearchProfile(const Domain & domain,VectorRef y,const Vector & initial_b,
    const EvaluationContext & context)
{
    const auto start=std::chrono::steady_clock::now();
    Profile profile{domain,y,context.scale,context,{}, {},0,0,{},0,0,{}};
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
    out.stop_reason=profile.failure.empty() ? "native-lm-stop" : profile.failure;
    out.evaluations=profile.evaluations; out.derivatives=profile.derivatives; out.accepted=accepted;
    out.references=profile.references; out.reference_seconds=profile.reference_seconds;
    out.stopped=status==Eigen::LevenbergMarquardtSpace::UserAsked || status==Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation || !profile.failure.empty();
    out.seconds=Seconds(start); return out;
}
}
