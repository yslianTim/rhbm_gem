#include "Numerics.hpp"
#include "CompactSvd.hpp"
#include "ResourceWork.hpp"
#include "InstrumentedLM.hpp"
#include "TiledDerivative.hpp"
#include "ActiveSetBoundedLM.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>

namespace rhbm_gem::core::joint_component {
WeakWidthStepAssessment AssessWeakWidthStep(const Matrix & jacobian,VectorRef step,const RankPolicy & policy)
{
    WeakWidthStepAssessment out;
    if(jacobian.rows()!=jacobian.cols() || jacobian.cols()!=step.size() ||
        jacobian.rows()==0 || !jacobian.allFinite() || !step.allFinite()) return out;
    const auto svd=EvaluateRank(jacobian,{policy,jacobian.cols(),-1,RankBoundary::StrictGreater},
        nullptr,CompactSvdVectors::Right);
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
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
double CpuSeconds(std::clock_t start)
{return static_cast<double>(std::clock()-start)/static_cast<double>(CLOCKS_PER_SEC);}
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
    ProfileSearchWork * telemetry{};
    ProfileEvaluationRole last_role{ProfileEvaluationRole::Unspecified};
    const JointProgressObserver & observer;
    const JointProgressComponent * progress_component{};
    std::chrono::steady_clock::time_point search_start;
    std::clock_t cpu_start{};
    int accepted_updates{};
    std::optional<double> accepted_objective;
    std::optional<double> accepted_gradient_inf_norm;
    const std::vector<bool> * fixed_width_mask{};
    std::optional<ActiveSetStep> last_bounded_subproblem;
    BoundedLmWork bounded_lm_work;
    bool has_bounded_step() const {return context.eta_upper_bound.has_value();}
    void RecordBoundedSubproblem(const ActiveSetStep & step)
    {
        ++bounded_lm_work.subproblem_solves;
        if(!step.valid) {++bounded_lm_work.failed_subproblems; return;}
        bounded_lm_work.active_set_iterations+=step.iterations;
        bounded_lm_work.releases+=step.releases;
        for(std::size_t i=0;i<step.active.size();++i)
            if(step.active[i] && (!fixed_width_mask || fixed_width_mask->size()!=step.active.size() ||
                !(*fixed_width_mask)[i])) ++bounded_lm_work.active_bound_occurrences;
        if(std::isfinite(step.kkt_inf_norm))
            bounded_lm_work.maximum_subproblem_kkt=std::isfinite(bounded_lm_work.maximum_subproblem_kkt) ?
                std::max(bounded_lm_work.maximum_subproblem_kkt,step.kkt_inf_norm) : step.kkt_inf_norm;
    }
    ActiveSetStep solve_bounded_step(const Matrix & jacobian,VectorRef response,VectorRef eta,
        VectorRef diagonal,double radius,double damping)
    {
        ActiveSetStep failed_step;
        if(!context.eta_upper_bound || !std::isfinite(*context.eta_upper_bound))
        {failed_step.status="invalid-width-upper-bound"; last_bounded_subproblem=failed_step;
            RecordBoundedSubproblem(failed_step); return failed_step;}
        for(Eigen::Index i=0;i<eta.size();++i)
            if(eta(i)>*context.eta_upper_bound)
            {failed_step.status="infeasible-current-width"; last_bounded_subproblem=failed_step;
                RecordBoundedSubproblem(failed_step); return failed_step;}
        const Vector upper=Vector::Constant(eta.size(),*context.eta_upper_bound)-eta;
        auto result=SolveBoundedTrustRegionStep(jacobian,response,diagonal,upper,radius,damping,
            fixed_width_mask ? *fixed_width_mask : std::vector<bool>{});
        last_bounded_subproblem=result;
        RecordBoundedSubproblem(result);
        return result;
    }
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
                if(!last_bounded_subproblem->active[i] ||
                    (fixed_width_mask && fixed_width_mask->size()==last_bounded_subproblem->active.size() &&
                        (*fixed_width_mask)[i])) continue;
                lm.active_bound_indices.push_back(i);
                lm.active_bound_multipliers.push_back(last_bounded_subproblem->multipliers(
                    static_cast<Eigen::Index>(i)));
            }
        }
#endif
        last_bounded_subproblem.reset();
        row.lm=std::move(lm);
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
            if(telemetry && evidence.reason=="independent-kkt-failed")
                ++telemetry->local_replay_kkt_rejections;
            else if(telemetry && !evidence.passed)
                ++telemetry->local_replay_disagreements;
            trusted=evidence.passed; row.trust=std::move(evidence);
        }
        if(proposed && trusted && step.lpNorm<Eigen::Infinity>()>2.302585092994046)
        {
            const auto prepared=PrepareDerivative(cached,scale,&context);
            const auto differential=ReduceDerivative(prepared,cached.residual,false);
            if(differential.valid)
            {
                Matrix jacobian=differential.jacobian;
                if(fixed_width_mask && fixed_width_mask->size()==static_cast<std::size_t>(jacobian.cols()))
                    for(Eigen::Index column=0;column<jacobian.cols();++column)
                        if((*fixed_width_mask)[static_cast<std::size_t>(column)]) jacobian.col(column).setZero();
                const auto assessment=AssessWeakWidthStep(jacobian,step,context.rank);
                if(assessment.available)
                {
                    row.weak_width_step_fraction=assessment.weak_step_fraction;
                    const double candidate_objective=cached.certificate.objective/(scale*scale);
                    const double old_objective=accepted_objective.value_or(candidate_objective);
                    if(RejectWeakWidthTrial(step.lpNorm<Eigen::Infinity>(),assessment.weak_step_fraction,
                        old_objective,candidate_objective))
                    {
                        if(telemetry) ++telemetry->width_safeguard_rejections;
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
        if(context.profile_cpu_budget_seconds>0 &&
            CpuSeconds(cpu_start)>=context.profile_cpu_budget_seconds)
        {failure="cpu-time-budget"; return false;}
        if(cached.valid && cached.eta.size()==eta.size() && (cached.eta.array()==eta.array()).all()) return true;
        if(evaluations>=context.profile_budget) {failure="profile-budget"; return false;}
        failure.clear();
        ResourcePhase evaluation(evaluations==0 ? "profile-evaluation" : "trial-evaluation",true,domain.rows,eta.size());
        const auto start=std::chrono::steady_clock::now();
        last_role=evaluations==0 ? ProfileEvaluationRole::InitialProfile : ProfileEvaluationRole::TrialProfile;
        cached=EvaluateProfile(domain,y,eta,false,&context,nullptr,last_role,telemetry); ++evaluations;
        joint_component::Trial row; row.endpoint=cached; row.evaluation=evaluations; row.seconds=Seconds(start);
        trace.push_back(std::move(row));
        Report();
        if(!cached.valid) failure="inner-"+cached.reason;
        return cached.valid;
    }
    int operator()(const Vector & eta,Vector & residual)
    {
        if(context.eta_upper_bound && (!std::isfinite(*context.eta_upper_bound) ||
            (eta.array()>*context.eta_upper_bound).any()))
        {failure="bounded-step-infeasible"; return -1;}
        if(!Get(eta)) return -1;
        residual=cached.residual/scale; return 0;
    }
    int linearize(const Vector & eta,const Vector &,Matrix & factor,Vector & response,Vector & norms)
    {
        if(context.profile_cpu_budget_seconds>0 &&
            CpuSeconds(cpu_start)>=context.profile_cpu_budget_seconds)
        {failure="cpu-time-budget"; return -1;}
        if(!Get(eta)) return -1;
        const auto prepare_started=std::chrono::steady_clock::now();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto derivative_before=DerivativeWorkForTesting();
#endif
        const auto prepared=PrepareDerivative(cached,scale,&context);
        const auto prepare_seconds=Seconds(prepare_started);
        const auto reduce_started=std::chrono::steady_clock::now();
        auto differential=ReduceDerivative(prepared,cached.residual,false); ++derivatives;
        const auto reduce_seconds=Seconds(reduce_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto derivative_after=DerivativeWorkForTesting();
#endif
        if(telemetry)
        {
            ProfileRoleWork work; work.derivative_preparations=1; work.derivative_prepare_seconds=prepare_seconds;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            work.derivative_raw_assembly_seconds=derivative_after.preparation.raw_assembly_seconds-
                derivative_before.preparation.raw_assembly_seconds;
            work.derivative_free_design_assembly_seconds=derivative_after.preparation.free_design_assembly_seconds-
                derivative_before.preparation.free_design_assembly_seconds;
            work.derivative_factor_match_seconds=derivative_after.preparation.factor_match_seconds-
                derivative_before.preparation.factor_match_seconds;
            work.derivative_factor_build_seconds=derivative_after.preparation.factor_build_seconds-
                derivative_before.preparation.factor_build_seconds;
            work.derivative_factor_compact_seconds=derivative_after.preparation.factor_compact_seconds-
                derivative_before.preparation.factor_compact_seconds;
            work.derivative_rank_seconds=derivative_after.preparation.rank_seconds-
                derivative_before.preparation.rank_seconds;
            work.derivative_least_squares_seconds=derivative_after.preparation.least_squares_seconds-
                derivative_before.preparation.least_squares_seconds;
            work.derivative_normal_solve_seconds=derivative_after.preparation.normal_solve_seconds-
                derivative_before.preparation.normal_solve_seconds;
            work.derivative_cancellation_check_seconds=derivative_after.preparation.cancellation_check_seconds-
                derivative_before.preparation.cancellation_check_seconds;
            work.derivative_cancellation_fallback_seconds=derivative_after.preparation.cancellation_fallback_seconds-
                derivative_before.preparation.cancellation_fallback_seconds;
#endif
            telemetry->Add(last_role,work);
            work={}; work.derivative_reductions=1; work.derivative_reduce_seconds=reduce_seconds;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            work.derivative_rows_seconds=derivative_after.rows_seconds-derivative_before.rows_seconds;
            work.derivative_jacobian_qr_seconds=derivative_after.jacobian_qr_seconds-
                derivative_before.jacobian_qr_seconds;
            work.derivative_norms_seconds=derivative_after.norms_seconds-derivative_before.norms_seconds;
            work.derivative_outer_overhead_seconds=derivative_after.outer_overhead_seconds-
                derivative_before.outer_overhead_seconds;
            work.tiled_qr_assembly_copy_seconds=derivative_after.jacobian_qr.assembly_copy_seconds-
                derivative_before.jacobian_qr.assembly_copy_seconds;
            work.tiled_qr_householder_seconds=derivative_after.jacobian_qr.householder_seconds-
                derivative_before.jacobian_qr.householder_seconds;
            work.tiled_qr_rhs_transform_seconds=derivative_after.jacobian_qr.rhs_transform_seconds-
                derivative_before.jacobian_qr.rhs_transform_seconds;
#endif
            telemetry->Add(last_role,work);
        }
        if(!differential.valid)
        {if(telemetry) ++telemetry->invalid_gradient_stops; failure=differential.reason; return -1;}
        if(fixed_width_mask && fixed_width_mask->size()==static_cast<std::size_t>(differential.jacobian.cols()))
            for(Eigen::Index column=0;column<differential.jacobian.cols();++column)
                if((*fixed_width_mask)[static_cast<std::size_t>(column)])
                {
                    differential.jacobian.col(column).setZero();
                    differential.jacobian_norms(column)=0.0;
                }
        factor=std::move(differential.jacobian); response=std::move(differential.response);
        norms=std::move(differential.jacobian_norms); return 0;
    }
    void Accept(const Vector & eta,int update)
    {
        const bool changed=update>accepted_updates;
        accepted_updates=update;
        UpdateAcceptedProfileObjective(accepted_objective,cached,context.scale,true);
        accepted_gradient_inf_norm=ProfileGradientInfinityNorm(cached,context);
        for(auto it=trace.rbegin();it!=trace.rend();++it)
            if((it->endpoint.eta.array()==eta.array()).all())
            {it->accepted=true; it->accepted_update=update; break;}
        if(changed || update==0) Report();
    }
};
}
namespace {
SearchResult SearchProfileImpl(const Domain & domain,VectorRef y,const Vector & initial_b,
    const EvaluationContext & context,const JointProgressObserver & observer,
    const JointProgressComponent * progress_component,ProfileSearchWork * telemetry,
    const std::vector<bool> * fixed_width_mask)
{
    ResourcePhase phase("search",true,domain.rows,initial_b.size());
    const auto start=std::chrono::steady_clock::now();
    ProfileSearchWork local_work;
    EvaluationContext search_context=context;
    const auto active_eta_max=search_context.eta_upper_bound;
    Profile profile{domain,y,search_context.scale,search_context,{}, {},0,0,{},telemetry ? &local_work : nullptr,
        ProfileEvaluationRole::Unspecified,observer,progress_component,start,std::clock(),0,{}, {},nullptr,{},{}
    };
    profile.fixed_width_mask=fixed_width_mask;
    Vector eta=initial_b.array().log(); int accepted{};
    auto search=[&](auto & lm) {
        lm.parameters.factor=.1; lm.parameters.ftol=1e-14; lm.parameters.xtol=1e-12;
        lm.parameters.gtol=1e-12; lm.parameters.maxfev=search_context.profile_budget;
        auto status=lm.minimizeInit(eta);
        if(profile.cached.valid)
        {
            const bool trusted=profile.Trial(eta,Vector::Zero(eta.size()),Vector::Ones(eta.size()),0,0,0,0,0,true);
            if(!trusted) return Eigen::LevenbergMarquardtSpace::UserAsked;
        }
        if(profile.cached.valid) profile.Accept(eta,0);
        while(status==Eigen::LevenbergMarquardtSpace::NotStarted || status==Eigen::LevenbergMarquardtSpace::Running)
        {
            if(accepted>=search_context.update_budget) {profile.failure="accepted-update-budget"; break;}
            const auto before=lm.iter; status=lm.minimizeOneStep(eta);
            if(lm.iter>before) {++accepted; profile.Accept(eta,accepted);}
        }
        return status;
    };
    Eigen::LevenbergMarquardtSpace::Status status;
    if(active_eta_max && (!std::isfinite(*active_eta_max) || (eta.array()>*active_eta_max).any()))
    {
        profile.failure=std::isfinite(*active_eta_max) ?
            "bounded-initial-width-infeasible" : "invalid-width-upper-bound";
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
    out.accepted_objective=profile.accepted_objective;
    out.accepted_gradient_inf_norm=profile.accepted_gradient_inf_norm;
    out.stop_reason=profile.failure.empty() ? "native-lm-stop" : profile.failure;
    out.evaluations=profile.evaluations; out.derivatives=profile.derivatives; out.accepted=accepted;
    out.stopped=status==Eigen::LevenbergMarquardtSpace::UserAsked || status==Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation || !profile.failure.empty();
    out.seconds=Seconds(start);
    if(telemetry)
    {
        if(out.lm_status>=0 && static_cast<std::size_t>(out.lm_status)<local_work.lm_status_counts.size())
            ++local_work.lm_status_counts[static_cast<std::size_t>(out.lm_status)];
        for(const auto & trial:out.trials)
        {
            if(trial.accepted) ++local_work.accepted_profile_trials;
            else ++local_work.rejected_profile_trials;
        }
        if(out.stop_reason=="accepted-update-budget") ++local_work.accepted_update_budget_stops;
        else if(out.stop_reason=="profile-budget") ++local_work.profile_budget_stops;
        else if(out.stop_reason=="cpu-time-budget") ++local_work.cpu_budget_stops;
        else if(out.stop_reason=="native-lm-stop") ++local_work.native_lm_stops;
        local_work.bounded_lm=profile.bounded_lm_work;
        local_work.total_seconds=out.seconds;
        local_work.lm_overhead_seconds=std::max(0.0,out.seconds-local_work.AttributedSeconds());
        out.profile_work=std::move(local_work);
    }
    else out.profile_work.bounded_lm=profile.bounded_lm_work;
    return out;
}
}

SearchResult SearchProfile(const Domain & domain,VectorRef y,const Vector & initial_b,
    const EvaluationContext & context,const JointProgressObserver & observer,
    const JointProgressComponent * progress_component,ProfileSearchWork * telemetry,
    const std::vector<bool> * fixed_width_mask)
{
    return SearchProfileImpl(domain,y,initial_b,context,observer,
        progress_component,telemetry,fixed_width_mask);
}
}
