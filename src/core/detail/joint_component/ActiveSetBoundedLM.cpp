#include "ActiveSetBoundedLM.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
BoundedLmTestWork & BoundedLmWorkForTesting() {static thread_local BoundedLmTestWork work; return work;}
#endif
namespace {
struct BoundedSubproblemContext
{
    double bound_tolerance{};
    std::size_t maximum_active_iterations{};
    double jacobian_norm{},diagonal_squared_norm{};
    bool kkt_norms_ready{};
};

struct BoundedLmWorkspace
{
    Eigen::MatrixXd augmented;
    Eigen::VectorXd rhs,adjusted,candidate,step,residual,gradient,multipliers,solution;
    std::vector<Eigen::Index> free_indices,hit_indices;
    std::vector<bool> active;
};

void ResizeWorkspaceVector(Eigen::VectorXd & value,Eigen::Index size)
{
    if(value.size()!=size) value.resize(size);
}

void ResizeWorkspaceMatrix(Eigen::MatrixXd & value,Eigen::Index rows,Eigen::Index columns)
{
    if(value.rows()!=rows || value.cols()!=columns) value.resize(rows,columns);
}

void PrepareWorkspace(BoundedLmWorkspace & workspace,Eigen::Index rows,Eigen::Index columns)
{
    if(workspace.active.size()!=static_cast<std::size_t>(columns))
        workspace.active.resize(static_cast<std::size_t>(columns));
    for(auto * value:{&workspace.candidate,&workspace.step,&workspace.gradient,&workspace.multipliers})
        ResizeWorkspaceVector(*value,columns);
    ResizeWorkspaceVector(workspace.adjusted,rows);
    ResizeWorkspaceVector(workspace.residual,rows);
    if(workspace.free_indices.capacity()<static_cast<std::size_t>(columns))
        workspace.free_indices.reserve(static_cast<std::size_t>(columns));
    if(workspace.hit_indices.capacity()<static_cast<std::size_t>(columns))
        workspace.hit_indices.reserve(static_cast<std::size_t>(columns));
}

void CopyActiveSetStep(ActiveSetStep & destination,const ActiveSetStep & source)
{
    destination=source;
}

bool ValidSubproblemInput(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,
    const Eigen::VectorXd & d,const Eigen::VectorXd & upper,double lambda)
{
    const Eigen::Index n=j.cols();
    return j.rows()==r.size() && n>0 && d.size()==n && upper.size()==n &&
        j.allFinite() && r.allFinite() && d.allFinite() && upper.allFinite() &&
        (d.array()>0.0).all() && !(upper.array()<0.0).any() &&
        std::isfinite(lambda) && lambda>=0.0;
}

BoundedSubproblemContext MakeBoundedSubproblemContext(const Eigen::MatrixXd & j,
    const Eigen::VectorXd & upper)
{
    BoundedSubproblemContext context;
    context.bound_tolerance=64.0*std::numeric_limits<double>::epsilon()*
        std::max(1.0,upper.cwiseAbs().maxCoeff());
    context.maximum_active_iterations=std::max<std::size_t>(40,
        20*static_cast<std::size_t>(j.cols()));
    return context;
}

double KktScale(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,const Eigen::VectorXd & d,
    const Eigen::VectorXd & p,double lambda,BoundedSubproblemContext & context)
{
    const double residual=(j*p+r).norm();
    if(!context.kkt_norms_ready)
    {
        context.jacobian_norm=j.norm();
        context.diagonal_squared_norm=d.squaredNorm();
        context.kkt_norms_ready=true;
    }
    return std::max(1.0,context.jacobian_norm*(residual+context.jacobian_norm*p.norm())+
        lambda*context.diagonal_squared_norm*std::max(1.0,p.norm()));
}

ActiveSetStep Failure(std::string status)
{
    ActiveSetStep result; result.status=std::move(status); return result;
}
void SolveAtDamping(BoundedSubproblemContext & context,BoundedLmWorkspace & workspace,
    const Eigen::MatrixXd & j,const Eigen::VectorXd & r,const Eigen::VectorXd & d,
    const Eigen::VectorXd & upper,double lambda,ActiveSetStep & result
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ,BoundedLmDampingTrace * trace
#endif
)
{
    const Eigen::Index n=j.cols();
    result.valid=false; result.status.clear(); result.damping=lambda;
    result.scaled_step_norm=0.0; result.kkt_inf_norm=0.0;
    result.iterations=0; result.releases=0;
    result.damping_solves=0; result.bracket_solves=0; result.bisection_solves=0;
    result.damping_qr_factorizations=0; result.damping_active_set_iterations=0; result.damping_releases=0;
    const auto fail=[&](const char * status) {result.status=status;};
    if(!std::isfinite(lambda) || lambda<0.0)
    {fail("invalid-subproblem-input"); return;}

    PrepareWorkspace(workspace,j.rows(),n);
    const double bound_tolerance=context.bound_tolerance;
    std::fill(workspace.active.begin(),workspace.active.end(),false);
    workspace.step.setZero();
    for(Eigen::Index i=0;i<n;++i)
    {
        workspace.active[static_cast<std::size_t>(i)]=upper(i)<=bound_tolerance;
        if(workspace.active[static_cast<std::size_t>(i)]) workspace.step(i)=upper(i);
    }

    std::size_t releases{};
    for(std::size_t iteration=0;iteration<context.maximum_active_iterations;++iteration)
    {
        ++result.damping_active_set_iterations;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(trace)
        {
            trace->active_set_sequence.push_back(workspace.active);
            trace->step_sequence.push_back(workspace.step);
            trace->iterations=iteration+1;
        }
#endif
        workspace.free_indices.clear();
        for(Eigen::Index i=0;i<n;++i)
            if(!workspace.active[static_cast<std::size_t>(i)])
                workspace.free_indices.push_back(i);

        workspace.candidate=workspace.step;
        if(!workspace.free_indices.empty())
        {
            const Eigen::Index free_count=static_cast<Eigen::Index>(workspace.free_indices.size());
            const Eigen::Index augmented_rows=j.rows()+(lambda>0.0 ? free_count : 0);
            ResizeWorkspaceMatrix(workspace.augmented,augmented_rows,free_count);
            ResizeWorkspaceVector(workspace.rhs,augmented_rows);
            ResizeWorkspaceVector(workspace.solution,free_count);
            workspace.augmented.setZero();
            workspace.rhs.setZero();
            workspace.adjusted=r;
            for(Eigen::Index i=0;i<n;++i)
                if(workspace.active[static_cast<std::size_t>(i)])
                    workspace.adjusted.noalias()+=j.col(i)*workspace.step(i);
            workspace.rhs.head(j.rows())=-workspace.adjusted;
            for(Eigen::Index k=0;k<free_count;++k)
                workspace.augmented.col(k).head(j.rows())=
                    j.col(workspace.free_indices[static_cast<std::size_t>(k)]);
            if(lambda>0.0)
            {
                const double root=std::sqrt(lambda);
                for(Eigen::Index k=0;k<free_count;++k)
                {
                    const auto index=workspace.free_indices[static_cast<std::size_t>(k)];
                    workspace.augmented(j.rows()+k,k)=root*d(index);
                }
            }
            ++result.damping_qr_factorizations;
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(workspace.augmented);
            qr.setThreshold(64.0*std::numeric_limits<double>::epsilon()*
                static_cast<double>(std::max(workspace.augmented.rows(),workspace.augmented.cols())));
            const auto rank=qr.rank();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(trace) trace->free_ranks.push_back(rank);
#endif
            if(rank!=free_count) {fail("rank-deficient-free-subproblem"); return;}
            workspace.solution=qr.solve(workspace.rhs);
            if(!workspace.solution.allFinite()) {fail("non-finite-free-solution"); return;}
            for(Eigen::Index k=0;k<free_count;++k)
                workspace.candidate(workspace.free_indices[static_cast<std::size_t>(k)])=
                    workspace.solution(k);
        }

        if(!workspace.candidate.allFinite()) {fail("non-finite-candidate-step"); return;}
        double first_hit=1.0;
        workspace.hit_indices.clear();
        for(const auto index:workspace.free_indices)
        {
            if(workspace.candidate(index)<=upper(index)) continue;
            const double direction=workspace.candidate(index)-workspace.step(index);
            if(!(direction>0.0)) {fail("inconsistent-active-set-direction"); return;}
            const double alpha=std::clamp((upper(index)-workspace.step(index))/direction,0.0,1.0);
            if(alpha<first_hit-bound_tolerance)
            {
                first_hit=alpha;
                workspace.hit_indices.assign(1,index);
            }
            else if(std::abs(alpha-first_hit)<=bound_tolerance)
                workspace.hit_indices.push_back(index);
        }
        if(!workspace.hit_indices.empty())
        {
            workspace.step.noalias()+=first_hit*(workspace.candidate-workspace.step);
            for(const auto index:workspace.hit_indices)
            {
                workspace.step(index)=upper(index);
                workspace.active[static_cast<std::size_t>(index)]=true;
            }
            continue;
        }
        workspace.step.swap(workspace.candidate);
        workspace.residual=j*workspace.step+r;
        workspace.gradient=j.transpose()*workspace.residual+
            lambda*d.cwiseProduct(d).cwiseProduct(workspace.step);
        if(!workspace.gradient.allFinite()) {fail("non-finite-subproblem-gradient"); return;}
        const double kkt_tolerance=2e-11*KktScale(j,r,d,workspace.step,lambda,context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(trace)
        {
            trace->gradient_sequence.push_back(workspace.gradient);
            trace->kkt_tolerance_sequence.push_back(kkt_tolerance);
        }
#endif
        Eigen::Index release=-1;
        double largest_wrong_multiplier=kkt_tolerance;
        for(Eigen::Index i=0;i<n;++i)
        {
            if(!workspace.active[static_cast<std::size_t>(i)]) continue;
            if(workspace.gradient(i)>largest_wrong_multiplier)
            {
                largest_wrong_multiplier=workspace.gradient(i);
                release=i;
            }
        }
        if(release>=0)
        {
            workspace.active[static_cast<std::size_t>(release)]=false;
            ++releases;
            ++result.damping_releases;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(trace) ++trace->releases;
#endif
            continue;
        }

        double kkt_inf{};
        workspace.multipliers.setZero();
        for(Eigen::Index i=0;i<n;++i)
        {
            if(workspace.active[static_cast<std::size_t>(i)])
                workspace.multipliers(i)=std::max(0.0,-workspace.gradient(i));
            else
                kkt_inf=std::max(kkt_inf,std::abs(workspace.gradient(i)));
            if(workspace.active[static_cast<std::size_t>(i)])
                kkt_inf=std::max(kkt_inf,std::max(0.0,workspace.gradient(i)));
        }
        if(kkt_inf>kkt_tolerance) {fail("subproblem-kkt-failed"); return;}
        result.valid=true; result.status="solved"; result.step=workspace.step;
        result.gradient=workspace.gradient; result.multipliers=workspace.multipliers;
        result.active=workspace.active; result.damping=lambda;
        result.scaled_step_norm=d.cwiseProduct(result.step).norm();
        result.kkt_inf_norm=kkt_inf; result.iterations=iteration+1; result.releases=releases;
        return;
    }
    fail("active-set-iteration-budget");
}

ActiveSetStep WithDamping(ActiveSetStep result,double damping)
{
    if(result.valid) result.damping=damping;
    return result;
}
}

namespace {
ActiveSetStep SolveBoundedTrustRegionStepImpl(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,
    const Eigen::VectorXd & d,const Eigen::VectorXd & upper,double trust_radius,double damping_seed,
    BoundedLmWorkspace & workspace)
{
    BoundedSubproblemContext context;
    ActiveSetStep candidate;
    std::size_t damping_solves{},bracket_solves{},bisection_solves{},damping_qr_factorizations{},
        damping_active_set_iterations{},damping_releases{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=BoundedLmWorkForTesting();
    const auto mark_search_accepted=[&]() {
        if(work.capture_trace && !work.damping_trace.empty())
            work.damping_trace.back().search_accepted=true;
    };
#endif
    const auto finish=[&](ActiveSetStep result) {
        result.damping_solves=damping_solves;
        result.bracket_solves=bracket_solves;
        result.bisection_solves=bisection_solves;
        result.damping_qr_factorizations=damping_qr_factorizations;
        result.damping_active_set_iterations=damping_active_set_iterations;
        result.damping_releases=damping_releases;
        return result;
    };
    const auto solve_candidate=[&](double lambda,int phase,ActiveSetStep & result) -> void {
        ++damping_solves;
        if(phase==1) ++bracket_solves;
        else if(phase==2) ++bisection_solves;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        BoundedLmDampingTrace trace;
        SolveAtDamping(context,workspace,j,r,d,upper,lambda,result,
            work.capture_trace ? &trace : nullptr);
        if(work.capture_trace)
        {
            trace.damping=lambda; trace.phase=phase; trace.status=result.status;
            if(result.valid)
            {
                trace.step=result.step; trace.gradient=result.gradient; trace.multipliers=result.multipliers;
                trace.active=result.active;
            }
            if(result.iterations) trace.iterations=result.iterations;
            if(result.valid) trace.releases=result.releases;
            trace.scaled_step_norm=result.scaled_step_norm; trace.kkt_inf_norm=result.kkt_inf_norm;
            trace.valid=result.valid;
            trace.within_trust_region=result.valid && result.scaled_step_norm<=trust_radius;
            work.damping_trace.push_back(std::move(trace));
        }
#else
        SolveAtDamping(context,workspace,j,r,d,upper,lambda,result);
#endif
        damping_qr_factorizations+=result.damping_qr_factorizations;
        damping_active_set_iterations+=result.damping_active_set_iterations;
        damping_releases+=result.damping_releases;
    };
    if(!std::isfinite(trust_radius) || !(trust_radius>0.0) || !std::isfinite(damping_seed) || damping_seed<0.0)
        return Failure("invalid-trust-region-input");
    if(!ValidSubproblemInput(j,r,d,upper,0.0)) return Failure("invalid-subproblem-input");
    context=MakeBoundedSubproblemContext(j,upper);
    solve_candidate(0.0,0,candidate);
    const double trust_tolerance=64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,trust_radius);
    constexpr double boundary_relative_tolerance=1e-13;
    if(candidate.valid && candidate.scaled_step_norm<=trust_radius+trust_tolerance)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        mark_search_accepted();
#endif
        return finish(WithDamping(std::move(candidate),0.0));
    }

    double low=0.0;
    double high=std::max(damping_seed,1e-12);
    ActiveSetStep feasible;
    for(int bracket=0;bracket<80;++bracket)
    {
        solve_candidate(high,1,candidate);
        if(candidate.valid && candidate.scaled_step_norm<=trust_radius)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            mark_search_accepted();
#endif
            CopyActiveSetStep(feasible,candidate); break;
        }
        if(!candidate.valid && candidate.status!="rank-deficient-free-subproblem")
            return finish(Failure(candidate.status));
        low=high;
        high*=4.0;
        if(!std::isfinite(high)) return finish(Failure("trust-region-damping-overflow"));
    }
    if(!feasible.valid) return finish(Failure("trust-region-damping-bracket-failed"));

    for(int iteration=0;iteration<64;++iteration)
    {
        const double middle=low+(high-low)*0.5;
        if(!(middle>low) || !(middle<high)) break;
        solve_candidate(middle,2,candidate);
        if(candidate.valid && candidate.scaled_step_norm<=trust_radius)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            mark_search_accepted();
#endif
            high=middle; CopyActiveSetStep(feasible,candidate);
            if(trust_radius-candidate.scaled_step_norm<=boundary_relative_tolerance*trust_radius)
                break;
        }
        else if(candidate.valid || candidate.status=="rank-deficient-free-subproblem") low=middle;
        else return finish(Failure(candidate.status));
    }
    feasible.damping=high;
    if(feasible.scaled_step_norm>trust_radius+trust_tolerance)
        return finish(Failure("trust-region-feasibility-failed"));
    return finish(std::move(feasible));
}
}

ActiveSetStep SolveBoundedTrustRegionStep(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,
    const Eigen::VectorXd & d,const Eigen::VectorXd & upper,double trust_radius,double damping_seed)
{
    BoundedLmWorkspace workspace;
    return SolveBoundedTrustRegionStepImpl(j,r,d,upper,trust_radius,damping_seed,workspace);
}

} // namespace rhbm_gem::core::joint_component
