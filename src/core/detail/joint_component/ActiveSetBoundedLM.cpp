#include "ActiveSetBoundedLM.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <chrono>
#endif

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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ++BoundedLmWorkForTesting().bound_tolerance_evaluations;
#endif
    context.bound_tolerance=64.0*std::numeric_limits<double>::epsilon()*
        std::max(1.0,upper.cwiseAbs().maxCoeff());
    context.maximum_active_iterations=std::max<std::size_t>(40,
        20*static_cast<std::size_t>(j.cols()));
    return context;
}

double KktScale(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,const Eigen::VectorXd & d,
    const Eigen::VectorXd & p,double lambda,BoundedSubproblemContext & context)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=BoundedLmWorkForTesting();
    ++work.kkt_scale_evaluations;
#endif
    const double residual=(j*p+r).norm();
    if(!context.kkt_norms_ready)
    {
        context.jacobian_norm=j.norm();
        context.diagonal_squared_norm=d.squaredNorm();
        context.kkt_norms_ready=true;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.jacobian_norm_evaluations;
        ++work.diagonal_squared_norm_evaluations;
#endif
    }
    return std::max(1.0,context.jacobian_norm*(residual+context.jacobian_norm*p.norm())+
        lambda*context.diagonal_squared_norm*std::max(1.0,p.norm()));
}

ActiveSetStep Failure(std::string status)
{
    ActiveSetStep result; result.status=std::move(status); return result;
}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct BoundedLmTimer
{
    BoundedLmTestWork & work;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
    ~BoundedLmTimer()
    {work.damping_search_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
};
#endif

ActiveSetStep SolveAtDamping(BoundedSubproblemContext & context,const Eigen::MatrixXd & j,
    const Eigen::VectorXd & r,const Eigen::VectorXd & d,const Eigen::VectorXd & upper,double lambda
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ,BoundedLmDampingTrace * trace
#endif
)
{
    const Eigen::Index n=j.cols();
    if(!std::isfinite(lambda) || lambda<0.0)
        return Failure("invalid-subproblem-input");

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=BoundedLmWorkForTesting();
    ++work.workspace_allocation_sites; // active mask
#endif
    const double bound_tolerance=context.bound_tolerance;
    std::vector<bool> active(static_cast<std::size_t>(n),false);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ++work.workspace_allocation_sites; // step
#endif
    Eigen::VectorXd step=Eigen::VectorXd::Zero(n);
    for(Eigen::Index i=0;i<n;++i)
    {
        active[static_cast<std::size_t>(i)]=upper(i)<=bound_tolerance;
        if(active[static_cast<std::size_t>(i)]) step(i)=upper(i);
    }

    std::size_t releases{};
    for(std::size_t iteration=0;iteration<context.maximum_active_iterations;++iteration)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(trace)
        {
            trace->active_set_sequence.push_back(active);
            trace->step_sequence.push_back(step);
            trace->iterations=iteration+1;
        }
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.active_set_iterations;
#endif
        std::vector<Eigen::Index> free_indices;
        for(Eigen::Index i=0;i<n;++i)
            if(!active[static_cast<std::size_t>(i)])
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(free_indices.size()==free_indices.capacity()) ++work.workspace_allocation_sites;
#endif
                free_indices.push_back(i);
            }

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.workspace_allocation_sites; // candidate
#endif
        Eigen::VectorXd candidate=step;
        if(!free_indices.empty())
        {
            const Eigen::Index free_count=static_cast<Eigen::Index>(free_indices.size());
            Eigen::MatrixXd augmented(j.rows()+(lambda>0.0 ? free_count : 0),free_count);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(augmented.size()) ++work.workspace_allocation_sites;
#endif
            Eigen::VectorXd rhs(augmented.rows());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(rhs.size()) ++work.workspace_allocation_sites;
#endif
            augmented.setZero();
            rhs.setZero();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(r.size()) ++work.workspace_allocation_sites;
#endif
            Eigen::VectorXd adjusted=r;
            for(Eigen::Index i=0;i<n;++i)
                if(active[static_cast<std::size_t>(i)]) adjusted.noalias()+=j.col(i)*step(i);
            rhs.head(j.rows())=-adjusted;
            for(Eigen::Index k=0;k<free_count;++k)
                augmented.col(k).head(j.rows())=j.col(free_indices[static_cast<std::size_t>(k)]);
            if(lambda>0.0)
            {
                const double root=std::sqrt(lambda);
                for(Eigen::Index k=0;k<free_count;++k)
                {
                    const auto index=free_indices[static_cast<std::size_t>(k)];
                    augmented(j.rows()+k,k)=root*d(index);
                }
            }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto bytes=static_cast<std::size_t>(n)*(2*sizeof(double))+
                ((static_cast<std::size_t>(n)+8*sizeof(std::size_t)-1)/(8*sizeof(std::size_t)))*sizeof(std::size_t)+
                static_cast<std::size_t>(n)*2*sizeof(Eigen::Index)+
                static_cast<std::size_t>(augmented.size()+rhs.size()+adjusted.size()+free_count)*sizeof(double);
            work.workspace_peak_bytes=std::max(work.workspace_peak_bytes,bytes);
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto qr_started=std::chrono::steady_clock::now();
#endif
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(augmented);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++work.dense_qr_factorizations;
            work.dense_qr_factorization_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-qr_started).count();
#endif
            qr.setThreshold(64.0*std::numeric_limits<double>::epsilon()*
                static_cast<double>(std::max(augmented.rows(),augmented.cols())));
            const auto rank=qr.rank();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(trace) trace->free_ranks.push_back(rank);
#endif
            if(rank!=free_count) return Failure("rank-deficient-free-subproblem");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++work.workspace_allocation_sites; // QR solution
#endif
            const Eigen::VectorXd solution=qr.solve(rhs);
            if(!solution.allFinite()) return Failure("non-finite-free-solution");
            for(Eigen::Index k=0;k<free_count;++k)
                candidate(free_indices[static_cast<std::size_t>(k)])=solution(k);
        }

        if(!candidate.allFinite()) return Failure("non-finite-candidate-step");
        double first_hit=1.0;
        std::vector<Eigen::Index> hit_indices;
        for(const auto index:free_indices)
        {
            if(candidate(index)<=upper(index)) continue;
            const double direction=candidate(index)-step(index);
            if(!(direction>0.0)) return Failure("inconsistent-active-set-direction");
            const double alpha=std::clamp((upper(index)-step(index))/direction,0.0,1.0);
            if(alpha<first_hit-bound_tolerance)
            {
                first_hit=alpha;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(hit_indices.empty() && !hit_indices.capacity()) ++work.workspace_allocation_sites;
#endif
                hit_indices.assign(1,index);
            }
            else if(std::abs(alpha-first_hit)<=bound_tolerance)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(hit_indices.size()==hit_indices.capacity()) ++work.workspace_allocation_sites;
#endif
                hit_indices.push_back(index);
            }
        }
        if(!hit_indices.empty())
        {
            step.noalias()+=first_hit*(candidate-step);
            for(const auto index:hit_indices)
            {
                step(index)=upper(index);
                active[static_cast<std::size_t>(index)]=true;
            }
            continue;
        }
        step=std::move(candidate);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.workspace_allocation_sites; // residual
#endif
        const Eigen::VectorXd residual=j*step+r;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.workspace_allocation_sites; // gradient
#endif
        Eigen::VectorXd gradient=j.transpose()*residual+
            lambda*d.cwiseProduct(d).cwiseProduct(step);
        if(!gradient.allFinite()) return Failure("non-finite-subproblem-gradient");
        const double kkt_tolerance=2e-11*KktScale(j,r,d,step,lambda,context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(trace)
        {
            trace->gradient_sequence.push_back(gradient);
            trace->kkt_tolerance_sequence.push_back(kkt_tolerance);
        }
#endif
        Eigen::Index release=-1;
        double largest_wrong_multiplier=kkt_tolerance;
        for(Eigen::Index i=0;i<n;++i)
        {
            if(!active[static_cast<std::size_t>(i)]) continue;
            if(gradient(i)>largest_wrong_multiplier)
            {
                largest_wrong_multiplier=gradient(i);
                release=i;
            }
        }
        if(release>=0)
        {
            active[static_cast<std::size_t>(release)]=false;
            ++releases;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++work.releases;
            if(trace) ++trace->releases;
#endif
            continue;
        }

        double kkt_inf{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.workspace_allocation_sites; // multipliers
#endif
        Eigen::VectorXd multipliers=Eigen::VectorXd::Zero(n);
        for(Eigen::Index i=0;i<n;++i)
        {
            if(active[static_cast<std::size_t>(i)])
                multipliers(i)=std::max(0.0,-gradient(i));
            else
                kkt_inf=std::max(kkt_inf,std::abs(gradient(i)));
            if(active[static_cast<std::size_t>(i)])
                kkt_inf=std::max(kkt_inf,std::max(0.0,gradient(i)));
        }
        if(kkt_inf>kkt_tolerance) return Failure("subproblem-kkt-failed");
        ActiveSetStep result;
        result.valid=true; result.status="solved"; result.step=std::move(step);
        result.gradient=std::move(gradient); result.multipliers=std::move(multipliers);
        result.active=std::move(active); result.damping=lambda;
        result.scaled_step_norm=d.cwiseProduct(result.step).norm();
        result.kkt_inf_norm=kkt_inf; result.iterations=iteration+1; result.releases=releases;
        return result;
    }
    return Failure("active-set-iteration-budget");
}

ActiveSetStep WithDamping(ActiveSetStep result,double damping)
{
    if(result.valid) result.damping=damping;
    return result;
}
}

ActiveSetStep SolveBoundedTrustRegionStep(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,
    const Eigen::VectorXd & d,const Eigen::VectorXd & upper,double trust_radius,double damping_seed)
{
    BoundedSubproblemContext context;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=BoundedLmWorkForTesting();
    BoundedLmTimer timer{work};
    const auto solve_candidate=[&](double lambda,int phase) {
        ++work.damping_solves;
        if(phase==1) ++work.bracket_solves;
        else if(phase==2) ++work.bisection_solves;
        BoundedLmDampingTrace trace;
        auto result=SolveAtDamping(context,j,r,d,upper,lambda,work.capture_trace ? &trace : nullptr);
        if(work.capture_trace)
        {
            trace.damping=lambda; trace.phase=phase; trace.status=result.status;
            trace.step=result.step; trace.gradient=result.gradient; trace.multipliers=result.multipliers;
            trace.active=result.active;
            if(result.iterations) trace.iterations=result.iterations;
            if(result.valid) trace.releases=result.releases;
            trace.scaled_step_norm=result.scaled_step_norm; trace.kkt_inf_norm=result.kkt_inf_norm;
            trace.valid=result.valid;
            trace.within_trust_region=result.valid && result.scaled_step_norm<=trust_radius;
            work.damping_trace.push_back(std::move(trace));
        }
        return result;
    };
    const auto mark_search_accepted=[&]() {
        if(work.capture_trace && !work.damping_trace.empty())
            work.damping_trace.back().search_accepted=true;
    };
#endif
    if(!std::isfinite(trust_radius) || !(trust_radius>0.0) || !std::isfinite(damping_seed) || damping_seed<0.0)
        return Failure("invalid-trust-region-input");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ++work.subproblem_validation_calls;
#endif
    if(!ValidSubproblemInput(j,r,d,upper,0.0)) return Failure("invalid-subproblem-input");
    context=MakeBoundedSubproblemContext(j,upper);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto unconstrained=solve_candidate(0.0,0);
#else
    auto unconstrained=SolveAtDamping(context,j,r,d,upper,0.0);
#endif
    const double trust_tolerance=64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,trust_radius);
    if(unconstrained.valid && unconstrained.scaled_step_norm<=trust_radius+trust_tolerance)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        mark_search_accepted();
#endif
        return WithDamping(std::move(unconstrained),0.0);
    }

    double low=0.0;
    double high=std::max(damping_seed,1e-12);
    ActiveSetStep feasible;
    for(int bracket=0;bracket<80;++bracket)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        auto candidate=solve_candidate(high,1);
#else
        auto candidate=SolveAtDamping(context,j,r,d,upper,high);
#endif
        if(candidate.valid && candidate.scaled_step_norm<=trust_radius)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            mark_search_accepted();
#endif
            feasible=std::move(candidate); break;
        }
        if(!candidate.valid && candidate.status!="rank-deficient-free-subproblem") return candidate;
        low=high;
        high*=4.0;
        if(!std::isfinite(high)) return Failure("trust-region-damping-overflow");
    }
    if(!feasible.valid) return Failure("trust-region-damping-bracket-failed");

    for(int iteration=0;iteration<64;++iteration)
    {
        const double middle=low+(high-low)*0.5;
        if(!(middle>low) || !(middle<high)) break;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        auto candidate=solve_candidate(middle,2);
#else
        auto candidate=SolveAtDamping(context,j,r,d,upper,middle);
#endif
        if(candidate.valid && candidate.scaled_step_norm<=trust_radius)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            mark_search_accepted();
#endif
            high=middle; feasible=std::move(candidate);
        }
        else if(candidate.valid || candidate.status=="rank-deficient-free-subproblem") low=middle;
        else return candidate;
    }
    feasible.damping=high;
    if(feasible.scaled_step_norm>trust_radius+trust_tolerance)
        return Failure("trust-region-feasibility-failed");
    return feasible;
}

} // namespace rhbm_gem::core::joint_component
