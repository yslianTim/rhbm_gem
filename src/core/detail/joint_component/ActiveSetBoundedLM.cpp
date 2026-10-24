#include "ActiveSetBoundedLM.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
double KktScale(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,const Eigen::VectorXd & d,
    const Eigen::VectorXd & p,double lambda)
{
    const double residual=(j*p+r).norm();
    return std::max(1.0,j.norm()*(residual+j.norm()*p.norm())+
        lambda*d.squaredNorm()*std::max(1.0,p.norm()));
}

ActiveSetStep Failure(std::string status)
{
    ActiveSetStep result; result.status=std::move(status); return result;
}

ActiveSetStep SolveAtDamping(const Eigen::MatrixXd & j,const Eigen::VectorXd & r,
    const Eigen::VectorXd & d,const Eigen::VectorXd & upper,double lambda)
{
    const Eigen::Index n=j.cols();
    if(j.rows()!=r.size() || n==0 || d.size()!=n || upper.size()!=n ||
        !j.allFinite() || !r.allFinite() || !d.allFinite() || !upper.allFinite() ||
        !(d.array()>0.0).all() || (upper.array()<0.0).any() ||
        !std::isfinite(lambda) || lambda<0.0)
        return Failure("invalid-subproblem-input");

    const double bound_tolerance=64.0*std::numeric_limits<double>::epsilon()*
        std::max(1.0,upper.cwiseAbs().maxCoeff());
    std::vector<bool> active(static_cast<std::size_t>(n),false);
    Eigen::VectorXd step=Eigen::VectorXd::Zero(n);
    for(Eigen::Index i=0;i<n;++i)
    {
        active[static_cast<std::size_t>(i)]=upper(i)<=bound_tolerance;
        if(active[static_cast<std::size_t>(i)]) step(i)=upper(i);
    }

    const std::size_t maximum_iterations=std::max<std::size_t>(40,20*static_cast<std::size_t>(n));
    std::size_t releases{};
    for(std::size_t iteration=0;iteration<maximum_iterations;++iteration)
    {
        std::vector<Eigen::Index> free_indices;
        for(Eigen::Index i=0;i<n;++i)
            if(!active[static_cast<std::size_t>(i)]) free_indices.push_back(i);

        Eigen::VectorXd candidate=step;
        if(!free_indices.empty())
        {
            const Eigen::Index free_count=static_cast<Eigen::Index>(free_indices.size());
            Eigen::MatrixXd augmented(j.rows()+(lambda>0.0 ? free_count : 0),free_count);
            Eigen::VectorXd rhs(augmented.rows());
            augmented.setZero();
            rhs.setZero();
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
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(augmented);
            qr.setThreshold(64.0*std::numeric_limits<double>::epsilon()*
                static_cast<double>(std::max(augmented.rows(),augmented.cols())));
            if(qr.rank()!=free_count) return Failure("rank-deficient-free-subproblem");
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
                first_hit=alpha; hit_indices.assign(1,index);
            }
            else if(std::abs(alpha-first_hit)<=bound_tolerance)
                hit_indices.push_back(index);
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
        const Eigen::VectorXd residual=j*step+r;
        Eigen::VectorXd gradient=j.transpose()*residual+
            lambda*d.cwiseProduct(d).cwiseProduct(step);
        if(!gradient.allFinite()) return Failure("non-finite-subproblem-gradient");
        const double kkt_tolerance=2e-11*KktScale(j,r,d,step,lambda);
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
            continue;
        }

        double kkt_inf{};
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
    if(!std::isfinite(trust_radius) || !(trust_radius>0.0) || !std::isfinite(damping_seed) || damping_seed<0.0)
        return Failure("invalid-trust-region-input");
    auto unconstrained=SolveAtDamping(j,r,d,upper,0.0);
    const double trust_tolerance=64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,trust_radius);
    if(unconstrained.valid && unconstrained.scaled_step_norm<=trust_radius+trust_tolerance)
        return WithDamping(std::move(unconstrained),0.0);

    double low=0.0;
    double high=std::max(damping_seed,1e-12);
    ActiveSetStep feasible;
    for(int bracket=0;bracket<80;++bracket)
    {
        auto candidate=SolveAtDamping(j,r,d,upper,high);
        if(candidate.valid && candidate.scaled_step_norm<=trust_radius)
        {feasible=std::move(candidate); break;}
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
        auto candidate=SolveAtDamping(j,r,d,upper,middle);
        if(candidate.valid && candidate.scaled_step_norm<=trust_radius)
        {high=middle; feasible=std::move(candidate);}
        else if(candidate.valid || candidate.status=="rank-deficient-free-subproblem") low=middle;
        else return candidate;
    }
    feasible.damping=high;
    if(feasible.scaled_step_norm>trust_radius+trust_tolerance)
        return Failure("trust-region-feasibility-failed");
    return feasible;
}

} // namespace rhbm_gem::core::joint_component
