#include "Numerics.hpp"
#include "ResourceWork.hpp"
#include "TiledDerivative.hpp"
#include "CompactSvd.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <numbers>
#include <utility>

namespace rhbm_gem::core::joint_component {
namespace {
void AddCompensated(double term,double & sum,double & compensation)
{
    const double next=sum+term;
    compensation+=std::abs(sum)>=std::abs(term) ? (sum-next)+term : (term-next)+sum;
    sum=next;
}
void AddCompensatedProduct(double left,double right,double & sum,double & compensation)
{
    const double product=left*right;
    AddCompensated(product,sum,compensation);
    AddCompensated(std::fma(left,right,-product),sum,compensation);
}
void AddLinearSolveWork(LinearSolveWork & destination,const LinearSolveWork & source)
{
#define ADD_COUNT(name) destination.name+=source.name
    ADD_COUNT(solve_calls); ADD_COUNT(valid_solves); ADD_COUNT(invalid_solves); ADD_COUNT(invalid_input_solves);
    ADD_COUNT(rank_deficient_solves); ADD_COUNT(nonfinite_solves); ADD_COUNT(active_set_budget_exhausted_solves);
    ADD_COUNT(other_invalid_solves); ADD_COUNT(factorization_solves); ADD_COUNT(active_set_iterations);
    ADD_COUNT(active_set_releases); ADD_COUNT(initial_kkt_count); ADD_COUNT(initial_kkt_le_1e10);
    ADD_COUNT(initial_kkt_le_1e8); ADD_COUNT(initial_kkt_le_1e7); ADD_COUNT(final_kkt_count);
    ADD_COUNT(final_kkt_le_1e10); ADD_COUNT(final_kkt_le_1e8); ADD_COUNT(final_kkt_le_1e7);
    ADD_COUNT(refinement_eligibility_count); ADD_COUNT(refinement_trigger_count); ADD_COUNT(refinement_iterations);
    ADD_COUNT(refinement_success_count); ADD_COUNT(refinement_failure_count); ADD_COUNT(refinement_factorizations);
#undef ADD_COUNT
    destination.initial_kkt_sum+=source.initial_kkt_sum;
    destination.initial_kkt_max=std::max(destination.initial_kkt_max,source.initial_kkt_max);
    destination.final_kkt_sum+=source.final_kkt_sum;
    destination.final_kkt_max=std::max(destination.final_kkt_max,source.final_kkt_max);
    destination.refinement_seconds+=source.refinement_seconds;
    destination.refinement_before_kkt_sum+=source.refinement_before_kkt_sum;
    destination.refinement_before_kkt_max=std::max(destination.refinement_before_kkt_max,source.refinement_before_kkt_max);
    destination.refinement_after_kkt_sum+=source.refinement_after_kkt_sum;
    destination.refinement_after_kkt_max=std::max(destination.refinement_after_kkt_max,source.refinement_after_kkt_max);
    destination.refinement_before_objective_sum+=source.refinement_before_objective_sum;
    destination.refinement_after_objective_sum+=source.refinement_after_objective_sum;
}
void AddKktSample(double value,std::size_t & count,std::size_t & le_1e10,std::size_t & le_1e8,
    std::size_t & le_1e7,double & sum,double & maximum)
{
    if(!std::isfinite(value) || value<0) return;
    ++count; sum+=value; maximum=std::max(maximum,value);
    if(value<=1e-10) ++le_1e10;
    if(value<=1e-8) ++le_1e8;
    if(value<=1e-7) ++le_1e7;
}
void AddRoleWork(ProfileRoleWork & destination,const ProfileRoleWork & source)
{
    destination.evaluations+=source.evaluations;
    destination.derivative_preparations+=source.derivative_preparations;
    destination.derivative_reductions+=source.derivative_reductions;
    destination.replay_checks+=source.replay_checks;
    destination.local_kkt_rejections+=source.local_kkt_rejections;
    AddLinearSolveWork(destination.linear_solve,source.linear_solve);
    destination.evaluation_seconds+=source.evaluation_seconds;
    destination.basis_seconds+=source.basis_seconds;
    destination.linear_matrix_preparation_seconds+=source.linear_matrix_preparation_seconds;
    destination.linear_rhs_solve_seconds+=source.linear_rhs_solve_seconds;
    destination.linear_certificate_seconds+=source.linear_certificate_seconds;
    destination.derivative_prepare_seconds+=source.derivative_prepare_seconds;
    destination.derivative_reduce_seconds+=source.derivative_reduce_seconds;
    destination.derivative_raw_assembly_seconds+=source.derivative_raw_assembly_seconds;
    destination.derivative_free_design_assembly_seconds+=source.derivative_free_design_assembly_seconds;
    destination.derivative_factor_match_seconds+=source.derivative_factor_match_seconds;
    destination.derivative_factor_build_seconds+=source.derivative_factor_build_seconds;
    destination.derivative_factor_compact_seconds+=source.derivative_factor_compact_seconds;
    destination.derivative_rank_seconds+=source.derivative_rank_seconds;
    destination.derivative_least_squares_seconds+=source.derivative_least_squares_seconds;
    destination.derivative_normal_solve_seconds+=source.derivative_normal_solve_seconds;
    destination.derivative_cancellation_check_seconds+=source.derivative_cancellation_check_seconds;
    destination.derivative_cancellation_fallback_seconds+=source.derivative_cancellation_fallback_seconds;
    destination.derivative_rows_seconds+=source.derivative_rows_seconds;
    destination.derivative_jacobian_qr_seconds+=source.derivative_jacobian_qr_seconds;
    destination.derivative_norms_seconds+=source.derivative_norms_seconds;
    destination.derivative_outer_overhead_seconds+=source.derivative_outer_overhead_seconds;
    destination.tiled_qr_assembly_copy_seconds+=source.tiled_qr_assembly_copy_seconds;
    destination.tiled_qr_householder_seconds+=source.tiled_qr_householder_seconds;
    destination.tiled_qr_rhs_transform_seconds+=source.tiled_qr_rhs_transform_seconds;
    destination.replay_trust_seconds+=source.replay_trust_seconds;
}
}
WidthBoundKkt CalculateWidthBoundKkt(const Vector & eta,const Vector & gradient,double upper)
{
    WidthBoundKkt out;
    if(!std::isfinite(upper) || eta.size()==0 || eta.size()!=gradient.size() ||
        !eta.allFinite() || !gradient.allFinite()) return out;
    const double tolerance=64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(upper));
    out.feasible=true;
    out.projected=0.0;
    for(Eigen::Index i=0;i<eta.size();++i)
    {
        const double slack=upper-eta(i);
        out.violation=std::max(out.violation,std::max(0.0,-slack));
        out.feasible &= slack>=-tolerance;
        const double projected=eta(i)-std::min(eta(i)-gradient(i),upper);
        out.projected=std::max(out.projected,std::abs(projected));
        if(slack<=tolerance)
        {
            out.boundary_gradient=out.active==0 ? gradient(i) : std::max(out.boundary_gradient,gradient(i));
            ++out.active;
            const double multiplier=std::max(0.0,-gradient(i));
            out.multiplier=std::max(out.multiplier,multiplier);
            out.complementarity=std::max(out.complementarity,std::abs(multiplier*slack));
        }
    }
    return out;
}
namespace {
void CaptureWidthBoundKkt(Certificate & certificate,const Vector & eta,const Vector & gradient,
    const EvaluationContext & context)
{
    if(!context.eta_upper_bound) return;
    const auto kkt=CalculateWidthBoundKkt(eta,gradient,*context.eta_upper_bound);
    certificate.width_bound_applied=true;
    certificate.width_bound_feasible=kkt.feasible;
    certificate.width_bound_projected_kkt=kkt.projected;
    certificate.maximum_width_bound_violation=kkt.violation;
    certificate.boundary_gradient_max=kkt.boundary_gradient;
    certificate.boundary_multiplier_max=kkt.multiplier;
    certificate.width_bound_complementarity=kkt.complementarity;
    certificate.active_width_bound_count=kkt.active;
}
}
bool IsValidJointConvergenceCriteria(const JointConvergenceCriteria & criteria)
{
    return std::isfinite(criteria.global_ac_kkt) && criteria.global_ac_kkt>=0.0 &&
        std::isfinite(criteria.global_width_gradient) && criteria.global_width_gradient>=0.0 &&
        std::isfinite(criteria.confirmed_log_width_change) && criteria.confirmed_log_width_change>=0.0;
}
namespace {
bool MeetsJointThreshold(double value,double threshold)
{return std::isfinite(value) && value>=0.0 && value<=threshold;}
}
ProfileRoleWork & ProfileSearchWork::Role(ProfileEvaluationRole role)
{
    switch(role)
    {
    case ProfileEvaluationRole::InitialProfile: return initial_profile;
    case ProfileEvaluationRole::TrialProfile: return trial_profile;
    case ProfileEvaluationRole::AcceptedEndpoint: return accepted_endpoint;
    case ProfileEvaluationRole::Unspecified: return total;
    }
    return total;
}
void ProfileSearchWork::Add(ProfileEvaluationRole role,const ProfileRoleWork & work)
{
    if(role!=ProfileEvaluationRole::Unspecified) AddRoleWork(Role(role),work);
    AddRoleWork(total,work);
}
void BoundedLmWork::Merge(const BoundedLmWork & other)
{
    subproblem_solves+=other.subproblem_solves;
    failed_subproblems+=other.failed_subproblems;
    active_set_iterations+=other.active_set_iterations;
    releases+=other.releases;
    active_bound_occurrences+=other.active_bound_occurrences;
    if(std::isfinite(other.maximum_subproblem_kkt))
        maximum_subproblem_kkt=std::isfinite(maximum_subproblem_kkt) ?
            std::max(maximum_subproblem_kkt,other.maximum_subproblem_kkt) : other.maximum_subproblem_kkt;
}
void ProfileSearchWork::Merge(const ProfileSearchWork & other)
{
    AddRoleWork(total,other.total); AddRoleWork(initial_profile,other.initial_profile);
    AddRoleWork(trial_profile,other.trial_profile); AddRoleWork(accepted_endpoint,other.accepted_endpoint);
    bounded_lm.Merge(other.bounded_lm);
    accepted_profile_trials+=other.accepted_profile_trials;
    rejected_profile_trials+=other.rejected_profile_trials;
    local_replay_kkt_rejections+=other.local_replay_kkt_rejections;
    local_replay_disagreements+=other.local_replay_disagreements;
    invalid_gradient_stops+=other.invalid_gradient_stops;
    width_safeguard_rejections+=other.width_safeguard_rejections;
    accepted_update_budget_stops+=other.accepted_update_budget_stops;
    profile_budget_stops+=other.profile_budget_stops;
    cpu_budget_stops+=other.cpu_budget_stops;
    native_lm_stops+=other.native_lm_stops;
    for(std::size_t i=0;i<lm_status_counts.size();++i) lm_status_counts[i]+=other.lm_status_counts[i];
    total_seconds+=other.total_seconds; lm_overhead_seconds+=other.lm_overhead_seconds;
}
double ProfileSearchWork::AttributedSeconds() const
{
    return total.basis_seconds+total.linear_matrix_preparation_seconds+total.linear_rhs_solve_seconds+
        total.linear_certificate_seconds+
        total.derivative_prepare_seconds+total.derivative_reduce_seconds+total.replay_trust_seconds;
}
namespace {
Vector StablePrediction(const Sparse & x,const Vector & beta)
{
    Vector prediction=Vector::Zero(x.rows()),compensation=prediction;
    for(Eigen::Index column=0;column<x.cols();++column)
        for(Sparse::InnerIterator entry(x,column);entry;++entry)
        {
            const Eigen::Index row=entry.row();
            double sum=prediction(row),correction=compensation(row);
            AddCompensatedProduct(entry.value(),beta(column),sum,correction);
            prediction(row)=sum; compensation(row)=correction;
        }
    prediction+=compensation;
    return prediction;
}
Vector StableResidual(const Sparse & x,VectorRef y,const Vector & beta)
{
    return StablePrediction(x,beta)-y;
}
Vector StableTransposeProduct(const Sparse & x,VectorRef values)
{
    Vector product=Vector::Zero(x.cols());
    for(Eigen::Index column=0;column<x.cols();++column)
    {
        double sum{},compensation{};
        for(Sparse::InnerIterator entry(x,column);entry;++entry)
            AddCompensatedProduct(entry.value(),values(entry.row()),sum,compensation);
        product(column)=sum+compensation;
    }
    return product;
}
double StableColumnDot(const Sparse & x,Eigen::Index column,VectorRef values)
{
    double sum{},compensation{};
    for(Sparse::InnerIterator entry(x,column);entry;++entry)
        AddCompensatedProduct(entry.value(),values(entry.row()),sum,compensation);
    return sum+compensation;
}
// Full-column TSQR from original rows. RHS columns undergo the same orthogonal
// transformations. Only the small R and transformed RHS survive each tile.
template<class Design>
std::pair<Matrix,Matrix> Reduce(const Design & x,const Matrix & rhs)
{
    TiledQR reduced(x.cols(),rhs.cols());
    for(Eigen::Index first=0;first<x.rows();first+=derivative_tile_rows)
    {
        const auto count=std::min(derivative_tile_rows,x.rows()-first);
        reduced.Append(Matrix(x.middleRows(first,count)),rhs.middleRows(first,count));
    }
    return {std::move(reduced.r),std::move(reduced.target)};
}
}
double RankPolicy::Relative(Eigen::Index columns) const
{return std::numeric_limits<double>::epsilon()*static_cast<double>(std::max(rows,columns));}
double RankPolicy::Absolute(Eigen::Index columns,double maximum) const
{return Relative(columns)*maximum;}
EvaluationContext CreateContext(std::shared_ptr<const JointProblemInput> input,const std::string & hash,const AuditPlan & plan)
{
    EvaluationContext c; c.snapshot_hash=hash; c.observations=Observe(input);
    const auto & y=*c.observations; const auto atoms=static_cast<Eigen::Index>(input->atom_ids.size());
    c.scale=std::max(1.0,y.norm()); c.rank={y.size(),2*atoms,atoms};
    c.linear.rank_relative=c.rank.Relative(2*atoms); c.linear.release_response_norm=y.norm(); c.audit=plan;
    c.atom_ids=Identities(std::shared_ptr<const std::vector<std::string>>(input,&input->atom_ids));
    c.row_ids=Identities(std::shared_ptr<const std::vector<std::string>>(input,&input->row_ids)); return c;
}
EvaluationContext CreateContext(VectorRef y,Eigen::Index atoms,const std::string & hash,const AuditPlan & plan)
{
    auto input=std::make_shared<JointProblemInput>(); input->observations.assign(y.data(),y.data()+y.size());
    for(Eigen::Index a=0;a<atoms;++a) input->atom_ids.push_back(std::to_string(a));
    for(Eigen::Index r=0;r<y.size();++r) input->row_ids.push_back(std::to_string(r));
    return CreateContext(std::move(input),hash,plan);
}
BasisValues EvaluateKernel(double square, double width, double cutoff)
{
    if (!(width>0.0) || !std::isfinite(width)) throw std::invalid_argument("Invalid matched width.");
    if (square > cutoff*cutoff) return {};
    const double r{std::sqrt(square)}, r2{square}, exponent{std::exp(-r2/(2.0*width*width))};
    BasisValues out;
    out.gaussian = std::pow(2.0*std::numbers::pi*width*width,-1.5)*exponent;
    out.gaussian_log_width = out.gaussian*(r2/(width*width)-3.0);
    const double center{std::sqrt(2.0/std::numbers::pi)/width};
    if (r < 1e-5) {out.charge=center; out.charge_log_width=-center;}
    else if (r <= 2.5) {out.charge=std::erf(r/width/std::sqrt(2.0))/r; out.charge_log_width=-center*exponent;}
    return out;
}
Certificate CertifyLinear(const Sparse & x,VectorRef y,const Eigen::VectorXd & beta,double observation_scale)
{
    Certificate out; out.evaluated=true;
    if (x.rows()!=y.size() || x.cols()!=beta.size() || !beta.allFinite() || !y.allFinite()) return out;
    Eigen::VectorXd norms(x.cols());
    for (Eigen::Index k=0;k<x.cols();++k) norms(k)=x.col(k).norm();
    if (!norms.allFinite() || (norms.array()<=0).any()) return out;
    const Eigen::VectorXd residual=StableResidual(x,y,beta);
    if (!residual.allFinite()) return out;
    const double s=observation_scale>0 ? observation_scale : std::max(1.0,y.norm());
    const Eigen::VectorXd u=norms.array()*beta.array()/s;
    const Eigen::VectorXd gradient=StableTransposeProduct(x,residual).array()/norms.array()/s;
    Eigen::VectorXd projected=u-gradient; std::vector<Eigen::Index> active; bool feasible=true;
    for (Eigen::Index k=0;k<beta.size();k+=2)
    {
        feasible &= beta(k)>=0; projected(k)=std::max(0.0,projected(k));
        if (beta(k)==0) active.push_back(k/2);
    }
    const double kkt=(u-projected).lpNorm<Eigen::Infinity>();
    out.available=true; out.feasible=feasible; out.projected_kkt=kkt; out.kkt_passed=feasible && kkt<=1e-10;
    out.rss=residual.squaredNorm(); out.objective=.5*residual.squaredNorm();
    out.residual_scale=residual.squaredNorm()/static_cast<double>(y.size());
    out.residual_rmse=residual.norm()/std::sqrt(static_cast<double>(y.size()));
    out.residual_max=residual.cwiseAbs().maxCoeff(); out.relative_residual=residual.norm()/s;
    out.active_atoms=active; return out;
}

std::optional<double> NormalizedProfileObjective(const Evaluation & evaluation,double scale)
{
    if(!evaluation.valid || !evaluation.certificate.available || !(scale>0) || !std::isfinite(scale))
        return std::nullopt;
    const double objective=evaluation.certificate.objective/(scale*scale);
    if(!std::isfinite(objective)) return std::nullopt;
    return objective;
}

std::optional<double> ProfileGradientInfinityNorm(const Evaluation & evaluation)
{
    if(!evaluation.valid || evaluation.gradient.size()==0 || !evaluation.gradient.allFinite())
        return std::nullopt;
    const double value=evaluation.gradient.lpNorm<Eigen::Infinity>();
    if(!std::isfinite(value)) return std::nullopt;
    return value;
}

std::optional<double> ProfileGradientInfinityNorm(const Evaluation & evaluation,const EvaluationContext & context)
{
    if(context.eta_upper_bound)
    {
        if(!evaluation.valid || evaluation.eta.size()!=evaluation.gradient.size()) return std::nullopt;
        const auto kkt=CalculateWidthBoundKkt(evaluation.eta,evaluation.gradient,
            *context.eta_upper_bound);
        if(!std::isfinite(kkt.projected)) return std::nullopt;
        return kkt.projected;
    }
    return ProfileGradientInfinityNorm(evaluation);
}

void UpdateAcceptedProfileObjective(std::optional<double> & accepted_objective,
    const Evaluation & evaluation,double scale,bool accepted)
{
    if(accepted) accepted_objective=NormalizedProfileObjective(evaluation,scale);
}

namespace {
struct EvaluationTelemetry
{
    ProfileSearchWork * target{};
    ProfileEvaluationRole role{ProfileEvaluationRole::Unspecified};
    ProfileRoleWork work;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
    EvaluationTelemetry(ProfileSearchWork * target,ProfileEvaluationRole role):target(target),role(role)
    {if(target) ++work.evaluations;}
    ~EvaluationTelemetry()
    {
        if(target)
        {
            work.evaluation_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            target->Add(role,work);
        }
    }
};
Evaluation Basis(const Domain & domain,VectorRef y,const Vector & eta)
{
    ResourcePhase phase("profile-basis-build",true,domain.rows,2*eta.size());
    Evaluation out; out.eta=eta;
    if (domain.rows!=y.size() || domain.rows==0 || eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) ||
        eta.size()==0 || !eta.allFinite() || !y.allFinite()) {out.reason="invalid-input"; return out;}
    const Vector widths=eta.array().exp();
    if (!widths.allFinite() || (widths.array()<=0).any()) {out.reason="invalid-b"; return out;}
    const Eigen::Index m=eta.size(); out.x.resize(y.size(),2*m); out.derivative.resize(y.size(),2*m);
    std::vector<Eigen::Triplet<double>> x,dx;
    std::size_t memberships{}; for (const auto & atom:domain.atoms) memberships+=atom.size();
    x.reserve(2*memberships); dx.reserve(2*memberships);
    for (Eigen::Index a=0;a<m;++a) for (const auto & p:domain.atoms[static_cast<std::size_t>(a)])
    {
        const auto b=EvaluateKernel(p.square,widths(a),2.5);
        const std::array<double,4> v{b.gaussian,b.charge,b.gaussian_log_width,b.charge_log_width};
        if (!std::all_of(v.begin(),v.end(),[](double n){return std::isfinite(n);}))
        {out.reason="nonfinite-basis"; return out;}
        for (Eigen::Index k=0;k<2;++k)
        {
            if (v[static_cast<std::size_t>(k)]!=0) x.emplace_back(p.row,2*a+k,v[static_cast<std::size_t>(k)]);
            if (v[static_cast<std::size_t>(k+2)]!=0) dx.emplace_back(p.row,2*a+k,v[static_cast<std::size_t>(k+2)]);
        }
    }
    out.x.setFromTriplets(x.begin(),x.end()); out.derivative.setFromTriplets(dx.begin(),dx.end());
    RecordSparseShape("profile-design",out.x.rows(),out.x.cols(),static_cast<std::size_t>(out.x.nonZeros()));
    out.valid=true; return out;
}
}
MemoryEstimate EndpointVerificationMemoryEstimate(const Domain & domain,VectorRef y,const EvaluationContext &,
    std::size_t resident_override)
{
    MemoryShape shape; shape.component_atoms=domain.atoms.size();
    shape.component_observations=static_cast<std::size_t>(domain.rows);
    for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        if(!CheckedAddBytes(shape.component_support_entries,domain.atoms[atom].size(),shape.component_support_entries))
        {
            auto estimate=PlanMemory(shape,MemoryStage::EndpointVerification);
            estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
        }
    const auto * input=domain.atoms.Storage();
    const auto snapshot=input ? EstimateSnapshotBytes(*input) : std::optional<std::size_t>{0};
    std::size_t resident{};
    std::size_t support_bytes{},observation_bytes{};
    if(!snapshot || !CheckedMultiplyBytes(shape.component_support_entries,sizeof(JointSupport),support_bytes) ||
        !CheckedMultiplyBytes(static_cast<std::size_t>(std::max<Eigen::Index>(0,y.size())),2*sizeof(double),observation_bytes) ||
        !CheckedAddBytes(*snapshot,domain.atoms.LocalSupportIndexBytes(),resident) ||
        !CheckedAddBytes(resident,support_bytes,resident) || !CheckedAddBytes(resident,observation_bytes,resident))
    {
        auto estimate=PlanMemory(shape,MemoryStage::EndpointVerification);
        estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
    }
    shape.resident_bytes=std::max(resident,resident_override);
    return PlanMemory(shape,MemoryStage::EndpointVerification);
}
Evaluation EvaluateProfile(const Domain & domain,VectorRef y,const Vector & eta,const EvaluationContext * context,
    const std::vector<LinearBlock> * blocks,
    ProfileEvaluationRole role,ProfileSearchWork * profile_work)
{
    ResourcePhase phase("ac-profile");
    ProfileEvaluationRoleScopeForTesting role_scope(role);
    EvaluationTelemetry telemetry(profile_work,role);
    const auto basis_started=std::chrono::steady_clock::now();
    auto out=Basis(domain,y,eta);
    if(profile_work) telemetry.work.basis_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-basis_started).count();
    if(!out.valid) return out; out.valid=false;
    const Eigen::Index m=eta.size();
    ResourcePhase linear_solve("linear-solve",true,out.x.rows(),out.x.cols(),static_cast<std::size_t>(out.x.nonZeros()));
    const auto solve_started=std::chrono::steady_clock::now();
    const auto solved=SolveLinear(out.x,y,Vector::Ones(y.size()),
        context ? &context->linear : nullptr,blocks,context ? context->scale : 0);
    if(profile_work)
    {
        auto & stats=telemetry.work.linear_solve;
        ++stats.solve_calls;
        stats.factorization_solves+=static_cast<std::size_t>(std::max(0,solved.solves));
        stats.active_set_iterations+=solved.active_set_iterations;
        stats.active_set_releases+=static_cast<std::size_t>(std::max(0,solved.releases));
        stats.refinement_eligibility_count+=solved.refinement_eligibility_count;
        stats.refinement_trigger_count+=solved.refinement_trigger_count;
        stats.refinement_iterations+=solved.refinement_iterations;
        stats.refinement_success_count+=solved.refinement_success_count;
        stats.refinement_failure_count+=solved.refinement_failure_count;
        stats.refinement_factorizations+=solved.refinement_factorizations;
        stats.refinement_seconds+=solved.refinement_seconds;
        AddKktSample(solved.initial_projected_kkt,stats.initial_kkt_count,stats.initial_kkt_le_1e10,
            stats.initial_kkt_le_1e8,stats.initial_kkt_le_1e7,stats.initial_kkt_sum,stats.initial_kkt_max);
        if(solved.refinement_trigger_count)
        {
            stats.refinement_before_kkt_sum+=solved.initial_projected_kkt;
            stats.refinement_before_kkt_max=std::max(stats.refinement_before_kkt_max,solved.initial_projected_kkt);
            stats.refinement_before_objective_sum+=solved.refinement_before_objective;
            stats.refinement_after_objective_sum+=solved.refinement_after_objective;
        }
        if(solved.valid) ++stats.valid_solves;
        else
        {
            ++stats.invalid_solves;
            if(solved.reason=="invalid-input") ++stats.invalid_input_solves;
            else if(solved.reason.starts_with("rank-deficient")) ++stats.rank_deficient_solves;
            else if(solved.reason=="nonfinite") ++stats.nonfinite_solves;
            else if(solved.reason=="active-set-budget-exhausted") ++stats.active_set_budget_exhausted_solves;
            else ++stats.other_invalid_solves;
        }
    }
    if(profile_work)
    {
        telemetry.work.linear_rhs_solve_seconds+=
            std::chrono::duration<double>(std::chrono::steady_clock::now()-solve_started).count();
    }
    out.beta=solved.beta;
    const auto certificate_started=std::chrono::steady_clock::now();
    out.certificate=CertifyLinear(out.x,y,out.beta,context ? context->scale : 0);
    if(profile_work) telemetry.work.linear_certificate_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-certificate_started).count();
    out.certificate.linear_solver_valid=solved.valid; out.certificate.linear_solver_reason=solved.reason;
    out.certificate.linear_solves=solved.solves; out.certificate.linear_releases=solved.releases;
    out.certificate.free_rank=solved.rank;
    if(blocks) out.certificate.block_factorizations=solved.block_factorizations;
    if(profile_work && out.certificate.available)
    {
        auto & stats=telemetry.work.linear_solve;
        AddKktSample(out.certificate.projected_kkt,stats.final_kkt_count,stats.final_kkt_le_1e10,
            stats.final_kkt_le_1e8,stats.final_kkt_le_1e7,stats.final_kkt_sum,stats.final_kkt_max);
        if(solved.refinement_trigger_count)
        {
            stats.refinement_after_kkt_sum+=out.certificate.projected_kkt;
            stats.refinement_after_kkt_max=std::max(stats.refinement_after_kkt_max,out.certificate.projected_kkt);
        }
    }
    if (!solved.valid) {out.reason=solved.reason; return out;}
    if (!out.certificate.kkt_passed)
    {if(profile_work) ++telemetry.work.local_kkt_rejections; out.reason="kkt-failed"; return out;}
    out.residual=StableResidual(out.x,y,out.beta); out.gradient=Vector::Zero(m);
    const double scale=context ? context->scale : std::max(1.0,y.norm());
    for (Eigen::Index k=0;k<2*m;++k)
        out.gradient(k/2)+=out.beta(k)*StableColumnDot(out.derivative,k,out.residual)/scale/scale;
    if(context) CaptureWidthBoundKkt(out.certificate,out.eta,out.gradient,*context);
    out.valid=out.residual.allFinite() && out.gradient.allFinite();
    out.reason=out.valid ? "qualified-inner" : "nonfinite-residual"; return out;
}

Evaluation EvaluateState(const Domain & domain,VectorRef y,const Vector & eta,const Vector & beta,const EvaluationContext & context,
    ProfileEvaluationRole role,ProfileSearchWork * profile_work)
{
    ProfileEvaluationRoleScopeForTesting role_scope(role);
    EvaluationTelemetry telemetry(profile_work,role);
    const auto basis_started=std::chrono::steady_clock::now();
    auto out=Basis(domain,y,eta);
    if(profile_work) telemetry.work.basis_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-basis_started).count();
    if(!out.valid) return out; out.valid=false;
    if(beta.size()!=out.x.cols() || !beta.allFinite()) {out.reason="invalid-coefficients"; return out;}
    out.beta=beta; out.residual=StableResidual(out.x,y,beta); out.gradient=Vector::Zero(eta.size());
    const auto certificate_started=std::chrono::steady_clock::now();
    out.certificate=CertifyLinear(out.x,y,beta,context.scale);
    if(profile_work) telemetry.work.linear_certificate_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-certificate_started).count();
    for(Eigen::Index k=0;k<beta.size();++k)
        out.gradient(k/2)+=beta(k)*StableColumnDot(out.derivative,k,out.residual)/context.scale/context.scale;
    CaptureWidthBoundKkt(out.certificate,out.eta,out.gradient,context);
    out.valid=out.residual.allFinite() && out.gradient.allFinite(); out.reason=out.valid ? "raw-state" : "nonfinite-state";
    return out;
}
Evaluation EvaluateStateStreaming(const Domain & domain,VectorRef y,const Vector & eta,const Vector & beta,
    const EvaluationContext & context,ProfileEvaluationRole role,ProfileSearchWork * profile_work)
{
    ResourcePhase phase("global-state-evaluation",true,domain.rows,2*eta.size());
    ProfileEvaluationRoleScopeForTesting role_scope(role);
    EvaluationTelemetry telemetry(profile_work,role);
    Evaluation out; out.eta=eta;
    if(domain.rows!=y.size() || domain.rows==0 || eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) ||
        eta.size()==0 || !eta.allFinite() || !y.allFinite())
    {out.reason="invalid-input"; return out;}
    if(beta.size()!=2*eta.size() || !beta.allFinite()) {out.reason="invalid-coefficients"; return out;}
    const Vector widths=eta.array().exp();
    if(!widths.allFinite() || (widths.array()<=0).any()) {out.reason="invalid-b"; return out;}

    Vector prediction=Vector::Zero(y.size()),norms_squared=Vector::Zero(beta.size());
    for(Eigen::Index atom=0;atom<eta.size();++atom)
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto basis=EvaluateKernel(support.square,widths(atom),2.5);
            if(!std::isfinite(basis.gaussian) || !std::isfinite(basis.charge) ||
                !std::isfinite(basis.gaussian_log_width) || !std::isfinite(basis.charge_log_width))
            {out.reason="nonfinite-basis"; return out;}
            prediction(support.row)+=beta(2*atom)*basis.gaussian+beta(2*atom+1)*basis.charge;
            norms_squared(2*atom)+=basis.gaussian*basis.gaussian;
            norms_squared(2*atom+1)+=basis.charge*basis.charge;
        }
    out.beta=beta; out.residual=std::move(prediction); out.residual-=y;
    Vector ac_gradient=Vector::Zero(beta.size());
    Vector gaussian_width_dot=Vector::Zero(eta.size()),charge_width_dot=Vector::Zero(eta.size());
    for(Eigen::Index atom=0;atom<eta.size();++atom)
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto basis=EvaluateKernel(support.square,widths(atom),2.5);
            const double residual=out.residual(support.row);
            ac_gradient(2*atom)+=basis.gaussian*residual;
            ac_gradient(2*atom+1)+=basis.charge*residual;
            gaussian_width_dot(atom)+=basis.gaussian_log_width*residual;
            charge_width_dot(atom)+=basis.charge_log_width*residual;
        }
    out.gradient.resize(eta.size());
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {
        out.gradient(atom)=beta(2*atom)*gaussian_width_dot(atom)/context.scale/context.scale;
        out.gradient(atom)+=beta(2*atom+1)*charge_width_dot(atom)/context.scale/context.scale;
    }
    if(!norms_squared.allFinite() || (norms_squared.array()<=0).any())
    {out.reason="rank-deficient"; return out;}
    const Vector norms=norms_squared.array().sqrt();
    const double certificate_scale=context.scale>0 ? context.scale : std::max(1.0,y.norm());
    const Vector scaled_beta=norms.array()*beta.array()/certificate_scale;
    const Vector scaled_gradient=ac_gradient.array()/norms.array()/certificate_scale;
    Vector projected=scaled_beta-scaled_gradient;
    Certificate certificate; certificate.evaluated=true;
    bool feasible=true;
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {
        const Eigen::Index amplitude=2*atom;
        feasible &= beta(amplitude)>=0;
        projected(amplitude)=std::max(0.0,projected(amplitude));
        if(beta(amplitude)==0) certificate.active_atoms.push_back(atom);
    }
    certificate.projected_kkt=(scaled_beta-projected).lpNorm<Eigen::Infinity>();
    certificate.available=true; certificate.feasible=feasible;
    certificate.kkt_passed=feasible && certificate.projected_kkt<=1e-10;
    certificate.rss=out.residual.squaredNorm(); certificate.objective=.5*certificate.rss;
    certificate.residual_scale=certificate.rss/static_cast<double>(y.size());
    certificate.residual_rmse=out.residual.norm()/std::sqrt(static_cast<double>(y.size()));
    certificate.residual_max=out.residual.cwiseAbs().maxCoeff();
    certificate.relative_residual=out.residual.norm()/certificate_scale;
    CaptureWidthBoundKkt(certificate,eta,out.gradient,context);
    out.certificate=std::move(certificate);
    out.valid=out.residual.allFinite() && out.gradient.allFinite();
    out.reason=out.valid ? "raw-state" : "nonfinite-state";
    return out;
}
TrustEvidence CheckTrust(const Domain & domain,VectorRef y,const Evaluation & e,const EvaluationContext & policy)
{
    return CheckTrust(domain,y,e,policy,CheckIndependentReplay(domain,y,e,policy));
}
TrustEvidence CheckTrust(const Domain &,VectorRef,const Evaluation & e,const EvaluationContext &,
    const IndependentReplayEvidence & replay)
{
    TrustEvidence out; out.primary_valid=e.valid; out.independent=replay;
    out.prediction_passed=replay.prediction_passed && replay.residual_passed && replay.objective_passed;
    out.gradient_passed=replay.kkt_agrees && replay.gradient_agrees;
    out.kkt_difference=replay.kkt_difference;
    out.prediction_difference=replay.prediction_difference;
    out.gradient_difference=replay.gradient_difference;
    out.cancellation_ratio=replay.cancellation_ratio;
    out.passed=replay.passed;
    out.reason=replay.reason;
    out.width_bound_feasible=replay.width_bound_feasible;
    return out;
}
TrustEvidence CheckReplay(const Domain & domain,VectorRef y,const Evaluation & e,const EvaluationContext & policy,
    ReplayDetailsForDiagnostics * diagnostic)
{
    const auto * context=&policy;
    TrustEvidence out; out.primary_valid=e.valid;
    if (!e.valid)
    {
        out.reason="invalid-primary";
        return out;
    }
    if(!e.certificate.available || !e.certificate.feasible || !e.certificate.kkt_passed)
    {out.reason="kkt-failed"; return out;}
    const double scale=context->scale;
    const Vector widths=e.eta.array().exp();
    Vector prediction=Vector::Zero(y.size()),compensation=prediction,absolute=prediction;
    // Independent scalar forward and compensated summation, using frozen support.
    auto add=[&](Eigen::Index row,double value) {
        const double sum=prediction(row),next=sum+value;
        compensation(row)+=std::abs(sum)>=std::abs(value) ? (sum-next)+value : (value-next)+sum;
        prediction(row)=next; absolute(row)+=std::abs(value);
    };
    auto add_product=[&](Eigen::Index row,double coefficient,double kernel) {
        const double product=coefficient*kernel;
        const double product_error=std::fma(coefficient,kernel,-product);
        add(row,product);
        add(row,product_error); absolute(row)-=std::abs(product_error);
    };
    for (Eigen::Index a=0;a<e.eta.size();++a)
    {
        const double b=widths(a);
        for (const auto & p:domain.atoms[static_cast<std::size_t>(a)])
        {
            if(p.square>2.5*2.5) continue;
            const double r=std::sqrt(p.square),g=std::pow(2*M_PI*b*b,-1.5)*std::exp(-p.square/(2*b*b));
            const double k=r<1e-5 ? std::sqrt(2/M_PI)/b : std::erf(r/b/std::sqrt(2.0))/r;
            add_product(p.row,e.beta(2*a),g); add_product(p.row,e.beta(2*a+1),k);
        }
    }
    prediction+=compensation;
    const Vector replay_residual=prediction-y,original_prediction=StablePrediction(e.x,e.beta);
    Vector norms(e.x.cols());
    for (Eigen::Index k=0;k<e.x.cols();++k) norms(k)=e.x.col(k).norm();
    const Vector u=norms.array()*e.beta.array()/scale;
    const Vector gradient=StableTransposeProduct(e.x,replay_residual).array()/norms.array()/scale;
    Vector projected=u-gradient;
    for (Eigen::Index k=0;k<e.beta.size();k+=2) projected(k)=std::max(0.0,projected(k));
    const double kkt=(u-projected).lpNorm<Eigen::Infinity>();
    Vector width_gradient=Vector::Zero(e.eta.size());
    for (Eigen::Index a=0;a<e.eta.size();++a)
    {
        const double b=widths(a);
        for (const auto & p:domain.atoms[static_cast<std::size_t>(a)])
        {
            if(p.square>2.5*2.5) continue;
            const double g=std::pow(2*M_PI*b*b,-1.5)*std::exp(-p.square/(2*b*b));
            const double dk=-std::sqrt(2/M_PI)/b*(p.square<1e-10 ? 1 : std::exp(-p.square/(2*b*b)));
            width_gradient(a)+=(e.beta(2*a)*g*(p.square/(b*b)-3)+e.beta(2*a+1)*dk)*replay_residual(p.row)/scale/scale;
        }
    }
    if(context->eta_upper_bound)
    {
        const auto bound=CalculateWidthBoundKkt(e.eta,width_gradient,*context->eta_upper_bound);
        out.width_bound_feasible=bound.feasible;
    }
    const bool prediction_ok=((prediction-original_prediction).array().abs()<=
        2e-12+2e-13*original_prediction.array().abs()).all();
    const bool gradient_ok=((width_gradient-e.gradient).array().abs()<=1e-13+2e-9*e.gradient.array().abs()).all();
    const double kkt_difference=std::abs(kkt-e.certificate.projected_kkt);
    if(diagnostic)
    {
        diagnostic->prediction=prediction;
        diagnostic->residual=replay_residual;
        diagnostic->width_gradient=width_gradient;
        diagnostic->objective=.5*replay_residual.squaredNorm();
        diagnostic->scaled_kkt=kkt;
    }
    out.prediction_passed=prediction_ok;
    out.gradient_passed=gradient_ok;
    out.kkt_difference=kkt_difference;
    out.prediction_difference=(prediction-original_prediction).lpNorm<Eigen::Infinity>();
    out.gradient_difference=(width_gradient-e.gradient).lpNorm<Eigen::Infinity>();
    out.cancellation_ratio=(absolute.array()/prediction.array().abs().max(1.0)).maxCoeff();
    const bool independent_kkt_ok=kkt<=1e-10;
    out.passed=prediction_ok && gradient_ok && independent_kkt_ok && kkt_difference<=1e-13 &&
        out.width_bound_feasible;
    if(out.passed) out.reason="trusted";
    else if(!independent_kkt_ok) out.reason="independent-kkt-failed";
    else if(!out.width_bound_feasible) out.reason="width-bound-infeasible";
    else out.reason="replay-disagreement";
    return out;
}

IndependentReplayEvidence CheckIndependentReplay(const Domain & domain,VectorRef y,const Evaluation & e,
    const EvaluationContext & policy,const Vector * assembled_prediction,const std::vector<std::size_t> * assembled_rows,
    const JointConvergenceCriteria & criteria)
{
    constexpr double pi=std::numbers::pi;
    IndependentReplayEvidence out; out.primary_valid=e.valid;
    if(!IsValidJointConvergenceCriteria(criteria))
    {out.reason="invalid-convergence-criteria"; return out;}
    if(!e.valid || domain.rows!=y.size() || e.residual.size()!=y.size() ||
        e.eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) || e.beta.size()!=2*e.eta.size())
    {out.reason="invalid-primary"; return out;}
    if(!e.certificate.available)
    {out.reason="primary-kkt-unavailable"; return out;}
    if(assembled_prediction && (!assembled_rows || assembled_rows->size()!=static_cast<std::size_t>(y.size()) ||
        std::any_of(assembled_rows->begin(),assembled_rows->end(),[&](std::size_t row) {
            return row>=static_cast<std::size_t>(assembled_prediction->size());
        })))
    {out.reason="invalid-assembled-prediction-mapping"; return out;}
    out.available=true;

    const auto atoms=e.eta.size();
    const double scale=policy.scale>0 ? policy.scale : std::max(1.0,y.norm());
    const Vector widths=e.eta.array().exp();
    if(!widths.allFinite() || (widths.array()<=0).any())
    {out.reason="invalid-widths"; return out;}
    Vector prediction=Vector::Zero(y.size()),prediction_error=prediction,absolute_sum=prediction;
    Vector column_norms_squared=Vector::Zero(e.beta.size());
    const auto add=[&](Eigen::Index row,double value) {
        const double sum=prediction(row),next=sum+value;
        prediction_error(row)+=std::abs(sum)>=std::abs(value) ? (sum-next)+value : (value-next)+sum;
        prediction(row)=next; absolute_sum(row)+=std::abs(value);
    };
    const auto add_product=[&](Eigen::Index row,double coefficient,double kernel_value) {
        const double product=coefficient*kernel_value;
        const double product_error=std::fma(coefficient,kernel_value,-product);
        add(row,product); add(row,product_error); absolute_sum(row)-=std::abs(product_error);
    };
    const auto kernel=[&](double square,double width) {
        struct Values {double gaussian{},charge{},gaussian_log_width{},charge_log_width{};};
        Values value;
        if(square>2.5*2.5) return value;
        const double radius=std::sqrt(square),width_squared=width*width;
        const double exponential=std::exp(-square/(2*width_squared));
        value.gaussian=exponential/std::pow(2*pi*width_squared,1.5);
        value.gaussian_log_width=value.gaussian*(square/width_squared-3.0);
        const double center=std::sqrt(2/(pi*width_squared));
        if(radius<1e-5)
        {value.charge=center; value.charge_log_width=-center;}
        else
        {value.charge=std::erf(radius/width/std::sqrt(2.0))/radius; value.charge_log_width=-center*exponential;}
        return value;
    };
    for(Eigen::Index atom=0;atom<atoms;++atom)
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto basis=kernel(support.square,widths(atom));
            add_product(support.row,e.beta(2*atom),basis.gaussian);
            add_product(support.row,e.beta(2*atom+1),basis.charge);
            column_norms_squared(2*atom)+=basis.gaussian*basis.gaussian;
            column_norms_squared(2*atom+1)+=basis.charge*basis.charge;
        }
    prediction+=prediction_error;
    if(!prediction.allFinite() || !column_norms_squared.allFinite())
    {out.reason="nonfinite-replay"; return out;}

    double rss{},rss_compensation{},prediction_difference{},residual_difference{},cancellation_ratio{},
        assembled_rss{},assembled_rss_compensation{},assembled_prediction_difference{},assembled_prediction_error_ratio{};
    bool prediction_ok=true,residual_ok=true,assembled_prediction_ok=true;
    for(Eigen::Index row=0;row<y.size();++row)
    {
        const double independent_residual=prediction(row)-y(row),primary_prediction=e.residual(row)+y(row);
        const double prediction_error_value=std::abs(prediction(row)-primary_prediction);
        const double residual_error_value=std::abs(independent_residual-e.residual(row));
        prediction_difference=std::max(prediction_difference,prediction_error_value);
        residual_difference=std::max(residual_difference,residual_error_value);
        const double prediction_tolerance=2e-12+2e-13*std::abs(primary_prediction);
        prediction_ok &= prediction_error_value<=prediction_tolerance;
        residual_ok &= residual_error_value<=prediction_tolerance;
        if(assembled_prediction)
        {
            const double assembled_value=(*assembled_prediction)(static_cast<Eigen::Index>(assembled_rows->at(
                static_cast<std::size_t>(row))));
            const double assembled_error=std::abs(prediction(row)-assembled_value);
            const double assembled_tolerance=2e-12+2e-13*std::abs(assembled_value);
            assembled_prediction_difference=std::max(assembled_prediction_difference,assembled_error);
            assembled_prediction_error_ratio=std::max(assembled_prediction_error_ratio,assembled_error/assembled_tolerance);
            assembled_prediction_ok &= std::isfinite(assembled_value) && assembled_error<=assembled_tolerance;
            const double assembled_residual=assembled_value-y(row);
            const double assembled_square=assembled_residual*assembled_residual;
            const double assembled_adjusted=assembled_square-assembled_rss_compensation;
            const double assembled_next=assembled_rss+assembled_adjusted;
            assembled_rss_compensation=(assembled_next-assembled_rss)-assembled_adjusted;
            assembled_rss=assembled_next;
        }
        const double square=independent_residual*independent_residual;
        const double adjusted=square-rss_compensation,next=rss+adjusted;
        rss_compensation=(next-rss)-adjusted; rss=next;
        cancellation_ratio=std::max(cancellation_ratio,absolute_sum(row)/std::max(1.0,std::abs(prediction(row))));
    }

    Vector ac_gradient=Vector::Zero(e.beta.size()),ac_compensation=ac_gradient;
    Vector width_gradient=Vector::Zero(atoms),width_compensation=width_gradient;
    const auto add_compensated=[](Vector & sum,Vector & compensation,Eigen::Index index,double value) {
        const double adjusted=value-compensation(index),next=sum(index)+adjusted;
        compensation(index)=(next-sum(index))-adjusted; sum(index)=next;
    };
    for(Eigen::Index atom=0;atom<atoms;++atom)
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto basis=kernel(support.square,widths(atom));
            const double residual=prediction(support.row)-y(support.row);
            add_compensated(ac_gradient,ac_compensation,2*atom,basis.gaussian*residual);
            add_compensated(ac_gradient,ac_compensation,2*atom+1,basis.charge*residual);
            add_compensated(width_gradient,width_compensation,atom,
                (e.beta(2*atom)*basis.gaussian_log_width+e.beta(2*atom+1)*basis.charge_log_width)*
                    residual/scale/scale);
        }
    const Vector column_norms=column_norms_squared.array().sqrt();
    if(!ac_gradient.allFinite() || !width_gradient.allFinite() || !column_norms.allFinite() ||
        (column_norms.array()<=0).any())
    {out.available=true; out.reason="rank-deficient-column"; return out;}

    const Vector scaled_beta=column_norms.array()*e.beta.array()/scale;
    const Vector scaled_gradient=ac_gradient.array()/column_norms.array()/scale;
    Vector projected=scaled_beta-scaled_gradient;
    bool feasible=true;
    for(Eigen::Index atom=0;atom<atoms;++atom)
    {feasible &= e.beta(2*atom)>=0; projected(2*atom)=std::max(0.0,projected(2*atom));}
    out.available=true; out.valid=true; out.feasible=feasible;
    out.primary_ac_kkt=e.certificate.projected_kkt;
    out.independent_ac_kkt=(scaled_beta-projected).lpNorm<Eigen::Infinity>();
    out.primary_width_gradient=e.gradient.lpNorm<Eigen::Infinity>();
    out.independent_width_gradient=width_gradient.lpNorm<Eigen::Infinity>();
    if(policy.eta_upper_bound)
    {
        const auto primary_width=CalculateWidthBoundKkt(e.eta,e.gradient,*policy.eta_upper_bound);
        const auto independent_width=CalculateWidthBoundKkt(e.eta,width_gradient,*policy.eta_upper_bound);
        out.width_bound_applied=true;
        out.width_bound_feasible=primary_width.feasible && independent_width.feasible;
        out.primary_width_projected_kkt=primary_width.projected;
        out.independent_width_projected_kkt=independent_width.projected;
        out.primary_width_kkt_passed=primary_width.feasible &&
            MeetsJointThreshold(primary_width.projected,criteria.global_width_gradient);
        out.independent_width_kkt_passed=independent_width.feasible &&
            MeetsJointThreshold(independent_width.projected,criteria.global_width_gradient);
        out.maximum_width_bound_violation=std::max(primary_width.violation,independent_width.violation);
        out.boundary_gradient_max=independent_width.boundary_gradient;
        out.boundary_multiplier_max=independent_width.multiplier;
        out.complementarity_residual=independent_width.complementarity;
        out.active_width_bounds=independent_width.active;
    }
    out.prediction_difference=prediction_difference; out.residual_difference=residual_difference;
    out.independent_objective=.5*rss;
    out.primary_objective=e.certificate.objective;
    out.primary_residual_norm=e.residual.norm();
    out.independent_residual_norm=std::sqrt(rss);
    out.objective_difference=std::abs(out.independent_objective-e.certificate.objective);
    out.kkt_difference=std::abs(out.independent_ac_kkt-out.primary_ac_kkt);
    out.cancellation_ratio=cancellation_ratio;
    if(assembled_prediction)
    {
        out.assembled_prediction_available=true;
        out.assembled_prediction_difference=assembled_prediction_difference;
        out.assembled_prediction_error_ratio=assembled_prediction_error_ratio;
        out.assembled_prediction_passed=assembled_prediction_ok;
        out.assembled_objective_difference=std::abs(.5*assembled_rss-.5*rss);
        out.assembled_objective_passed=std::isfinite(out.assembled_objective_difference) &&
            out.assembled_objective_difference<=1e-12+2e-12*std::abs(.5*rss);
    }
    out.primary_kkt_passed=e.certificate.feasible &&
        MeetsJointThreshold(out.primary_ac_kkt,criteria.global_ac_kkt);
    out.independent_kkt_passed=feasible &&
        MeetsJointThreshold(out.independent_ac_kkt,criteria.global_ac_kkt);
    out.prediction_passed=prediction_ok;
    out.residual_passed=residual_ok;
    out.objective_passed=out.objective_difference<=1e-12+2e-12*std::abs(e.certificate.objective);
    out.kkt_agrees=out.kkt_difference<=1e-13;
    out.gradient_agrees=true;
    double maximum_gradient_difference{};
    for(Eigen::Index atom=0;atom<atoms;++atom)
    {
        const double primary=e.gradient(atom),independent=width_gradient(atom);
        const double difference=std::abs(independent-primary);
        const double enclosure=1e-13+2e-9*std::abs(primary);
        if(atom==0 || difference>maximum_gradient_difference)
        {
            maximum_gradient_difference=difference;
            out.maximum_gradient_difference_primary=primary;
            out.maximum_gradient_difference_independent=independent;
            out.maximum_gradient_difference_bound=enclosure;
            out.maximum_gradient_difference_ratio=difference/enclosure;
            if(static_cast<std::size_t>(atom)<policy.atom_ids.size())
                out.maximum_gradient_difference_atom_id=policy.atom_ids[static_cast<std::size_t>(atom)];
        }
        if(difference>enclosure)
        {
            out.gradient_agrees=false;
            const double ratio=difference/enclosure;
            if(out.gradient_enclosure_failure_atom_id.empty() || ratio>out.gradient_enclosure_failure_ratio)
            {
                out.gradient_enclosure_failure_primary=primary;
                out.gradient_enclosure_failure_independent=independent;
                out.gradient_enclosure_failure_difference=difference;
                out.gradient_enclosure_failure_bound=enclosure;
                out.gradient_enclosure_failure_ratio=ratio;
                if(static_cast<std::size_t>(atom)<policy.atom_ids.size())
                    out.gradient_enclosure_failure_atom_id=policy.atom_ids[static_cast<std::size_t>(atom)];
            }
        }
    }
    out.gradient_difference=maximum_gradient_difference;
    out.passed=out.primary_kkt_passed && out.independent_kkt_passed && out.prediction_passed &&
        out.residual_passed && out.objective_passed && out.kkt_agrees && out.gradient_agrees &&
        (!out.width_bound_applied || (out.width_bound_feasible && out.primary_width_kkt_passed &&
            out.independent_width_kkt_passed)) &&
        (!assembled_prediction || (out.assembled_prediction_passed && out.assembled_objective_passed));
    out.reason=out.passed ? "verified" : !out.primary_kkt_passed ? "primary-kkt-failed" :
        !out.independent_kkt_passed ? "independent-kkt-failed" : !out.prediction_passed ?
        "prediction-replay-disagreement" : !out.residual_passed ? "residual-replay-disagreement" :
        !out.objective_passed ? "objective-replay-disagreement" : !out.kkt_agrees ?
        "independent-kkt-disagreement" : !out.gradient_agrees ? "per-coordinate-replay-enclosure" :
        out.width_bound_applied && !out.width_bound_feasible ? "width-bound-infeasible" :
        out.width_bound_applied && !out.primary_width_kkt_passed ? "primary-bounded-width-kkt-failed" :
        out.width_bound_applied && !out.independent_width_kkt_passed ? "independent-bounded-width-kkt-failed" :
        assembled_prediction && !out.assembled_prediction_passed ? "assembled-prediction-disagreement" :
        assembled_prediction && !out.assembled_objective_passed ? "assembled-objective-disagreement" :
        "replay-disagreement";
    out.valid=std::isfinite(out.independent_ac_kkt) && std::isfinite(out.independent_width_gradient) &&
        std::isfinite(out.objective_difference);
    if(!out.valid) {out.passed=false; out.reason="nonfinite-replay-evidence";}
    return out;
}

namespace {
template<class Design> Spectrum SpectrumRecord(const Design & input,const RankPolicy & policy,Eigen::Index columns,bool normalize)
{
    Spectrum out; out.rows=policy.rows; out.columns=columns; out.full_design_columns=columns;
    if(input.cols()==0) {out.available=false; out.reason="empty-matrix"; return out;}
    if(input.rows()==0) {out.reason="unobserved-columns"; return out;}
    Design x=input; Vector norms(input.cols());
    for(Eigen::Index k=0;k<x.cols();++k) {norms(k)=x.col(k).norm(); if(normalize && norms(k)>0) x.col(k)/=norms(k);}
    const auto reduced=Reduce(x,Matrix(x.rows(),0));
    const auto svd=EvaluateRank(reduced.first,{policy,columns,-1,RankBoundary::StrictGreater}); const auto & values=svd.singular_values;
    if(!svd.valid) {out.available=false; out.reason="spectrum-factorization-failed"; return out;}
    const double threshold=svd.threshold;
    out.singular_values=values; out.rank=svd.rank; out.threshold=threshold;
    out.column_norms=norms; out.minimum=values(values.size()-1); out.condition=values(0)/values(values.size()-1); return out;
}
}
Spectrum ComputeSpectrum(const Sparse & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return SpectrumRecord(x,p,columns,normalize);}
Spectrum ComputeSpectrum(const Matrix & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return SpectrumRecord(x,p,columns,normalize);}
namespace {
Assessment AssessProfileImpl(const Domain & domain,VectorRef y,const Vector & eta,const EvaluationContext & policy,
    const Vector * supplied_beta,const Vector * assembled_prediction,const std::vector<std::size_t> * assembled_rows)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting total_stage("total-assessment",y.size(),eta.size());
#endif
    const auto verification_memory=EndpointVerificationMemoryEstimate(domain,y,policy);
    if(verification_memory.status!=MemoryStatus::Pass)
    {
        Assessment assessment; assessment.primary.eta=eta;
        assessment.primary.reason="resource-limit:endpoint-verification:"+verification_memory.reason;
        assessment.failure=assessment.primary.reason;
        assessment.verification_memory_estimate=verification_memory;
        return assessment;
    }
    if(verification_memory.estimated_peak_bytes>verification_memory.budget_bytes/2)
        LogMemoryEstimate(verification_memory);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting primary_stage("endpoint-primary-evaluation",y.size(),eta.size(),&verification_memory);
#endif
    const auto primary_started=std::chrono::steady_clock::now();
    const auto endpoint=supplied_beta ? EvaluateStateStreaming(domain,y,eta,*supplied_beta,policy) :
        EvaluateProfile(domain,y,eta,&policy);
    const double primary_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-primary_started).count();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    primary_stage.Finish();
    AssessmentStageTimerForTesting replay_stage("independent-streaming-replay",y.size(),eta.size(),&verification_memory);
#endif
    const auto replay_started=std::chrono::steady_clock::now();
    const auto replay=CheckIndependentReplay(domain,y,endpoint,policy,assembled_prediction,assembled_rows);
    const double replay_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-replay_started).count();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    replay_stage.Finish();
#endif
    auto assessment=VerifyJointEndpoint(domain,y,endpoint,replay,policy);
    assessment.verification_memory_estimate=verification_memory;
    assessment.primary_endpoint_seconds=primary_seconds;
    assessment.independent_replay_seconds=replay_seconds;
    return assessment;
}
}
Assessment AssessProfile(const Domain & domain,VectorRef y,const Vector & eta,const EvaluationContext & policy,
    const Vector * supplied_beta,const Vector * assembled_prediction,const std::vector<std::size_t> * assembled_rows)
{
    return AssessProfileImpl(domain,y,eta,policy,supplied_beta,assembled_prediction,assembled_rows);
}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
AssessmentWork & AssessmentWorkForTesting() {static thread_local AssessmentWork work; return work;}
namespace {
AssessmentStageObserverForTesting & AssessmentStageObserver()
{static thread_local AssessmentStageObserverForTesting observer{}; return observer;}
void * & AssessmentStageObserverContext()
{static thread_local void * context{}; return context;}
AssessmentStageWork & AssessmentStage(const std::string & name)
{
    auto & stages=AssessmentWorkForTesting().stages;
    const auto found=std::find_if(stages.begin(),stages.end(),[&](const auto & stage){return stage.name==name;});
    if(found!=stages.end()) return *found;
    AssessmentStageWork stage; stage.name=name; stages.push_back(std::move(stage)); return stages.back();
}
void NotifyAssessmentStageObserver()
{
    if(const auto observer=AssessmentStageObserver()) observer(AssessmentWorkForTesting(),AssessmentStageObserverContext());
}
}
void SetAssessmentStageObserverForTesting(AssessmentStageObserverForTesting observer,void * context)
{AssessmentStageObserver()=observer; AssessmentStageObserverContext()=context;}
void BeginAssessmentStageForTesting(const std::string & name,Eigen::Index rows,Eigen::Index columns,
    const MemoryEstimate * estimate)
{
    auto & work=AssessmentWorkForTesting(); auto & stage=AssessmentStage(name);
    ++stage.calls; stage.rows=rows; stage.columns=columns;
    if(estimate)
    {
        stage.memory_estimate_available=true; stage.memory_stage=estimate->stage; stage.memory_status=estimate->status;
        stage.estimated_resident_bytes=estimate->resident_bytes; stage.estimated_temporary_bytes=estimate->temporary_bytes;
        stage.estimated_peak_bytes=estimate->estimated_peak_bytes; stage.budget_bytes=estimate->budget_bytes;
        stage.memory_preflight_reason=estimate->reason;
        stage.sparse_design_bytes=estimate->sparse_design_bytes; stage.dense_factor_bytes=estimate->dense_factor_bytes;
        stage.factor_workspace_bytes=estimate->factor_workspace_bytes;
        stage.factor_rows=estimate->factor_rows; stage.factor_columns=estimate->factor_columns;
        stage.maximum_tile_rows=estimate->shape.maximum_tile_rows;
    }
    stage.current_rss_before_bytes=CurrentRssBytes();
    stage.process_peak_rss_bytes=ProcessPeakRssBytes();
    work.stage_stack.push_back(name); work.active_stage=name; work.last_stage=name;
    NotifyAssessmentStageObserver();
}
void FinishAssessmentStageForTesting(const std::string & name,double seconds,bool completed)
{
    auto & work=AssessmentWorkForTesting(); auto & stage=AssessmentStage(name);
    stage.seconds+=seconds; if(completed) ++stage.completed_calls;
    stage.current_rss_after_bytes=CurrentRssBytes();
    const auto peak=ProcessPeakRssBytes();
    if(peak && (!stage.process_peak_rss_bytes || *peak>*stage.process_peak_rss_bytes)) stage.process_peak_rss_bytes=peak;
    if(!work.stage_stack.empty()) work.stage_stack.pop_back();
    work.active_stage=work.stage_stack.empty() ? std::string{} : work.stage_stack.back();
    work.last_stage=name;
    NotifyAssessmentStageObserver();
}
AssessmentStageTimerForTesting::AssessmentStageTimerForTesting(std::string name,Eigen::Index rows,Eigen::Index columns,
    const MemoryEstimate * estimate)
    :name_(std::move(name)),uncaught_exceptions_(std::uncaught_exceptions())
{
    BeginAssessmentStageForTesting(name_,rows,columns,estimate);
    started_=std::chrono::steady_clock::now();
}
AssessmentStageTimerForTesting::~AssessmentStageTimerForTesting()
{
    if(finished_) return;
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
    FinishAssessmentStageForTesting(name_,seconds,std::uncaught_exceptions()==uncaught_exceptions_);
}
void AssessmentStageTimerForTesting::Finish()
{
    if(finished_) return;
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
    FinishAssessmentStageForTesting(name_,seconds,true); finished_=true;
}
#endif
Assessment VerifyJointEndpoint(const Domain & domain,VectorRef y,const Evaluation & endpoint,
    const IndependentReplayEvidence & replay,const EvaluationContext & context,
    const JointConvergenceCriteria & criteria)
{
    Assessment out; out.primary=endpoint; out.independent_replay=replay; out.convergence_criteria=criteria;
    if(!IsValidJointConvergenceCriteria(criteria))
    {out.failure="invalid-convergence-criteria"; return out;}
    const auto atom_count=static_cast<Eigen::Index>(domain.atoms.size());
    bool valid=domain.rows==y.size() && y.allFinite() && endpoint.valid &&
        endpoint.eta.size()==atom_count && endpoint.beta.size()==2*atom_count &&
        endpoint.gradient.size()==atom_count && endpoint.eta.allFinite() && endpoint.beta.allFinite() &&
        endpoint.gradient.allFinite();
    if(valid)
        for(Eigen::Index atom=0;atom<atom_count;++atom)
        {
            const double width=std::exp(endpoint.eta(atom));
            valid &= endpoint.beta(2*atom)>=0 && std::isfinite(width) && width>0;
        }
    if(!valid)
    {
        out.primary.valid=false;
        out.primary.reason="invalid-accepted-parameter-state";
        out.failure=out.primary.reason;
        return out;
    }
    if(!replay.available || !replay.valid)
    {
        out.failure=replay.reason.empty() ? "endpoint-replay-unavailable" : replay.reason;
        return out;
    }
    out.inner=replay.feasible && replay.primary_kkt_passed &&
        MeetsJointThreshold(replay.primary_ac_kkt,criteria.global_ac_kkt) &&
        replay.independent_kkt_passed && MeetsJointThreshold(replay.independent_ac_kkt,criteria.global_ac_kkt);
    out.gradient=MeetsJointThreshold(endpoint.gradient.lpNorm<Eigen::Infinity>(),criteria.global_width_gradient) &&
        MeetsJointThreshold(replay.primary_width_gradient,criteria.global_width_gradient) &&
        MeetsJointThreshold(replay.independent_width_gradient,criteria.global_width_gradient) && replay.gradient_agrees;
    if(context.eta_upper_bound)
        out.gradient=replay.width_bound_applied && replay.width_bound_feasible && replay.primary_width_kkt_passed &&
            replay.independent_width_kkt_passed &&
            MeetsJointThreshold(replay.primary_width_projected_kkt,criteria.global_width_gradient) &&
            MeetsJointThreshold(replay.independent_width_projected_kkt,criteria.global_width_gradient) && replay.gradient_agrees;
    const bool replay_consistent=replay.valid && replay.prediction_passed && replay.residual_passed &&
        replay.objective_passed && replay.kkt_agrees && replay.gradient_agrees &&
        (!replay.assembled_prediction_available ||
            (replay.assembled_prediction_passed && replay.assembled_objective_passed));
    out.endpoint_verified=replay.passed && replay_consistent && out.inner;
    out.numerical_qualified=out.endpoint_verified && out.gradient;
    out.failure=out.numerical_qualified ? "none" : !out.inner ? "global-ac-kkt-failed" :
        !replay_consistent ? (replay.reason.empty() ? "endpoint-verification-failed" : replay.reason) :
        !out.gradient ? "global-width-gradient-failed" :
        replay.reason.empty() ? "endpoint-verification-failed" : replay.reason;
    return out;
}
Assessment AssessEvaluated(const Domain & domain,VectorRef y,const Evaluation & endpoint,
    const IndependentReplayEvidence & replay,const EvaluationContext & policy,bool,std::size_t)
{return VerifyJointEndpoint(domain,y,endpoint,replay,policy);}

}
