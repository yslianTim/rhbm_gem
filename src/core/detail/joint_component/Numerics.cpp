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
void AddRoleWork(ProfileRoleWork & destination,const ProfileRoleWork & source)
{
    destination.evaluations+=source.evaluations;
    destination.derivative_preparations+=source.derivative_preparations;
    destination.derivative_reductions+=source.derivative_reductions;
    destination.replay_checks+=source.replay_checks;
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
const char * ProfileEvaluationRoleName(ProfileEvaluationRole role)
{
    switch(role)
    {
    case ProfileEvaluationRole::Unspecified: return "unspecified";
    case ProfileEvaluationRole::InitialProfile: return "initial-profile";
    case ProfileEvaluationRole::TrialProfile: return "trial-profile";
    case ProfileEvaluationRole::AcceptedEndpoint: return "accepted-endpoint";
    case ProfileEvaluationRole::Reference: return "reference";
    }
    return "unspecified";
}
ProfileRoleWork & ProfileSearchWork::Role(ProfileEvaluationRole role)
{
    switch(role)
    {
    case ProfileEvaluationRole::InitialProfile: return initial_profile;
    case ProfileEvaluationRole::TrialProfile: return trial_profile;
    case ProfileEvaluationRole::AcceptedEndpoint: return accepted_endpoint;
    case ProfileEvaluationRole::Reference: return reference_evaluation;
    case ProfileEvaluationRole::Unspecified: return total;
    }
    return total;
}
void ProfileSearchWork::Add(ProfileEvaluationRole role,const ProfileRoleWork & work)
{
    if(role!=ProfileEvaluationRole::Unspecified) AddRoleWork(Role(role),work);
    AddRoleWork(total,work);
}
void ProfileSearchWork::Merge(const ProfileSearchWork & other)
{
    AddRoleWork(total,other.total); AddRoleWork(initial_profile,other.initial_profile);
    AddRoleWork(trial_profile,other.trial_profile); AddRoleWork(accepted_endpoint,other.accepted_endpoint);
    AddRoleWork(reference_evaluation,other.reference_evaluation);
    total_seconds+=other.total_seconds; lm_overhead_seconds+=other.lm_overhead_seconds;
}
double ProfileSearchWork::AttributedSeconds() const
{
    return total.basis_seconds+total.linear_matrix_preparation_seconds+total.linear_rhs_solve_seconds+
        total.linear_certificate_seconds+
        total.derivative_prepare_seconds+total.derivative_reduce_seconds+total.replay_trust_seconds;
}
namespace {
double Difference(const Vector & a,const Vector & b) {return ((a-b).array().abs()/(1+a.array().abs().max(b.array().abs()))).maxCoeff();}
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
CompactSvdResult Decompose(const Matrix & r,Eigen::Index rows,CompactSvdVectors vectors=CompactSvdVectors::None,const Vector * rhs=nullptr)
{
    return EvaluateRank(r,{{rows,0,0},r.cols()},rhs,vectors);
}
}
JointSolverConfiguration ResolveJointSolverConfiguration(const FixedNeighborSearchPolicy & policy)
{return {policy.core_atoms};}
std::string_view FixedNeighborBlockOrderName(FixedNeighborBlockOrder order)
{
    switch(order)
    {
    case FixedNeighborBlockOrder::Forward: return "forward";
    case FixedNeighborBlockOrder::Reverse: return "reverse";
    }
    return {};
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
    if(input->selection_domain && !CheckedMultiplyBytes(input->selection_domain->target_indices.size(),3,c.target_columns))
        c.target_columns=std::numeric_limits<std::size_t>::max();
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
    const double r{std::sqrt(square)}, r2{r*r}, b2{width*width}, exponent{std::exp(-r2/(2.0*b2))};
    BasisValues out;
    out.gaussian = std::pow(2.0*std::numbers::pi*b2,-1.5)*exponent;
    out.gaussian_log_width = out.gaussian*(r2/b2-3.0);
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
    const Eigen::VectorXd residual=x*beta-y;
    if (!residual.allFinite()) return out;
    const double s=observation_scale>0 ? observation_scale : std::max(1.0,y.norm());
    const Eigen::VectorXd u=norms.array()*beta.array()/s;
    const Eigen::VectorXd gradient=(x.transpose()*residual).array()/norms.array()/s;
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
MemoryEstimate ReferenceSolveMemoryEstimate(const Domain & domain,VectorRef,const EvaluationContext & context,
    std::size_t resident_override)
{
    MemoryShape shape; shape.component_atoms=domain.atoms.size(); shape.component_observations=static_cast<std::size_t>(domain.rows);
    shape.component_target_columns=context.target_columns;
    const auto overflow=[&] {
        auto estimate=PlanMemory(shape,MemoryStage::ReferenceSolve);
        estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
    };
    const auto add_product=[](std::size_t & total,std::size_t count,std::size_t width) {
        std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
    };
    if(!CheckedMultiplyBytes(shape.component_atoms,2,shape.component_free_columns)) return overflow();
    for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        if(!CheckedAddBytes(shape.component_support_entries,domain.atoms[atom].size(),shape.component_support_entries))
            return overflow();
    const auto * input=domain.atoms.Storage();
    const auto snapshot=input ? EstimateSnapshotBytes(*input) : std::nullopt;
    std::size_t view_bytes{};
    if(!snapshot || !CheckedAddBytes(*snapshot,domain.atoms.LocalSupportIndexBytes(),view_bytes) ||
        !add_product(view_bytes,shape.component_observations,2*sizeof(double)) ||
        !add_product(view_bytes,shape.component_observations,2*sizeof(Eigen::Index)) ||
        !add_product(view_bytes,shape.component_atoms,2*sizeof(Eigen::Index)+4*sizeof(double)) ||
        !add_product(view_bytes,shape.component_support_entries,4*(sizeof(double)+sizeof(Eigen::Index)))) return overflow();
    shape.resident_bytes=std::max(view_bytes,resident_override);
    return PlanMemory(shape,MemoryStage::ReferenceSolve);
}
MemoryEstimate EndpointVerificationMemoryEstimate(const Domain & domain,VectorRef y,const EvaluationContext & context,
    std::size_t resident_override)
{
    MemoryShape shape; shape.component_atoms=domain.atoms.size();
    shape.component_observations=static_cast<std::size_t>(domain.rows);
    shape.component_target_columns=context.target_columns;
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
Evaluation EvaluateProfile(const Domain & domain,VectorRef y,const Vector & eta,bool reference,const EvaluationContext * context,
    const std::vector<LinearBlock> * blocks,
    ProfileEvaluationRole role,ProfileSearchWork * profile_work)
{
    ResourcePhase phase(reference ? "reference" : "ac-profile");
    ProfileEvaluationRoleScopeForTesting role_scope(role);
    EvaluationTelemetry telemetry(profile_work,role);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(reference) ++AssessmentWorkForTesting().reference_evaluations;
#endif
    if(reference)
    {
        const EvaluationContext fallback_context;
        const auto estimate=ReferenceSolveMemoryEstimate(domain,y,context ? *context : fallback_context);
        if(estimate.status!=MemoryStatus::Pass || estimate.estimated_peak_bytes>estimate.budget_bytes/2)
            LogMemoryEstimate(estimate);
        if(estimate.status!=MemoryStatus::Pass)
        {Evaluation unavailable_evaluation; unavailable_evaluation.eta=eta;
            unavailable_evaluation.reason="resource-limit:reference-solve"; return unavailable_evaluation;}
    }
    const auto basis_started=std::chrono::steady_clock::now();
    auto out=Basis(domain,y,eta);
    if(profile_work) telemetry.work.basis_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-basis_started).count();
    if(!out.valid) return out; out.valid=false;
    const Eigen::Index m=eta.size();
    ResourcePhase linear_solve("linear-solve",true,out.x.rows(),out.x.cols(),static_cast<std::size_t>(out.x.nonZeros()));
    const auto solve_started=std::chrono::steady_clock::now();
    const auto solved=SolveLinear(out.x,y,Vector::Ones(y.size()),reference,true,nullptr,context ? &context->linear : nullptr,blocks);
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
    out.certificate.linear_solves=solved.solves; out.certificate.free_rank=solved.rank;
    if(blocks) out.certificate.block_factorizations=solved.block_factorizations;
    if (!solved.valid) {out.reason=solved.reason; return out;}
    if (!out.certificate.kkt_passed) {out.reason="kkt-failed"; return out;}
    out.residual=out.x*out.beta-y; out.gradient=Vector::Zero(m);
    const double scale=context ? context->scale : std::max(1.0,y.norm());
    for (Eigen::Index k=0;k<2*m;++k)
        out.gradient(k/2)+=out.beta(k)*out.derivative.col(k).dot(out.residual)/scale/scale;
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
    out.beta=beta; out.residual=out.x*beta-y; out.gradient=Vector::Zero(eta.size());
    const auto certificate_started=std::chrono::steady_clock::now();
    out.certificate=CertifyLinear(out.x,y,beta,context.scale);
    if(profile_work) telemetry.work.linear_certificate_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-certificate_started).count();
    for(Eigen::Index k=0;k<beta.size();++k)
        out.gradient(k/2)+=beta(k)*out.derivative.col(k).dot(out.residual)/context.scale/context.scale;
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
    return out;
}
TrustEvidence CheckTrust(const Domain & domain,VectorRef y,const Evaluation & e,const EvaluationContext & policy,const Evaluation & reference)
{
    if(!reference.valid && reference.reason.starts_with("resource-limit:reference-solve"))
    {TrustEvidence unavailable_trust; unavailable_trust.primary_valid=e.valid;
        unavailable_trust.reason=reference.reason; return unavailable_trust;}
    auto out=CheckReplay(domain,y,e,policy); out.reference=reference;
    if(!out.primary_valid) return out;
    out.coefficient_difference=reference.beta.size()==e.beta.size() ? Difference(e.beta,reference.beta) :
        std::numeric_limits<double>::infinity();
    const bool reference_gradient_ok=reference.valid && reference.gradient.size()==e.gradient.size() &&
        ((reference.gradient-e.gradient).array().abs()<=1e-13+2e-9*e.gradient.array().abs()).all();
    out.gradient_passed &= reference_gradient_ok;
    out.passed=out.passed && reference.valid && reference.certificate.kkt_passed &&
        out.coefficient_difference<=1e-10 && reference_gradient_ok;
    if(!reference.valid) out.reason="invalid-reference";
    else if(out.reason=="trusted" && !out.passed) out.reason="reference-disagreement";
    return out;
}
TrustEvidence CheckReplay(const Domain & domain,VectorRef y,const Evaluation & e,const EvaluationContext & policy)
{
    const auto * context=&policy;
    TrustEvidence out; out.primary_valid=e.valid;
    if (!e.valid)
    {
        out.reason="invalid-primary";
        if(e.x.cols()>0 && e.x.nonZeros()>0 && context->audit.trial_details)
            out.design=DesignSpectrum(e.x,Vector::Ones(y.size()),&policy.rank);
        return out;
    }
    if(!e.certificate.available || !e.certificate.feasible || !e.certificate.kkt_passed)
    {out.reason="kkt-failed"; return out;}
    const double scale=context->scale;
    Vector prediction=Vector::Zero(y.size()),compensation=prediction,absolute=prediction;
    // Independent scalar forward and compensated summation, using frozen support.
    auto add=[&](Eigen::Index row,double value) {
        const double increment=value-compensation(row),sum=prediction(row)+increment;
        compensation(row)=(sum-prediction(row))-increment; prediction(row)=sum; absolute(row)+=std::abs(value);
    };
    for (Eigen::Index a=0;a<e.eta.size();++a)
    {
        const double b=std::exp(e.eta(a));
        for (const auto & p:domain.atoms[static_cast<std::size_t>(a)])
        {
            const double r=std::sqrt(p.square),g=std::pow(2*M_PI*b*b,-1.5)*std::exp(-p.square/(2*b*b));
            const double k=r<1e-5 ? std::sqrt(2/M_PI)/b : std::erf(r/b/std::sqrt(2.0))/r;
            add(p.row,e.beta(2*a)*g); add(p.row,e.beta(2*a+1)*k);
        }
    }
    const Vector replay_residual=prediction-y,original_prediction=e.x*e.beta;
    Vector norms(e.x.cols());
    for (Eigen::Index k=0;k<e.x.cols();++k) norms(k)=e.x.col(k).norm();
    const Vector u=norms.array()*e.beta.array()/scale;
    const Vector gradient=(e.x.transpose()*replay_residual).array()/norms.array()/scale;
    Vector projected=u-gradient;
    for (Eigen::Index k=0;k<e.beta.size();k+=2) projected(k)=std::max(0.0,projected(k));
    const double kkt=(u-projected).lpNorm<Eigen::Infinity>();
    Vector width_gradient=Vector::Zero(e.eta.size());
    for (Eigen::Index a=0;a<e.eta.size();++a)
    {
        const double b=std::exp(e.eta(a));
        for (const auto & p:domain.atoms[static_cast<std::size_t>(a)])
        {
            const double g=std::pow(2*M_PI*b*b,-1.5)*std::exp(-p.square/(2*b*b));
            const double dk=-std::sqrt(2/M_PI)/b*(p.square<1e-10 ? 1 : std::exp(-p.square/(2*b*b)));
            width_gradient(a)+=(e.beta(2*a)*g*(p.square/(b*b)-3)+e.beta(2*a+1)*dk)*replay_residual(p.row)/scale/scale;
        }
    }
    const bool prediction_ok=((prediction-original_prediction).array().abs()<=
        2e-12+2e-13*original_prediction.array().abs()).all();
    const bool gradient_ok=((width_gradient-e.gradient).array().abs()<=1e-13+2e-9*e.gradient.array().abs()).all();
    const double kkt_difference=std::abs(kkt-e.certificate.projected_kkt);
    out.prediction_passed=prediction_ok;
    out.gradient_passed=gradient_ok;
    out.kkt_difference=kkt_difference;
    out.prediction_difference=(prediction-original_prediction).lpNorm<Eigen::Infinity>();
    out.gradient_difference=(width_gradient-e.gradient).lpNorm<Eigen::Infinity>();
    out.cancellation_ratio=(absolute.array()/prediction.array().abs().max(1.0)).maxCoeff();
    if (context->audit.trial_details) out.design=DesignSpectrum(e.x,Vector::Ones(y.size()),&policy.rank);
    out.passed=prediction_ok && gradient_ok && kkt_difference<=1e-13;
    out.reason=out.passed ? "trusted" : "replay-disagreement";
    return out;
}

IndependentReplayEvidence CheckIndependentReplay(const Domain & domain,VectorRef y,const Evaluation & e,
    const EvaluationContext & policy,const Vector * assembled_prediction,const std::vector<std::size_t> * assembled_rows)
{
    constexpr double pi=std::numbers::pi;
    IndependentReplayEvidence out; out.primary_valid=e.valid;
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
        const double adjusted=value-prediction_error(row);
        const double next=prediction(row)+adjusted;
        prediction_error(row)=(next-prediction(row))-adjusted;
        prediction(row)=next; absolute_sum(row)+=std::abs(value);
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
            add(support.row,e.beta(2*atom)*basis.gaussian);
            add(support.row,e.beta(2*atom+1)*basis.charge);
            column_norms_squared(2*atom)+=basis.gaussian*basis.gaussian;
            column_norms_squared(2*atom+1)+=basis.charge*basis.charge;
        }
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
    out.prediction_difference=prediction_difference; out.residual_difference=residual_difference;
    out.independent_objective=.5*rss;
    out.objective_difference=std::abs(out.independent_objective-e.certificate.objective);
    out.kkt_difference=std::abs(out.independent_ac_kkt-out.primary_ac_kkt);
    out.gradient_difference=(width_gradient-e.gradient).lpNorm<Eigen::Infinity>();
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
    out.primary_kkt_passed=e.certificate.feasible && e.certificate.kkt_passed && out.primary_ac_kkt<=1e-10;
    out.independent_kkt_passed=feasible && out.independent_ac_kkt<=1e-10;
    out.prediction_passed=prediction_ok;
    out.residual_passed=residual_ok;
    out.objective_passed=out.objective_difference<=1e-12+2e-12*std::abs(e.certificate.objective);
    out.kkt_agrees=out.kkt_difference<=1e-13;
    out.gradient_agrees=((width_gradient-e.gradient).array().abs()<=
        1e-13+2e-9*e.gradient.array().abs()).all();
    out.passed=out.primary_kkt_passed && out.independent_kkt_passed && out.prediction_passed &&
        out.residual_passed && out.objective_passed && out.kkt_agrees && out.gradient_agrees &&
        (!assembled_prediction || (out.assembled_prediction_passed && out.assembled_objective_passed));
    out.reason=out.passed ? "verified" : !out.primary_kkt_passed ? "primary-kkt-failed" :
        !out.independent_kkt_passed ? "independent-kkt-failed" :
        assembled_prediction && !out.assembled_prediction_passed ? "assembled-prediction-disagreement" :
        assembled_prediction && !out.assembled_objective_passed ? "assembled-objective-disagreement" :
        "replay-disagreement";
    out.valid=std::isfinite(out.independent_ac_kkt) && std::isfinite(out.independent_width_gradient) &&
        std::isfinite(out.objective_difference);
    if(!out.valid) {out.passed=false; out.reason="nonfinite-replay-evidence";}
    return out;
}

namespace {
Spectrum CompactSpectrum(const CompactSvdResult & svd,Eigen::Index rows)
{
    const auto & v=svd.singular_values; Spectrum out;
    if(!svd.valid) {out.available=false; out.reason="spectrum-factorization-failed"; return out;}
    out.rank=svd.rank; out.rows=rows; out.columns=v.size(); out.singular_values=v;
    out.minimum=v(v.size()-1); out.condition=v(0)/v(v.size()-1);
    out.threshold=svd.threshold; return out;
}
template<class Design> Spectrum SpectrumRecord(const Design & input,const RankPolicy & policy,Eigen::Index columns,bool normalize)
{
    Spectrum out; out.rows=policy.rows; out.columns=columns;
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
Spectrum DesignSpectrum(const Sparse & x,const Vector & weights,const RankPolicy * policy)
{
    Vector scales(x.cols());
    for(Eigen::Index k=0;k<x.cols();++k)
    {
        double sum{}; for(Sparse::InnerIterator e(x,k);e;++e) sum+=weights(e.row())*e.value()*e.value();
        scales(k)=std::sqrt(sum);
    }
    for(Eigen::Index k=0;k<x.cols();++k) if(scales(k)==0) scales(k)=1;
    const auto reduced=ReferenceQR(x,weights,scales,Vector::Zero(x.rows()));
    const auto svd=EvaluateRank(reduced.first,{{policy ? policy->rows : x.rows(),x.cols(),0},x.cols(),-1,RankBoundary::StrictGreater}); const auto & values=svd.singular_values;
    if(!svd.valid) {Spectrum out; out.available=false; out.reason="spectrum-factorization-failed"; return out;}
    Spectrum out; out.rows=x.rows(); out.columns=x.cols(); out.rank=svd.rank; out.minimum=values.tail(1)(0);
    out.condition=values(0)/values.tail(1)(0); out.singular_values=values; out.threshold=svd.threshold; return out;
}
Assessment AssessProfile(const Domain & domain,VectorRef y,const Vector & eta,const EvaluationContext & policy,
    const Vector * supplied_beta,const Vector * assembled_prediction,const std::vector<std::size_t> * assembled_rows)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting total_stage("total-assessment",y.size(),eta.size());
#endif
    const auto verification_memory=EndpointVerificationMemoryEstimate(domain,y,policy);
    if(verification_memory.status!=MemoryStatus::Pass)
    {
        Evaluation primary; primary.eta=eta;
        primary.reason="resource-limit:endpoint-verification:"+verification_memory.reason;
        IndependentReplayEvidence replay; replay.reason=primary.reason;
        auto assessment=AssessEvaluated(domain,y,primary,replay,policy,supplied_beta!=nullptr);
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
        EvaluateProfile(domain,y,eta,false,&policy);
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
    auto assessment=AssessEvaluated(domain,y,endpoint,replay,policy,supplied_beta!=nullptr);
    assessment.verification_memory_estimate=verification_memory;
    assessment.primary_endpoint_seconds=primary_seconds;
    assessment.independent_replay_seconds=replay_seconds;
    return assessment;
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
        stage.right_singular_vector_bytes=estimate->right_singular_vector_bytes;
        stage.right_null_space_bytes=estimate->right_null_space_bytes;
        stage.factor_rows=estimate->factor_rows; stage.factor_columns=estimate->factor_columns;
        stage.maximum_tile_rows=estimate->shape.maximum_tile_rows;
        stage.attributed_memory_stages=estimate->attributed_stages;
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
bool SameAssessmentPolicy(const EvaluationContext & a,const EvaluationContext & b)
{
    return a.scale==b.scale && a.rank.rows==b.rank.rows && a.rank.design_columns==b.rank.design_columns &&
        a.rank.width_columns==b.rank.width_columns && a.linear.rank_relative==b.linear.rank_relative &&
        a.linear.active_set_iteration_factor==b.linear.active_set_iteration_factor && a.linear.release_factor==b.linear.release_factor &&
        a.linear.release_response_norm==b.linear.release_response_norm && a.atom_ids==b.atom_ids && a.row_ids==b.row_ids;
}
namespace {
MemoryEstimate RankAssessmentMemoryEstimate(const Domain & domain,VectorRef,const Evaluation & endpoint,
    const Evaluation * reference,const EvaluationContext & context,std::size_t resident_override)
{
    MemoryShape shape; shape.component_atoms=domain.atoms.size(); shape.component_observations=static_cast<std::size_t>(domain.rows);
    shape.component_target_columns=context.target_columns;
    if(context.target_columns%3==0 && context.target_columns/3<=shape.component_atoms)
    {
        shape.component_target_atoms=context.target_columns/3;
        shape.component_halo_atoms=shape.component_atoms-shape.component_target_atoms;
    }
    shape.component_width_columns=shape.component_atoms;
    shape.maximum_tile_rows=8192;
    if(!CheckedMultiplyBytes(shape.component_atoms,2,shape.component_design_columns) ||
        !CheckedMultiplyBytes(shape.component_atoms,3,shape.component_jacobian_columns))
    {
        auto estimate=PlanMemory(shape,MemoryStage::RankAssessment);
        estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
    }
    const auto overflow=[&] {
        auto estimate=PlanMemory(shape,MemoryStage::RankAssessment);
        estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
    };
    for(Eigen::Index k=0;k<endpoint.beta.size();++k)
        if(k%2 || endpoint.beta(k)>0)
            if(!CheckedAddBytes(shape.component_free_columns,1,shape.component_free_columns)) return overflow();
    for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        if(!CheckedAddBytes(shape.component_support_entries,domain.atoms[atom].size(),shape.component_support_entries))
            return overflow();
    const auto * input=domain.atoms.Storage();
    const auto snapshot=input ? EstimateSnapshotBytes(*input) : std::nullopt;
    std::size_t resident{};
    const auto add_product=[](std::size_t & total,std::size_t count,std::size_t width) {
        std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
    };
    const auto add_sparse=[&](const Sparse & matrix) {
        std::size_t outer{};
        return matrix.outerSize()>=0 && matrix.nonZeros()>=0 &&
            CheckedAddBytes(static_cast<std::size_t>(matrix.outerSize()),1,outer) &&
            add_product(resident,static_cast<std::size_t>(matrix.nonZeros()),sizeof(double)+sizeof(Eigen::Index)) &&
            add_product(resident,outer,sizeof(Eigen::Index));
    };
    if(!snapshot || !CheckedAddBytes(*snapshot,domain.atoms.LocalSupportIndexBytes(),resident) ||
        !add_product(resident,static_cast<std::size_t>(std::max<Eigen::Index>(0,domain.rows)),3*sizeof(double)) ||
        !add_product(resident,shape.component_atoms,16*sizeof(double)+4*sizeof(Eigen::Index)) ||
        !add_product(resident,shape.component_observations,2*sizeof(Eigen::Index)) ||
        !add_product(resident,shape.component_atoms,2*sizeof(Eigen::Index)) ||
        !add_sparse(endpoint.x) || !add_sparse(endpoint.derivative) ||
        (reference && (!add_sparse(reference->x) || !add_sparse(reference->derivative)))) return overflow();
    shape.resident_bytes=std::max(resident,resident_override);
    return PlanMemory(shape,MemoryStage::RankAssessment);
}
}
namespace {
Assessment AssessEvaluatedImpl(const Domain & domain,VectorRef y,const Evaluation & endpoint,const Evaluation * reference,
    const IndependentReplayEvidence * replay,const EvaluationContext & policy,bool supplied,std::size_t resident_override)
{
    ResourcePhase phase("assessment");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ++AssessmentWorkForTesting().assessments;
#endif
    const auto * context=&policy; const double scale=context->scale; const auto & eta=endpoint.eta; Assessment out;
    out.primary=endpoint;
    out.reference_free=reference==nullptr;
    if(reference) out.reference=*reference;
    if(replay) out.independent_replay=*replay;
    if(supplied && endpoint.valid && !endpoint.certificate.kkt_passed && !replay)
    {out.primary.valid=false; out.primary.reason="assembled-kkt-failed";}
    if(!out.primary.valid || (reference && !reference->valid) ||
        (replay && (!replay->available || !replay->valid)))
    {
        out.failure=reference && !reference->valid ? reference->reason : replay && !replay->valid ? replay->reason : "inner-solve";
        return out;
    }
    if(reference) out.coefficient_difference=Difference(endpoint.beta,reference->beta);
    out.endpoint_verified=replay && replay->passed && replay->primary_kkt_passed && replay->independent_kkt_passed &&
        replay->primary_width_gradient<=1e-12;
    struct RankAssessmentTimer
    {
        double & value;
        std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
        ~RankAssessmentTimer()
        {value=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
    } rank_timer{out.rank_assessment_seconds};
    out.rank_memory_estimate=RankAssessmentMemoryEstimate(domain,y,endpoint,reference,policy,resident_override);
    const auto design_memory=PlanMemory(out.rank_memory_estimate.shape,MemoryStage::DesignSpectrum,
        out.rank_memory_estimate.budget_bytes);
    const auto projected_memory=PlanMemory(out.rank_memory_estimate.shape,MemoryStage::ProjectedWidthAssessment,
        out.rank_memory_estimate.budget_bytes);
    const auto correction_memory=PlanMemory(out.rank_memory_estimate.shape,MemoryStage::LocalCorrection,
        out.rank_memory_estimate.budget_bytes);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting identifiability_stage("numerical-identifiability",y.size(),
        static_cast<Eigen::Index>(out.rank_memory_estimate.shape.component_design_columns),&out.rank_memory_estimate);
#endif
    if(out.rank_memory_estimate.status!=MemoryStatus::Pass ||
        out.rank_memory_estimate.estimated_peak_bytes>out.rank_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.rank_memory_estimate);
    if(out.rank_memory_estimate.status!=MemoryStatus::Pass)
    {out.failure="resource-limit:rank-assessment"; return out;}
    Evaluation rebuilt;
    const Evaluation * evaluated=&endpoint;
    if(endpoint.x.cols()==0 || endpoint.derivative.cols()==0)
    {
        rebuilt=EvaluateState(domain,y,endpoint.eta,endpoint.beta,policy);
        if(!rebuilt.valid) {out.failure=rebuilt.reason; return out;}
        evaluated=&rebuilt;
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting design_stage("design-spectrum",evaluated->x.rows(),evaluated->x.cols(),&design_memory);
#endif
    out.design=DesignSpectrum(evaluated->x,Vector::Ones(y.size()),&policy.rank);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    design_stage.Finish();
#endif
    out.inner=reference ? out.coefficient_difference<=1e-10 && out.design->rank==evaluated->x.cols() :
        replay->primary_kkt_passed && replay->independent_kkt_passed;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting projected_overall("projected-width-assessment",evaluated->x.rows(),eta.size(),
        &projected_memory);
    AssessmentStageTimerForTesting prepare_stage("derivative-preparation",evaluated->derivative.rows(),
        evaluated->derivative.cols(),&projected_memory);
    ++AssessmentWorkForTesting().compact_attempts;
#endif
    const auto prepared=PrepareDerivativeCompact(*evaluated,scale,context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    prepare_stage.Finish();
    AssessmentStageTimerForTesting reduce_stage("derivative-reduction",prepared.raw.rows(),prepared.raw.cols(),
        &projected_memory);
#endif
    bool compact_differential=!prepared.reference_order;
    ReducedDifferential differential;
    if(!compact_differential)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++AssessmentWorkForTesting().compact_other_fallbacks;
#endif
        differential=ReduceDerivative(prepared,evaluated->residual);
    }
    else
    {
        differential=ReduceDerivativeCompact(prepared,evaluated->residual);
        if(!differential.valid)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++AssessmentWorkForTesting().compact_other_fallbacks;
#endif
            differential=ReduceDerivative(prepared,evaluated->residual);
            compact_differential=false;
        }
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    reduce_stage.Finish();
#endif
    if(!differential.valid) {out.failure=differential.reason; return out;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    CompactSvdResult widths,normalized_svd;
    const auto width_spectra=[&]() {
        {
            AssessmentStageTimerForTesting width_stage("projected-width-spectrum",differential.projected.rows(),
                differential.projected.cols(),&projected_memory);
            widths=Decompose(differential.projected,context->rank.rows,CompactSvdVectors::Right);
        }
        if(!widths.valid) return false;
        {
            AssessmentStageTimerForTesting normalized_stage("normalized-width-spectrum",differential.projected.rows(),
                differential.projected.cols(),&projected_memory);
            Matrix normalized=differential.projected;
            for(Eigen::Index k=0;k<differential.projected_norms.size();++k)
                if(differential.projected_norms(k)>0) normalized.col(k)/=differential.projected_norms(k);
            normalized_svd=Decompose(normalized,context->rank.rows);
        }
        return normalized_svd.valid;
    };
#else
    auto widths=Decompose(differential.projected,context->rank.rows,CompactSvdVectors::Right);
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(!width_spectra()) {out.failure="spectrum-factorization-failed"; return out;}
    out.widths=CompactSpectrum(widths,context->rank.rows);
    out.widths->column_norms=differential.projected_norms;
    out.normalized_widths=CompactSpectrum(normalized_svd,context->rank.rows);
#else
    if(!widths.valid) {out.failure="spectrum-factorization-failed"; return out;}
    out.widths=CompactSpectrum(widths,context->rank.rows);
    out.widths->column_norms=differential.projected_norms;
    Matrix normalized=differential.projected;
    for(Eigen::Index k=0;k<differential.projected_norms.size();++k)
        if(differential.projected_norms(k)>0) normalized.col(k)/=differential.projected_norms(k);
    const auto normalized_svd=Decompose(normalized,context->rank.rows);
    if(!normalized_svd.valid) {out.failure="spectrum-factorization-failed"; return out;}
    out.normalized_widths=CompactSpectrum(normalized_svd,context->rank.rows);
#endif
    const auto projected_width_rank=widths.rank;
    out.weak_directions.resize(eta.size(),std::min<Eigen::Index>(3,widths.right_vectors.cols()));
    for(Eigen::Index k=0;k<out.weak_directions.cols();++k)
    {
        Vector direction=widths.right_vectors.col(widths.right_vectors.cols()-1-k);
        Eigen::Index sign{}; direction.cwiseAbs().maxCoeff(&sign); if(direction(sign)<0) direction=-direction;
        out.weak_directions.col(k)=direction;
    }
    Matrix{}.swap(widths.right_vectors);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    projected_overall.Finish();
#endif
    Vector response=-differential.response;
    CompactSvdResult correction_svd;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting correction_overall("local-correction",differential.jacobian.rows(),
        differential.jacobian.cols(),&correction_memory);
#endif
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting correction_stage("correction-jacobian-spectrum",
            differential.jacobian.rows(),differential.jacobian.cols(),&correction_memory);
#endif
        correction_svd=Decompose(differential.jacobian,context->rank.rows,CompactSvdVectors::None,&response);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        correction_stage.Finish();
#endif
    }
    if(compact_differential && !CompactRankDecisionSafe(correction_svd))
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++AssessmentWorkForTesting().compact_boundary_fallbacks;
#endif
        differential=ReduceDerivative(prepared,evaluated->residual);
        if(!differential.valid) {out.failure=differential.reason; return out;}
        response=-differential.response;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting fallback_stage("correction-jacobian-spectrum",
            differential.jacobian.rows(),differential.jacobian.cols(),&correction_memory);
#endif
        correction_svd=Decompose(differential.jacobian,context->rank.rows,CompactSvdVectors::None,&response);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        fallback_stage.Finish();
#endif
    }
    else if(compact_differential)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++AssessmentWorkForTesting().compact_accepted;
#endif
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    correction_overall.Finish();
#endif
    if(!correction_svd.valid) {out.failure="spectrum-factorization-failed"; return out;}
    const Vector & correction=correction_svd.solution; out.correction=correction;
    out.jacobian=CompactSpectrum(correction_svd,context->rank.rows);
    out.gradient=endpoint.gradient.lpNorm<Eigen::Infinity>()<=1e-12 &&
        (reference ? reference->gradient.lpNorm<Eigen::Infinity>()<=1e-12 : replay->gradient_agrees);
    out.local=correction.allFinite() && correction.lpNorm<Eigen::Infinity>()<=1e-10;
    out.identified=projected_width_rank==eta.size() && correction_svd.rank==eta.size();
    const bool design_full_rank=out.design->rank==evaluated->x.cols();
    out.qualified=out.endpoint_verified && design_full_rank && out.identified && out.gradient && out.local;
    out.failure=reference && !out.inner ? "inner-solve" : !design_full_rank ? "design-unidentified" :
        !out.identified ? "width-unidentified" : !out.gradient || !out.local ? "b-not-stationary" :
        !out.endpoint_verified ? "endpoint-verification-failed" : "none";
    return out;
}
}
Assessment AssessEvaluated(const Domain & domain,VectorRef y,const Evaluation & endpoint,const Evaluation & reference,
    const EvaluationContext & policy,bool supplied)
{return AssessEvaluatedImpl(domain,y,endpoint,&reference,nullptr,policy,supplied,0);}
Assessment AssessEvaluated(const Domain & domain,VectorRef y,const Evaluation & endpoint,
    const IndependentReplayEvidence & replay,const EvaluationContext & policy,bool supplied,std::size_t resident_override)
{return AssessEvaluatedImpl(domain,y,endpoint,nullptr,&replay,policy,supplied,resident_override);}

}
