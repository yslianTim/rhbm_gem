#include "Numerics.hpp"
#include "SparseFactor.hpp"
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
    destination.linear_symbolic_seconds+=source.linear_symbolic_seconds;
    destination.linear_numeric_seconds+=source.linear_numeric_seconds;
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
    return total.basis_seconds+total.linear_matrix_preparation_seconds+total.linear_symbolic_seconds+
        total.linear_numeric_seconds+total.linear_rhs_solve_seconds+total.linear_certificate_seconds+
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
{return {ActiveSparseBackend(),policy.core_atoms};}
std::string_view SparseBackendName(SparseBackend backend)
{
    switch(backend)
    {
    case SparseBackend::Eigen: return "EIGEN";
    case SparseBackend::Spqr: return "SPQR";
    }
    return {};
}
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
Evaluation EvaluateProfile(const Domain & domain,VectorRef y,const Vector & eta,bool reference,const EvaluationContext * context,
    const std::vector<LinearBlock> * blocks,LinearWorkspace * workspace,const void * workspace_identity,
    ProfileEvaluationRole role,ProfileSearchWork * profile_work)
{
    ResourcePhase phase(reference ? "reference" : "ac-profile");
    ProfileEvaluationRoleScopeForTesting role_scope(role);
    EvaluationTelemetry telemetry(profile_work,role);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const auto search_stage=ResourceWorkForTesting().active_search_stage;
    FactorCreationRoleScopeForTesting factor_role(reference ? "reference" :
        search_stage=="trial-evaluation" ? "trial" : "profile");
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(reference) ++AssessmentWorkForTesting().reference_evaluations;
#endif
    if(workspace) workspace->Bind(&domain,workspace_identity ? workspace_identity : y.data(),context ? &context->linear : nullptr);
    const auto matrix_started=std::chrono::steady_clock::now();
    const auto basis_started=std::chrono::steady_clock::now();
    auto out=Basis(domain,y,eta);
    if(profile_work) telemetry.work.basis_seconds+=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-basis_started).count();
    SparseWorkForTesting().matrix_preparation_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-matrix_started).count();
    if(!out.valid) return out; out.valid=false;
    const Eigen::Index m=eta.size();
    ResourcePhase linear_solve("linear-solve",true,out.x.rows(),out.x.cols(),static_cast<std::size_t>(out.x.nonZeros()));
    const auto sparse_before=SparseWorkForTesting();
    const auto solved=SolveLinear(out.x,y,Vector::Ones(y.size()),reference,true,nullptr,context ? &context->linear : nullptr,blocks,workspace);
    const auto & sparse_after=SparseWorkForTesting();
    if(profile_work)
    {
        telemetry.work.linear_matrix_preparation_seconds+=
            sparse_after.matrix_preparation_seconds-sparse_before.matrix_preparation_seconds;
        telemetry.work.linear_symbolic_seconds+=sparse_after.symbolic_seconds-sparse_before.symbolic_seconds;
        telemetry.work.linear_numeric_seconds+=sparse_after.numeric_seconds-sparse_before.numeric_seconds;
        telemetry.work.linear_rhs_solve_seconds+=(sparse_after.q_seconds-sparse_before.q_seconds)+
            (sparse_after.triangular_seconds-sparse_before.triangular_seconds)+
            (sparse_after.least_squares_seconds-sparse_before.least_squares_seconds);
    }
    if(!reference) out.factor=solved.factor;
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
TrustEvidence CheckTrust(const Domain & domain,VectorRef y,const Evaluation & e,const EvaluationContext & policy)
{
    const auto reference=EvaluateProfile(domain,y,e.eta,true,&policy);
    return CheckTrust(domain,y,e,policy,reference);
}
TrustEvidence CheckTrust(const Domain & domain,VectorRef y,const Evaluation & e,const EvaluationContext & policy,const Evaluation & reference)
{
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
    Spectrum out; out.rank=svd.rank; out.minimum=values.tail(1)(0);
    out.condition=values(0)/values.tail(1)(0); out.singular_values=values; out.threshold=svd.threshold; return out;
}
Assessment AssessProfile(const Domain & domain,VectorRef y,const Vector & eta,const EvaluationContext & policy,const Vector * supplied_beta)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting total_stage("total-assessment",y.size(),eta.size());
    AssessmentStageTimerForTesting primary_stage("endpoint-primary-evaluation",y.size(),eta.size());
#endif
    const auto endpoint=supplied_beta ? EvaluateState(domain,y,eta,*supplied_beta,policy) : EvaluateProfile(domain,y,eta,false,&policy);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    primary_stage.Finish();
    AssessmentStageTimerForTesting reference_stage("reference-evaluation",y.size(),eta.size());
#endif
    const auto reference=EvaluateProfile(domain,y,eta,true,&policy);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    reference_stage.Finish();
#endif
    return AssessEvaluated(domain,y,endpoint,reference,policy,supplied_beta!=nullptr);
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
    stages.push_back({name}); return stages.back();
}
void NotifyAssessmentStageObserver()
{
    if(const auto observer=AssessmentStageObserver()) observer(AssessmentWorkForTesting(),AssessmentStageObserverContext());
}
}
void SetAssessmentStageObserverForTesting(AssessmentStageObserverForTesting observer,void * context)
{AssessmentStageObserver()=observer; AssessmentStageObserverContext()=context;}
void BeginAssessmentStageForTesting(const std::string & name,Eigen::Index rows,Eigen::Index columns)
{
    auto & work=AssessmentWorkForTesting(); auto & stage=AssessmentStage(name);
    ++stage.calls; stage.rows=rows; stage.columns=columns;
    work.stage_stack.push_back(name); work.active_stage=name; work.last_stage=name;
    NotifyAssessmentStageObserver();
}
void FinishAssessmentStageForTesting(const std::string & name,double seconds,bool completed)
{
    auto & work=AssessmentWorkForTesting(); auto & stage=AssessmentStage(name);
    stage.seconds+=seconds; if(completed) ++stage.completed_calls;
    if(!work.stage_stack.empty()) work.stage_stack.pop_back();
    work.active_stage=work.stage_stack.empty() ? std::string{} : work.stage_stack.back();
    work.last_stage=name;
    NotifyAssessmentStageObserver();
}
AssessmentStageTimerForTesting::AssessmentStageTimerForTesting(std::string name,Eigen::Index rows,Eigen::Index columns)
    :name_(std::move(name)),uncaught_exceptions_(std::uncaught_exceptions())
{
    BeginAssessmentStageForTesting(name_,rows,columns);
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
Assessment AssessEvaluated(const Domain &,VectorRef y,const Evaluation & endpoint,const Evaluation & reference,
    const EvaluationContext & policy,bool supplied)
{
    ResourcePhase phase("assessment");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ++AssessmentWorkForTesting().assessments;
#endif
    const auto * context=&policy; const double scale=context->scale; const auto & eta=endpoint.eta; Assessment out;
    out.primary=endpoint; out.reference=reference;
    if(supplied && endpoint.valid && !endpoint.certificate.kkt_passed)
    {out.primary.valid=false; out.primary.reason="assembled-kkt-failed";}
    if(!out.primary.valid || !reference.valid) {out.failure="inner-solve"; return out;}
    const double difference=Difference(endpoint.beta,reference.beta); out.coefficient_difference=difference;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting design_stage("design-spectrum",endpoint.x.rows(),endpoint.x.cols());
#endif
    out.design=DesignSpectrum(endpoint.x,Vector::Ones(y.size()),&policy.rank);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    design_stage.Finish();
#endif
    out.inner=difference<=1e-10 && out.design->rank==endpoint.x.cols();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting prepare_stage("derivative-preparation",endpoint.derivative.rows(),endpoint.derivative.cols());
    const auto reduction_kind=JacobianReductionForTesting();
    const bool compact_route=reduction_kind==JacobianReductionKindForTesting::CompactStackQr;
    if(compact_route) ++AssessmentWorkForTesting().compact_attempts;
#else
    constexpr bool compact_route=true;
#endif
    const auto prepared=compact_route ?
        PrepareDerivativeCompact(endpoint,scale,context) : PrepareDerivative(endpoint,scale,context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    prepare_stage.Finish();
    AssessmentStageTimerForTesting reduce_stage("derivative-reduction",prepared.raw.rows(),prepared.raw.cols());
#endif
    bool compact_differential=compact_route;
    ReducedDifferential differential;
    if(compact_differential && prepared.reference_order)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++AssessmentWorkForTesting().compact_other_fallbacks;
#endif
        differential=ReduceDerivative(prepared,endpoint.residual);
        compact_differential=false;
    }
    else differential=compact_differential ? ReduceDerivativeCompact(prepared,endpoint.residual) :
        ReduceDerivative(prepared,endpoint.residual);
    if(compact_differential && !differential.valid)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++AssessmentWorkForTesting().compact_other_fallbacks;
#endif
        differential=ReduceDerivative(prepared,endpoint.residual);
        compact_differential=false;
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    reduce_stage.Finish();
#endif
    if(!differential.valid) {out.failure=differential.reason; return out;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    CompactSvdResult widths,normalized_svd;
    const auto width_spectra=[&]() {
        {
            AssessmentStageTimerForTesting width_stage("projected-width-spectrum",differential.projected.rows(),differential.projected.cols());
            widths=Decompose(differential.projected,context->rank.rows,CompactSvdVectors::Right);
        }
        if(!widths.valid) return false;
        {
            AssessmentStageTimerForTesting normalized_stage("normalized-width-spectrum",differential.projected.rows(),differential.projected.cols());
            Matrix normalized=differential.projected;
            for(Eigen::Index k=0;k<differential.projected_norms.size();++k)
                if(differential.projected_norms(k)>0) normalized.col(k)/=differential.projected_norms(k);
            normalized_svd=Decompose(normalized,context->rank.rows);
        }
        return normalized_svd.valid;
    };
#else
    const auto widths=Decompose(differential.projected,context->rank.rows,CompactSvdVectors::Right);
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(!width_spectra()) {out.failure="spectrum-factorization-failed"; return out;}
    if(differential.projected_candidate &&
        (!CompactRankDecisionSafe(widths) || !CompactRankDecisionSafe(normalized_svd)))
    {
        auto & projected_work=DerivativeWorkForTesting().projected_reduction;
        ++projected_work.fallbacks;
        projected_work.fallback_reason="rank-decision-boundary";
        const auto saved=ProjectedReductionForTesting();
        const auto candidate_kind=projected_work.kind;
        ProjectedReductionForTesting()=ProjectedReductionKindForTesting::ObservationTiledQr;
        differential=ReduceDerivativeCompact(prepared,endpoint.residual);
        ProjectedReductionForTesting()=saved;
        projected_work.kind=candidate_kind;
        if(!differential.valid || !width_spectra()) {out.failure="spectrum-factorization-failed"; return out;}
    }
    else if(differential.projected_candidate)
        ++DerivativeWorkForTesting().projected_reduction.accepted;
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
    out.weak_directions.resize(eta.size(),std::min<Eigen::Index>(3,widths.right_vectors.cols()));
    for(Eigen::Index k=0;k<out.weak_directions.cols();++k)
    {
        Vector direction=widths.right_vectors.col(widths.right_vectors.cols()-1-k);
        Eigen::Index sign{}; direction.cwiseAbs().maxCoeff(&sign); if(direction(sign)<0) direction=-direction;
        out.weak_directions.col(k)=direction;
    }
    Vector response=-differential.response;
    CompactSvdResult correction_svd;
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting correction_stage("correction-jacobian-spectrum",differential.jacobian.rows(),differential.jacobian.cols());
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
        differential=ReduceDerivative(prepared,endpoint.residual);
        if(!differential.valid) {out.failure=differential.reason; return out;}
        response=-differential.response;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting fallback_stage("correction-jacobian-spectrum",differential.jacobian.rows(),differential.jacobian.cols());
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
    if(!correction_svd.valid) {out.failure="spectrum-factorization-failed"; return out;}
    const Vector & correction=correction_svd.solution; out.correction=correction;
    out.jacobian=CompactSpectrum(correction_svd,context->rank.rows);
    out.gradient=endpoint.gradient.lpNorm<Eigen::Infinity>()<=1e-12 && reference.gradient.lpNorm<Eigen::Infinity>()<=1e-12;
    out.local=correction.allFinite() && correction.lpNorm<Eigen::Infinity>()<=1e-10;
    out.identified=widths.rank==eta.size() && correction_svd.rank==eta.size();
    out.failure=!out.inner ? "inner-solve" : !out.identified ? "width-unidentified" :
        !out.gradient || !out.local ? "b-not-stationary" : "none";
    return out;
}

}
