#include "Numerics.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
double Difference(const Vector & a,const Vector & b) {return ((a-b).array().abs()/(1+a.array().abs().max(b.array().abs()))).maxCoeff();}
// Full-column TSQR from original rows. RHS columns undergo the same orthogonal
// transformations. Only the small R and transformed RHS survive each tile.
template<class Design>
std::pair<Matrix,Matrix> Reduce(const Design & x,const Matrix & rhs)
{
    Matrix r(0,x.cols()),target(0,rhs.cols());
    constexpr Eigen::Index tile=8192;
    for (Eigen::Index first=0;first<x.rows();first+=tile)
    {
        const Eigen::Index n=std::min(tile,x.rows()-first),prior=r.rows();
        Matrix a(prior+n,x.cols()),b(prior+n,rhs.cols());
        a.topRows(prior)=r; b.topRows(prior)=target;
        a.bottomRows(n)=Matrix(x.middleRows(first,n)); b.bottomRows(n)=rhs.middleRows(first,n);
        const Eigen::HouseholderQR<Matrix> qr(a);
        const Matrix transformed=qr.householderQ().adjoint()*b;
        const Eigen::Index keep=std::min(a.rows(),a.cols());
        r=qr.matrixQR().topRows(keep).template triangularView<Eigen::Upper>(); target=transformed.topRows(keep);
    }
    return {std::move(r),std::move(target)};
}
Eigen::JacobiSVD<Matrix> Decompose(const Matrix & r,Eigen::Index rows)
{
    Eigen::JacobiSVD<Matrix> out(r,Eigen::ComputeThinU|Eigen::ComputeThinV);
    out.setThreshold(eps*static_cast<double>(std::max(rows,r.cols()))); return out;
}
}
double RankPolicy::Relative(Eigen::Index columns) const
{return std::numeric_limits<double>::epsilon()*static_cast<double>(std::max(rows,columns));}
double RankPolicy::Absolute(Eigen::Index columns,double maximum) const
{return Relative(columns)*maximum;}
EvaluationContext CreateContext(const Vector & y,Eigen::Index atoms,const std::string & hash,const AuditPlan & plan)
{
    EvaluationContext c; c.snapshot_hash=hash; c.observations=std::make_shared<const Eigen::VectorXd>(y);
    c.scale=std::max(1.0,y.norm()); c.rank={y.size(),2*atoms,atoms};
    c.linear.rank_relative=c.rank.Relative(2*atoms);
    c.linear.release_response_norm=y.norm();
    c.audit=plan;
    for(Eigen::Index a=0;a<atoms;++a) c.atom_ids.push_back(std::to_string(a));
    for(Eigen::Index r=0;r<y.size();++r) c.row_ids.push_back(std::to_string(r));
    return c;
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
Certificate CertifyLinear(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & beta,double observation_scale)
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

namespace {
Evaluation Basis(const Domain & domain,const Vector & y,const Vector & eta)
{
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
    out.valid=true; return out;
}
}
Evaluation EvaluateProfile(const Domain & domain,const Vector & y,const Vector & eta,bool reference,const EvaluationContext * context,
    const std::vector<LinearBlock> * blocks)
{
    auto out=Basis(domain,y,eta); if(!out.valid) return out; out.valid=false;
    const Eigen::Index m=eta.size();
    const auto solved=SolveLinear(out.x,y,Vector::Ones(y.size()),reference,true,nullptr,context ? &context->linear : nullptr,blocks);
    out.beta=solved.beta; out.certificate=CertifyLinear(out.x,y,out.beta,context ? context->scale : 0);
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

Evaluation EvaluateState(const Domain & domain,const Vector & y,const Vector & eta,const Vector & beta,const EvaluationContext & context)
{
    auto out=Basis(domain,y,eta); if(!out.valid) return out; out.valid=false;
    if(beta.size()!=out.x.cols() || !beta.allFinite()) {out.reason="invalid-coefficients"; return out;}
    out.beta=beta; out.residual=out.x*beta-y; out.gradient=Vector::Zero(eta.size());
    out.certificate=CertifyLinear(out.x,y,beta,context.scale);
    for(Eigen::Index k=0;k<beta.size();++k)
        out.gradient(k/2)+=beta(k)*out.derivative.col(k).dot(out.residual)/context.scale/context.scale;
    out.valid=out.residual.allFinite() && out.gradient.allFinite(); out.reason=out.valid ? "raw-state" : "nonfinite-state";
    return out;
}
Vector ComputeLocalCorrection(const Evaluation & e,const Differential & d,const EvaluationContext & context,double absolute)
{
    if(!e.valid || !d.valid) return {};
    const auto reduced=Reduce(d.jacobian,Matrix(-e.residual/context.scale));
    auto svd=Decompose(reduced.first,context.rank.rows);
    if(absolute>=0 && svd.singularValues()(0)>0) svd.setThreshold(absolute/svd.singularValues()(0));
    if(svd.rank()!=d.jacobian.cols()) return {};
    return svd.solve(reduced.second.col(0));
}

Differential DifferentiateProfile(const Evaluation & e,double scale,const EvaluationContext * context,double free_design_threshold)
{
    Differential out;
    if (!e.valid || !(scale>0) || !std::isfinite(scale)) {out.reason="invalid-inner"; return out;}
    std::vector<Eigen::Index> free;
    for (Eigen::Index k=0;k<e.beta.size();++k) if (k%2 || e.beta(k)>0) free.push_back(k);
    const Eigen::Index n=e.x.rows(),m=e.beta.size()/2,p=static_cast<Eigen::Index>(free.size());
    Sparse z(n,p); std::vector<Eigen::Triplet<double>> entries;
    Matrix t=Matrix::Zero(p,m),raw=Matrix::Zero(n,m);
    for (Eigen::Index k=0;k<2*m;++k) for (Sparse::InnerIterator entry(e.derivative,k);entry;++entry)
        raw(entry.row(),k/2)+=entry.value()*e.beta(k);
    for (Eigen::Index col=0;col<p;++col)
    {
        const auto k=free[static_cast<std::size_t>(col)]; const double norm=e.x.col(k).norm();
        if (!(norm>0)) {out.reason="zero-free-column"; return out;}
        for (Sparse::InnerIterator entry(e.x,k);entry;++entry) entries.emplace_back(entry.row(),col,entry.value()/norm);
        t(col,k/2)=e.derivative.col(k).dot(e.residual)/norm;
    }
    z.setFromTriplets(entries.begin(),entries.end());
    const auto reduced=Reduce(z,raw); auto svd=Decompose(reduced.first,context ? context->rank.rows : n);
    if(free_design_threshold>=0 && svd.singularValues()(0)>0) svd.setThreshold(free_design_threshold/svd.singularValues()(0));
    if (svd.rank()!=p) {out.reason="rank-deficient-free-design"; return out;}
    const Matrix coefficients=reduced.first.triangularView<Eigen::Upper>().solve(reduced.second);
    // R^{-1} R^{-T} T via orthogonal factors; never form X^T X.
    const Matrix adjoint=reduced.first.transpose().triangularView<Eigen::Lower>().solve(t);
    const Matrix correction=reduced.first.triangularView<Eigen::Upper>().solve(adjoint);
    out.projected=(raw-z*coefficients)/scale;
    out.jacobian=out.projected-z*correction/scale;
    out.valid=out.jacobian.allFinite() && out.projected.allFinite();
    out.reason=out.valid ? "full-profile-derivative" : "nonfinite-derivative"; return out;
}

TrustEvidence CheckTrust(const Domain & domain,const Vector & y,const Evaluation & e,const EvaluationContext & policy)
{
    const auto reference=EvaluateProfile(domain,y,e.eta,true,&policy);
    return CheckTrust(domain,y,e,policy,reference);
}
TrustEvidence CheckTrust(const Domain & domain,const Vector & y,const Evaluation & e,const EvaluationContext & policy,const Evaluation & reference)
{
    const auto * context=&policy;
    TrustEvidence out; out.reference=reference; out.primary_valid=e.valid;
    if (!e.valid)
    {
        out.reason="invalid-primary";
        if(e.x.cols()>0 && e.x.nonZeros()>0 && context->audit.trial_details)
            out.design=DesignSpectrum(e.x,Vector::Ones(y.size()));
        return out;
    }
    const double difference=reference.beta.size()==e.beta.size() ? Difference(e.beta,reference.beta) :
        std::numeric_limits<double>::infinity(),scale=context ? context->scale : std::max(1.0,y.norm());
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
    const bool reference_gradient_ok=reference.valid && ((reference.gradient-e.gradient).array().abs()<=1e-13+2e-9*e.gradient.array().abs()).all();
    const double kkt_difference=std::abs(kkt-e.certificate.projected_kkt);
    out.coefficient_difference=difference; out.prediction_passed=prediction_ok;
    out.gradient_passed=gradient_ok && reference_gradient_ok;
    out.kkt_difference=kkt_difference;
    out.prediction_difference=(prediction-original_prediction).lpNorm<Eigen::Infinity>();
    out.gradient_difference=(width_gradient-e.gradient).lpNorm<Eigen::Infinity>();
    out.cancellation_ratio=(absolute.array()/prediction.array().abs().max(1.0)).maxCoeff();
    if (context->audit.trial_details) out.design=DesignSpectrum(e.x,Vector::Ones(y.size()));
    out.passed=difference<=1e-10 && prediction_ok && gradient_ok && reference_gradient_ok && kkt_difference<=1e-13;
    out.reason=out.passed ? "trusted" : !reference.valid ? "invalid-reference" : "replay-disagreement";
    return out;
}

namespace {
Spectrum CompactSpectrum(const Eigen::JacobiSVD<Matrix> & svd,Eigen::Index rows)
{
    const auto & v=svd.singularValues(); Spectrum out;
    out.rank=svd.rank(); out.rows=rows; out.columns=svd.cols(); out.singular_values=v;
    out.minimum=v(v.size()-1); out.condition=v(0)/v(v.size()-1);
    out.threshold=eps*static_cast<double>(std::max(rows,svd.cols()))*v(0); return out;
}
template<class Design> Spectrum SpectrumRecord(const Design & input,const RankPolicy & policy,Eigen::Index columns,bool normalize)
{
    Spectrum out; out.rows=policy.rows; out.columns=columns;
    if(input.cols()==0) {out.available=false; out.reason="empty-matrix"; return out;}
    if(input.rows()==0) {out.reason="unobserved-columns"; return out;}
    Design x=input; Vector norms(input.cols());
    for(Eigen::Index k=0;k<x.cols();++k) {norms(k)=x.col(k).norm(); if(normalize && norms(k)>0) x.col(k)/=norms(k);}
    const auto reduced=Reduce(x,Matrix(x.rows(),0));
    const Eigen::JacobiSVD<Matrix> svd(reduced.first,Eigen::ComputeThinV); const auto & values=svd.singularValues();
    const double threshold=policy.Absolute(columns,values(0));
    out.singular_values=values; out.rank=(values.array()>threshold).count(); out.threshold=threshold;
    out.column_norms=norms; out.minimum=values(values.size()-1); out.condition=values(0)/values(values.size()-1); return out;
}
}
Spectrum ComputeSpectrum(const Sparse & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return SpectrumRecord(x,p,columns,normalize);}
Spectrum ComputeSpectrum(const Matrix & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return SpectrumRecord(x,p,columns,normalize);}
Spectrum DesignSpectrum(const Sparse & x,const Vector & weights)
{
    Vector scales(x.cols());
    for(Eigen::Index k=0;k<x.cols();++k)
    {
        double sum{}; for(Sparse::InnerIterator e(x,k);e;++e) sum+=weights(e.row())*e.value()*e.value();
        scales(k)=std::sqrt(sum);
    }
    for(Eigen::Index k=0;k<x.cols();++k) if(scales(k)==0) scales(k)=1;
    const auto reduced=ReferenceQR(x,weights,scales,Vector::Zero(x.rows()));
    const Eigen::JacobiSVD<Matrix> svd(reduced.first); const auto values=svd.singularValues();
    const double threshold=eps*static_cast<double>(std::max(x.rows(),x.cols()));
    Spectrum out; out.rank=(values.array()>threshold*values(0)).count(); out.minimum=values.tail(1)(0);
    out.condition=values(0)/values.tail(1)(0); out.singular_values=values; out.threshold=threshold*values(0); return out;
}
Assessment AssessProfile(const Domain & domain,const Vector & y,const Vector & eta,const EvaluationContext & policy,const Vector * supplied_beta)
{
    const auto endpoint=supplied_beta ? EvaluateState(domain,y,eta,*supplied_beta,policy) : EvaluateProfile(domain,y,eta,false,&policy);
    const auto reference=EvaluateProfile(domain,y,eta,true,&policy);
    return AssessEvaluated(domain,y,endpoint,reference,policy,supplied_beta!=nullptr);
}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
AssessmentWork & AssessmentWorkForTesting() {static thread_local AssessmentWork work; return work;}
#endif
bool SameAssessmentPolicy(const EvaluationContext & a,const EvaluationContext & b)
{
    return a.scale==b.scale && a.rank.rows==b.rank.rows && a.rank.design_columns==b.rank.design_columns &&
        a.rank.width_columns==b.rank.width_columns && a.linear.rank_relative==b.linear.rank_relative &&
        a.linear.active_set_iteration_factor==b.linear.active_set_iteration_factor && a.linear.release_factor==b.linear.release_factor &&
        a.linear.release_response_norm==b.linear.release_response_norm && a.atom_ids==b.atom_ids && a.row_ids==b.row_ids &&
        a.audit.directions.rows()==b.audit.directions.rows() && a.audit.directions.cols()==b.audit.directions.cols() &&
        (a.audit.directions.array()==b.audit.directions.array()).all();
}
Assessment AssessEvaluated(const Domain & domain,const Vector & y,const Evaluation & endpoint,const Evaluation & reference,
    const EvaluationContext & policy,bool supplied)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ++AssessmentWorkForTesting().assessments;
#endif
    const auto * context=&policy; const double scale=context->scale; const auto & eta=endpoint.eta; Assessment out;
    out.primary=endpoint; out.reference=reference;
    if(supplied && endpoint.valid && !endpoint.certificate.kkt_passed)
    {out.primary.valid=false; out.primary.reason="assembled-kkt-failed";}
    if(!out.primary.valid || !reference.valid) {out.failure="inner-solve"; return out;}
    const double difference=Difference(endpoint.beta,reference.beta); out.coefficient_difference=difference;
    out.design=DesignSpectrum(endpoint.x,Vector::Ones(y.size()));
    const auto differential=DifferentiateProfile(endpoint,scale,context);
    if(!differential.valid) {out.failure=differential.reason; return out;}
    const auto width_reduced=Reduce(differential.projected,Matrix(y.size(),0));
    const auto widths=Decompose(width_reduced.first,context->rank.rows);
    out.widths=CompactSpectrum(widths,context->rank.rows);
    const Vector norms=differential.projected.colwise().norm(); out.widths->column_norms=norms;
    Matrix normalized=width_reduced.first;
    for(Eigen::Index k=0;k<norms.size();++k) if(norms(k)>0) normalized.col(k)/=norms(k);
    out.normalized_widths=CompactSpectrum(Decompose(normalized,context->rank.rows),context->rank.rows);
    out.weak_directions.resize(eta.size(),std::min<Eigen::Index>(3,widths.matrixV().cols()));
    for(Eigen::Index k=0;k<out.weak_directions.cols();++k)
    {
        Vector direction=widths.matrixV().col(widths.matrixV().cols()-1-k);
        Eigen::Index sign{}; direction.cwiseAbs().maxCoeff(&sign); if(direction(sign)<0) direction=-direction;
        out.weak_directions.col(k)=direction;
    }
    const auto correction_reduced=Reduce(differential.jacobian,Matrix(-endpoint.residual/scale));
    const auto correction_svd=Decompose(correction_reduced.first,context->rank.rows);
    const Vector correction=correction_svd.solve(correction_reduced.second.col(0)); out.correction=correction;
    out.jacobian=CompactSpectrum(correction_svd,context->rank.rows);
    std::vector<Vector> directions{Vector::Ones(eta.size()).normalized(),Vector(eta.size()),out.weak_directions.col(0)};
    for(Eigen::Index k=0;k<eta.size();++k) directions[1](k)=k%2 ? -1 : 1;
    directions[1].normalize();
    if(context->audit.directions.size())
    {directions.clear(); for(Eigen::Index k=0;k<context->audit.directions.cols();++k) directions.push_back(context->audit.directions.col(k));}
    bool verified=true;
    for(std::size_t k=0;k<directions.size();++k) for(double h:{1e-4,5e-5})
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentWorkForTesting().directional_evaluations+=2;
#endif
        const auto plus=EvaluateProfile(domain,y,eta+h*directions[k],true,context),minus=EvaluateProfile(domain,y,eta-h*directions[k],true,context);
        bool passed=false,same_face=false; double error=std::numeric_limits<double>::infinity();
        if(plus.valid && minus.valid)
        {
            same_face=plus.certificate.active_atoms==endpoint.certificate.active_atoms && minus.certificate.active_atoms==endpoint.certificate.active_atoms;
            const Vector analytic=differential.jacobian*directions[k];
            const Vector finite=(plus.residual-minus.residual)/(2*h*scale);
            error=(finite-analytic).norm()/std::max({1e-12,finite.norm(),analytic.norm()});
            passed=same_face && error<=1e-6;
        }
        verified &= passed; out.derivatives.push_back({k,h,error,same_face,passed,plus.valid,minus.valid});
    }
    out.derivative_verified=verified;
    out.inner=difference<=1e-10 && out.design->rank==endpoint.x.cols();
    out.gradient=endpoint.gradient.lpNorm<Eigen::Infinity>()<=1e-12 && reference.gradient.lpNorm<Eigen::Infinity>()<=1e-12;
    out.local=correction.allFinite() && correction.lpNorm<Eigen::Infinity>()<=1e-10;
    out.identified=widths.rank()==eta.size() && correction_svd.rank()==eta.size();
    out.qualified=out.inner && out.gradient && out.local && out.identified && verified;
    out.failure=!out.inner ? "inner-solve" : !verified ? "derivative-unverified" : !out.identified ? "width-unidentified" :
        !out.gradient || !out.local ? "b-not-stationary" : "none";
    return out;
}

}
