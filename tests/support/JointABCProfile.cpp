#include "support/JointABCProfile.hpp"
#include "support/FixedBOracle.hpp"
#include "support/InstrumentedLM.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <unsupported/Eigen/NonLinearOptimization>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sys/resource.h>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;
constexpr double eps=std::numeric_limits<double>::epsilon();
j::value Number(double x) {return std::isfinite(x) ? j::value(x) : j::value(nullptr);}
j::array Values(const Vector & v) {j::array out; for (double x:v) out.push_back(Number(x)); return out;}
Vector Parse(const j::value & v)
{
    Vector out(static_cast<Eigen::Index>(v.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k)));
    return out;
}
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
double Difference(const Vector & a,const Vector & b)
{return ((a-b).array().abs()/(1+a.array().abs().max(b.array().abs()))).maxCoeff();}
void Write(const fs::path & path,const j::value & value)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';
}
std::ofstream CSV(const fs::path & path,const std::string & header)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<header<<'\n'<<std::setprecision(17); return out;
}
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
j::object Spectrum(const Eigen::JacobiSVD<Matrix> & svd,Eigen::Index rows)
{
    const auto & v=svd.singularValues();
    return {{"rank",svd.rank()},{"singular_values",Values(v)},{"minimum_singular",Number(v(v.size()-1))},
        {"condition",Number(v(0)/v(v.size()-1))},
        {"rank_threshold",Number(eps*static_cast<double>(std::max(rows,svd.cols()))*v(0))}};
}
j::object Endpoint(const Evaluation & e)
{
    auto out=e.certificate;
    out["valid"]=e.valid; out["reason"]=e.reason; out["beta"]=Values(e.beta);
    out["eta"]=Values(e.eta); out["b"]=Values(e.eta.array().exp()); out["b_gradient"]=Values(e.gradient);
    out["b_gradient_inf"]=e.gradient.size() ? Number(e.gradient.lpNorm<Eigen::Infinity>()) : j::value(nullptr);
    return out;
}
struct Profile
{
    const Domain & domain;
    const Vector & y;
    double scale;
    const EvaluationContext & context;
    Evaluation cached;
    j::array trace;
    int evaluations{},derivatives{};
    std::string failure;
    std::string variant;
    int references{};
    double reference_seconds{};
    bool guarded() const {return variant=="guarded" || variant=="guarded-log";}
    bool retry() const {return guarded() && evaluations<context.profile_budget && failure!="unrepresentable-step";}
    bool Trial(const Vector & accepted,const Vector & step,const Vector & diagonal,
        double radius,double damping,double actual,double predicted,double ratio,bool proposed)
    {
        if (trace.empty()) return false;
        auto & row=trace.back().as_object();
        row["lm"]=j::object{{"accepted_eta",Values(accepted)},{"step",Values(step)},
            {"diagonal",Values(diagonal)},{"radius",Number(radius)},{"damping",Number(damping)},
            {"actual_decrease",cached.valid ? Number(actual) : j::value(nullptr)},
            {"predicted_decrease",Number(predicted)},{"ratio",cached.valid ? Number(ratio) : j::value(nullptr)},
            {"proposed_acceptance",proposed}};
        bool trusted=cached.valid;
        if (proposed || context.audit.trial_details)
        {
            const auto start=std::chrono::steady_clock::now();
            auto evidence=Trust(domain,y,cached,&context); ++references;
            reference_seconds+=Seconds(start); trusted=evidence.at("passed").as_bool();
            row["trust"]=std::move(evidence);
        }
        if (guarded() && !trusted) failure="untrusted-trial";
        return trusted;
    }
    int values() const {return static_cast<int>(domain.rows);}
    bool Get(const Vector & eta)
    {
        if (cached.valid && cached.eta.size()==eta.size() && (cached.eta.array()==eta.array()).all()) return true;
        if (evaluations>=context.profile_budget) {failure="profile-budget"; return false;}
        if (guarded()) failure.clear();
        const auto start=std::chrono::steady_clock::now();
        cached=Evaluate(domain,y,eta,false,&context); ++evaluations;
        auto row=Endpoint(cached); row["evaluation"]=evaluations; row["accepted"]=false;
        row["seconds"]=Seconds(start); trace.push_back(row);
        if (!cached.valid) failure="inner-"+cached.reason;
        return cached.valid;
    }
    int operator()(const Vector & eta,Vector & residual)
    {
        if (!Get(eta)) return -1;
        residual=cached.residual/scale; return 0;
    }
    int df(const Vector & eta,Matrix & jacobian)
    {
        if (!Get(eta)) return -1;
        auto differential=Differentiate(cached,scale,&context); ++derivatives;
        if (!differential.valid) {failure=differential.reason; return -1;}
        jacobian=std::move(differential.jacobian); return 0;
    }
    void Accept(const Vector & eta,int update)
    {
        for (auto it=trace.rbegin();it!=trace.rend();++it)
            if ((Parse(it->at("eta")).array()==eta.array()).all())
            {it->as_object()["accepted"]=true; it->as_object()["accepted_update"]=update; break;}
    }
};
}

Domain::Domain(const unique_grid::Grid & grid,const std::vector<Atom> & geometry)
    :rows(static_cast<Eigen::Index>(grid.voxels.size())),atoms(geometry.size())
{
    for (std::size_t a=0;a<geometry.size();++a) for (Eigen::Index p=0;p<rows;++p)
    {
        const double square=SquareDistance(grid.voxels[static_cast<std::size_t>(p)].position,geometry[a].position);
        if (square<=2.5*2.5) atoms[a].push_back({p,square});
    }
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
        const auto b=EvaluateBasis(p.square,widths(a),2.5);
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
Evaluation Evaluate(const Domain & domain,const Vector & y,const Vector & eta,bool reference,const EvaluationContext * context,
    const std::vector<joint_ac::LinearBlock> * blocks)
{
    auto out=Basis(domain,y,eta); if(!out.valid) return out; out.valid=false;
    const Eigen::Index m=eta.size();
    const auto solved=joint_ac::WeightedSolve(out.x,y,Vector::Ones(y.size()),reference,true,nullptr,context ? &context->linear : nullptr,blocks);
    out.beta=solved.beta; out.certificate=fixed_b::Certificate(out.x,y,out.beta,context ? context->scale : 0);
    out.certificate["linear_solves"]=solved.solves; out.certificate["free_rank"]=solved.rank;
    if(blocks) out.certificate["block_factorizations"]=solved.block_factorizations;
    if (!solved.valid) {out.reason=solved.reason; return out;}
    if (!out.certificate.at("kkt_passed").as_bool()) {out.reason="kkt-failed"; return out;}
    out.residual=out.x*out.beta-y; out.gradient=Vector::Zero(m);
    const double scale=context ? context->scale : std::max(1.0,y.norm());
    for (Eigen::Index k=0;k<2*m;++k)
        out.gradient(k/2)+=out.beta(k)*out.derivative.col(k).dot(out.residual)/scale/scale;
    out.valid=out.residual.allFinite() && out.gradient.allFinite();
    out.reason=out.valid ? "qualified-inner" : "nonfinite-residual"; return out;
}

Evaluation AtState(const Domain & domain,const Vector & y,const Vector & eta,const Vector & beta,const EvaluationContext & context)
{
    auto out=Basis(domain,y,eta); if(!out.valid) return out; out.valid=false;
    if(beta.size()!=out.x.cols() || !beta.allFinite()) {out.reason="invalid-coefficients"; return out;}
    out.beta=beta; out.residual=out.x*beta-y; out.gradient=Vector::Zero(eta.size());
    out.certificate=fixed_b::Certificate(out.x,y,beta,context.scale);
    for(Eigen::Index k=0;k<beta.size();++k)
        out.gradient(k/2)+=beta(k)*out.derivative.col(k).dot(out.residual)/context.scale/context.scale;
    out.valid=out.residual.allFinite() && out.gradient.allFinite(); out.reason=out.valid ? "raw-state" : "nonfinite-state";
    return out;
}
namespace {
template<class Design> j::object SpectrumRecord(const Design & input,const RankPolicy & policy,Eigen::Index columns,bool normalize)
{
    if(input.cols()==0) return {{"available",false},{"reason","empty-matrix"}};
    if(input.rows()==0) return {{"available",true},{"singular_values",j::array{}},{"rank",0},
        {"rank_threshold",0},{"global_rows",policy.rows},{"global_columns",columns},{"reason","unobserved-columns"}};
    Design x=input; Vector norms(input.cols());
    for(Eigen::Index k=0;k<x.cols();++k) {norms(k)=x.col(k).norm(); if(normalize && norms(k)>0) x.col(k)/=norms(k);}
    const auto reduced=Reduce(x,Matrix(x.rows(),0));
    const Eigen::JacobiSVD<Matrix> svd(reduced.first,Eigen::ComputeThinV); const auto & values=svd.singularValues();
    const double threshold=policy.Absolute(columns,values(0));
    return {{"available",true},{"singular_values",Values(values)},{"rank",(values.array()>threshold).count()},
        {"rank_threshold",threshold},{"global_rows",policy.rows},{"global_columns",columns},
        {"column_norms",Values(norms)},{"minimum_singular",Number(values(values.size()-1))},
        {"condition",Number(values(0)/values(values.size()-1))}};
}
}
j::object MatrixSpectrum(const Sparse & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return SpectrumRecord(x,p,columns,normalize);}
j::object MatrixSpectrum(const Matrix & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return SpectrumRecord(x,p,columns,normalize);}
Vector LocalCorrection(const Evaluation & e,const Differential & d,const EvaluationContext & context,double absolute)
{
    if(!e.valid || !d.valid) return {};
    const auto reduced=Reduce(d.jacobian,Matrix(-e.residual/context.scale));
    auto svd=Decompose(reduced.first,context.rank.rows);
    if(absolute>=0 && svd.singularValues()(0)>0) svd.setThreshold(absolute/svd.singularValues()(0));
    if(svd.rank()!=d.jacobian.cols()) return {};
    return svd.solve(reduced.second.col(0));
}

Differential Differentiate(const Evaluation & e,double scale,const EvaluationContext * context,double free_design_threshold)
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

j::object Trust(const Domain & domain,const Vector & y,const Evaluation & e,const EvaluationContext * context)
{
    const auto reference=Evaluate(domain,y,e.eta,true,context);
    j::object out{{"reference",Endpoint(reference)},{"passed",false}};
    if (!e.valid)
    {
        out["reason"]="invalid-primary";
        if(e.x.cols()>0 && e.x.nonZeros()>0 && (context ? context->audit.trial_details : RegisteredAudit(static_cast<Eigen::Index>(domain.atoms.size())).trial_details))
            out["design_spectrum"]=joint_ac::SparseSpectrum(e.x,Vector::Ones(y.size()));
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
    const double kkt_difference=std::abs(kkt-j::value_to<double>(e.certificate.at("projected_kkt")));
    out["scaled_coefficient_difference"]=Number(difference); out["prediction_passed"]=prediction_ok;
    out["gradient_passed"]=gradient_ok && reference_gradient_ok;
    out["kkt_replay_difference"]=Number(kkt_difference);
    out["maximum_prediction_difference"]=Number((prediction-original_prediction).lpNorm<Eigen::Infinity>());
    out["maximum_gradient_difference"]=Number((width_gradient-e.gradient).lpNorm<Eigen::Infinity>());
    out["maximum_cancellation_ratio"]=Number((absolute.array()/prediction.array().abs().max(1.0)).maxCoeff());
    if ((context ? context->audit.trial_details : RegisteredAudit(static_cast<Eigen::Index>(domain.atoms.size())).trial_details)) out["design_spectrum"]=joint_ac::SparseSpectrum(e.x,Vector::Ones(y.size()));
    out["passed"]=difference<=1e-10 && prediction_ok && gradient_ok && reference_gradient_ok && kkt_difference<=1e-13;
    out["reason"]=out.at("passed").as_bool() ? "trusted" : !reference.valid ? "invalid-reference" : "replay-disagreement";
    return out;
}

j::object Fit(const Domain & domain,const Vector & y,const Vector & initial_b,j::object * resources,const std::string & variant,const EvaluationContext * provided)
{
    const auto fallback=provided ? EvaluationContext{} : MakeContext(y,static_cast<Eigen::Index>(domain.atoms.size()));
    const auto * context=provided ? provided : &fallback;
    const auto start=std::chrono::steady_clock::now();
    auto measure=[](auto since) {
        struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        const auto bytes=usage.ru_maxrss;
#else
        const auto bytes=usage.ru_maxrss*1024;
#endif
        return j::object{{"seconds",Seconds(since)},{"process_peak_rss_bytes",bytes}};
    };
    const double scale=context ? context->scale : std::max(1.0,y.norm());
    Profile profile{domain,y,scale,*context,{}, {},0,0,{},variant};
    Vector eta=initial_b.array().log(); int accepted{};
    auto search=[&](auto & lm) {
    lm.parameters.factor=.1; lm.parameters.ftol=1e-14; lm.parameters.xtol=1e-12;
    lm.parameters.gtol=1e-12; lm.parameters.maxfev=context->profile_budget;
    if (variant=="guarded-log") {lm.useExternalScaling=true; lm.diag=Vector::Ones(eta.size());}
    auto status=lm.minimizeInit(eta);
    if (!variant.empty() && profile.cached.valid)
    {
        const bool trusted=profile.Trial(eta,Vector::Zero(eta.size()),Vector::Ones(eta.size()),0,0,0,0,0,true);
        if (profile.guarded() && !trusted) return Eigen::LevenbergMarquardtSpace::UserAsked;
    }
    if (profile.cached.valid) profile.Accept(eta,0);
    while (status==Eigen::LevenbergMarquardtSpace::NotStarted || status==Eigen::LevenbergMarquardtSpace::Running)
    {
        if (accepted>=context->update_budget) {profile.failure="accepted-update-budget"; break;}
        const auto before=lm.iter;
        status=lm.minimizeOneStep(eta);
        if (lm.iter>before) {++accepted; profile.Accept(eta,accepted);}
    }
    return status;
    };
    Eigen::LevenbergMarquardtSpace::Status status;
    if (variant.empty()) {Eigen::LevenbergMarquardt<Profile> lm(profile); status=search(lm);}
    else {InstrumentedLM<Profile> lm(profile); status=search(lm);}
    auto initial=profile.trace.empty() ? Endpoint(profile.cached) : profile.trace.front().as_object();
    for (const char * key:{"evaluation","accepted","accepted_update","seconds","trust","lm"}) initial.erase(key);
    if (resources) (*resources)["search"]=measure(start);
    const auto audit_start=std::chrono::steady_clock::now();
    auto finish=[&](j::object & result) {
        if (resources) (*resources)["endpoint_audit"]=measure(audit_start);
        result["seconds"]=Seconds(start);
        return result;
    };
    j::object out{{"schema_version",1},{"experiment","joint-abc-profile"},{"alpha",0},
        {"execution_complete",true},{"joint_qualified",false},{"initial",initial},
        {"lm_status",static_cast<int>(status)},{"stop_reason",profile.failure.empty() ? "native-lm-stop" : profile.failure},
        {"profile_evaluations",profile.evaluations},{"jacobian_evaluations",profile.derivatives},
        {"accepted_updates",accepted},{"trials",profile.trace},{"row_count",y.size()},
        {"settings",j::object{{"factor",.1},{"ftol",1e-14},{"xtol",1e-12},{"gtol",1e-12},
            {"profile_budget",context->profile_budget},{"accepted_update_budget",context->update_budget}}},
        {"linear_solver","sparse-qr-householder-1024"},{"reference_solver","independent-tsqr-8192-svd"},
        {"residual_scale",scale},{"variance_semantics","descriptive RSS/N; zero permitted"}};
    if (!variant.empty())
    {
        out["variant"]=variant; out["search_reference_evaluations"]=profile.references;
        out["initial_accepted"]=!profile.trace.empty() && profile.trace.front().at("accepted").as_bool();
        out["search_reference_seconds"]=profile.reference_seconds;
        out["search_stopped_without_convergence"]=status==Eigen::LevenbergMarquardtSpace::UserAsked ||
            status==Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation || !profile.failure.empty();
    }
    auto assessment=Assess(domain,y,eta,*context);
    for(auto & field:assessment) out[field.key()]=std::move(field.value());
    return finish(out);
}

j::object Assess(const Domain & domain,const Vector & y,const Vector & eta,const EvaluationContext & policy,const Vector * supplied_beta)
{
    const auto * context=&policy; const double scale=context->scale;
    j::object out{{"joint_qualified",false}};
    // Fresh endpoint evaluation, independent of the optimizer's cached trial.
    auto endpoint=supplied_beta ? AtState(domain,y,eta,*supplied_beta,*context) : Evaluate(domain,y,eta,false,context);
    const auto reference=Evaluate(domain,y,eta,true,context);
    if(supplied_beta && endpoint.valid && !endpoint.certificate.at("kkt_passed").as_bool())
    {endpoint.valid=false; endpoint.reason="assembled-kkt-failed";}
    out["primary"]=Endpoint(endpoint); out["reference"]=Endpoint(reference);
    out["endpoint_evaluations"]=2; out["directional_evaluations"]=0;
    if (!endpoint.valid || !reference.valid)
    {
        out["qualification_failure"]="inner-solve"; return out;
    }
    const double difference=Difference(endpoint.beta,reference.beta);
    out["scaled_reference_difference"]=Number(difference);
    out["design_spectrum"]=joint_ac::SparseSpectrum(endpoint.x,Vector::Ones(y.size()));
    const auto differential=Differentiate(endpoint,scale,context);
    if (!differential.valid)
    {
        out["qualification_failure"]=differential.reason; return out;
    }
    const auto width_reduced=Reduce(differential.projected,Matrix(y.size(),0));
    const auto widths=Decompose(width_reduced.first,context->rank.rows);
    auto spectrum=Spectrum(widths,context->rank.rows);
    const Vector norms=differential.projected.colwise().norm(); spectrum["column_norms"]=Values(norms);
    Matrix normalized=width_reduced.first;
    for (Eigen::Index k=0;k<norms.size();++k) if (norms(k)>0) normalized.col(k)/=norms(k);
    spectrum["column_normalized"]=Spectrum(Decompose(normalized,context->rank.rows),context->rank.rows);
    spectrum["active_face_only"]=!endpoint.certificate.at("active_atoms").as_array().empty();
    j::array weak;
    for (Eigen::Index k=0;k<std::min<Eigen::Index>(3,widths.matrixV().cols());++k)
    {
        Vector direction=widths.matrixV().col(widths.matrixV().cols()-1-k);
        Eigen::Index sign{}; direction.cwiseAbs().maxCoeff(&sign); if (direction(sign)<0) direction=-direction;
        weak.push_back(Values(direction));
    }
    spectrum["weak_directions"]=weak; out["width_spectrum"]=spectrum;
    const auto correction_reduced=Reduce(differential.jacobian,Matrix(-endpoint.residual/scale));
    const auto correction_svd=Decompose(correction_reduced.first,context->rank.rows);
    const Vector correction=correction_svd.solve(correction_reduced.second.col(0));
    out["local_correction"]=Values(correction); out["local_correction_inf"]=Number(correction.lpNorm<Eigen::Infinity>());
    out["profile_jacobian_spectrum"]=Spectrum(correction_svd,context->rank.rows);
    std::vector<Vector> directions{Vector::Ones(eta.size()).normalized(),Vector(eta.size()),Parse(weak[0])};
    for (Eigen::Index k=0;k<eta.size();++k) directions[1](k)=k%2 ? -1 : 1;
    directions[1].normalize();
    if(context->audit.directions.size())
    {directions.clear(); for(Eigen::Index k=0;k<context->audit.directions.cols();++k) directions.push_back(context->audit.directions.col(k));}
    j::array checks; bool verified=true; int validation_evaluations{};
    for (std::size_t k=0;k<directions.size();++k) for (double h:{1e-4,5e-5})
    {
        const auto plus=Evaluate(domain,y,eta+h*directions[k],true,context),minus=Evaluate(domain,y,eta-h*directions[k],true,context);
        validation_evaluations+=2;
        bool passed=false,same_face=false; double error=std::numeric_limits<double>::infinity();
        if (plus.valid && minus.valid)
        {
            same_face=plus.certificate.at("active_atoms")==endpoint.certificate.at("active_atoms") &&
                minus.certificate.at("active_atoms")==endpoint.certificate.at("active_atoms");
            const Vector analytic=differential.jacobian*directions[k];
            const Vector finite=(plus.residual-minus.residual)/(2*h*scale);
            error=(finite-analytic).norm()/std::max({1e-12,finite.norm(),analytic.norm()});
            passed=same_face && error<=1e-6;
        }
        verified &= passed;
        checks.emplace_back(j::object{{"direction",k},{"h",h},{"relative_l2_difference",Number(error)},
            {"same_active_face",same_face},{"passed",passed},{"plus_valid",plus.valid},{"minus_valid",minus.valid}});
    }
    out["directional_evaluations"]=validation_evaluations; out["derivative_checks"]=checks;
    out["derivative_verified"]=verified;
    const bool inner=difference<=1e-10 && j::value_to<int>(out.at("design_spectrum").at("rank"))==endpoint.x.cols();
    const bool gradient=endpoint.gradient.lpNorm<Eigen::Infinity>()<=1e-12 && reference.gradient.lpNorm<Eigen::Infinity>()<=1e-12;
    const bool local=correction.allFinite() && correction.lpNorm<Eigen::Infinity>()<=1e-10;
    const bool identified=widths.rank()==eta.size() && correction_svd.rank()==eta.size();
    out["joint_qualified"]=inner && gradient && local && identified && verified;
    out["qualification_checks"]=j::object{{"inner",inner},{"b_gradient",gradient},
        {"local_correction",local},{"identified",identified},{"derivative",verified}};
    out["qualification_failure"]=!inner ? "inner-solve" : !verified ? "derivative-unverified" :
        !identified ? "width-unidentified" : !gradient || !local ? "b-not-stationary" : "none";
    return out;
}

void Run(const std::string & manifest,const std::string & map,const std::string & checkpoint,const std::string & output_path)
{
    const auto data=fixed_b::Prepare(manifest,map,checkpoint,output_path,"joint-abc-profile");
    const Domain domain(data.grid,data.atoms); const fs::path output(output_path);
    fs::create_directories(output/"weak-directions"); j::array completed;
    for (bool quantized:{false,true}) for (const std::string start:{"checkpoint","narrower","wider","mixed"})
    {
        const std::string name=start+(quantized ? "-float32" : "-double");
        Vector initial=data.checkpoint_b;
        for (Eigen::Index a=0;a<initial.size();++a)
        {
            const int serial=j::value_to<int>(data.identities.at(static_cast<std::size_t>(a)).at("serial_id"));
            if (start=="narrower" || (start=="mixed" && serial%2)) initial(a)*=.8;
            else if (start=="wider" || start=="mixed") initial(a)*=1.2;
        }
        const auto & y=quantized ? data.y32 : data.y64;
        std::cout<<"Starting "<<name<<std::endl;
        auto fit=Fit(domain,y,initial); fit["case"]=name; fit["initial_b"]=Values(initial);
        struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        fit["process_peak_rss_bytes"]=usage.ru_maxrss;
#else
        fit["process_peak_rss_bytes"]=usage.ru_maxrss*1024;
#endif
        const auto & endpoint=fit.at("primary");
        if (endpoint.at("valid").as_bool())
        {
            auto model=data.atoms; const auto beta=Parse(endpoint.at("beta")),b=Parse(endpoint.at("b"));
            // Output replay is not an optimizer evaluation or an endpoint certificate.
            for (std::size_t a=0;a<model.size();++a)
            {model[a].amplitude=beta(static_cast<Eigen::Index>(2*a)); model[a].charge=beta(static_cast<Eigen::Index>(2*a+1)); model[a].width=b(static_cast<Eigen::Index>(a));}
            auto residual=CSV(output/"residuals"/(name+".csv"),"row,prediction,residual");
            for (Eigen::Index p=0;p<y.size();++p)
            {
                // Match sparse column accumulation order on the original rows.
                double prediction{};
                for (const auto & atom:model)
                {
                    const auto basis=EvaluateBasis(SquareDistance(data.grid.voxels[static_cast<std::size_t>(p)].position,atom.position),atom.width,2.5);
                    prediction+=atom.amplitude*basis.gaussian; prediction+=atom.charge*basis.charge;
                }
                residual<<p<<','<<prediction<<','<<prediction-y(p)<<'\n';
            }
            residual.close();
            fit["residual_sha256"]=rhbm_gem::core::simulation::FileSha256((output/"residuals"/(name+".csv")).string());
        }
        if (fit.contains("width_spectrum"))
        {
            auto weak=CSV(output/"weak-directions"/(name+".csv"),"serial_id,weakest,second,third");
            const auto & directions=fit.at("width_spectrum").at("weak_directions");
            for (std::size_t a=0;a<data.atoms.size();++a)
            {
                weak<<j::value_to<int>(data.identities[a].at("serial_id"));
                for (std::size_t k=0;k<3;++k) weak<<','<<j::value_to<double>(directions.at(k).at(a));
                weak<<'\n';
            }
        }
        Write(output/"fits"/(name+".json"),fit); completed.emplace_back(name);
        Write(output/"completion.json",j::object{{"complete",completed.size()==8},{"execution_complete",completed.size()==8},{"cases",completed}});
        std::cout<<name<<" qualified="<<fit.at("joint_qualified")<<" seconds="<<fit.at("seconds")<<std::endl;
    }
}
} // namespace second_stage_test::matched::joint_abc
