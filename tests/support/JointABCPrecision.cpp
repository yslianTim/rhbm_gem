#include "support/JointABCPrecision.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <Eigen/QR>
#include <chrono>
#include <set>
#include <map>
#include <sys/resource.h>

namespace second_stage_test::matched::certification {
namespace {
namespace j=boost::json;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
template<class T> using M=Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template<class T> using V=Eigen::Matrix<T,Eigen::Dynamic,1>;
using P50=boost::multiprecision::cpp_dec_float_50;
using P100=boost::multiprecision::cpp_dec_float_100;
template<class T> std::string String(const T & v) {return v.str(std::numeric_limits<T>::max_digits10,std::ios_base::scientific);}
j::value Number(double v) {return std::isfinite(v) ? j::value(v) : j::value(nullptr);}
j::array Values(const Vector & v) {j::array a; for(double x:v) a.push_back(Number(x)); return a;}
long PeakRSS()
{
    struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    return usage.ru_maxrss;
#else
    return usage.ru_maxrss*1024;
#endif
}
template<class T> V<T> Promote(const Vector & v) {V<T> r(v.size()); for(Eigen::Index i=0;i<v.size();++i) r(i)=T(v(i)); return r;}
template<class T> Vector Demote(const V<T> & v) {Vector r(v.size()); for(Eigen::Index i=0;i<v.size();++i) r(i)=v(i).template convert_to<double>(); return r;}

template<class T> void Column(const joint_abc::Domain & domain,Eigen::Index a,const T & eta,M<T> & x,M<T> * derivative=nullptr)
{
    using boost::multiprecision::sqrt; using boost::multiprecision::exp;
    const T b=exp(eta),pi=boost::math::constants::pi<T>(),b2=b*b,normal=T(1)/(sqrt(2*pi*b2)*(2*pi*b2));
    const T center=sqrt(2/pi)/b;
    std::map<double,std::array<T,4>> cache;
    for(const auto & s:domain.atoms[static_cast<std::size_t>(a)])
    {
        const auto found=cache.find(s.square);
        if(found!=cache.end())
        {
            x(s.row,2*a)=found->second[0]; x(s.row,2*a+1)=found->second[1];
            if(derivative) {(*derivative)(s.row,2*a)=found->second[2]; (*derivative)(s.row,2*a+1)=found->second[3];}
            continue;
        }
        const T square(s.square),r=sqrt(square),exponent=exp(-square/(2*b2));
        x(s.row,2*a)=normal*exponent;
        x(s.row,2*a+1)=s.square<1e-10 ? center : T(boost::math::erf(r/(b*sqrt(T(2))))/r);
        if(derivative)
        {
            (*derivative)(s.row,2*a)=x(s.row,2*a)*(square/b2-3);
            (*derivative)(s.row,2*a+1)=s.square<1e-10 ? T(-center) : T(-center*exponent);
        }
        cache.emplace(s.square,std::array<T,4>{x(s.row,2*a),x(s.row,2*a+1),
            derivative ? (*derivative)(s.row,2*a) : T(0),derivative ? (*derivative)(s.row,2*a+1) : T(0)});
    }
}
template<class T> M<T> Design(const joint_abc::Domain & domain,const V<T> & eta,M<T> * derivative=nullptr)
{
    M<T> x=M<T>::Zero(domain.rows,2*eta.size());
    if(derivative) *derivative=M<T>::Zero(domain.rows,2*eta.size());
    for(Eigen::Index a=0;a<eta.size();++a) Column(domain,a,eta(a),x,derivative);
    return x;
}

template<class T> struct Linear
{
    V<T> beta,residual,norms;
    std::vector<Eigen::Index> free;
    M<T> z,r;
    Eigen::HouseholderQR<M<T>> qr;
    T kkt{};
    bool valid{};
};
template<class T> Linear<T> Solve(const M<T> & x,const V<T> & y,const std::set<Eigen::Index> & active)
{
    using boost::multiprecision::abs;
    Linear<T> out; out.beta=V<T>::Zero(x.cols()); out.norms=x.colwise().norm();
    for(Eigen::Index k=0;k<x.cols();++k) if(!active.contains(k)) out.free.push_back(k);
    out.z.resize(x.rows(),static_cast<Eigen::Index>(out.free.size()));
    for(std::size_t i=0;i<out.free.size();++i)
    {
        const auto k=out.free[i]; if(out.norms(k)==0) return out;
        out.z.col(static_cast<Eigen::Index>(i))=x.col(k)/out.norms(k);
    }
    out.qr.compute(out.z); const Eigen::Index p=out.z.cols();
    out.r=out.qr.matrixQR().topRows(p).template triangularView<Eigen::Upper>();
    const T tolerance=pow(T(10),-std::numeric_limits<T>::digits10+10);
    for(Eigen::Index k=0;k<p;++k) if(abs(out.r(k,k))<tolerance) return out;
    const V<T> coefficients=out.qr.solve(y);
    for(std::size_t i=0;i<out.free.size();++i) out.beta(out.free[i])=coefficients(static_cast<Eigen::Index>(i))/out.norms(out.free[i]);
    out.residual=x*out.beta-y; const T scale=std::max(T(1),T(y.norm()));
    const V<T> g=x.transpose()*out.residual;
    for(Eigen::Index k=0;k<x.cols();++k)
    {
        const T u=out.norms(k)*out.beta(k)/scale,gradient=g(k)/out.norms(k)/scale;
        const T projected=k%2 ? T(u-gradient) : std::max(T(0),T(u-gradient));
        out.kkt=std::max(out.kkt,T(abs(u-projected)));
    }
    out.valid=true; return out;
}
template<class T> Linear<T> Constrained(const M<T> & x,const V<T> & y,std::set<Eigen::Index> active)
{
    using boost::multiprecision::abs;
    std::set<std::set<Eigen::Index>> visited;
    for(Eigen::Index iteration=0;iteration<4*x.cols()*x.cols();++iteration)
    {
        if(!visited.insert(active).second) break;
        auto out=Solve(x,y,active); if(!out.valid) return out;
        Eigen::Index negative=-1;
        for(Eigen::Index k=0;k<x.cols();k+=2) if(out.beta(k)<0 && (negative<0 || out.beta(k)*out.norms(k)<out.beta(negative)*out.norms(negative))) negative=k;
        if(negative>=0) {active.insert(negative); continue;}
        const V<T> gradient=x.transpose()*out.residual; Eigen::Index release=-1;
        const T tolerance=pow(T(10),-std::numeric_limits<T>::digits10+10)*std::max(T(1),T(y.norm()));
        for(auto k:active) if(gradient(k)/out.norms(k)<-tolerance && (release<0 || gradient(k)/out.norms(k)<gradient(release)/out.norms(release))) release=k;
        if(release<0) return out;
        active.erase(release);
    }
    return {};
}

template<class T> struct Precision
{
    bool valid{},feasible{};
    V<T> beta,gradient,correction;
    M<T> directional;
    T kkt{};
};
template<class T> Precision<T> Calculate(const joint_abc::Domain & domain,const Vector & y,const joint_abc::Evaluation & e,const Matrix & directions)
{
    Precision<T> out; const V<T> eta=Promote<T>(e.eta),response=Promote<T>(y);
    M<T> dx; const auto x=Design(domain,eta,&dx);
    std::set<Eigen::Index> active; for(Eigen::Index k=0;k<e.beta.size();k+=2) if(e.beta(k)==0) active.insert(k);
    auto inner=Solve(x,response,active); if(!inner.valid) return out;
    out.beta=inner.beta; out.kkt=inner.kkt; out.feasible=true;
    for(Eigen::Index k=0;k<out.beta.size();k+=2) out.feasible &= out.beta(k)>=0;
    const T scale=std::max(T(1),T(response.norm()));
    M<T> raw=M<T>::Zero(domain.rows,eta.size());
    for(Eigen::Index k=0;k<x.cols();++k) raw.col(k/2)+=dx.col(k)*inner.beta(k);
    // Dense QR differentiated stationarity; no rounded double basis/solution reuse.
    M<T> rhs(static_cast<Eigen::Index>(inner.free.size()),eta.size()); rhs.setZero();
    for(std::size_t i=0;i<inner.free.size();++i)
    {
        const auto k=inner.free[i]; rhs(static_cast<Eigen::Index>(i),k/2)=dx.col(k).dot(inner.residual)/inner.norms(k);
    }
    const M<T> adjoint=inner.r.transpose().template triangularView<Eigen::Lower>().solve(rhs);
    const M<T> correction=inner.r.template triangularView<Eigen::Upper>().solve(adjoint);
    const M<T> jacobian=(raw-inner.z*inner.qr.solve(raw)-inner.z*correction)/scale;
    out.gradient=raw.transpose()*inner.residual/(scale*scale);
    const Eigen::HouseholderQR<M<T>> qr(jacobian);
    const M<T> r=qr.matrixQR().topRows(eta.size()).template triangularView<Eigen::Upper>();
    const T tolerance=pow(T(10),-std::numeric_limits<T>::digits10+10);
    for(Eigen::Index k=0;k<eta.size();++k) if(abs(r(k,k))<tolerance) return out;
    out.correction=qr.solve(-inner.residual/scale);
    out.directional=jacobian*directions.template cast<T>(); out.valid=true; return out;
}
template<class A,class B> P100 Compare(const A & a,const B & b)
{
    P100 difference{};
    for(Eigen::Index k=0;k<a.size();++k)
    {
        const P100 x(String(a.data()[k])),y(String(b.data()[k]));
        difference=std::max(difference,P100(abs(x-y)/(1+std::max(P100(abs(x)),P100(abs(y))))));
    }
    return difference;
}

template<class T> j::object Scan(const joint_abc::Domain & domain,const Vector & y,const Vector & eta0,const Vector & beta0)
{
    const auto eta=Promote<T>(eta0),beta=Promote<T>(beta0),response=Promote<T>(y);
    M<T> derivative;
    const auto base=Design(domain,eta,&derivative);
    // A single dense QR of [X0,y] preserves the entire constant subspace.
    // At each scan point only two columns change. Orthogonally reducing those
    // columns gives an equivalent small dense constrained LS problem, avoiding
    // repeated large QR factorizations without forming normal equations.
    const Eigen::Index p=base.cols()+1;
    M<T> augmented(base.rows(),p); augmented.leftCols(base.cols())=base; augmented.col(p-1)=response;
    const Eigen::HouseholderQR<M<T>> common(augmented);
    const M<T> r=common.matrixQR().topRows(p).template triangularView<Eigen::Upper>();
    V<T> reduced_y=V<T>::Zero(p+2); reduced_y.head(p)=r.col(p-1);
    V<T> prediction=V<T>::Zero(p+2); prediction.head(p)=r.leftCols(base.cols())*beta;
    const T scale=std::max(T(1),T(response.norm())),pi=boost::math::constants::pi<T>();
    std::set<Eigen::Index> active; for(Eigen::Index k=0;k<beta.size();k+=2) if(beta(k)==0) active.insert(k);
    j::array rows;
    for(Eigen::Index atom:{1,5,9})
    {
        const V<T> first_order=derivative.col(2*atom)*beta(2*atom)+derivative.col(2*atom+1)*beta(2*atom+1)+
            base.col(2*atom)*(4*pi*exp(2*eta(atom))*beta(2*atom+1));
        T near_zero_max{}; std::size_t near_zero_rows{};
        for(const auto & s:domain.atoms[static_cast<std::size_t>(atom)]) if(s.square<1e-10)
        {++near_zero_rows; near_zero_max=std::max(near_zero_max,T(abs(first_order(s.row))));}
        T previous{};
        for(int k=4;k<=24;++k)
        {
            const T step=pow(T(2),-k); V<T> changed_eta=eta,path_beta=beta;
            changed_eta(atom)+=step; path_beta(2*atom)+=4*pi*exp(2*eta(atom))*beta(2*atom+1)*step;
            auto changed=base; Column(domain,atom,changed_eta(atom),changed);
            const M<T> rotated=common.householderQ().adjoint()*M<T>(changed.middleCols(2*atom,2));
            const Eigen::HouseholderQR<M<T>> tail(rotated.bottomRows(rotated.rows()-p));
            M<T> design=M<T>::Zero(p+2,base.cols()); design.topRows(p)=r.leftCols(base.cols());
            design.block(0,2*atom,p,2)=rotated.topRows(p);
            design.block(p,2*atom,2,2)=tail.matrixQR().topRows(2).template triangularView<Eigen::Upper>();
            const auto fixed=Solve(design,reduced_y,active),profile=Constrained(design,reduced_y,active);
            const V<T> delta=design*path_beta-prediction,residual=design*path_beta-reduced_y;
            const T change=delta.norm(),structure_loss=change*change/2;
            rows.push_back(j::object{{"atom",atom},{"k",k},{"step",String(step)},
                {"first_order_prediction_norm",String(T(first_order.norm()))},
                {"near_zero_rows",near_zero_rows},{"near_zero_first_order_max",String(near_zero_max)},
                {"prediction_change",String(change)},{"structure_loss",String(structure_loss)},
                {"prediction_slope",previous>0 ? Number((log(previous/change)/log(T(2))).template convert_to<double>()) : j::value(nullptr)},
                {"path_objective",String(T(residual.squaredNorm()/2))},
                {"fixed_face_valid",fixed.valid},{"profile_valid",profile.valid},
                {"fixed_face_objective",fixed.valid ? j::value(String(T(fixed.residual.squaredNorm()/2))) : j::value(nullptr)},
                {"profile_objective",profile.valid ? j::value(String(T(profile.residual.squaredNorm()/2))) : j::value(nullptr)},
                {"profile_kkt",profile.valid ? j::value(String(profile.kkt)) : j::value(nullptr)},
                {"path_feasible",path_beta(2*atom)>=0},{"scale",String(scale)}});
            previous=change;
        }
    }
    return {{"digits",std::numeric_limits<T>::digits10},{"rows",rows},{"active_a",static_cast<std::int64_t>(active.size())}};
}
}

j::object PrecisionAudit(const joint_abc::Domain & domain,const Vector & y,const joint_abc::Evaluation & e,const Matrix & directions)
{
    const auto start=std::chrono::steady_clock::now();
    const auto a=Calculate<P50>(domain,y,e,directions);
    const auto b=Calculate<P100>(domain,y,e,directions);
    j::object out{{"precisions",j::array{50,100}},{"solver","independent-dense-householder-qr"},
        {"valid50",a.valid},{"valid100",b.valid},{"agreement_passed",false},{"derivative_passed",false}};
    if(a.valid && b.valid)
    {
        const P100 difference=std::max({Compare(a.beta,b.beta),Compare(a.gradient,b.gradient),Compare(a.correction,b.correction),Compare(a.directional,b.directional)});
        const auto analytic=joint_abc::Differentiate(e,std::max(1.0,y.norm())); j::array errors; bool derivative=true;
        for(Eigen::Index k=0;k<directions.cols();++k)
        {
            const Vector reference=Demote<P100>(b.directional.col(k)),observed=analytic.jacobian*directions.col(k);
            const double error=(reference-observed).norm()/std::max({1e-12,reference.norm(),observed.norm()}); errors.push_back(Number(error)); derivative &= error<=1e-6;
        }
        out["maximum_scaled_precision_difference"]=String(difference); out["agreement_passed"]=difference<=P100("1e-20");
        out["derivative_relative_errors"]=errors; out["derivative_passed"]=derivative && a.feasible && b.feasible;
        out["fixed_face_feasible"]=a.feasible && b.feasible; out["projected_kkt"]=String(b.kkt);
        out["local_correction"]=Values(Demote(b.correction)); out["local_correction_inf"]=Number(Demote(b.correction).lpNorm<Eigen::Infinity>());
        out["local_correction_passed"]=a.feasible && b.feasible && b.correction.cwiseAbs().maxCoeff()<=P100("1e-10");
        out["b_gradient"]=Values(Demote(b.gradient));
    }
    out["process_peak_rss_bytes"]=PeakRSS();
    out["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); return out;
}

j::object BoundaryAudit(const joint_abc::Domain & domain,const Vector & y,const Vector & eta,const Vector & beta)
{
    const auto start=std::chrono::steady_clock::now();
    auto a=Scan<P50>(domain,y,eta,beta),b=Scan<P100>(domain,y,eta,beta);
    bool agreement=true; P100 maximum{};
    for(std::size_t k=0;k<a.at("rows").as_array().size();++k)
    {
        const auto & x=a.at("rows").at(k); const auto & z=b.at("rows").at(k);
        agreement &= x.at("fixed_face_valid")==z.at("fixed_face_valid") && x.at("profile_valid")==z.at("profile_valid");
        for(const char * key:{"prediction_change","structure_loss","path_objective","fixed_face_objective","profile_objective"})
        {
            if(x.at(key).is_null() || z.at(key).is_null()) {agreement &= x.at(key).is_null() && z.at(key).is_null(); continue;}
            const P100 u(j::value_to<std::string>(x.at(key))),v(j::value_to<std::string>(z.at(key)));
            maximum=std::max(maximum,P100(abs(u-v)/(1+std::max(P100(abs(u)),P100(abs(v))))));
        }
    }
    return {{"schema_version",1},{"interpretation","positive feasible compensation directions; no regular parameter recovery claim"},
        {"precision50",a},{"precision100",b},{"agreement_passed",agreement && maximum<=P100("1e-20")},
        {"maximum_scaled_precision_difference",String(maximum)},
        {"process_peak_rss_bytes",PeakRSS()},
        {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
}
} // namespace second_stage_test::matched::certification
