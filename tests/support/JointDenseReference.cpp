#include "JointDenseReference.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"

namespace second_stage_test::matched::joint_abc {
using namespace runtime;
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
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
Vector DenseLocalCorrection(const Evaluation & e,const DenseDifferential & d,const EvaluationContext & context,double absolute)
{
    if(!e.valid || !d.valid) return {};
    const auto reduced=Reduce(d.jacobian,Matrix(-e.residual/context.scale));
    auto svd=Decompose(reduced.first,context.rank.rows);
    if(absolute>=0 && svd.singularValues()(0)>0) svd.setThreshold(absolute/svd.singularValues()(0));
    if(svd.rank()!=d.jacobian.cols()) return {};
    return svd.solve(reduced.second.col(0));
}

DenseDifferential DenseDifferentiate(const Evaluation & e,double scale,const EvaluationContext * context,double free_design_threshold)
{
    DenseDifferential out;
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

DenseDifferential MaterializeDerivative(const runtime::Evaluation & e,double scale,const EvaluationContext * context,double absolute,Eigen::Index tile)
{
    const auto prepared=PrepareDerivative(e,scale,context,absolute,tile); DenseDifferential out;
    out.valid=prepared.valid; out.reason=prepared.reason; if(!out.valid) return out;
    out.projected.resize(e.x.rows(),e.eta.size()); out.jacobian.resize(e.x.rows(),e.eta.size());
    Matrix projected,jacobian;
    for(Eigen::Index first=0;first<e.x.rows();first+=tile)
    {
        const auto count=std::min(tile,e.x.rows()-first); prepared.Rows(first,count,projected,jacobian);
        out.projected.middleRows(first,count)=projected; out.jacobian.middleRows(first,count)=jacobian;
    }
    out.valid=out.projected.allFinite() && out.jacobian.allFinite();
    if(!out.valid) out.reason="nonfinite-derivative"; return out;
}

}
