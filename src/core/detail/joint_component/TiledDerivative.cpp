#include "TiledDerivative.hpp"
#include "SparseFactor.hpp"
#include "CompactSvd.hpp"
#include <cmath>
#include <optional>

namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
DerivativeWork & DerivativeWorkForTesting() {static thread_local DerivativeWork work; return work;}
#endif
TiledDifferential PrepareDerivative(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile)
{
    auto & work=SparseWorkForTesting(); ++work.derivative_preparations;
    ResourcePhase phase("derivative-prepare");
    WorkTimer preparation_timer(work.derivative_seconds);
    TiledDifferential out;
    if(!e.valid || !(scale>0) || !std::isfinite(scale)) {out.reason="invalid-inner"; return out;}
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
    std::vector<Eigen::Index> free;
    for(Eigen::Index k=0;k<e.beta.size();++k) if(k%2 || e.beta(k)>0) free.push_back(k);
    const Eigen::Index n=e.x.rows(),m=e.beta.size()/2,p=static_cast<Eigen::Index>(free.size());
    out.scale=scale; out.free_design.resize(n,p); out.raw.resize(n,m);
    std::vector<Eigen::Triplet<double>> entries,raw;
    Matrix t=Matrix::Zero(p,m);
    RecordDenseShape("derivative-t",p,m);
    RecordDenseShape("derivative-coefficients",p,m);
    RecordDenseShape("derivative-correction",p,m);
    for(Eigen::Index k=0;k<2*m;++k) for(Sparse::InnerIterator entry(e.derivative,k);entry;++entry)
        raw.emplace_back(entry.row(),k/2,entry.value()*e.beta(k));
    for(Eigen::Index col=0;col<p;++col)
    {
        const auto k=free[static_cast<std::size_t>(col)]; const double norm=e.x.col(k).norm();
        if(!(norm>0)) {out.reason="zero-free-column"; return out;}
        for(Sparse::InnerIterator entry(e.x,k);entry;++entry) entries.emplace_back(entry.row(),col,entry.value()/norm);
        t(col,k/2)=e.derivative.col(k).dot(e.residual)/norm;
    }
    out.free_design.setFromTriplets(entries.begin(),entries.end()); out.raw.setFromTriplets(raw.begin(),raw.end());
    const auto cancelled=[&]() {
        constexpr double relative_roundoff=64*std::numeric_limits<double>::epsilon()/1e-10;
        for(Eigen::Index k=0;k<m;++k)
        {
            const Vector column=Vector(out.raw.col(k));
            if(column.norm()>0 && (column-out.free_design*out.coefficients.col(k)).norm()<=relative_roundoff*column.norm()) return true;
        }
        return false;
    };
    const auto reference_raw=[&]() {
        ++SparseWorkForTesting().cancellation_reductions;
        out.reference_order=true;
        // Match the independent reference's multiply-add accumulation before
        // normalization amplifies a nearly cancelled projected column.
        raw.clear();
        for(Eigen::Index a=0;a<m;++a)
        {
            Sparse::InnerIterator ga(e.derivative,2*a),ch(e.derivative,2*a+1);
            while(ga || ch)
            {
                const auto row=!ch || (ga && ga.row()<ch.row()) ? ga.row() : ch.row();
                double value=0;
                if(ga && ga.row()==row) {value=std::fma(ga.value(),e.beta(2*a),value); ++ga;}
                if(ch && ch.row()==row) {value=std::fma(ch.value(),e.beta(2*a+1),value); ++ch;}
                raw.emplace_back(row,a,value);
            }
        }
        out.raw.setFromTriplets(raw.begin(),raw.end());
    };
    if(SparseBackendEnabled())
    {
        LinearWorkspace workspace;
        const Sparse design=out.free_design;
        auto factor=e.factor;
        try {
            if(!factor || !factor->Matches(design,free)) factor=workspace.Factor(design,free,0);
            else ++SparseWorkForTesting().factor_reuses;
            Matrix compact;
            {++work.derivative_compacts; WorkTimer timer(work.derivative_compact_seconds); compact=factor->Compact();}
            const auto svd=EvaluateRank(compact,{{context ? context->rank.rows : n,2*m,m},p,absolute});
            if(!svd.valid) {out.reason="nonfinite-derivative"; return out;}
            if(svd.rank!=p) {out.reason="rank-deficient-free-design"; return out;}
            out.coefficients.resize(p,m); out.correction.resize(p,m);
            for(Eigen::Index first=0;first<m;first+=16)
            {
                const auto count=std::min<Eigen::Index>(16,m-first);
                out.coefficients.middleCols(first,count)=factor->LeastSquares(Matrix(out.raw.middleCols(first,count)));
                out.correction.middleCols(first,count)=factor->NormalSolve(t.middleCols(first,count));
            }
        } catch(const std::runtime_error &) {out.reason="sparse-derivative-failed"; return out;}
        out.valid=out.coefficients.allFinite() && out.correction.allFinite();
        if(!out.valid) {out.reason="nonfinite-derivative"; return out;}
        // Near-complete cancellation amplifies QR ordering roundoff in the
        // normalized width spectrum. Use reference-compatible tiled arithmetic in
        // that regime; this changes computation, never acceptance thresholds.
        if(!cancelled()) {out.reason="full-profile-derivative"; return out;}
        reference_raw();
        out.valid=false;
    }
    std::optional<WorkTimer> cancellation_timer;
    if(SparseBackendEnabled()) cancellation_timer.emplace(work.cancellation_seconds);
    const auto tiled_solve=[&]() {
        TiledQR qr(p,m);
        {
            ++work.derivative_compacts; WorkTimer compact_timer(work.derivative_compact_seconds);
            for(Eigen::Index first=0;first<n;first+=tile)
            {
                const auto count=std::min(tile,n-first);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                auto & tiles=DerivativeWorkForTesting();
                tiles.maximum_generated_rows=std::max(tiles.maximum_generated_rows,count);
                tiles.maximum_reduction_rows=std::max(tiles.maximum_reduction_rows,count+qr.r.rows());
#endif
                qr.Append(Matrix(out.free_design.middleRows(first,count)),Matrix(out.raw.middleRows(first,count)),out.reference_order);
            }
        }
        const auto svd=EvaluateRank(qr.r,{{context ? context->rank.rows : n,2*m,m},p,absolute});
        if(!svd.valid) {out.reason="nonfinite-derivative"; return false;}
        if(svd.rank!=p) {out.reason="rank-deficient-free-design"; return false;}
        out.coefficients=qr.r.triangularView<Eigen::Upper>().solve(qr.target);
        const Matrix adjoint=qr.r.transpose().triangularView<Eigen::Lower>().solve(t);
        out.correction=qr.r.triangularView<Eigen::Upper>().solve(adjoint);
        return true;
    };
    if(!tiled_solve()) return out;
    if(!out.reference_order && cancelled())
    {
        reference_raw();
        if(!tiled_solve()) return out;
    }
    out.valid=out.coefficients.allFinite() && out.correction.allFinite();
    out.reason=out.valid ? "full-profile-derivative" : "nonfinite-derivative"; return out;
}
void TiledDifferential::Rows(Eigen::Index first,Eigen::Index count,Matrix & projected,Matrix & jacobian) const
{
    if(reference_order)
    {
        const Sparse design=free_design.middleRows(first,count);
        projected=(Matrix(raw.middleRows(first,count))-design*coefficients)/scale;
        jacobian=projected-design*correction/scale;
        return;
    }
    projected=(Matrix(raw.middleRows(first,count))-free_design.middleRows(first,count)*coefficients)/scale;
    jacobian=projected-free_design.middleRows(first,count)*correction/scale;
}
ReducedDifferential ReduceDerivative(const TiledDifferential & d,VectorRef residual,bool widths,Eigen::Index tile)
{
    ResourcePhase phase("derivative-reduce");
    ReducedDifferential out; out.reason=d.reason; if(!d.valid) return out;
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
    const auto n=d.raw.rows(),m=d.raw.cols();
    TiledQR projected(m,0),jacobian(m,1);
    out.projected_norms=Vector::Zero(m); out.jacobian_norms=Vector::Zero(m);
    Matrix p,j;
    for(Eigen::Index first=0;first<n;first+=tile)
    {
        const auto count=std::min(tile,n-first); d.Rows(first,count,p,j);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        auto & work=DerivativeWorkForTesting(); work.maximum_generated_rows=std::max(work.maximum_generated_rows,count);
        work.maximum_reduction_rows=std::max(work.maximum_reduction_rows,count+jacobian.r.rows());
#endif
        if(!p.allFinite() || !j.allFinite()) {out.reason="nonfinite-derivative"; return out;}
        if(widths) projected.Append(p,Matrix(count,0));
        jacobian.Append(j,Matrix(residual.segment(first,count)/d.scale));
        for(Eigen::Index k=0;k<m;++k)
        {
            if(widths) out.projected_norms(k)=std::hypot(out.projected_norms(k),p.col(k).norm());
            out.jacobian_norms(k)=std::hypot(out.jacobian_norms(k),j.col(k).blueNorm());
        }
    }
    out.projected=std::move(projected.r); out.jacobian=std::move(jacobian.r);
    out.response=jacobian.target.col(0); out.valid=true; return out;
}
}
