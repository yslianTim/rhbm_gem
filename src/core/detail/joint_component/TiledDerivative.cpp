#include "TiledDerivative.hpp"
#include "CompactSvd.hpp"
#include <cmath>
#include <optional>
#include <utility>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <chrono>
#endif

namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
DerivativeWork & DerivativeWorkForTesting() {static thread_local DerivativeWork work; return work;}
namespace {
struct DerivativeTimer
{
    double & seconds;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
    ~DerivativeTimer()
    {seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
};
}
#endif
namespace {
Spectrum FreeDesignSpectrum(const CompactSvdResult & svd,Eigen::Index rows,Eigen::Index columns,
    Eigen::Index full_columns,Vector norms)
{
    Spectrum out; out.rows=rows; out.columns=columns; out.full_design_columns=full_columns;
    out.column_norms=std::move(norms);
    if(!svd.valid)
    {out.available=false; out.reason="spectrum-factorization-failed"; return out;}
    out.rank=svd.rank; out.singular_values=svd.singular_values; out.threshold=svd.threshold;
    out.minimum=out.singular_values(out.singular_values.size()-1);
    out.condition=out.singular_values(0)/out.minimum;
    out.threshold_sensitive=!CompactRankDecisionSafe(svd);
    return out;
}
TiledDifferential PrepareDerivativeImpl(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile,bool projected_only)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=DerivativeWorkForTesting();
    ++work.prepare_calls;
    DerivativeTimer timer{work.prepare_seconds};
#endif
    TiledDifferential out;
    out.projected_only=projected_only;
    if(!e.valid || !(scale>0) || !std::isfinite(scale)) {out.reason="invalid-inner"; return out;}
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
    std::vector<Eigen::Index> free;
    std::size_t free_design_nonzeros{};
    const auto * x_outer=e.x.outerIndexPtr();
    const auto * x_inner_nonzeros=e.x.innerNonZeroPtr();
    for(Eigen::Index k=0;k<e.beta.size();++k) if(k%2 || e.beta(k)>0)
    {
        free.push_back(k);
        const auto count=x_inner_nonzeros ? x_inner_nonzeros[k] : x_outer[k+1]-x_outer[k];
        free_design_nonzeros+=static_cast<std::size_t>(count);
    }
    const Eigen::Index n=e.x.rows(),m=e.beta.size()/2,p=static_cast<Eigen::Index>(free.size());
    out.scale=scale; out.free_design.resize(n,p); out.raw.resize(n,m);
    Vector free_column_norms(projected_only ? p : 0);
    std::vector<Eigen::Triplet<double>> entries,raw;
    entries.reserve(free_design_nonzeros);
    raw.reserve(static_cast<std::size_t>(e.derivative.nonZeros()));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    work.triplet_reserved_entries+=entries.capacity()+raw.capacity();
    work.triplet_reserved_bytes+=(entries.capacity()+raw.capacity())*sizeof(Eigen::Triplet<double>);
#endif
    Matrix t;
    if(!projected_only)
    {
        t=Matrix::Zero(p,m);
    }
    for(Eigen::Index k=0;k<2*m;++k) for(Sparse::InnerIterator entry(e.derivative,k);entry;++entry)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(raw.size()==raw.capacity()) ++work.triplet_capacity_growths;
#endif
        raw.emplace_back(entry.row(),k/2,entry.value()*e.beta(k));
    }
    for(Eigen::Index col=0;col<p;++col)
    {
        const auto k=free[static_cast<std::size_t>(col)]; const double norm=e.x.col(k).norm();
        if(!(norm>0))
        {
            if(!projected_only)
            {
                out.reason="zero-free-column"; return out;
            }
        }
        if(projected_only) free_column_norms(col)=norm;
        if(norm>0) for(Sparse::InnerIterator entry(e.x,k);entry;++entry)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(entries.size()==entries.capacity()) ++work.triplet_capacity_growths;
#endif
            entries.emplace_back(entry.row(),col,entry.value()/norm);
        }
        if(!projected_only) t(col,k/2)=e.derivative.col(k).dot(e.residual)/norm;
    }
    out.free_design.setFromTriplets(entries.begin(),entries.end());
    out.raw.setFromTriplets(raw.begin(),raw.end());
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
        out.reference_order=true;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.reference_order_recomputations;
#endif
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(raw.size()==raw.capacity()) ++work.triplet_capacity_growths;
#endif
                raw.emplace_back(row,a,value);
            }
        }
        out.raw.setFromTriplets(raw.begin(),raw.end());
    };
    // EIGEN's qualified derivative path keeps the observation-sized work
    // bounded by tiled Householder reductions rather than materializing a
    // sparse Q action.
    const auto tiled_solve=[&]() {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.tiled_qr_invocations;
#endif
        TiledQR qr(p,m);
        {
            for(Eigen::Index first=0;first<n;first+=tile)
            {
                const auto count=std::min(tile,n-first);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                work.maximum_generated_rows=std::max(work.maximum_generated_rows,count);
                work.maximum_reduction_rows=std::max(work.maximum_reduction_rows,count+qr.r.rows());
                work.maximum_tile_workspace_bytes=std::max(work.maximum_tile_workspace_bytes,
                    static_cast<std::size_t>((count+qr.r.rows())*(p+m))*sizeof(double));
                ++work.tile_count;
#endif
                const Matrix design(out.free_design.middleRows(first,count));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                const auto qr_started=std::chrono::steady_clock::now();
#endif
                qr.Append(design,Matrix(out.raw.middleRows(first,count)),out.reference_order);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                work.tiled_qr_seconds+=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-qr_started).count();
#endif
            }
        }
        const auto boundary=projected_only ? RankBoundary::StrictGreater : RankBoundary::SvdNative;
        const auto svd=EvaluateRank(qr.r,{{context ? context->rank.rows : n,2*m,m},p,absolute,boundary},
            projected_only ? CompactSvdVectors::Right : CompactSvdVectors::None);
        if(projected_only) out.free_design_rank=FreeDesignSpectrum(svd,
            context ? context->rank.rows : n,p,2*m,free_column_norms);
        if(!svd.valid) {out.reason="nonfinite-derivative"; return false;}
        if(svd.rank!=p && !projected_only)
        {out.reason="rank-deficient-free-design"; return false;}
        if(projected_only)
        {
            if(svd.rank==p)
            {
                qr.r.triangularView<Eigen::Upper>().solveInPlace(qr.target);
                out.coefficients=std::move(qr.target);
            }
            else
            {
                out.coefficients=Matrix::Zero(p,m);
                if(svd.rank)
                {
                    const Matrix right=svd.right_vectors.leftCols(svd.rank);
                    Matrix left=qr.r*right;
                    for(Eigen::Index k=0;k<svd.rank;++k) left.col(k)/=svd.singular_values(k);
                    Matrix response=left.transpose()*qr.target;
                    for(Eigen::Index k=0;k<svd.rank;++k) response.row(k)/=svd.singular_values(k);
                    out.coefficients=right*response;
                }
            }
        }
        else
        {
            out.coefficients=qr.r.triangularView<Eigen::Upper>().solve(qr.target.leftCols(m));
            const Matrix adjoint=qr.r.transpose().triangularView<Eigen::Lower>().solve(t);
            out.correction=qr.r.triangularView<Eigen::Upper>().solve(adjoint);
        }
        return true;
    };
    if(!tiled_solve()) return out;
    const bool cancellation_detected=cancelled();
    if(!out.reference_order && cancellation_detected)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.cancellation_fallbacks;
#endif
        reference_raw();
        if(!tiled_solve()) return out;
    }
    out.valid=out.coefficients.allFinite() && (projected_only || out.correction.allFinite());
    out.reason=out.valid ? projected_only ? "projected-width-derivative" : "full-profile-derivative" :
        "nonfinite-derivative";
    return out;
}
}
TiledDifferential PrepareDerivative(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile)
{return PrepareDerivativeImpl(e,scale,context,absolute,tile,false);}
TiledDifferential PrepareProjectedWidthDerivative(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile)
{return PrepareDerivativeImpl(e,scale,context,absolute,tile,true);}
void TiledDifferential::Rows(Eigen::Index first,Eigen::Index count,Matrix & projected,Matrix & jacobian) const
{
    if(reference_order)
    {
        const Sparse design=free_design.middleRows(first,count);
        projected=(Matrix(raw.middleRows(first,count))-design*coefficients)/scale;
        if(projected_only) {jacobian.resize(count,0); return;}
        jacobian=projected-design*correction/scale;
        return;
    }
    projected=(Matrix(raw.middleRows(first,count))-free_design.middleRows(first,count)*coefficients)/scale;
    if(projected_only) {jacobian.resize(count,0); return;}
    jacobian=projected-free_design.middleRows(first,count)*correction/scale;
}
namespace {
ReducedDifferential ReduceDerivativeImpl(const TiledDifferential & d,VectorRef residual,Eigen::Index tile)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=DerivativeWorkForTesting();
    ++work.reduction_calls;
    DerivativeTimer timer{work.reduction_seconds};
#endif
    ReducedDifferential out; out.reason=d.reason; if(!d.valid) return out;
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
    const auto n=d.raw.rows(),m=d.raw.cols();
    TiledQR jacobian(m,1);
    out.jacobian_norms=Vector::Zero(m);
    Matrix projected,j;
    for(Eigen::Index first=0;first<n;first+=tile)
    {
        const auto count=std::min(tile,n-first);
        d.Rows(first,count,projected,j);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        work.maximum_generated_rows=std::max(work.maximum_generated_rows,count);
        work.maximum_reduction_rows=std::max(work.maximum_reduction_rows,count+jacobian.r.rows());
        work.maximum_tile_workspace_bytes=std::max(work.maximum_tile_workspace_bytes,
            static_cast<std::size_t>((count+jacobian.r.rows())*(m+1))*sizeof(double));
        ++work.tile_count; ++work.tiled_qr_invocations;
#endif
        if(!projected.allFinite() || !j.allFinite()) {out.reason="nonfinite-derivative"; return out;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto qr_started=std::chrono::steady_clock::now();
#endif
        jacobian.Append(j,Matrix(residual.segment(first,count)/d.scale));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        work.tiled_qr_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-qr_started).count();
#endif
        for(Eigen::Index k=0;k<m;++k)
            out.jacobian_norms(k)=std::hypot(out.jacobian_norms(k),j.col(k).blueNorm());
    }
    out.jacobian=std::move(jacobian.r); out.response=jacobian.target.col(0);
    out.valid=true; return out;
}
}
ReducedDifferential ReduceDerivative(const TiledDifferential & d,VectorRef residual,Eigen::Index tile)
{return ReduceDerivativeImpl(d,residual,tile);}
ReducedProjectedWidthDifferential ReduceProjectedWidthDerivative(const TiledDifferential & d,Eigen::Index tile)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=DerivativeWorkForTesting();
    ++work.reduction_calls;
    DerivativeTimer timer{work.reduction_seconds};
#endif
    ReducedProjectedWidthDifferential out; out.reason=d.reason; if(!d.valid || !d.projected_only) return out;
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
    const auto n=d.raw.rows(),m=d.raw.cols();
    TiledQR projected(m,0);
    out.projected_norms=Vector::Zero(m); Matrix p,j;
    for(Eigen::Index first=0;first<n;first+=tile)
    {
        const auto count=std::min(tile,n-first);
        d.Rows(first,count,p,j);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        work.maximum_generated_rows=std::max(work.maximum_generated_rows,count);
        work.maximum_reduction_rows=std::max(work.maximum_reduction_rows,count+projected.r.rows());
        work.maximum_tile_workspace_bytes=std::max(work.maximum_tile_workspace_bytes,
            static_cast<std::size_t>((count+projected.r.rows())*m)*sizeof(double));
        ++work.tile_count; ++work.tiled_qr_invocations;
#endif
        if(!p.allFinite()) {out.reason="nonfinite-projected-width-derivative"; return out;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto qr_started=std::chrono::steady_clock::now();
#endif
        projected.Append(p,Matrix(count,0));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        work.tiled_qr_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-qr_started).count();
#endif
        for(Eigen::Index k=0;k<m;++k)
            out.projected_norms(k)=std::hypot(out.projected_norms(k),p.col(k).norm());
    }
    out.projected=std::move(projected.r); out.valid=true; return out;
}
}
