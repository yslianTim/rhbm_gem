#include "TiledDerivative.hpp"
#include "CompactSvd.hpp"
#include <chrono>
#include <cmath>
#include <optional>
#include <utility>

namespace rhbm_gem::core::joint_component {
namespace {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
#endif
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
    TiledDifferential out;
    out.projected_only=projected_only;
    if(!e.valid || !(scale>0) || !std::isfinite(scale)) {out.reason="invalid-inner"; return out;}
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool collect_breakdown=HotspotAttributionEnabled();
    const auto raw_construction_started=collect_breakdown ? std::chrono::steady_clock::now() :
        std::chrono::steady_clock::time_point{};
#endif
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(collect_breakdown) out.breakdown.workspace_resize_events+=2;
#endif
    Vector free_column_norms(projected_only ? p : 0);
    std::vector<Eigen::Triplet<double>> entries,raw;
    entries.reserve(free_design_nonzeros);
    raw.reserve(static_cast<std::size_t>(e.derivative.nonZeros()));
    Matrix t;
    if(!projected_only)
    {
        t=Matrix::Zero(p,m);
    }
    for(Eigen::Index k=0;k<2*m;++k) for(Sparse::InnerIterator entry(e.derivative,k);entry;++entry)
    {
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
            entries.emplace_back(entry.row(),col,entry.value()/norm);
        }
        if(!projected_only) t(col,k/2)=e.derivative.col(k).dot(e.residual)/norm;
    }
    out.free_design.setFromTriplets(entries.begin(),entries.end());
    out.raw.setFromTriplets(raw.begin(),raw.end());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(collect_breakdown)
    {
        out.breakdown.workspace_resize_events+=2;
        out.breakdown.raw_construction_seconds+=Seconds(raw_construction_started);
    }
#endif
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
    // EIGEN's qualified derivative path keeps the observation-sized work
    // bounded by tiled Householder reductions rather than materializing a
    // sparse Q action.
    const auto tiled_solve=[&]() {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto tiled_qr_started=collect_breakdown ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
        double row_census_seconds{};
        if(collect_breakdown) ++out.breakdown.tiled_qr_calls;
#endif
        TiledQR qr(p,m);
        {
            for(Eigen::Index first=0;first<n;first+=tile)
            {
                const auto count=std::min(tile,n-first);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                TiledQrTileCensus census;
                census.stage=TiledQrStage::DerivativePrepare;
                census.input_rows=static_cast<std::size_t>(count);
                census.factor_rows=static_cast<std::size_t>(qr.r.rows()+count);
                census.columns=static_cast<std::size_t>(p);
                census.rhs_columns=static_cast<std::size_t>(m);
                census.retained_rows=static_cast<std::size_t>(std::min<Eigen::Index>(
                    qr.r.rows()+count,p));
                if(collect_breakdown)
                {
                    const auto census_started=std::chrono::steady_clock::now();
                    for(Eigen::Index row=first;row<first+count;++row)
                    {
                        bool design_structural{},derivative_structural{},design_nonzero{},derivative_nonzero{};
                        for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator entry(out.free_design,row);entry;++entry)
                        {
                            ++census.design_structural_nonzeros;
                            ++census.structural_nonzeros;
                            design_structural=true;
                            if(entry.value()!=0.0)
                            {
                                ++census.design_exact_nonzeros;
                                ++census.exact_nonzeros;
                                design_nonzero=true;
                            }
                        }
                        for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator entry(out.raw,row);entry;++entry)
                        {
                            ++census.derivative_structural_nonzeros;
                            ++census.structural_nonzeros;
                            derivative_structural=true;
                            if(entry.value()!=0.0)
                            {
                                ++census.derivative_exact_nonzeros;
                                ++census.exact_nonzeros;
                                derivative_nonzero=true;
                            }
                        }
                        census.structural_rows+=design_structural || derivative_structural;
                        census.design_nonzero_rows+=design_nonzero;
                        census.derivative_nonzero_rows+=derivative_nonzero;
                        census.both_zero_rows+=!design_nonzero && !derivative_nonzero;
                        const bool active=design_nonzero || derivative_nonzero;
                        census.nonzero_design_rows+=active;
                        census.exact_zero_rows+=!active;
                        if(!active && e.residual(row)!=0.0) ++census.zero_design_nonzero_rhs_rows;
                    }
                    census.row_census_seconds=std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-census_started).count();
                    row_census_seconds+=census.row_census_seconds;
                }
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                const auto materialization_started=collect_breakdown ? std::chrono::steady_clock::now() :
                    std::chrono::steady_clock::time_point{};
#endif
                const Matrix design(out.free_design.middleRows(first,count));
                const Matrix raw_tile(out.raw.middleRows(first,count));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(collect_breakdown)
                {
                    census.dense_tile_materialization_seconds=std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-materialization_started).count();
                    const auto element_count=static_cast<std::size_t>(count)*static_cast<std::size_t>(p+m)+
                        census.factor_rows*static_cast<std::size_t>(p+m)+census.retained_rows+
                        (out.reference_order ? census.factor_rows*static_cast<std::size_t>(p+m) : 0);
                    census.estimated_peak_temporary_bytes=element_count*sizeof(double);
                }
                TiledQRAppendTiming append_timing;
                qr.Append(design,raw_tile,out.reference_order,collect_breakdown ? &append_timing : nullptr);
                if(collect_breakdown)
                {
                    census.dense_assembly_seconds=append_timing.dense_assembly_seconds;
                    census.householder_factorization_seconds=append_timing.householder_factorization_seconds;
                    census.qtrhs_transform_seconds=append_timing.qtrhs_transform_seconds;
                    census.r_extraction_seconds=append_timing.r_extraction_seconds;
                    AddTiledQrTileCensus(out.breakdown.tiled_qr_tiles,census);
                }
#else
                qr.Append(design,raw_tile,out.reference_order);
#endif
            }
        }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(collect_breakdown)
        {
            out.breakdown.tiled_qr_seconds+=Seconds(tiled_qr_started)-row_census_seconds;
        }
        const auto svd_started=collect_breakdown ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
        if(collect_breakdown) ++out.breakdown.rank_svd_calls;
#endif
        const auto boundary=projected_only ? RankBoundary::StrictGreater : RankBoundary::SvdNative;
        const auto svd=EvaluateRank(qr.r,{{context ? context->rank.rows : n,2*m,m},p,absolute,boundary},
            projected_only ? CompactSvdVectors::Right : CompactSvdVectors::None);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(collect_breakdown) out.breakdown.rank_svd_seconds+=Seconds(svd_started);
#endif
        if(projected_only) out.free_design_rank=FreeDesignSpectrum(svd,
            context ? context->rank.rows : n,p,2*m,free_column_norms);
        if(!svd.valid) {out.reason="nonfinite-derivative"; return false;}
        if(svd.rank!=p && !projected_only)
        {out.reason="rank-deficient-free-design"; return false;}
        if(projected_only)
        {
            if(svd.rank==p)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                const auto solve_started=collect_breakdown ? std::chrono::steady_clock::now() :
                    std::chrono::steady_clock::time_point{};
                if(collect_breakdown) ++out.breakdown.triangular_solve_calls;
#endif
                qr.r.triangularView<Eigen::Upper>().solveInPlace(qr.target);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(collect_breakdown) out.breakdown.triangular_solve_seconds+=Seconds(solve_started);
#endif
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto coefficient_started=collect_breakdown ? std::chrono::steady_clock::now() :
                std::chrono::steady_clock::time_point{};
            if(collect_breakdown) ++out.breakdown.triangular_solve_calls;
#endif
            out.coefficients=qr.r.triangularView<Eigen::Upper>().solve(qr.target.leftCols(m));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(collect_breakdown) out.breakdown.triangular_solve_seconds+=Seconds(coefficient_started);
            const auto correction_started=collect_breakdown ? std::chrono::steady_clock::now() :
                std::chrono::steady_clock::time_point{};
            if(collect_breakdown) ++out.breakdown.triangular_solve_calls;
#endif
            const Matrix adjoint=qr.r.transpose().triangularView<Eigen::Lower>().solve(t);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(collect_breakdown) out.breakdown.triangular_solve_seconds+=Seconds(correction_started);
            const auto correction_solve_started=collect_breakdown ? std::chrono::steady_clock::now() :
                std::chrono::steady_clock::time_point{};
            if(collect_breakdown) ++out.breakdown.triangular_solve_calls;
#endif
            out.correction=qr.r.triangularView<Eigen::Upper>().solve(adjoint);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(collect_breakdown)
            {
                out.breakdown.triangular_solve_seconds+=Seconds(correction_solve_started);
                out.breakdown.correction_matrix_seconds+=Seconds(correction_started);
            }
#endif
        }
        return true;
    };
    if(!tiled_solve()) return out;
    const bool cancellation_detected=cancelled();
    if(!out.reference_order && cancellation_detected)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto fallback_started=collect_breakdown ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
        if(collect_breakdown) ++out.breakdown.cancellation_fallbacks;
#endif
        reference_raw();
        const bool fallback_valid=tiled_solve();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(collect_breakdown) out.breakdown.cancellation_fallback_seconds+=Seconds(fallback_started);
#endif
        if(!fallback_valid) return out;
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
void TiledDifferential::Rows(Eigen::Index first,Eigen::Index count,Matrix & raw_tile,
    Matrix & projected,Matrix & jacobian) const
{
    raw_tile=raw.middleRows(first,count);
    if(reference_order)
    {
        const Sparse design=free_design.middleRows(first,count);
        projected=(raw_tile-design*coefficients)/scale;
        if(projected_only) {jacobian.resize(count,0); return;}
        jacobian=projected-design*correction/scale;
        return;
    }
    projected=(raw_tile-free_design.middleRows(first,count)*coefficients)/scale;
    if(projected_only) {jacobian.resize(count,0); return;}
    jacobian=projected-free_design.middleRows(first,count)*correction/scale;
}
void TiledDifferential::Rows(Eigen::Index first,Eigen::Index count,Matrix & projected,Matrix & jacobian) const
{
    Matrix raw_tile;
    Rows(first,count,raw_tile,projected,jacobian);
}
namespace {
ReducedDifferential ReduceDerivativeImpl(const TiledDifferential & d,VectorRef residual,Eigen::Index tile)
{
    ReducedDifferential out; out.reason=d.reason; if(!d.valid) return out;
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool collect_breakdown=HotspotAttributionEnabled();
#endif
    const auto n=d.raw.rows(),m=d.raw.cols();
    TiledQR jacobian(m,1);
    out.jacobian_norms=Vector::Zero(m);
    Matrix raw_tile,projected,j,residual_tile;
    for(Eigen::Index first=0;first<n;first+=tile)
    {
        const auto count=std::min(tile,n-first);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto corrected_started=collect_breakdown ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
#endif
        d.Rows(first,count,raw_tile,projected,j);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        TiledQrTileCensus census;
        census.stage=TiledQrStage::DerivativeReduce;
        census.input_rows=static_cast<std::size_t>(count);
        census.factor_rows=static_cast<std::size_t>(jacobian.r.rows()+count);
        census.columns=static_cast<std::size_t>(m);
        census.rhs_columns=1;
        census.retained_rows=static_cast<std::size_t>(std::min<Eigen::Index>(jacobian.r.rows()+count,m));
        if(collect_breakdown)
        {
            census.corrected_jacobian_seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-corrected_started).count();
            out.breakdown.corrected_jacobian_seconds+=census.corrected_jacobian_seconds;
            for(Eigen::Index row=0;row<count;++row)
            {
                ++census.structural_rows;
                bool corrected_nonzero{},projection_nonzero{};
                for(Eigen::Index col=0;col<m;++col)
                {
                    ++census.structural_nonzeros;
                    corrected_nonzero|=j(row,col)!=0.0;
                    projection_nonzero|=projected(row,col)!=0.0;
                    census.corrected_jacobian_nonzeros+=j(row,col)!=0.0;
                    census.exact_nonzeros+=j(row,col)!=0.0;
                }
                census.corrected_jacobian_nonzero_rows+=corrected_nonzero;
                census.nonzero_design_rows+=corrected_nonzero;
                census.exact_zero_rows+=!corrected_nonzero;
                if(!corrected_nonzero && residual(first+row)!=0.0) ++census.zero_design_nonzero_rhs_rows;
                census.projection_nonzero_rows+=projection_nonzero;
                census.rhs_nonzero_rows+=residual(first+row)!=0.0;
            }
        }
#endif
        if(residual_tile.rows()!=count || residual_tile.cols()!=1)
        {
            residual_tile.resize(count,1);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(collect_breakdown) ++out.breakdown.workspace_resize_events;
#endif
        }
        residual_tile=residual.segment(first,count)/d.scale;
        if(!projected.allFinite() || !j.allFinite()) {out.reason="nonfinite-derivative"; return out;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto append_started=collect_breakdown ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        TiledQRAppendTiming append_timing;
        jacobian.Append(j,residual_tile,false,collect_breakdown ? &append_timing : nullptr);
        if(collect_breakdown)
        {
            census.dense_assembly_seconds=append_timing.dense_assembly_seconds;
            census.householder_factorization_seconds=append_timing.householder_factorization_seconds;
            census.qtrhs_transform_seconds=append_timing.qtrhs_transform_seconds;
            census.r_extraction_seconds=append_timing.r_extraction_seconds;
            census.estimated_peak_temporary_bytes=(static_cast<std::size_t>(count)*
                static_cast<std::size_t>(m+1)+census.factor_rows*static_cast<std::size_t>(m+1)+
                census.retained_rows)*sizeof(double);
            AddTiledQrTileCensus(out.breakdown.tiled_qr_tiles,census);
        }
#else
        jacobian.Append(j,residual_tile);
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(collect_breakdown) out.breakdown.tiled_qr_seconds+=Seconds(append_started);
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
    ReducedProjectedWidthDifferential out; out.reason=d.reason; if(!d.valid || !d.projected_only) return out;
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool collect_breakdown=HotspotAttributionEnabled();
#endif
    const auto n=d.raw.rows(),m=d.raw.cols();
    TiledQR projected(m,0);
    out.projected_norms=Vector::Zero(m); Matrix raw_tile,p,j;
    for(Eigen::Index first=0;first<n;first+=tile)
    {
        const auto count=std::min(tile,n-first);
        d.Rows(first,count,raw_tile,p,j);
        if(!p.allFinite()) {out.reason="nonfinite-projected-width-derivative"; return out;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        TiledQrTileCensus census;
        census.stage=TiledQrStage::ProjectedWidth;
        census.input_rows=static_cast<std::size_t>(count);
        census.factor_rows=static_cast<std::size_t>(projected.r.rows()+count);
        census.columns=static_cast<std::size_t>(m);
        census.retained_rows=static_cast<std::size_t>(std::min<Eigen::Index>(
            projected.r.rows()+count,m));
        if(collect_breakdown)
        {
            census.structural_rows=static_cast<std::size_t>(count);
            for(Eigen::Index row=0;row<count;++row)
            {
                bool active{};
                for(Eigen::Index col=0;col<m;++col)
                {
                    const bool nonzero=p(row,col)!=0.0;
                    ++census.structural_nonzeros;
                    active|=nonzero;
                    census.projection_nonzeros+=nonzero;
                    census.exact_nonzeros+=nonzero;
                }
                census.nonzero_design_rows+=active;
                census.exact_zero_rows+=!active;
                census.projection_nonzero_rows+=active;
            }
        }
        TiledQRAppendTiming append_timing;
        projected.Append(p,Matrix(count,0),false,collect_breakdown ? &append_timing : nullptr);
        if(collect_breakdown)
        {
            census.dense_assembly_seconds=append_timing.dense_assembly_seconds;
            census.householder_factorization_seconds=append_timing.householder_factorization_seconds;
            census.qtrhs_transform_seconds=append_timing.qtrhs_transform_seconds;
            census.r_extraction_seconds=append_timing.r_extraction_seconds;
            census.estimated_peak_temporary_bytes=(static_cast<std::size_t>(count)*
                static_cast<std::size_t>(m)+census.factor_rows*static_cast<std::size_t>(m)+
                census.retained_rows)*sizeof(double);
            AddTiledQrTileCensus(out.breakdown.tiled_qr_tiles,census);
            RecordProjectedWidthTileCensus(census);
        }
#else
        projected.Append(p,Matrix(count,0));
#endif
        for(Eigen::Index k=0;k<m;++k)
            out.projected_norms(k)=std::hypot(out.projected_norms(k),p.col(k).norm());
    }
    out.projected=std::move(projected.r); out.valid=true; return out;
}
}
