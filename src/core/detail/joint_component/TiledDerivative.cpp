#include "TiledDerivative.hpp"
#include "SparseFactor.hpp"
#include "CompactSvd.hpp"
#include <cmath>
#include <optional>

namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
DerivativeWork & DerivativeWorkForTesting() {static thread_local DerivativeWork work; return work;}
JacobianReductionKindForTesting & JacobianReductionForTesting()
{static thread_local auto kind=JacobianReductionKindForTesting::ObservationTsqr; return kind;}
ProjectedReductionKindForTesting & ProjectedReductionForTesting()
{static thread_local auto kind=ProjectedReductionKindForTesting::ObservationTiledQr; return kind;}
namespace {
TiledQrTelemetry AccumulateTiledQr(const TiledQrTelemetry & accumulated,const TiledQrTelemetry & current)
{
    auto result=accumulated;
    if(!current.role.empty()) result.role=current.role;
    result.append_calls+=current.append_calls; result.rows_processed+=current.rows_processed;
    result.columns=current.columns; result.responses=current.responses;
    result.maximum_assembled_rows=std::max(result.maximum_assembled_rows,current.maximum_assembled_rows);
    result.maximum_dense_design_bytes=std::max(result.maximum_dense_design_bytes,current.maximum_dense_design_bytes);
    result.maximum_dense_response_bytes=std::max(result.maximum_dense_response_bytes,current.maximum_dense_response_bytes);
    result.qr_seconds+=current.qr_seconds;
    result.assembly_copy_seconds+=current.assembly_copy_seconds;
    result.householder_seconds+=current.householder_seconds;
    result.rhs_transform_seconds+=current.rhs_transform_seconds;
    return result;
}
}
#endif
namespace {
TiledDifferential PrepareDerivativeImpl(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile,bool compact_jacobian)
{
    auto & work=SparseWorkForTesting(); ++work.derivative_preparations;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & derivative_work=DerivativeWorkForTesting();
#endif
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
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        WorkTimer timer(derivative_work.preparation.raw_assembly_seconds);
#endif
        for(Eigen::Index k=0;k<2*m;++k) for(Sparse::InnerIterator entry(e.derivative,k);entry;++entry)
            raw.emplace_back(entry.row(),k/2,entry.value()*e.beta(k));
    }
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        WorkTimer timer(derivative_work.preparation.free_design_assembly_seconds);
#endif
        for(Eigen::Index col=0;col<p;++col)
        {
            const auto k=free[static_cast<std::size_t>(col)]; const double norm=e.x.col(k).norm();
            if(!(norm>0)) {out.reason="zero-free-column"; return out;}
            for(Sparse::InnerIterator entry(e.x,k);entry;++entry) entries.emplace_back(entry.row(),col,entry.value()/norm);
            t(col,k/2)=e.derivative.col(k).dot(e.residual)/norm;
        }
    }
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        WorkTimer timer(derivative_work.preparation.free_design_assembly_seconds);
#endif
        out.free_design.setFromTriplets(entries.begin(),entries.end());
    }
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        WorkTimer timer(derivative_work.preparation.raw_assembly_seconds);
#endif
        out.raw.setFromTriplets(raw.begin(),raw.end());
    }
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        WorkTimer timer(derivative_work.preparation.cancellation_fallback_seconds);
#endif
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
            bool matches{};
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                WorkTimer timer(derivative_work.preparation.factor_match_seconds);
#endif
                matches=factor && factor->Matches(design,free);
            }
            if(!matches)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                WorkTimer timer(derivative_work.preparation.factor_build_seconds);
#endif
                factor=workspace.Factor(design,free,0);
            }
            else ++SparseWorkForTesting().factor_reuses;
            Matrix compact;
            {
                ++work.derivative_compacts; WorkTimer timer(work.derivative_compact_seconds);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                WorkTimer stage(derivative_work.preparation.factor_compact_seconds);
#endif
                compact=factor->Compact();
            }
            CompactSvdResult svd;
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                WorkTimer timer(derivative_work.preparation.rank_seconds);
#endif
                svd=EvaluateRank(compact,{{context ? context->rank.rows : n,2*m,m},p,absolute});
            }
            if(!svd.valid) {out.reason="nonfinite-derivative"; return out;}
            if(svd.rank!=p) {out.reason="rank-deficient-free-design"; return out;}
            if(compact_jacobian)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                AssessmentStageTimerForTesting projection_stage("compact-design-response",n,p);
#endif
                out.free_design_factor=compact;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                out.free_design_factor_for_testing=factor;
#endif
                const Matrix response=e.residual/scale;
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    WorkTimer timer(derivative_work.preparation.least_squares_seconds);
#endif
                    out.free_design_response=(compact*factor->LeastSquares(response)).col(0);
                }
            }
            out.coefficients.resize(p,m); out.correction.resize(p,m);
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                WorkTimer timer(derivative_work.preparation.least_squares_seconds);
#endif
                for(Eigen::Index first=0;first<m;first+=16)
                {
                    const auto count=std::min<Eigen::Index>(16,m-first);
                    out.coefficients.middleCols(first,count)=factor->LeastSquares(Matrix(out.raw.middleCols(first,count)));
                }
            }
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                WorkTimer timer(derivative_work.preparation.normal_solve_seconds);
#endif
                for(Eigen::Index first=0;first<m;first+=16)
                {
                    const auto count=std::min<Eigen::Index>(16,m-first);
                    out.correction.middleCols(first,count)=factor->NormalSolve(t.middleCols(first,count));
                }
            }
        } catch(const std::runtime_error &) {out.reason="sparse-derivative-failed"; return out;}
        out.valid=out.coefficients.allFinite() && out.correction.allFinite();
        if(!out.valid) {out.reason="nonfinite-derivative"; return out;}
        // Near-complete cancellation amplifies QR ordering roundoff in the
        // normalized width spectrum. Use reference-compatible tiled arithmetic in
        // that regime; this changes computation, never acceptance thresholds.
        bool cancellation_detected{};
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            WorkTimer timer(derivative_work.preparation.cancellation_check_seconds);
#endif
            cancellation_detected=cancelled();
        }
        if(!cancellation_detected) {out.reason="full-profile-derivative"; return out;}
        reference_raw();
        out.valid=false;
    }
    std::optional<WorkTimer> cancellation_timer;
    if(SparseBackendEnabled()) cancellation_timer.emplace(work.cancellation_seconds);
    const auto tiled_solve=[&]() {
        const Eigen::Index response_count=m+(compact_jacobian ? 1 : 0);
        TiledQR qr(p,response_count);
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
                const Matrix design(out.free_design.middleRows(first,count));
                if(compact_jacobian)
                {
                    Matrix response(count,m+1); response.leftCols(m)=Matrix(out.raw.middleRows(first,count));
                    response.col(m)=e.residual.segment(first,count)/scale;
                    qr.Append(design,response,out.reference_order);
                }
                else
                    qr.Append(design,Matrix(out.raw.middleRows(first,count)),out.reference_order);
            }
        }
        const auto svd=EvaluateRank(qr.r,{{context ? context->rank.rows : n,2*m,m},p,absolute});
        if(!svd.valid) {out.reason="nonfinite-derivative"; return false;}
        if(svd.rank!=p) {out.reason="rank-deficient-free-design"; return false;}
        if(compact_jacobian)
        {
            out.free_design_factor=qr.r;
            out.free_design_response=qr.target.col(m);
        }
        out.coefficients=qr.r.triangularView<Eigen::Upper>().solve(qr.target.leftCols(m));
        const Matrix adjoint=qr.r.transpose().triangularView<Eigen::Lower>().solve(t);
        out.correction=qr.r.triangularView<Eigen::Upper>().solve(adjoint);
        return true;
    };
    if(!tiled_solve()) return out;
    bool cancellation_detected{};
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        WorkTimer timer(derivative_work.preparation.cancellation_check_seconds);
#endif
        cancellation_detected=cancelled();
    }
    if(!out.reference_order && cancellation_detected)
    {
        reference_raw();
        if(!tiled_solve()) return out;
    }
    out.valid=out.coefficients.allFinite() && out.correction.allFinite();
    out.reason=out.valid ? "full-profile-derivative" : "nonfinite-derivative"; return out;
}
}
TiledDifferential PrepareDerivative(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile)
{return PrepareDerivativeImpl(e,scale,context,absolute,tile,false);}
TiledDifferential PrepareDerivativeCompact(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile)
{return PrepareDerivativeImpl(e,scale,context,absolute,tile,true);}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
TiledDifferential PrepareDerivativeForTesting(const Evaluation & e,double scale,const EvaluationContext * context,
    double absolute,Eigen::Index tile,bool compact_jacobian)
{return compact_jacobian ? PrepareDerivativeCompact(e,scale,context,absolute,tile) :
    PrepareDerivative(e,scale,context,absolute,tile);}
#endif
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
namespace {
ReducedDifferential ReduceDerivativeImpl(const TiledDifferential & d,VectorRef residual,bool widths,Eigen::Index tile,bool compact)
{
    ResourcePhase phase("derivative-reduce");
    ReducedDifferential out; out.reason=d.reason; if(!d.valid) return out;
    if(compact && !widths) {out.reason="compact-jacobian-requires-widths"; return out;}
    if(tile<=0) throw std::invalid_argument("Invalid derivative tile size.");
    const auto n=d.raw.rows(),m=d.raw.cols();
    TiledQR projected(m,compact ? 1 : 0,"derivative-projected-qr"),jacobian(m,1,"derivative-jacobian-qr");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto & work=DerivativeWorkForTesting(); const auto projected_base=work.projected_qr;
    const auto jacobian_base=work.jacobian_qr,compact_base=work.compact_jacobian_qr;
    struct ReductionAttribution
    {
        DerivativeWork & work;
        std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
        double rows_seconds{},jacobian_qr_seconds{},norms_seconds{};
        ~ReductionAttribution()
        {
            work.rows_seconds+=rows_seconds;
            work.jacobian_qr_seconds+=jacobian_qr_seconds;
            work.norms_seconds+=norms_seconds;
            const double total=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            work.outer_overhead_seconds+=std::max(0.0,total-rows_seconds-jacobian_qr_seconds-norms_seconds);
        }
    } attribution{work};
    const bool structured=compact && widths &&
        ProjectedReductionForTesting()==ProjectedReductionKindForTesting::StructuredCompactQr;
    const bool tail_candidate=compact && widths &&
        ProjectedReductionForTesting()==ProjectedReductionKindForTesting::ProjectedTailQr;
    const bool census=compact && widths &&
        ProjectedReductionForTesting()==ProjectedReductionKindForTesting::ProjectedTailCensus;
    auto & projected_work=work.projected_reduction;
    projected_work.kind=structured ? "structured-compact-qr" : tail_candidate ? "projected-tail-qr" :
        census ? "projected-tail-census" : "observation-tiled-qr";
    projected_work.ordering=structured ? "SuiteSparseQR_FIXED" : tail_candidate ? "SuiteSparseQR_COLAMD" :
        census ? "sparse-qmult-free-design-factor" : "none";
    projected_work.factor_rows=m; projected_work.factor_columns=m;
    if(structured || tail_candidate || census) ++projected_work.attempts;
#endif
    out.projected_norms=Vector::Zero(m); out.jacobian_norms=Vector::Zero(m);
    Matrix p,j;
    for(Eigen::Index first=0;first<n;first+=tile)
    {
        const auto count=std::min(tile,n-first);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.tile_count;
        {
            AssessmentStageTimerForTesting rows_stage("derivative-rows",count,m);
            WorkTimer attribution_timer(attribution.rows_seconds);
            d.Rows(first,count,p,j);
            work.maximum_generated_rows=std::max(work.maximum_generated_rows,count);
            work.maximum_reduction_rows=std::max(work.maximum_reduction_rows,count+jacobian.r.rows());
        }
#else
        d.Rows(first,count,p,j);
#endif
        if(!p.allFinite() || !j.allFinite()) {out.reason="nonfinite-derivative"; return out;}
        if(widths)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(!structured && !tail_candidate)
            {
                AssessmentStageTimerForTesting projected_stage("derivative-projected-qr",count,m);
                if(compact)
                {
                    Matrix response(count,1); response.col(0)=residual.segment(first,count)/d.scale;
                    projected.Append(p,response);
                }
                else projected.Append(p,Matrix(count,0));
                work.projected_qr=AccumulateTiledQr(projected_base,projected.telemetry);
                projected_work.observation_projected_rows_processed+=static_cast<std::size_t>(count);
                projected_work.maximum_dense_bytes=std::max(projected_work.maximum_dense_bytes,
                    projected.telemetry.maximum_dense_design_bytes);
            }
#else
            if(compact)
            {
                Matrix response(count,1); response.col(0)=residual.segment(first,count)/d.scale;
                projected.Append(p,response);
            }
            else projected.Append(p,Matrix(count,0));
#endif
        }
        if(!compact)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            AssessmentStageTimerForTesting jacobian_stage("derivative-jacobian-qr",count,m);
            WorkTimer attribution_timer(attribution.jacobian_qr_seconds);
#endif
            jacobian.Append(j,Matrix(residual.segment(first,count)/d.scale));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            work.jacobian_qr=AccumulateTiledQr(jacobian_base,jacobian.telemetry);
#endif
        }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        {
            AssessmentStageTimerForTesting norms_stage("derivative-norms",count,m);
            WorkTimer attribution_timer(attribution.norms_seconds);
#endif
            for(Eigen::Index k=0;k<m;++k)
            {
                if(widths) out.projected_norms(k)=std::hypot(out.projected_norms(k),p.col(k).norm());
                out.jacobian_norms(k)=std::hypot(out.jacobian_norms(k),j.col(k).blueNorm());
            }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        }
#endif
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(structured)
    {
        StructuredProjectedQrResultForTesting candidate;
        {
            AssessmentStageTimerForTesting projected_stage("derivative-projected-qr",n,m);
            const auto started=std::chrono::steady_clock::now();
            candidate=StructuredProjectedQrForTesting(d.free_design,d.raw,residual,d.scale);
            projected_work.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            projected_stage.Finish();
        }
        projected_work.symbolic_seconds+=candidate.symbolic_seconds;
        projected_work.numeric_seconds+=candidate.numeric_seconds;
        projected_work.sparse_rows=candidate.sparse_rows;
        projected_work.sparse_columns=candidate.sparse_columns;
        projected_work.sparse_nonzeros+=candidate.sparse_nonzeros;
        projected_work.factor_rows=candidate.factor_rows;
        projected_work.factor_columns=candidate.factor_columns;
        projected_work.compact_projected_rows_processed+=static_cast<std::size_t>(candidate.factor_rows);
        projected_work.maximum_dense_bytes=std::max(projected_work.maximum_dense_bytes,candidate.maximum_dense_bytes);
        if(candidate.valid)
        {
            projected.r=std::move(candidate.factor);
            projected.target=candidate.response;
            out.projected_candidate=true;
        }
        else
        {
            projected_work.fallback_reason=candidate.reason;
            ++projected_work.fallbacks;
            for(Eigen::Index first=0;first<n;first+=tile)
            {
                const auto count=std::min(tile,n-first);
                d.Rows(first,count,p,j);
                Matrix response(count,1); response.col(0)=residual.segment(first,count)/d.scale;
                projected.Append(p,response);
                projected_work.observation_projected_rows_processed+=static_cast<std::size_t>(count);
            }
        }
    }
    else if(tail_candidate)
    {
        ProjectedTailQrResultForTesting candidate;
        {
            AssessmentStageTimerForTesting projected_stage("derivative-projected-qr",n,m);
            const auto started=std::chrono::steady_clock::now();
            if(!d.free_design_factor_for_testing)
                candidate.reason="free-design-factor-unavailable";
            else candidate=d.free_design_factor_for_testing->ProjectedTailQrForTesting(d.raw,residual,d.scale);
            projected_work.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            projected_stage.Finish();
        }
        projected_work.q_transform_seconds+=candidate.q_transform_seconds;
        projected_work.tail_extract_seconds+=candidate.tail_extract_seconds;
        projected_work.symbolic_seconds+=candidate.tail_symbolic_seconds;
        projected_work.numeric_seconds+=candidate.tail_numeric_seconds;
        projected_work.tail_qmult_seconds+=candidate.tail_qmult_seconds;
        projected_work.tail_compact_seconds+=candidate.tail_compact_seconds;
        if(!candidate.ordering.empty()) projected_work.ordering="SuiteSparseQR_"+candidate.ordering;
        projected_work.sparse_rows=n; projected_work.sparse_columns=m;
        projected_work.q_transformed_nonzeros=candidate.q_transformed_nonzeros;
        projected_work.tail_nonzeros=candidate.tail_nonzeros;
        projected_work.tail_factor_nonzeros=candidate.tail_factor_nonzeros;
        projected_work.tail_factor_storage_bytes=candidate.tail_factor_storage_bytes;
        projected_work.q_transformed_storage_bytes=candidate.q_transformed_storage_bytes;
        projected_work.tail_storage_bytes=candidate.tail_storage_bytes;
        projected_work.tail_rows=candidate.tail_rows;
        projected_work.sparse_nonzeros=candidate.tail_nonzeros;
        projected_work.observations=n;
        projected_work.free_design_columns=d.free_design.cols();
        projected_work.width_columns=m;
        projected_work.raw_nonzeros=static_cast<std::size_t>(d.raw.nonZeros());
        projected_work.raw_density=n>0 && m>0 ?
            static_cast<double>(projected_work.raw_nonzeros)/(static_cast<double>(n)*static_cast<double>(m)) : 0.;
        projected_work.tail_density=candidate.tail_rows>0 && m>0 ?
            static_cast<double>(candidate.tail_nonzeros)/(static_cast<double>(candidate.tail_rows)*static_cast<double>(m)) : 0.;
        projected_work.q_transformed_density=n>0 && m>0 ?
            static_cast<double>(candidate.q_transformed_nonzeros)/(static_cast<double>(n)*static_cast<double>(m)) : 0.;
        projected_work.factor_rows=m; projected_work.factor_columns=m;
        projected_work.compact_projected_rows_processed+=static_cast<std::size_t>(m);
        projected_work.maximum_dense_bytes=std::max(projected_work.maximum_dense_bytes,
            static_cast<std::size_t>(m)*static_cast<std::size_t>(m)*sizeof(double));
        if(candidate.valid)
        {
            projected.r=std::move(candidate.factor);
            projected.target=candidate.response;
            out.projected_candidate=true;
        }
        else
        {
            projected_work.fallback_reason=candidate.reason;
            ++projected_work.fallbacks;
            for(Eigen::Index first=0;first<n;first+=tile)
            {
                const auto count=std::min(tile,n-first);
                d.Rows(first,count,p,j);
                Matrix response(count,1); response.col(0)=residual.segment(first,count)/d.scale;
                projected.Append(p,response);
                projected_work.observation_projected_rows_processed+=static_cast<std::size_t>(count);
            }
        }
    }
    else if(census)
    {
        projected_work.observations=n;
        projected_work.free_design_columns=d.free_design.cols();
        projected_work.width_columns=m;
        projected_work.raw_nonzeros=static_cast<std::size_t>(d.raw.nonZeros());
        projected_work.raw_density=n>0 && m>0 ?
            static_cast<double>(projected_work.raw_nonzeros)/(static_cast<double>(n)*static_cast<double>(m)) : 0.;
        try
        {
            if(!d.free_design_factor_for_testing) throw std::runtime_error("free-design-factor-unavailable");
            const Sparse raw=d.raw;
            const auto transformed=d.free_design_factor_for_testing->OrthogonalTransposeTailSparseForTesting(raw);
            projected_work.sparse_rows=transformed.rows;
            projected_work.sparse_columns=transformed.columns;
            projected_work.sparse_nonzeros=transformed.transformed_nonzeros;
            projected_work.q_transformed_nonzeros=transformed.transformed_nonzeros;
            projected_work.tail_rows=transformed.tail_rows;
            projected_work.tail_nonzeros=transformed.tail_nonzeros;
            projected_work.q_transformed_storage_bytes=transformed.transformed_storage_bytes;
            projected_work.tail_storage_bytes=transformed.tail_storage_bytes;
            projected_work.q_transformed_density=n>0 && m>0 ?
                static_cast<double>(transformed.transformed_nonzeros)/(static_cast<double>(n)*static_cast<double>(m)) : 0.;
            projected_work.tail_density=transformed.tail_rows>0 && m>0 ?
                static_cast<double>(transformed.tail_nonzeros)/(static_cast<double>(transformed.tail_rows)*static_cast<double>(m)) : 0.;
            projected_work.q_transform_seconds=transformed.q_transform_seconds;
            projected_work.seconds+=transformed.q_transform_seconds;
        }
        catch(const std::exception & error)
        {
            projected_work.fallback_reason=error.what();
        }
        ++projected_work.fallbacks;
        if(projected_work.fallback_reason.empty())
            projected_work.fallback_reason="census-only-observation-tiled-fallback";
    }
    work.projected_qr=AccumulateTiledQr(projected_base,projected.telemetry);
    work.jacobian_qr=AccumulateTiledQr(jacobian_base,jacobian.telemetry);
    if(compact) out.projected_response=projected.target.col(0);
#endif
    out.projected=std::move(projected.r);
    if(compact)
    {
        const auto design_columns=d.free_design_factor.rows();
        if(d.free_design_factor.cols()!=design_columns || d.free_design_response.size()!=design_columns)
        {out.reason="compact-design-factor-unavailable"; return out;}
        Matrix stack(out.projected.rows()+design_columns,m);
        stack.topRows(out.projected.rows())=out.projected;
        stack.bottomRows(design_columns)=-(d.free_design_factor*(d.correction/d.scale));
        Vector response(stack.rows()); response.head(projected.target.rows())=projected.target.col(0);
        response.tail(design_columns)=d.free_design_response;
        TiledQR compact_jacobian(m,1,"derivative-compact-jacobian-qr");
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        {
            AssessmentStageTimerForTesting compact_stage("derivative-compact-jacobian-qr",stack.rows(),m);
            compact_jacobian.Append(stack,response);
        }
        work.compact_jacobian_qr=AccumulateTiledQr(compact_base,compact_jacobian.telemetry);
#else
        compact_jacobian.Append(stack,response);
#endif
        out.jacobian=std::move(compact_jacobian.r); out.response=compact_jacobian.target.col(0);
    }
    else
    {
        out.jacobian=std::move(jacobian.r); out.response=jacobian.target.col(0);
    }
    out.valid=true; return out;
}
}
ReducedDifferential ReduceDerivative(const TiledDifferential & d,VectorRef residual,bool widths,Eigen::Index tile)
{return ReduceDerivativeImpl(d,residual,widths,tile,false);}
ReducedDifferential ReduceDerivativeCompact(const TiledDifferential & d,VectorRef residual,bool widths,Eigen::Index tile)
{return ReduceDerivativeImpl(d,residual,widths,tile,true);}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
ReducedDifferential ReduceDerivativeForTesting(const TiledDifferential & d,VectorRef residual,bool widths,
    JacobianReductionKindForTesting kind,Eigen::Index tile)
{return kind==JacobianReductionKindForTesting::CompactStackQr ?
    ReduceDerivativeCompact(d,residual,widths,tile) : ReduceDerivative(d,residual,widths,tile);}
#endif
}
