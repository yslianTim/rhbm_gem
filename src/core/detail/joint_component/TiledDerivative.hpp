#pragma once
#include "Numerics.hpp"
#include "TiledQR.hpp"

namespace rhbm_gem::core::joint_component {
struct TiledDifferential
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> free_design,raw;
    Matrix coefficients,correction;
    std::optional<Spectrum> free_design_rank;
    double scale{};
    bool valid{},reference_order{},projected_only{};
    std::string reason;
    void Rows(Eigen::Index first,Eigen::Index count,Matrix & projected,Matrix & jacobian) const;
};
struct ReducedDifferential
{
    Matrix projected,jacobian; // Compact R factors, never observation-sized.
    Vector response,projected_norms,jacobian_norms;
    bool valid{};
    std::string reason;
};
TiledDifferential PrepareDerivative(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows);
TiledDifferential PrepareProjectedWidthDerivative(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows);
ReducedDifferential ReduceDerivative(const TiledDifferential &,VectorRef,bool=true,Eigen::Index=derivative_tile_rows);
ReducedDifferential ReduceProjectedWidthDerivative(const TiledDifferential &,Eigen::Index=derivative_tile_rows);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct DerivativeWork
{
    struct Preparation
    {
        double raw_assembly_seconds{},free_design_assembly_seconds{},factor_match_seconds{},factor_build_seconds{},
            factor_compact_seconds{},rank_seconds{},least_squares_seconds{},normal_solve_seconds{},
            cancellation_check_seconds{},cancellation_fallback_seconds{};
    } preparation;
    double rows_seconds{},jacobian_qr_seconds{},norms_seconds{},outer_overhead_seconds{};
    std::size_t free_design_qr_factorizations{},free_design_rank_assessments{};
    Eigen::Index maximum_generated_rows{},maximum_reduction_rows{};
    std::size_t tile_count{};
    TiledQrTelemetry projected_qr,jacobian_qr;
};
DerivativeWork & DerivativeWorkForTesting();
#endif
}
