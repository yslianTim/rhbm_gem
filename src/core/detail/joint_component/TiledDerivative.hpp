#pragma once
#include "Numerics.hpp"
#include "TiledQR.hpp"

namespace rhbm_gem::core::joint_component {
struct TiledDifferential
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> free_design,raw;
    Matrix coefficients,correction;
    std::optional<Spectrum> free_design_rank;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    DerivativeBreakdown breakdown;
#endif
    double scale{};
    bool valid{},reference_order{},projected_only{};
    std::string reason;
    void Rows(Eigen::Index first,Eigen::Index count,Matrix & projected,Matrix & jacobian) const;
    void Rows(Eigen::Index first,Eigen::Index count,Matrix & raw_tile,Matrix & projected,Matrix & jacobian) const;
};
struct ReducedDifferential
{
    Matrix jacobian; // Compact R factor, never observation-sized.
    Vector response,jacobian_norms;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    DerivativeBreakdown breakdown;
#endif
    bool valid{};
    std::string reason;
};
struct ReducedProjectedWidthDifferential
{
    Matrix projected; // Compact R factor, never observation-sized.
    Vector projected_norms;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    DerivativeBreakdown breakdown;
#endif
    bool valid{};
    std::string reason;
};
TiledDifferential PrepareDerivative(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows);
TiledDifferential PrepareProjectedWidthDerivative(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows);
ReducedDifferential ReduceDerivative(const TiledDifferential &,VectorRef,Eigen::Index=derivative_tile_rows);
ReducedProjectedWidthDifferential ReduceProjectedWidthDerivative(const TiledDifferential &,
    Eigen::Index=derivative_tile_rows);
}
