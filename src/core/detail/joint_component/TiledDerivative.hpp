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
    Matrix jacobian; // Compact R factor, never observation-sized.
    Vector response,jacobian_norms;
    bool valid{};
    std::string reason;
};
struct ReducedProjectedWidthDifferential
{
    Matrix projected; // Compact R factor, never observation-sized.
    Vector projected_norms;
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct DerivativeWork
{
    Eigen::Index maximum_generated_rows{},maximum_reduction_rows{};
};
DerivativeWork & DerivativeWorkForTesting();
#endif
}
