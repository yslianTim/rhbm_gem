#pragma once
#include "Numerics.hpp"
#include "TiledQR.hpp"

namespace rhbm_gem::core::joint_component {
struct TiledDifferential
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> free_design,raw;
    Matrix free_design_factor;
    Vector free_design_response;
    Matrix coefficients,correction;
    double scale{};
    bool valid{},reference_order{};
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
TiledDifferential PrepareDerivativeCompact(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows);
ReducedDifferential ReduceDerivative(const TiledDifferential &,VectorRef,bool=true,Eigen::Index=derivative_tile_rows);
ReducedDifferential ReduceDerivativeCompact(const TiledDifferential &,VectorRef,bool=true,
    Eigen::Index=derivative_tile_rows);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
enum class JacobianReductionKindForTesting {ObservationTsqr,CompactStackQr};
JacobianReductionKindForTesting & JacobianReductionForTesting();
TiledDifferential PrepareDerivativeForTesting(const Evaluation &,double,const EvaluationContext *,double,
    Eigen::Index,bool);
ReducedDifferential ReduceDerivativeForTesting(const TiledDifferential &,VectorRef,bool,
    JacobianReductionKindForTesting,Eigen::Index=derivative_tile_rows);
struct DerivativeWork
{
    Eigen::Index maximum_generated_rows{},maximum_reduction_rows{};
    std::size_t tile_count{};
    TiledQrTelemetry projected_qr,jacobian_qr,compact_jacobian_qr;
};
DerivativeWork & DerivativeWorkForTesting();
#endif
}
