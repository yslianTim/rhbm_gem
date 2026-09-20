#pragma once
#include "Numerics.hpp"
#include "TiledQR.hpp"

namespace rhbm_gem::core::joint_component {
struct TiledDifferential
{
    Eigen::SparseMatrix<double,Eigen::RowMajor> free_design,raw;
    Matrix coefficients,correction;
    double scale{};
    bool valid{};
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
ReducedDifferential ReduceDerivative(const TiledDifferential &,VectorRef,bool=true,Eigen::Index=derivative_tile_rows);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct DerivativeWork {Eigen::Index maximum_generated_rows{},maximum_reduction_rows{};};
DerivativeWork & DerivativeWorkForTesting();
#endif
}
