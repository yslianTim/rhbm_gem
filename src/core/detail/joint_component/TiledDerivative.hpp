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
    void Rows(Eigen::Index first,Eigen::Index count,Matrix & raw_tile,Matrix & projected,Matrix & jacobian) const;
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
struct DerivativeWorkspace
{
    std::vector<Eigen::Index> free_indices;
    std::vector<Eigen::Triplet<double>> free_design_entries,raw_entries;
    Matrix t,raw_tile,projected_tile,jacobian_tile,residual_tile;
    Vector free_column_norms;
    TiledQR qr;
};
TiledDifferential PrepareDerivative(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows,DerivativeWorkspace * =nullptr);
TiledDifferential PrepareProjectedWidthDerivative(const Evaluation &,double,const EvaluationContext *,double=-1,
    Eigen::Index=derivative_tile_rows);
ReducedDifferential ReduceDerivative(const TiledDifferential &,VectorRef,Eigen::Index=derivative_tile_rows,
    DerivativeWorkspace * =nullptr);
ReducedProjectedWidthDifferential ReduceProjectedWidthDerivative(const TiledDifferential &,
    Eigen::Index=derivative_tile_rows);
}
