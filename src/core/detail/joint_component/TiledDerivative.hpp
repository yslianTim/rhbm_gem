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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::shared_ptr<FreeDesignFactor> free_design_factor_for_testing;
#endif
    void Rows(Eigen::Index first,Eigen::Index count,Matrix & projected,Matrix & jacobian) const;
};
struct ReducedDifferential
{
    Matrix projected,jacobian; // Compact R factors, never observation-sized.
    Vector response,projected_norms,jacobian_norms;
    bool valid{};
    std::string reason;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    bool projected_candidate{};
    Vector projected_response;
#endif
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
enum class ProjectedReductionKindForTesting {ObservationTiledQr,StructuredCompactQr,ProjectedTailCensus,ProjectedTailQr};
ProjectedReductionKindForTesting & ProjectedReductionForTesting();
TiledDifferential PrepareDerivativeForTesting(const Evaluation &,double,const EvaluationContext *,double,
    Eigen::Index,bool);
ReducedDifferential ReduceDerivativeForTesting(const TiledDifferential &,VectorRef,bool,
    JacobianReductionKindForTesting,Eigen::Index=derivative_tile_rows);
struct ProjectedReductionWorkForTesting
{
    std::string kind{"observation-tiled-qr"},ordering{"none"},fallback_reason;
    std::size_t attempts{},accepted{},fallbacks{},sparse_nonzeros{};
    std::size_t raw_nonzeros{},q_transformed_nonzeros{},tail_nonzeros{},
        q_transformed_storage_bytes{},tail_storage_bytes{},tail_factor_nonzeros{},tail_factor_storage_bytes{};
    std::size_t observation_projected_rows_processed{},compact_projected_rows_processed{},maximum_dense_bytes{};
    Eigen::Index observations{},free_design_columns{},width_columns{},tail_rows{},
        sparse_rows{},sparse_columns{},factor_rows{},factor_columns{};
    double raw_density{},q_transformed_density{},tail_density{};
    double seconds{},symbolic_seconds{},numeric_seconds{},q_transform_seconds{},tail_extract_seconds{},
        tail_qmult_seconds{},tail_compact_seconds{};
};
struct DerivativeWork
{
    struct Preparation
    {
        double raw_assembly_seconds{},free_design_assembly_seconds{},factor_match_seconds{},factor_build_seconds{},
            factor_compact_seconds{},rank_seconds{},least_squares_seconds{},normal_solve_seconds{},
            cancellation_check_seconds{},cancellation_fallback_seconds{};
    } preparation;
    double rows_seconds{},jacobian_qr_seconds{},norms_seconds{},outer_overhead_seconds{};
    Eigen::Index maximum_generated_rows{},maximum_reduction_rows{};
    std::size_t tile_count{};
    TiledQrTelemetry projected_qr,jacobian_qr,compact_jacobian_qr;
    ProjectedReductionWorkForTesting projected_reduction;
};
DerivativeWork & DerivativeWorkForTesting();
#endif
}
