#pragma once
#include "Numerics.hpp"
#include "ResourceWork.hpp"
#include <chrono>
#include <span>
#include <cstdint>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif

namespace rhbm_gem::core::joint_component {
struct SparseFactorState;
// Read-only exported SPQR storage. The owning immutable factor must outlive this view.
struct RankFactorView
{
    const Sparse * design{};
    Eigen::Index rows{},columns{},reflectors{};
    std::span<const int64_t> r_outer,r_inner,h_outer,h_inner,permutation,row_permutation;
    std::span<const double> r_values,h_values,tau;
};
struct SparseWork
{
    std::size_t q_actions{},triangular_solves{},compact_extractions{},fixed_factorizations{},factor_storage_bytes{};
    double least_squares_seconds{},q_seconds{},triangular_seconds{},compact_seconds{},fixed_factor_seconds{};
    std::size_t symbolic{},numeric{},symbolic_reuses{},factor_reuses{},reference{},factor_nonzeros{},cancellation_reductions{};
    Eigen::Index symbolic_rows{},symbolic_columns{},numeric_rows{},numeric_columns{},fixed_factor_rows{},fixed_factor_columns{};
    std::size_t symbolic_input_nonzeros{},numeric_input_nonzeros{},fixed_factor_input_nonzeros{},fixed_factor_nonzeros{};
    std::size_t fixed_factor_storage_bytes{};
    double matrix_preparation_seconds{},symbolic_seconds{},numeric_seconds{},reference_seconds{},reference_svd_seconds{},derivative_seconds{};
    std::size_t derivative_preparations{},derivative_compacts{},reference_compacts{},free_design_svds{},reference_svds{},reference_solves{},bdc_svds{},jacobi_retries{};
    double derivative_compact_seconds{},reference_compact_seconds{},free_design_svd_seconds{},reference_solve_seconds{},cancellation_seconds{},jacobi_retry_seconds{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::size_t native_operator_factorizations{};
    double native_operator_factor_seconds{};
#endif
};
SparseWork & SparseWorkForTesting();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct StructuredProjectedQrResultForTesting
{
    Matrix factor;
    Vector response;
    Eigen::Index sparse_rows{},sparse_columns{},factor_rows{},factor_columns{};
    std::size_t sparse_nonzeros{};
    std::size_t maximum_dense_bytes{};
    double symbolic_seconds{},numeric_seconds{},orthogonal_seconds{};
    bool valid{},input_order_preserved{};
    std::string reason;
};
StructuredProjectedQrResultForTesting StructuredProjectedQrForTesting(
    const Sparse &,const Sparse &,VectorRef,double);
struct ProjectedTailTransformForTesting
{
    Sparse tail;
    Eigen::Index rows{},columns{},tail_rows{};
    std::size_t input_nonzeros{},transformed_nonzeros{},tail_nonzeros{},
        transformed_storage_bytes{},tail_storage_bytes{};
    double q_transform_seconds{},tail_extract_seconds{};
};
struct ProjectedTailQrResultForTesting
{
    Matrix factor;
    Vector response;
    Eigen::Index rows{},free_design_columns{},width_columns{},tail_rows{};
    std::size_t q_transformed_nonzeros{},tail_nonzeros{},q_transformed_storage_bytes{},tail_storage_bytes{},
        tail_factor_nonzeros{},tail_factor_storage_bytes{};
    double q_transform_seconds{},tail_extract_seconds{},tail_symbolic_seconds{},
        tail_numeric_seconds{},tail_qmult_seconds{},tail_compact_seconds{},seconds{};
    bool valid{},columns_restored{};
    std::string ordering,reason;
};
struct FactorResidencyRecord
{
    std::size_t factor_id{},generation{};
    std::string kind,role,created_at_stage,destroyed_at_stage,search_stage,destroyed_search_stage;
    Eigen::Index rows{},columns{};
    std::size_t nonzeros{},r_nonzeros{},h_nonzeros{},owned_factor_bytes{};
    double created_seconds{},destroyed_seconds{};
    bool alive{};
};
struct FactorConstructionRecord
{
    std::size_t factor_id{},generation{};
    std::string kind,role,stage;
    Eigen::Index rows{},columns{};
    std::size_t nonzeros{},owned_factor_bytes_estimate{};
    std::function<std::size_t()> estimate_owned_bytes;
};
struct FactorResidencyWork
{
    std::vector<FactorResidencyRecord> factors;
    std::vector<FactorConstructionRecord> constructing;
    std::size_t next_factor_id{1},maximum_concurrent_factor_count{},maximum_concurrent_owned_bytes{};
    std::string maximum_concurrent_factor_stage,maximum_concurrent_owned_bytes_stage;
    std::size_t peak_rss_bytes_observed{},peak_rss_factor_count{},peak_rss_owned_bytes{};
    std::vector<std::size_t> peak_rss_factor_ids;
    std::vector<std::string> peak_rss_factor_generations;
    std::string peak_rss_stage;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
};
FactorResidencyWork & FactorResidencyWorkForTesting();
void ResetFactorResidencyWorkForTesting();
std::string & FactorCreationRoleForTesting();
class FactorCreationRoleScopeForTesting
{
    std::string previous_;
public:
    explicit FactorCreationRoleScopeForTesting(std::string);
    ~FactorCreationRoleScopeForTesting();
    FactorCreationRoleScopeForTesting(const FactorCreationRoleScopeForTesting &)=delete;
    FactorCreationRoleScopeForTesting & operator=(const FactorCreationRoleScopeForTesting &)=delete;
};
enum class SpqrOrdering {Colamd,Default,Best,Metis};
SpqrOrdering & SpqrOrderingForTesting();
const char * SpqrOrderingName(SpqrOrdering);
bool SpqrOrderingAvailable(SpqrOrdering);
enum class OperatorFactorRepresentation {ExportedFixed,NativeQr};
OperatorFactorRepresentation & OperatorFactorRepresentationForTesting();
const char * OperatorFactorRepresentationName(OperatorFactorRepresentation);
class OperatorFactorRepresentationScopeForTesting
{
    OperatorFactorRepresentation previous_;
public:
    explicit OperatorFactorRepresentationScopeForTesting(OperatorFactorRepresentation);
    ~OperatorFactorRepresentationScopeForTesting();
    OperatorFactorRepresentationScopeForTesting(const OperatorFactorRepresentationScopeForTesting &)=delete;
    OperatorFactorRepresentationScopeForTesting & operator=(const OperatorFactorRepresentationScopeForTesting &)=delete;
};
#endif
// Inclusive elapsed time, including early returns and exception unwinding.
struct WorkTimer
{
    double & seconds;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
    explicit WorkTimer(double & value):seconds(value) {}
    ~WorkTimer() {seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
};
bool SparseBackendEnabled();
class FreeDesignFactor
{
    friend class LinearWorkspace;
    std::shared_ptr<SparseFactorState> state_;
    [[maybe_unused]] std::size_t generation_{};
    FreeDesignFactor(std::shared_ptr<SparseFactorState>,std::size_t);
    void Check() const;
public:
    static std::shared_ptr<FreeDesignFactor> Fixed(const Sparse &,const std::vector<Eigen::Index> &);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    static std::shared_ptr<FreeDesignFactor> NativeFixedForTesting(const Sparse &,const std::vector<Eigen::Index> &);
    const Sparse & DesignForTesting() const;
    const std::vector<Eigen::Index> & ColumnsForTesting() const;
    double ToleranceForTesting() const;
    // Census sparse Q^T rhs in place; materialize_tail is reserved for QR parity tests and candidates.
    ProjectedTailTransformForTesting OrthogonalTransposeTailSparseForTesting(const Sparse &,bool=false) const;
    Matrix OrthogonalTransposeForTesting(const Matrix &) const;
    ProjectedTailQrResultForTesting ProjectedTailQrForTesting(const Sparse &,VectorRef,double) const;
#endif
    bool Matches(const Sparse &,const std::vector<Eigen::Index> &) const;
    int Rank() const;
    std::optional<RankFactorView> RankView() const;
    Matrix Compact() const;
    Matrix LeastSquares(const Matrix &) const;
    Matrix PseudoInverseTranspose(const Matrix &) const;
    Matrix ProjectComplement(const Matrix &) const;
    Matrix NormalSolve(const Matrix &) const;
};
class LinearWorkspace
{
    std::shared_ptr<SparseFactorState> state_;
    bool copy_on_write_{};
public:
    LinearWorkspace();
    void Bind(const void * domain,const void * observations,const LinearPolicy *);
    void EnableCopyOnWrite() {copy_on_write_=true;}
    std::shared_ptr<FreeDesignFactor> Factor(const Sparse &,const std::vector<Eigen::Index> &,double);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    void HandoffForTesting();
#endif
};
std::pair<Matrix,Vector> SparseReferenceQR(const Sparse &,const Vector &,const Vector &,VectorRef);
}
