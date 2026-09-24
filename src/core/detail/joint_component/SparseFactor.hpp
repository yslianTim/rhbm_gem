#pragma once
#include "Numerics.hpp"
#include "ResourceWork.hpp"
#include <chrono>
#include <span>
#include <cstdint>

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
    double matrix_preparation_seconds{},symbolic_seconds{},numeric_seconds{},reference_seconds{},reference_svd_seconds{},derivative_seconds{};
    std::size_t derivative_preparations{},derivative_compacts{},reference_compacts{},free_design_svds{},reference_svds{},reference_solves{},bdc_svds{},jacobi_retries{};
    double derivative_compact_seconds{},reference_compact_seconds{},free_design_svd_seconds{},reference_solve_seconds{},cancellation_seconds{},jacobi_retry_seconds{};
};
SparseWork & SparseWorkForTesting();
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
public:
    LinearWorkspace();
    void Bind(const void * domain,const void * observations,const LinearPolicy *);
    std::shared_ptr<FreeDesignFactor> Factor(const Sparse &,const std::vector<Eigen::Index> &,double);
};
std::pair<Matrix,Vector> SparseReferenceQR(const Sparse &,const Vector &,const Vector &,VectorRef);
}
