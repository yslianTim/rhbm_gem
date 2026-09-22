#pragma once
#include "Numerics.hpp"

namespace rhbm_gem::core::joint_component {
struct SparseFactorState;
struct SparseWork
{
    std::size_t symbolic{},numeric{},symbolic_reuses{},factor_reuses{},reference{},factor_nonzeros{},cancellation_reductions{};
    double matrix_preparation_seconds{},symbolic_seconds{},numeric_seconds{},reference_seconds{},reference_svd_seconds{},derivative_seconds{};
};
SparseWork & SparseWorkForTesting();
bool SparseBackendEnabled();
class FreeDesignFactor
{
    friend class LinearWorkspace;
    std::shared_ptr<SparseFactorState> state_;
    [[maybe_unused]] std::size_t generation_{};
    FreeDesignFactor(std::shared_ptr<SparseFactorState>,std::size_t);
    void Check() const;
public:
    bool Matches(const Sparse &,const std::vector<Eigen::Index> &) const;
    int Rank() const;
    Matrix Compact() const;
    Matrix LeastSquares(const Matrix &) const;
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
