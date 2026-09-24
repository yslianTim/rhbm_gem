#pragma once
#include "SparseFactor.hpp"

namespace rhbm_gem::core::joint_component {
struct OperatorWork
{
    std::size_t preparations{},applications{},adjoints{},rank_checks{},normals{};
    double preparation_seconds{},rank_seconds{},apply_seconds{},adjoint_seconds{},normal_seconds{},design_seconds{},factor_seconds{},compact_seconds{},svd_seconds{};
};
OperatorWork & OperatorWorkForTesting();
// A unique immutable identity; equal dimensions do not imply equal states.
struct LinearizationIdentity {};
class ProfileJacobianOperator
{
    Sparse raw_;
    Vector contraction_;
    Indices owners_;
    std::shared_ptr<FreeDesignFactor> factor_;
    std::shared_ptr<const LinearizationIdentity> identity_;
    double scale_{};
    bool valid_{};
    std::string reason_;
    void Check(VectorRef,Eigen::Index) const;
public:
    ProfileJacobianOperator(const Evaluation &,const EvaluationContext &,double absolute=-1);
    bool Valid() const {return valid_;}
    const std::string & Reason() const {return reason_;}
    Eigen::Index Rows() const {return raw_.rows();}
    Eigen::Index Columns() const {return raw_.cols();}
    Eigen::Index FreeColumns() const {return contraction_.size();}
    Eigen::Index RawNonZeros() const {return raw_.nonZeros();}
    const std::shared_ptr<const LinearizationIdentity> & Identity() const {return identity_;}
    Vector Apply(VectorRef) const;
    Vector ApplyAdjoint(VectorRef) const;
    Vector ApplyNormal(VectorRef) const;
};
}
