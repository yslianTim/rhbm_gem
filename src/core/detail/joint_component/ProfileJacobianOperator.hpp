#pragma once
#include "SparseFactor.hpp"
#include "FreeDesignRank.hpp"

namespace rhbm_gem::core::joint_component {
struct OperatorWork
{
    std::size_t preparations{},applications{},adjoints{},rank_checks{},rank_rows{},rank_columns{},rank_entries{},rank_workspace_bytes{},
        rank_compact_extractions{},rank_free_design_svds{},rank_design_nonzeros{},rank_r_nonzeros{},rank_reflector_nonzeros{},rank_reflectors{},normals{};
    FreeDesignRankWorkStage rank_work_stage{FreeDesignRankWorkStage::None};
    std::optional<std::size_t> rank_estimated_total_entries,rank_estimated_remaining_entries,rank_estimated_reconstruction_entries;
    FreeDesignRankCertificate rank_certificate{FreeDesignRankCertificate::None};
    FreeDesignLocalWitness rank_local_witness;
    double preparation_seconds{},rank_seconds{},apply_seconds{},adjoint_seconds{},normal_seconds{},design_seconds{},factor_seconds{},compact_seconds{},svd_seconds{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::size_t native_factor_accepted{},native_factor_fallbacks{},accepted_factor_reuse_attempts{},
        accepted_factor_reuse_accepted{},accepted_factor_reuse_fallbacks{};
    std::string factor_ownership,accepted_factor_reuse_fallback_reason;
#endif
    std::string rank_status,rank_reason;
};
OperatorWork & OperatorWorkForTesting();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
enum class OperatorFactorOwnershipKindForTesting {DedicatedFixed,DedicatedNative,ReuseAcceptedCopyOnWrite,ReuseAcceptedHandoff};
OperatorFactorOwnershipKindForTesting & OperatorFactorOwnershipForTesting();
const char * OperatorFactorOwnershipName(OperatorFactorOwnershipKindForTesting);
class OperatorFactorOwnershipScopeForTesting
{
    OperatorFactorOwnershipKindForTesting previous_;
public:
    explicit OperatorFactorOwnershipScopeForTesting(OperatorFactorOwnershipKindForTesting);
    ~OperatorFactorOwnershipScopeForTesting();
    OperatorFactorOwnershipScopeForTesting(const OperatorFactorOwnershipScopeForTesting &)=delete;
    OperatorFactorOwnershipScopeForTesting & operator=(const OperatorFactorOwnershipScopeForTesting &)=delete;
};
#endif
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
    FreeDesignRankResult rank_evidence_;
    std::string reason_;
    void Check(VectorRef,Eigen::Index) const;
public:
    ProfileJacobianOperator(const Evaluation &,const EvaluationContext &,double absolute=-1,FreeDesignRankBackend=FreeDesignRankBackend::Dense);
    const FreeDesignRankResult & RankEvidence() const {return rank_evidence_;}
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
