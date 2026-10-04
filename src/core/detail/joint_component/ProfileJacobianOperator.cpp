#include "ProfileJacobianOperator.hpp"
#include "CompactSvd.hpp"
#include "OperatorSearch.hpp"
#include <cmath>

namespace rhbm_gem::core::joint_component {
OperatorWork & OperatorWorkForTesting() {static thread_local OperatorWork work; return work;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
OperatorFactorOwnershipKindForTesting & OperatorFactorOwnershipForTesting()
{static thread_local auto kind=OperatorFactorOwnershipKindForTesting::ReuseAcceptedCopyOnWrite; return kind;}
const char * OperatorFactorOwnershipName(OperatorFactorOwnershipKindForTesting kind)
{
    switch(kind)
    {
    case OperatorFactorOwnershipKindForTesting::DedicatedFixed: return "dedicated-fixed";
    case OperatorFactorOwnershipKindForTesting::DedicatedNative: return "dedicated-native";
    case OperatorFactorOwnershipKindForTesting::ReuseAcceptedCopyOnWrite: return "reuse-accepted-copy-on-write";
    case OperatorFactorOwnershipKindForTesting::ReuseAcceptedHandoff: return "reuse-accepted-handoff";
    case OperatorFactorOwnershipKindForTesting::ReuseAcceptedEvictBeforeTrial: return "evict-before-trial";
    }
    return "unknown";
}
OperatorFactorOwnershipScopeForTesting::OperatorFactorOwnershipScopeForTesting(OperatorFactorOwnershipKindForTesting kind)
    :previous_(OperatorFactorOwnershipForTesting())
{OperatorFactorOwnershipForTesting()=kind;}
OperatorFactorOwnershipScopeForTesting::~OperatorFactorOwnershipScopeForTesting()
{OperatorFactorOwnershipForTesting()=previous_;}
#endif
namespace {
struct RankWorkAudit
{
    OperatorWork & operator_work;
    SparseWork & sparse_work;
    std::size_t compact_extractions{sparse_work.compact_extractions};
    std::size_t free_design_svds{sparse_work.free_design_svds};
    ~RankWorkAudit()
    {
        operator_work.rank_compact_extractions+=sparse_work.compact_extractions-compact_extractions;
        operator_work.rank_free_design_svds+=sparse_work.free_design_svds-free_design_svds;
    }
};
const char * RankStatus(FreeDesignRankStatus status)
{
    switch(status)
    {
    case FreeDesignRankStatus::Unavailable: return "unavailable";
    case FreeDesignRankStatus::FullRank: return "full-rank";
    case FreeDesignRankStatus::Deficient: return "deficient";
    }
    return "unavailable";
}
void RecordRankEvidence(OperatorWork & work,const FreeDesignRankResult & evidence)
{
    work.rank_entries+=evidence.entries;
    work.rank_workspace_bytes=std::max(work.rank_workspace_bytes,evidence.workspace_bytes);
    work.rank_work_stage=evidence.work_stage;
    work.rank_estimated_total_entries=evidence.estimated_total_entries;
    work.rank_estimated_remaining_entries=evidence.estimated_remaining_entries;
    work.rank_estimated_reconstruction_entries=evidence.estimated_reconstruction_entries;
    work.rank_certificate=evidence.certificate;
    work.rank_local_witness=evidence.local_witness;
    work.rank_design_nonzeros=evidence.design_nonzeros;
    work.rank_r_nonzeros=evidence.factor_r_nonzeros;
    work.rank_reflector_nonzeros=evidence.reflector_nonzeros;
    work.rank_reflectors=evidence.reflector_count;
    work.rank_status=RankStatus(evidence.status);
    work.rank_reason=evidence.reason;
}
}
ProfileJacobianOperator::ProfileJacobianOperator(const Evaluation & e,const EvaluationContext & context,double absolute,FreeDesignRankBackend backend)
    :identity_(std::make_shared<const LinearizationIdentity>()),scale_(context.scale)
{
    const auto n=e.x.rows(),m=e.eta.size();
    bool reused_rank_evidence_ready{};
    ResourcePhase phase("operator-prepare",true,n,e.x.cols(),static_cast<std::size_t>(e.x.nonZeros()));
    auto & work=OperatorWorkForTesting(); ++work.preparations; WorkTimer timer(work.preparation_seconds);
    bool reuse_accepted_factor{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const auto ownership=OperatorFactorOwnershipForTesting();
    work.factor_ownership=OperatorFactorOwnershipName(ownership);
    reuse_accepted_factor=ownership==OperatorFactorOwnershipKindForTesting::ReuseAcceptedCopyOnWrite ||
        ownership==OperatorFactorOwnershipKindForTesting::ReuseAcceptedHandoff ||
        ownership==OperatorFactorOwnershipKindForTesting::ReuseAcceptedEvictBeforeTrial;
#else
    reuse_accepted_factor=true;
#endif
    reason_="invalid-inner";
    if(!e.valid || !(scale_>0) || !std::isfinite(scale_) || m<=0 || e.beta.size()!=2*m ||
        e.x.cols()!=2*m || e.derivative.rows()!=n || e.derivative.cols()!=2*m || e.residual.size()!=n ||
        !e.beta.allFinite() || !e.residual.allFinite()) return;
    const auto design_started=std::chrono::steady_clock::now();
    Indices free;
    for(Eigen::Index k=0;k<2*m;++k) if(k%2 || e.beta(k)>0) free.push_back(k);
    const auto p=static_cast<Eigen::Index>(free.size());
    work.rank_rows=static_cast<std::size_t>(n); work.rank_columns=static_cast<std::size_t>(p);
    if(n<p) {reason_="rank-deficient-free-design"; return;}
    Sparse design(n,p); raw_.resize(n,m); contraction_.resize(p); owners_.resize(free.size());
    std::vector<Eigen::Triplet<double>> entries,raw;
    for(Eigen::Index k=0;k<2*m;++k) for(Sparse::InnerIterator v(e.derivative,k);v;++v)
    {
        if(!std::isfinite(v.value())) {reason_="nonfinite-derivative"; return;}
        raw.emplace_back(v.row(),k/2,v.value()*e.beta(k));
    }
    for(Eigen::Index c=0;c<p;++c)
    {
        const auto k=free[static_cast<std::size_t>(c)]; const double norm=e.x.col(k).norm();
        if(!(norm>0) || !std::isfinite(norm)) {reason_="zero-free-column"; return;}
        for(Sparse::InnerIterator v(e.x,k);v;++v) entries.emplace_back(v.row(),c,v.value()/norm);
        owners_[static_cast<std::size_t>(c)]=k/2;
        contraction_(c)=e.derivative.col(k).dot(e.residual)/norm;
    }
    if(!contraction_.allFinite()) {reason_="nonfinite-derivative"; return;}
    design.setFromTriplets(entries.begin(),entries.end()); raw_.setFromTriplets(raw.begin(),raw.end());
    work.design_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-design_started).count();
    if(reuse_accepted_factor)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++work.accepted_factor_reuse_attempts;
#endif
        const auto reuse_fallback=[&](const char * reason) {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++work.accepted_factor_reuse_fallbacks; work.accepted_factor_reuse_fallback_reason=reason;
#else
            (void)reason;
#endif
        };
        if(!e.factor || !e.factor->Matches(design,free))
            reuse_fallback("accepted-factor-mismatch");
        else if(backend!=FreeDesignRankBackend::SpqrBounds)
        {
            factor_=e.factor;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++work.accepted_factor_reuse_accepted;
#endif
        }
        else
        {
            ResourcePhase rank_phase("operator-rank");
            ResourcePhase rank_stage("rank-certificate",true,n,p,static_cast<std::size_t>(design.nonZeros()));
            ++work.rank_checks; WorkTimer rank_timer(work.rank_seconds);
            RankWorkAudit rank_audit{work,SparseWorkForTesting()};
            auto evidence=EvaluateFreeDesignRank(design,nullptr,{context.rank,p,absolute},context.search.operator_rank.budget);
            if(evidence.status==FreeDesignRankStatus::FullRank)
            {
                RecordRankEvidence(work,evidence); rank_evidence_=std::move(evidence);
                reused_rank_evidence_ready=true; factor_=e.factor;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                ++work.accepted_factor_reuse_accepted;
#endif
            }
            else reuse_fallback("rank-certificate-needs-fixed-factor-view");
        }
    }
    try {
#if defined(RHBM_GEM_TEST_INSTRUMENTATION) && defined(RHBM_GEM_JOINT_SPQR)
        const bool native_requested=ownership==OperatorFactorOwnershipKindForTesting::DedicatedNative ||
            (ownership==OperatorFactorOwnershipKindForTesting::DedicatedFixed &&
                OperatorFactorRepresentationForTesting()==OperatorFactorRepresentation::NativeQr);
        bool native_eligible=false;
        FreeDesignRankResult native_evidence;
        if(native_requested && backend==FreeDesignRankBackend::SpqrBounds)
        {
            ResourcePhase rank_phase("operator-rank");
            ResourcePhase rank_stage("rank-certificate",true,n,p,static_cast<std::size_t>(design.nonZeros()));
            ++work.rank_checks; WorkTimer rank_timer(work.rank_seconds);
            RankWorkAudit rank_audit{work,SparseWorkForTesting()};
            native_evidence=EvaluateFreeDesignRank(design,nullptr,{context.rank,p,absolute},context.search.operator_rank.budget);
            RecordRankEvidence(work,native_evidence);
            native_eligible=native_evidence.status==FreeDesignRankStatus::FullRank &&
                native_evidence.certificate==FreeDesignRankCertificate::LocalSupport;
        }
#endif
        // Production reuses a matched accepted factor; dedicated factors remain the bounded-rank fallback.
        {ResourcePhase factor_stage("fixed-operator-factor",true,n,p,static_cast<std::size_t>(design.nonZeros()));
        WorkTimer factor_timer(work.factor_seconds);
#if defined(RHBM_GEM_TEST_INSTRUMENTATION) && defined(RHBM_GEM_JOINT_SPQR)
        if(!factor_ && native_eligible)
        {
            try {factor_=FreeDesignFactor::NativeFixedForTesting(design,free);}
            catch(const std::runtime_error &) {native_eligible=false;}
        }
        if(!factor_ && !native_eligible) factor_=FreeDesignFactor::Fixed(design,free);
        if(native_requested)
        {
            if(native_eligible) ++work.native_factor_accepted;
            else ++work.native_factor_fallbacks;
        }
#else
        if(!factor_) factor_=FreeDesignFactor::Fixed(design,free);
#endif
        }
        {
            ResourcePhase rank_phase("operator-rank");
            ResourcePhase rank_stage("rank-certificate",true,n,p,static_cast<std::size_t>(design.nonZeros()));
            if(backend==FreeDesignRankBackend::SpqrBounds)
            {
#if defined(RHBM_GEM_TEST_INSTRUMENTATION) && defined(RHBM_GEM_JOINT_SPQR)
                if(native_eligible)
                {rank_evidence_=std::move(native_evidence); reused_rank_evidence_ready=true;}
#endif
                if(!reused_rank_evidence_ready)
                {
                    ++work.rank_checks; WorkTimer rank_timer(work.rank_seconds);
                    RankWorkAudit rank_audit{work,SparseWorkForTesting()};
                    rank_evidence_=EvaluateFreeDesignRank(design,factor_.get(),{context.rank,p,absolute},context.search.operator_rank.budget);
                    RecordRankEvidence(work,rank_evidence_);
                }
                if(rank_evidence_.status!=FreeDesignRankStatus::FullRank)
                {
                    reason_=rank_evidence_.status==FreeDesignRankStatus::Deficient ? "rank-deficient-free-design" : rank_evidence_.reason;
                    factor_.reset(); return;
                }
            }
            else
            {
                ++work.rank_checks; WorkTimer rank_timer(work.rank_seconds);
                RankWorkAudit rank_audit{work,SparseWorkForTesting()};
                Matrix compact;
                {WorkTimer compact_timer(work.compact_seconds); compact=factor_->Compact();}
                CompactSvdResult svd;
                {WorkTimer svd_timer(work.svd_seconds); svd=EvaluateRank(compact,{context.rank,p,absolute});}
                rank_evidence_.status=!svd.valid ? FreeDesignRankStatus::Unavailable :
                    svd.rank!=p || factor_->Rank()!=p ? FreeDesignRankStatus::Deficient : FreeDesignRankStatus::FullRank;
                rank_evidence_.certificate=svd.valid ? FreeDesignRankCertificate::DenseOracle : FreeDesignRankCertificate::None;
                rank_evidence_.reason=!svd.valid ? "rank-factorization-failed" :
                    rank_evidence_.status==FreeDesignRankStatus::Deficient ? "rank-deficient-free-design" : "rank-dense-oracle";
                work.rank_certificate=rank_evidence_.certificate;
                work.rank_status=RankStatus(rank_evidence_.status); work.rank_reason=rank_evidence_.reason;
                if(!svd.valid) {reason_="nonfinite-derivative"; factor_.reset(); return;}
                if(svd.rank!=p || factor_->Rank()!=p)
                {reason_="rank-deficient-free-design"; factor_.reset(); return;}
            }
        }
        valid_=true; reason_="full-profile-operator";
    } catch(const std::runtime_error &) {factor_.reset(); reason_="operator-factorization-failed";}
}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
void ProfileJacobianOperator::EvictAcceptedFactorForTesting()
{
    if(!factor_) throw std::logic_error("Cannot evict an unavailable accepted operator factor");
    evicted_design_=factor_->DesignForTesting();
    evicted_columns_=factor_->ColumnsForTesting();
    evicted_tolerance_=factor_->ToleranceForTesting();
    factor_.reset();
}
void ProfileJacobianOperator::RebuildAcceptedFactorForTesting()
{
    if(factor_) return;
    if(evicted_design_.rows()==0 || evicted_columns_.empty())
        throw std::logic_error("Accepted operator factor has no rebuild snapshot");
    ResourcePhase phase("accepted-factor-rebuild",true,evicted_design_.rows(),evicted_design_.cols(),
        static_cast<std::size_t>(evicted_design_.nonZeros()));
    auto & work=SearchWorkForTesting(); ++work.accepted_factor_rebuilds;
    WorkTimer timer(work.accepted_factor_rebuild_seconds);
    FactorCreationRoleScopeForTesting role("accepted-rebuild");
    LinearWorkspace workspace;
    factor_=workspace.Factor(evicted_design_,evicted_columns_,evicted_tolerance_);
}
#endif
void ProfileJacobianOperator::Check(VectorRef rhs,Eigen::Index size) const
{
    if(!valid_) throw std::logic_error("Unavailable profile operator: "+reason_);
    if(rhs.size()!=size || !rhs.allFinite()) throw std::invalid_argument("Invalid profile operator RHS");
}
Vector ProfileJacobianOperator::Apply(VectorRef v) const
{
    ResourcePhase phase("operator-apply");
    Check(v,Columns()); auto & work=OperatorWorkForTesting(); ++work.applications; WorkTimer timer(work.apply_seconds);
    RecordDenseShape("operator-observation-vector",Rows(),1);
    RecordDenseShape("operator-free-vector",FreeColumns(),1);
    Vector t(contraction_.size());
    for(Eigen::Index k=0;k<t.size();++k) t(k)=contraction_(k)*v(owners_[static_cast<std::size_t>(k)]);
    const Vector projected=factor_->ProjectComplement(raw_*v);
    const Vector out=(projected-factor_->PseudoInverseTranspose(t).col(0))/scale_;
    if(!out.allFinite()) throw std::runtime_error("Nonfinite profile operator action");
    return out;
}
Vector ProfileJacobianOperator::ApplyAdjoint(VectorRef w) const
{
    ResourcePhase phase("operator-adjoint");
    Check(w,Rows()); auto & work=OperatorWorkForTesting(); ++work.adjoints; WorkTimer timer(work.adjoint_seconds);
    RecordDenseShape("operator-observation-vector",Rows(),1);
    RecordDenseShape("operator-free-vector",FreeColumns(),1);
    Vector out=raw_.transpose()*factor_->ProjectComplement(w).col(0);
    const Vector z=factor_->LeastSquares(w);
    for(Eigen::Index k=0;k<z.size();++k) out(owners_[static_cast<std::size_t>(k)])-=contraction_(k)*z(k);
    out/=scale_;
    if(!out.allFinite()) throw std::runtime_error("Nonfinite profile operator adjoint");
    return out;
}
Vector ProfileJacobianOperator::ApplyNormal(VectorRef v) const
{
    ResourcePhase phase("operator-normal");
    Check(v,Columns()); auto & work=OperatorWorkForTesting(); ++work.normals; WorkTimer timer(work.normal_seconds);
    RecordDenseShape("operator-observation-vector",Rows(),1);
    RecordDenseShape("operator-free-vector",FreeColumns(),1);
    Vector t(contraction_.size());
    for(Eigen::Index k=0;k<t.size();++k) t(k)=contraction_(k)*v(owners_[static_cast<std::size_t>(k)]);
    Vector out=raw_.transpose()*factor_->ProjectComplement(raw_*v).col(0);
    const Vector solved=factor_->NormalSolve(t);
    for(Eigen::Index k=0;k<solved.size();++k) out(owners_[static_cast<std::size_t>(k)])+=contraction_(k)*solved(k);
    out/=scale_; out/=scale_;
    if(!out.allFinite()) throw std::runtime_error("Nonfinite profile operator normal action");
    return out;
}

}
