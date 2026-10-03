#include "ProfileJacobianOperator.hpp"
#include "CompactSvd.hpp"
#include <cmath>

namespace rhbm_gem::core::joint_component {
OperatorWork & OperatorWorkForTesting() {static thread_local OperatorWork work; return work;}
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
}
ProfileJacobianOperator::ProfileJacobianOperator(const Evaluation & e,const EvaluationContext & context,double absolute,FreeDesignRankBackend backend)
    :identity_(std::make_shared<const LinearizationIdentity>()),scale_(context.scale)
{
    const auto n=e.x.rows(),m=e.eta.size();
    ResourcePhase phase("operator-prepare",true,n,e.x.cols(),static_cast<std::size_t>(e.x.nonZeros()));
    auto & work=OperatorWorkForTesting(); ++work.preparations; WorkTimer timer(work.preparation_seconds);
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
    try {
        // Dedicated ownership: a trial's mutable workspace cannot expire this factor.
        {ResourcePhase factor_stage("fixed-operator-factor",true,n,p,static_cast<std::size_t>(design.nonZeros()));
        WorkTimer factor_timer(work.factor_seconds);
        factor_=FreeDesignFactor::Fixed(design,free);}
        {
            ResourcePhase rank_phase("operator-rank");
            ResourcePhase rank_stage("rank-certificate",true,n,p,static_cast<std::size_t>(design.nonZeros()));
            ++work.rank_checks; WorkTimer rank_timer(work.rank_seconds);
            RankWorkAudit rank_audit{work,SparseWorkForTesting()};
            if(backend==FreeDesignRankBackend::SpqrBounds)
            {
                rank_evidence_=EvaluateFreeDesignRank(design,factor_.get(),{context.rank,p,absolute},context.search.operator_rank.budget);
                work.rank_entries+=rank_evidence_.entries;
                work.rank_workspace_bytes=std::max(work.rank_workspace_bytes,rank_evidence_.workspace_bytes);
                work.rank_work_stage=rank_evidence_.work_stage;
                work.rank_estimated_total_entries=rank_evidence_.estimated_total_entries;
                work.rank_estimated_remaining_entries=rank_evidence_.estimated_remaining_entries;
                work.rank_estimated_reconstruction_entries=rank_evidence_.estimated_reconstruction_entries;
                work.rank_certificate=rank_evidence_.certificate;
                work.rank_local_witness=rank_evidence_.local_witness;
                work.rank_design_nonzeros=rank_evidence_.design_nonzeros;
                work.rank_r_nonzeros=rank_evidence_.factor_r_nonzeros;
                work.rank_reflector_nonzeros=rank_evidence_.reflector_nonzeros;
                work.rank_reflectors=rank_evidence_.reflector_count;
                work.rank_status=RankStatus(rank_evidence_.status); work.rank_reason=rank_evidence_.reason;
                if(rank_evidence_.status!=FreeDesignRankStatus::FullRank)
                {
                    reason_=rank_evidence_.status==FreeDesignRankStatus::Deficient ? "rank-deficient-free-design" : rank_evidence_.reason;
                    factor_.reset(); return;
                }
            }
            else
            {
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
