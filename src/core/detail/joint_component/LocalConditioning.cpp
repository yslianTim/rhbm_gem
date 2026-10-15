#include "LocalConditioning.hpp"
#include "CompactSvd.hpp"
#include "ResourceWork.hpp"
#include "TiledDerivative.hpp"
#include "TiledQR.hpp"
#include <chrono>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}

struct LocalSpectrum
{
    Eigen::Index rank{};
    double minimum{unavailable},maximum{unavailable},condition{unavailable};
    bool available{},threshold_sensitive{};
    std::string status{"unavailable"},reason;
};

LocalSpectrum SpectrumOf(const Matrix & factor,Eigen::Index rows,Eigen::Index columns,const RankPolicy & policy)
{
    LocalSpectrum out;
    if(rows<=0 || columns<=0 || factor.cols()!=columns || factor.rows()==0)
    {out.reason="empty-design"; return out;}
    const auto svd=EvaluateRank(factor,{policy,columns,-1,RankBoundary::StrictGreater});
    if(!svd.valid || svd.singular_values.size()==0 || !svd.singular_values.allFinite())
    {out.reason="spectrum-factorization-failed"; return out;}
    out.available=true; out.rank=svd.rank; out.maximum=svd.singular_values(0);
    out.minimum=rows<columns ? 0.0 : svd.singular_values(svd.singular_values.size()-1);
    out.threshold_sensitive=!CompactRankDecisionSafe(svd) && out.maximum>0;
    if(out.rank==columns && out.minimum>0)
    {
        out.condition=out.maximum/out.minimum;
        out.status=out.threshold_sensitive ? "near-singular" : "full-rank";
    }
    else
    {
        out.status="rank-deficient";
        if(out.threshold_sensitive) out.reason="rank-threshold-sensitive";
    }
    if(!std::isfinite(out.condition)) out.condition=unavailable;
    return out;
}

void MarkUnavailable(JointBlockConditioning & evidence,const std::string & reason)
{
    evidence.status="unavailable"; evidence.reason=reason;
    evidence.free_ac_status=evidence.projected_width_status=evidence.profile_status="unavailable";
    evidence.free_ac_reason=evidence.projected_width_reason=evidence.profile_reason=reason;
}

void SetSpectrum(std::string & status,std::string & reason,std::optional<std::size_t> & rank,
    std::optional<double> & condition,const LocalSpectrum & spectrum)
{
    status=spectrum.status; reason=spectrum.reason;
    if(!spectrum.available) return;
    rank=static_cast<std::size_t>(spectrum.rank);
    if(std::isfinite(spectrum.condition)) condition=spectrum.condition;
}
}

LocalConditioningAssessment AssessFixedNeighborLocalConditioning(
    const Domain & domain,VectorRef y,const Vector & eta,const Vector & beta,const EvaluationContext & context,
    std::string block_id,std::size_t resident_bytes)
{
    const auto started=Clock::now(); LocalConditioningAssessment out;
    auto & evidence=out.evidence;
    evidence.block_id=std::move(block_id); evidence.core_atom_count=domain.atoms.size();
    evidence.affected_row_count=static_cast<std::size_t>(std::max<Eigen::Index>(domain.rows,0));
    if(domain.atoms.size()>std::numeric_limits<std::size_t>::max()/2)
    {MarkUnavailable(evidence,"size-overflow"); out.seconds=Seconds(started); return out;}
    MemoryShape shape; shape.component_atoms=domain.atoms.size();
    shape.component_observations=evidence.affected_row_count; shape.component_blocks=1;
    shape.component_free_columns=2*domain.atoms.size(); shape.component_width_columns=domain.atoms.size();
    shape.maximum_block_rows=evidence.affected_row_count; shape.maximum_neighbor_atoms=domain.atoms.size();
    shape.maximum_tile_rows=static_cast<std::size_t>(derivative_tile_rows); shape.resident_bytes=resident_bytes;
    for(const auto & support:domain.atoms)
        if(!CheckedAddBytes(shape.component_support_entries,support.size(),shape.component_support_entries))
        {shape.component_support_entries=std::numeric_limits<std::size_t>::max(); break;}
    shape.maximum_block_support_entries=shape.component_support_entries;
    out.memory_estimate=PlanMemory(shape,MemoryStage::LocalConditioning);
    if(out.memory_estimate.status!=MemoryStatus::Overflow)
        evidence.workspace_peak_memory_estimate_bytes=out.memory_estimate.estimated_peak_bytes;
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {
        if(out.memory_estimate.status==MemoryStatus::Overflow) out.memory_estimate.reason="size-overflow";
        MarkUnavailable(evidence,std::string("resource-limit:local-conditioning:")+out.memory_estimate.reason);
        out.seconds=Seconds(started); return out;
    }
    if(eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) ||
        beta.size()!=2*eta.size() || y.size()!=domain.rows)
    {MarkUnavailable(evidence,"invalid-local-dimensions"); out.seconds=Seconds(started); return out;}
    if(domain.rows<=0)
    {MarkUnavailable(evidence,"zero-support"); out.seconds=Seconds(started); return out;}
    if(!eta.allFinite() || !beta.allFinite() || !y.allFinite())
    {MarkUnavailable(evidence,"nonfinite-local-state"); out.seconds=Seconds(started); return out;}
    if(!(context.scale>0) || !std::isfinite(context.scale))
    {MarkUnavailable(evidence,"invalid-scale"); out.seconds=Seconds(started); return out;}

    std::vector<Eigen::Index> free_columns;
    free_columns.reserve(static_cast<std::size_t>(beta.size()));
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {
        if(beta(2*atom)==0) ++evidence.active_amplitude_count;
        else free_columns.push_back(2*atom);
        free_columns.push_back(2*atom+1);
    }
    evidence.free_ac_columns=free_columns.size();
    const auto state=EvaluateState(domain,y,eta,beta,context,ProfileEvaluationRole::AcceptedEndpoint);
    if(!state.valid)
    {MarkUnavailable(evidence,state.reason.empty() ? "local-state-evaluation-failed" : state.reason);
        out.seconds=Seconds(started); return out;}

    const auto p=static_cast<Eigen::Index>(free_columns.size());
    Vector norms(p);
    bool zero_column=false;
    for(Eigen::Index k=0;k<p;++k)
    {
        norms(k)=state.x.col(free_columns[static_cast<std::size_t>(k)]).norm();
        zero_column|=!(norms(k)>0);
    }
    TiledQR ac_factor(p,0,"local-conditioning-ac");
    for(Eigen::Index first=0;first<domain.rows;first+=derivative_tile_rows)
    {
        const auto count=std::min(derivative_tile_rows,domain.rows-first);
        Matrix tile=Matrix::Zero(count,p);
        for(Eigen::Index column=0;column<p;++column)
        {
            const auto source=free_columns[static_cast<std::size_t>(column)];
            if(!(norms(column)>0)) continue;
            for(Sparse::InnerIterator entry(state.x,source);entry;++entry)
                if(entry.row()>=first && entry.row()<first+count)
                    tile(entry.row()-first,column)=entry.value()/norms(column);
        }
        ac_factor.Append(tile,Matrix(count,0));
    }
    const RankPolicy ac_policy{domain.rows,p,0};
    const auto ac_spectrum=SpectrumOf(ac_factor.r,domain.rows,p,ac_policy);
    evidence.free_ac_status=zero_column ? "zero-column" : ac_spectrum.status;
    evidence.free_ac_reason=zero_column ? "zero-free-ac-column" : ac_spectrum.reason;
    if(ac_spectrum.available)
    {
        evidence.free_ac_rank=static_cast<std::size_t>(ac_spectrum.rank);
        if(std::isfinite(ac_spectrum.minimum)) evidence.free_ac_sigma_min=ac_spectrum.minimum;
        if(std::isfinite(ac_spectrum.maximum)) evidence.free_ac_sigma_max=ac_spectrum.maximum;
        if(std::isfinite(ac_spectrum.condition)) evidence.free_ac_condition=ac_spectrum.condition;
    }

    const auto projected_prepared=PrepareProjectedWidthDerivative(state,context.scale,&context,-1,
        derivative_tile_rows,true);
    if(projected_prepared.valid)
    {
        const auto projected=ReduceProjectedWidthDerivative(projected_prepared);
        if(projected.valid)
        {
            const auto spectrum=SpectrumOf(projected.projected,domain.rows,eta.size(),context.rank);
            SetSpectrum(evidence.projected_width_status,evidence.projected_width_reason,
                evidence.projected_width_rank,evidence.projected_width_condition,spectrum);
            if(evidence.active_amplitude_count && spectrum.available)
            {
                evidence.projected_width_status="active-face-"+evidence.projected_width_status;
                evidence.projected_width_reason="conditional-on-active-amplitude-face";
            }
        }
        else evidence.projected_width_reason=projected.reason;
    }
    else evidence.projected_width_reason=projected_prepared.reason;
    if(evidence.projected_width_status=="not-assessed") evidence.projected_width_status="unavailable";

    if(evidence.active_amplitude_count)
    {
        evidence.profile_status="active-set-boundary";
        evidence.profile_reason="profile-derivative-may-switch-active-face";
    }
    else if(!ac_spectrum.available || ac_spectrum.rank!=p)
    {
        evidence.profile_status="unavailable";
        evidence.profile_reason="free-ac-rank-deficient";
    }
    else
    {
        const auto prepared=PrepareDerivative(state,context.scale,&context);
        if(!prepared.valid)
        {
            evidence.profile_status="unavailable";
            evidence.profile_reason=prepared.reason;
        }
        else
        {
            const auto reduced=ReduceDerivative(prepared,state.residual,false);
            if(!reduced.valid)
            {
                evidence.profile_status="unavailable";
                evidence.profile_reason=reduced.reason;
            }
            else
            {
                const auto spectrum=SpectrumOf(reduced.jacobian,domain.rows,eta.size(),context.rank);
                SetSpectrum(evidence.profile_status,evidence.profile_reason,
                    evidence.profile_rank,evidence.profile_condition,spectrum);
                if(ac_spectrum.threshold_sensitive && spectrum.available)
                {
                    evidence.profile_status="conditional-free-ac-near-singular";
                    evidence.profile_reason="free-ac-rank-threshold-sensitive";
                }
            }
        }
    }
    evidence.status="assessed";
    if(evidence.free_ac_status=="unavailable" || evidence.projected_width_status=="unavailable" ||
        evidence.profile_status=="unavailable")
    {evidence.status="partial"; evidence.reason="one-or-more-local-diagnostics-unavailable";}
    out.seconds=Seconds(started); return out;
}
}
