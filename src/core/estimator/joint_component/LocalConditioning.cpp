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
    Eigen::Index rank{},columns{};
    double minimum{unavailable},maximum{unavailable},condition{unavailable},threshold{unavailable},
        threshold_margin_relative{unavailable};
    bool available{},threshold_sensitive{};
    std::string status{"unavailable"},reason;
};

LocalSpectrum SpectrumOf(const Matrix & factor,Eigen::Index rows,Eigen::Index columns,const RankPolicy & policy)
{
    LocalSpectrum out;
    if(rows<=0 || columns<=0 || factor.cols()!=columns || factor.rows()==0)
    {out.reason="empty-design"; return out;}
    Matrix normalized=factor;
    for(Eigen::Index column=0;column<columns;++column)
    {
        const double norm=normalized.col(column).norm();
        if(norm>0) normalized.col(column)/=norm;
    }
    const auto svd=EvaluateRank(normalized,{policy,columns,-1,RankBoundary::StrictGreater});
    if(!svd.valid || svd.singular_values.size()==0 || !svd.singular_values.allFinite())
    {out.reason="spectrum-factorization-failed"; return out;}
    out.available=true; out.rank=svd.rank; out.columns=columns; out.maximum=svd.singular_values(0);
    out.minimum=rows<columns ? 0.0 : svd.singular_values(svd.singular_values.size()-1);
    out.threshold=svd.threshold;
    if(out.maximum>0)
    {
        double distance=std::numeric_limits<double>::infinity();
        for(Eigen::Index k=0;k<svd.singular_values.size();++k)
            distance=std::min(distance,std::abs(svd.singular_values(k)-svd.threshold));
        out.threshold_margin_relative=distance/out.maximum;
        const double roundoff_margin=64*std::numeric_limits<double>::epsilon()*
            static_cast<double>(std::max(normalized.rows(),normalized.cols()))*out.maximum;
        out.threshold_sensitive=!CompactRankDecisionSafe(svd) || distance<=roundoff_margin;
    }
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

std::string ConditioningClass(const LocalSpectrum & spectrum,bool zero_column=false)
{
    if(zero_column) return "Singular";
    if(!spectrum.available) return "Unavailable";
    if(spectrum.rank<spectrum.columns || spectrum.minimum==0) return "Singular";
    if(!std::isfinite(spectrum.condition)) return "Unavailable";
    return spectrum.threshold_sensitive || spectrum.condition*std::numeric_limits<double>::epsilon()>=
        LocalConditioningRoundoffAmplificationWarning ? "IllConditioned" : "WellConditioned";
}

void MarkUnavailable(JointBlockConditioning & evidence,const std::string & reason)
{
    evidence.status="unavailable"; evidence.reason=reason;
    evidence.free_ac_status=evidence.projected_width_status=evidence.profile_status="unavailable";
    evidence.free_ac_reason=evidence.projected_width_reason=evidence.profile_reason=reason;
    evidence.free_ac_conditioning=evidence.projected_width_conditioning=evidence.profile_conditioning="Unavailable";
    evidence.conditioning_classification="Unavailable";
}

void SetSpectrum(std::string & status,std::string & reason,std::optional<std::size_t> & rank,
    std::optional<double> & minimum,std::optional<double> & maximum,std::optional<double> & condition,
    std::optional<double> & threshold,std::optional<double> & threshold_margin,
    std::optional<bool> & threshold_sensitive,const LocalSpectrum & spectrum)
{
    status=spectrum.status; reason=spectrum.reason;
    if(!spectrum.available) return;
    rank=static_cast<std::size_t>(spectrum.rank);
    if(std::isfinite(spectrum.minimum)) minimum=spectrum.minimum;
    if(std::isfinite(spectrum.maximum)) maximum=spectrum.maximum;
    if(std::isfinite(spectrum.condition)) condition=spectrum.condition;
    if(std::isfinite(spectrum.threshold)) threshold=spectrum.threshold;
    if(std::isfinite(spectrum.threshold_margin_relative)) threshold_margin=spectrum.threshold_margin_relative;
    threshold_sensitive=spectrum.threshold_sensitive;
}
}

LocalConditioningAssessment AssessBcdLocalConditioning(
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
    evidence.local_tau.resize(static_cast<std::size_t>(eta.size()));
    for(Eigen::Index k=0;k<eta.size();++k) evidence.local_tau[static_cast<std::size_t>(k)]=std::exp(eta(k));

    // Endpoint-face columns from constraint status, not the iterative solver's working set.
    std::vector<Eigen::Index> free_columns;
    free_columns.reserve(static_cast<std::size_t>(beta.size()));
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {
        if(beta(2*atom)==0) ++evidence.active_amplitude_count;
        else free_columns.push_back(2*atom);
        free_columns.push_back(2*atom+1);
    }
    evidence.free_ac_columns=free_columns.size();
    const auto state=EvaluateState(domain,y,eta,beta,context);
    if(!state.valid)
    {MarkUnavailable(evidence,state.reason.empty() ? "local-state-evaluation-failed" : state.reason);
        out.seconds=Seconds(started); return out;}
    if(state.certificate.available)
    {
        evidence.local_ac_projected_kkt=state.certificate.projected_kkt;
        if(context.width_bounds.upper)
            evidence.active_width_bound_targets=state.certificate.active_width_bound_count;
    }
    const auto local_width_kkt=ProfileGradientInfinityNorm(state,context);
    if(local_width_kkt) evidence.local_projected_width_kkt=*local_width_kkt;

    const auto p=static_cast<Eigen::Index>(free_columns.size());
    Vector norms(p);
    bool zero_column=false;
    for(Eigen::Index k=0;k<p;++k)
    {
        norms(k)=state.x.col(free_columns[static_cast<std::size_t>(k)]).norm();
        zero_column|=!(norms(k)>0);
    }
    TiledQR ac_factor(p,0);
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
        if(std::isfinite(ac_spectrum.threshold)) evidence.free_ac_rank_threshold=ac_spectrum.threshold;
        if(std::isfinite(ac_spectrum.threshold_margin_relative))
            evidence.free_ac_rank_threshold_margin_relative=ac_spectrum.threshold_margin_relative;
        evidence.free_ac_threshold_sensitive=ac_spectrum.threshold_sensitive;
    }
    evidence.free_ac_conditioning=ConditioningClass(ac_spectrum,zero_column);

    const auto projected_prepared=PrepareProjectedWidthDerivative(state,context.scale,&context);
    if(projected_prepared.valid)
    {
        const auto projected=ReduceProjectedWidthDerivative(projected_prepared);
        if(projected.valid)
        {
            evidence.projected_width_column_norms.assign(projected.projected_norms.data(),
                projected.projected_norms.data()+projected.projected_norms.size());
            const auto spectrum=SpectrumOf(projected.projected,domain.rows,eta.size(),context.rank);
            SetSpectrum(evidence.projected_width_status,evidence.projected_width_reason,
                evidence.projected_width_rank,evidence.projected_width_sigma_min,
                evidence.projected_width_sigma_max,evidence.projected_width_condition,
                evidence.projected_width_rank_threshold,evidence.projected_width_rank_threshold_margin_relative,
                evidence.projected_width_threshold_sensitive,spectrum);
            evidence.projected_width_conditioning=ConditioningClass(spectrum);
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
    if(evidence.projected_width_conditioning=="NotAssessed") evidence.projected_width_conditioning="Unavailable";

    if(evidence.active_amplitude_count)
    {
        evidence.profile_status="active-set-boundary";
        evidence.profile_reason="profile-derivative-may-switch-active-face";
        evidence.profile_conditioning="ActiveFace";
    }
    else if(!ac_spectrum.available || ac_spectrum.rank!=p)
    {
        evidence.profile_status="unavailable";
        evidence.profile_reason="free-ac-rank-deficient";
        evidence.profile_conditioning=ac_spectrum.available ? "Singular" : "Unavailable";
    }
    else
    {
        const auto prepared=PrepareDerivative(state,context.scale,&context);
        if(!prepared.valid)
        {
            evidence.profile_status="unavailable";
            evidence.profile_reason=prepared.reason;
            evidence.profile_conditioning="Unavailable";
        }
        else
        {
            const auto reduced=ReduceDerivative(prepared,state.residual);
            if(!reduced.valid)
            {
                evidence.profile_status="unavailable";
                evidence.profile_reason=reduced.reason;
                evidence.profile_conditioning="Unavailable";
            }
            else
            {
                const auto spectrum=SpectrumOf(reduced.jacobian,domain.rows,eta.size(),context.rank);
                SetSpectrum(evidence.profile_status,evidence.profile_reason,
                    evidence.profile_rank,evidence.profile_sigma_min,evidence.profile_sigma_max,
                    evidence.profile_condition,evidence.profile_rank_threshold,
                    evidence.profile_rank_threshold_margin_relative,evidence.profile_threshold_sensitive,spectrum);
                evidence.profile_conditioning=ConditioningClass(spectrum);
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
    if(evidence.free_ac_conditioning=="Singular" || evidence.projected_width_conditioning=="Singular" ||
        evidence.profile_conditioning=="Singular")
        evidence.conditioning_classification="Singular";
    else if(evidence.active_amplitude_count)
        evidence.conditioning_classification="ActiveFace";
    else if(evidence.free_ac_conditioning=="IllConditioned" ||
        evidence.projected_width_conditioning=="IllConditioned" || evidence.profile_conditioning=="IllConditioned")
        evidence.conditioning_classification="IllConditioned";
    else if(evidence.free_ac_conditioning=="Unavailable" ||
        evidence.projected_width_conditioning=="Unavailable" || evidence.profile_conditioning=="Unavailable")
        evidence.conditioning_classification="Unavailable";
    else evidence.conditioning_classification="WellConditioned";
    out.seconds=Seconds(started); return out;
}

AcCompensationAssessment AssessAcCompensation(
    const Domain & domain,VectorRef y,const Vector & eta,const Vector & beta,const EvaluationContext & context,
    Eigen::Index target,std::size_t resident_bytes)
{
    AcCompensationAssessment out;
    if(target<0 || target>=eta.size() || static_cast<std::size_t>(eta.size())!=domain.atoms.size() ||
        beta.size()!=2*eta.size() || domain.rows!=y.size())
    {out.reason="invalid-local-dimensions"; return out;}
    if(!domain.rows || domain.atoms[static_cast<std::size_t>(target)].empty())
    {out.reason="empty-target-support"; return out;}
    if(!eta.allFinite() || !beta.allFinite() || !y.allFinite() ||
        !std::isfinite(context.scale) || context.scale<=0)
    {out.reason="invalid-local-state"; return out;}
    for(Eigen::Index a=0;a<beta.size();a+=2)
    {
        if(beta(a)<0) {out.reason="infeasible-amplitude"; return out;}
        out.conditional_on_active_face|=beta(a)==0;
    }
    out.amplitude_active=beta(2*target)==0;
    if(context.width_bounds.upper)
    {
        const Vector zero=Vector::Zero(eta.size());
        const auto bound=CalculateWidthBoundKkt(eta,zero,context.width_bounds);
        if(!bound.feasible) {out.reason="infeasible-width-bound"; return out;}
        const auto face=CalculateWidthBoundKkt(eta.segment(target,1),zero.head(1),context.width_bounds);
        out.width_upper_bound_active=face.active_upper!=0;
        out.width_lower_bound_active=face.active_lower!=0;
    }

    MemoryShape shape;
    shape.component_atoms=domain.atoms.size(); shape.component_observations=static_cast<std::size_t>(domain.rows);
    shape.component_free_columns=static_cast<std::size_t>(beta.size());
    shape.maximum_tile_rows=static_cast<std::size_t>(std::min(domain.rows,derivative_tile_rows));
    shape.resident_bytes=resident_bytes;
    for(const auto & support:domain.atoms)
        if(!CheckedAddBytes(shape.component_support_entries,support.size(),shape.component_support_entries))
        {out.reason="size-overflow"; return out;}
    out.memory_estimate=PlanMemory(shape,MemoryStage::LocalConditioning);
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {out.reason="resource-limit:local-conditioning"; return out;}

    const auto state=EvaluateState(domain,y,eta,beta,context);
    if(!state.valid) {out.reason=state.reason; return out;}
    const Eigen::Index g=2*target,k=g+1;
    Vector gaussian=state.x.col(g),charge=state.x.col(k);
    const double gaussian_norm=gaussian.stableNorm(),charge_norm=charge.stableNorm();
    if(std::isfinite(gaussian_norm)) out.gaussian_norm=gaussian_norm;
    if(std::isfinite(charge_norm)) out.charge_norm=charge_norm;
    const auto normalizable=[](double norm) {
        return std::isfinite(norm) && norm>=std::numeric_limits<double>::min();
    };
    const bool valid_g=normalizable(gaussian_norm),valid_k=normalizable(charge_norm);
    if(valid_g) gaussian/=gaussian_norm;
    if(valid_k) charge/=charge_norm;
    if(valid_g && valid_k) out.raw_correlation=std::clamp(gaussian.dot(charge),-1.0,1.0);
    else out.raw_reason="zero-tiny-or-nonfinite-column-norm";

    const double a_signal=beta(g)==0 ? 0 : beta(g)*gaussian_norm;
    const double c_signal=beta(k)==0 ? 0 : beta(k)*charge_norm;
    const double signal_scale=std::max(std::abs(a_signal),std::abs(c_signal));
    const bool valid_signal=std::isfinite(a_signal) && std::isfinite(c_signal) &&
        (beta(g)==0 || valid_g || gaussian_norm==0) &&
        (beta(k)==0 || valid_k || charge_norm==0);
    if(!valid_signal) out.cancellation_reason="unrepresentable-contribution-norm";
    else if(signal_scale==0) out.cancellation_reason="zero-contributions";
    else if(a_signal==0 || c_signal==0) out.cancellation_ratio=1.0;
    else
    {
        const double a=a_signal/signal_scale,c=c_signal/signal_scale;
        const Vector net=a*gaussian+c*charge;
        const double denominator=net.stableNorm(),numerator=std::abs(a)+std::abs(c);
        if(denominator<=8*std::numeric_limits<double>::epsilon()*numerator)
            out.cancellation_reason="net-prediction-at-roundoff";
        else out.cancellation_ratio=std::max(1.0,numerator/denominator);
    }

    Indices columns;
    for(Eigen::Index j=0;j<beta.size();++j)
        if(j%2 || beta(j)>0) columns.push_back(j);
    const Eigen::Index p=static_cast<Eigen::Index>(columns.size());
    const auto position=[&](Eigen::Index j) {
        return static_cast<Eigen::Index>(std::lower_bound(columns.begin(),columns.end(),j)-columns.begin());
    };
    const Eigen::Index g_position=out.amplitude_active ? -1 : position(g),k_position=position(k);
    Vector norms(p);
    // Basis assembly stores compressed columns. Use stable norms without
    // densifying the nuisance design.
    bool usable_nuisance=true;
    for(Eigen::Index j=0;j<p;++j)
    {
        const auto col=columns[static_cast<std::size_t>(j)];
        const auto first=state.x.outerIndexPtr()[col], count=state.x.outerIndexPtr()[col+1]-first;
        norms(j)=count ? Eigen::Map<const Vector>(state.x.valuePtr()+first,count).stableNorm() : 0;
        if(j!=g_position && j!=k_position)
            usable_nuisance&=norms(j)==0 || normalizable(norms(j));
    }
    if(!out.amplitude_active) norms(g_position)=gaussian_norm;
    norms(k_position)=charge_norm;
    Matrix compact;
    {
        TiledQR qr(p,0);
        for(Eigen::Index first=0;first<domain.rows;first+=derivative_tile_rows)
        {
            const auto count=std::min(derivative_tile_rows,domain.rows-first);
            Matrix tile=Matrix::Zero(count,p);
            for(Eigen::Index j=0;j<p;++j)
            {
                if(!normalizable(norms(j))) continue;
                const auto source=columns[static_cast<std::size_t>(j)];
                for(Sparse::InnerIterator entry(state.x,source);entry && entry.row()<first+count;++entry)
                    if(entry.row()>=first) tile(entry.row()-first,j)=entry.value()/norms(j);
            }
            qr.Append(tile,Matrix(count,0));
        }
        compact=std::move(qr.r);
    }
    // With fewer rows than columns, keep the structural null directions too.
    if(compact.rows()<p)
    {
        const auto observed=compact.rows();
        compact.conservativeResize(p,Eigen::NoChange);
        compact.bottomRows(p-observed).setZero();
    }
    const double roundoff=64*std::numeric_limits<double>::epsilon()*static_cast<double>(p);
    const auto threshold_sensitive=[&](const CompactSvdResult & svd) {
        return !CompactRankDecisionSafe(svd) ||
            (svd.singular_values.array()-svd.threshold).abs().minCoeff()<=roundoff*svd.singular_values(0);
    };
    const Eigen::Index q=p-(out.amplitude_active ? 1 : 2);
    Matrix pair(p,out.amplitude_active ? 1 : 2);
    if(!out.amplitude_active) pair.col(0)=compact.col(g_position);
    pair.col(pair.cols()-1)=compact.col(k_position);
    if(!usable_nuisance) out.conditional_reason="tiny-or-nonfinite-nuisance-column-norm";
    else if(q>0)
    {
        Matrix nuisance_transposed(q,p);
        Eigen::Index row{};
        for(Eigen::Index j=0;j<p;++j)
            if(j!=g_position && j!=k_position) nuisance_transposed.row(row++)=compact.col(j).transpose();
        // Right vectors of the transpose span the nuisance column space.
        const auto svd=EvaluateRank(nuisance_transposed,{{domain.rows,q,0},q,-1,RankBoundary::StrictGreater},
            CompactSvdVectors::Right);
        if(svd.valid)
        {
            out.nuisance_rank=svd.rank;
            out.nuisance_rank_threshold_sensitive=threshold_sensitive(svd);
            const auto span=svd.right_vectors.leftCols(svd.rank);
            pair-=(span*(span.transpose()*pair)).eval();
        }
        else out.conditional_reason="nuisance-factorization-failed";
    }
    else out.nuisance_rank=0;
    if(out.nuisance_rank.has_value())
    {
        const double gn=out.amplitude_active ? 0 : pair.col(0).norm(),kn=pair.col(pair.cols()-1).norm();
        if(!out.amplitude_active && valid_g) out.residualized_gaussian_norm=gn*gaussian_norm;
        if(valid_k) out.residualized_charge_norm=kn*charge_norm;
        if(out.amplitude_active) out.conditional_reason="conditional-on-active-amplitude-face";
        else if(!valid_g || !valid_k) out.conditional_reason=out.raw_reason;
        else if(gn<=roundoff || kn<=roundoff) out.conditional_reason="residualized-column-at-roundoff";
        else out.conditional_correlation=std::clamp(pair.col(0).dot(pair.col(1))/gn/kn,-1.0,1.0);
    }

    const auto svd=EvaluateRank(compact,{{domain.rows,p,0},p,-1,RankBoundary::StrictGreater},
        CompactSvdVectors::Right);
    if(!svd.valid) {out.weak_reason="free-ac-factorization-failed"; return out;}
    out.free_rank=svd.rank; out.singular_values=svd.singular_values;
    out.free_rank_threshold_sensitive=threshold_sensitive(svd);
    if(!norms.unaryExpr(normalizable).all())
    {out.weak_reason="zero-tiny-or-nonfinite-free-column-norm"; return out;}
    const double cutoff=std::max(svd.threshold,svd.singular_values(p-1)+roundoff*svd.singular_values(0));
    out.weak_dimension=(svd.singular_values.array()<=cutoff).count();
    const auto weak=svd.right_vectors.rightCols(out.weak_dimension);
    out.weak_atom_shares=Vector::Zero(eta.size());
    for(Eigen::Index j=0;j<p;++j)
        out.weak_atom_shares(columns[static_cast<std::size_t>(j)]/2)+=
            weak.row(j).squaredNorm()/static_cast<double>(out.weak_dimension);
    if(a_signal>0 && c_signal<0 && valid_signal)
        out.weak_cancellation_alignment=(a_signal/signal_scale*weak.row(g_position)+
            c_signal/signal_scale*weak.row(k_position)).norm()/std::hypot(a_signal/signal_scale,c_signal/signal_scale);
    if(out.weak_dimension==1)
    {
        out.physical_weak_direction=Vector::Zero(beta.size());
        for(Eigen::Index j=0;j<p;++j)
            out.physical_weak_direction(columns[static_cast<std::size_t>(j)])=weak(j,0)/norms(j);
        if(!out.physical_weak_direction.allFinite())
        {
            out.physical_weak_direction.resize(0);
            out.weak_reason="unrepresentable-physical-direction";
        }
    }
    return out;
}
}
