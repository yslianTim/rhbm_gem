#include "TargetEvidence.hpp"
#include "TiledDerivative.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <algorithm>
#include <chrono>
#include <numeric>

namespace rhbm_gem::core::joint_component {
namespace {
using Status=JointCheckStatus;
std::vector<double> Values(const Vector & v) {return {v.data(),v.data()+v.size()};}
bool Sensitive(const CompactSvdResult & svd)
{
    return svd.threshold>0 && ((svd.singular_values.array()>=.5*svd.threshold) &&
        (svd.singular_values.array()<=2*svd.threshold)).any();
}
Matrix SquareCompact(const Matrix & r)
{
    Matrix out=Matrix::Zero(r.cols(),r.cols()); out.topRows(r.rows())=r; return out;
}
}
TargetGeometry BuildTargetGeometry(const Domain & domain,VectorRef y,const JointState & state,const EvaluationContext & context)
{
    TargetGeometry out;
    const auto m=static_cast<Eigen::Index>(state.b.size());
    if(!m || state.ac.size()!=static_cast<std::size_t>(2*m) || state.log_b.size()!=static_cast<std::size_t>(m))
    {out.reason="missing-full-parameter-state"; return out;}
    const VectorMap eta(state.log_b.data(),m),beta(state.ac.data(),2*m);
    const Vector widths=eta.array().exp();
    for(Eigen::Index a=0;a<m;++a)
        if(!std::isfinite(state.b[static_cast<std::size_t>(a)]) || !(widths(a)>0) ||
            std::abs(widths(a)-state.b[static_cast<std::size_t>(a)])>8*std::numeric_limits<double>::epsilon()*widths(a))
        {out.reason="inconsistent-saved-widths"; return out;}
    out.boundary=false;
    for(Eigen::Index a=0;a<m;++a) out.boundary |= beta(2*a)==0;
    if(out.boundary) {out.reason="active-amplitude-boundary"; return out;}
    out.endpoint=EvaluateState(domain,y,eta,beta,context);
    if(!out.endpoint.valid || !out.endpoint.certificate.kkt_passed)
    {out.reason="invalid-saved-state"; return out;}
    out.boundary=!out.endpoint.certificate.active_atoms.empty();
    if(out.boundary) {out.reason="active-amplitude-boundary"; return out;}
    const auto & e=out.endpoint;
    Sparse z(y.size(),3*m); std::vector<Eigen::Triplet<double>> entries;
    for(Eigen::Index a=0;a<m;++a)
    {
        for(int k=0;k<2;++k)
        {
            for(Sparse::InnerIterator v(e.x,2*a+k);v;++v) entries.emplace_back(v.row(),3*a+k,v.value());
            for(Sparse::InnerIterator v(e.derivative,2*a+k);v;++v)
                entries.emplace_back(v.row(),3*a+2,v.value()*beta(2*a+k));
        }
    }
    z.setFromTriplets(entries.begin(),entries.end()); out.scales.resize(z.cols());
    for(Eigen::Index k=0;k<z.cols();++k)
    {
        out.scales(k)=z.col(k).norm();
        if(!std::isfinite(out.scales(k))) {out.reason="nonfinite-jacobian"; return out;}
        // An exact zero column stays zero; its coordinate is explicitly unidentified.
        if(out.scales(k)==0) out.scales(k)=1;
        z.col(k)/=out.scales(k);
    }
    TiledQR qr(z.cols(),0);
    for(Eigen::Index first=0;first<z.rows();first+=derivative_tile_rows)
    {
        const auto count=std::min(derivative_tile_rows,z.rows()-first);
        qr.Append(Matrix(z.middleRows(first,count)),Matrix(count,0));
    }
    out.svd=EvaluateRank(SquareCompact(qr.r),{context.rank,z.cols()},nullptr,CompactSvdVectors::Right);
    if(!out.svd.valid) {out.reason="jacobian-factorization-failed"; return out;}
    out.threshold_sensitive=Sensitive(out.svd);
    out.valid=true; out.reason=out.threshold_sensitive ? "rank-threshold-sensitive" : ""; return out;
}
JointTargetEvidence AssessTargetState(const JointProblem & problem,const JointParameterLayout & layout,
    const std::optional<JointState> & state,const std::vector<JointCheck> & existing,JointEvidenceScope scope,std::size_t original_rows)
{
    const auto & data=JointProblemAccess::Get(problem); const auto & input=problem.Input();
    JointTargetEvidence out; out.original_rows=original_rows;
    auto check=[&](const char * name,Status status,std::string reason,std::optional<double> value={},std::optional<double> threshold={}) {
        out.checks.push_back({name,status,scope,value,threshold,std::move(reason)});
    };
    for(const auto & c:existing) if(c.name=="inner" || c.name=="kkt" || c.name=="width-stationarity" || c.name=="assembled-profile") out.checks.push_back(c);
    TargetGeometry geometry;
    auto context=ProfileContext(data.context,layout,static_cast<Eigen::Index>(original_rows));
    const auto domain=ProfileDomain(data.domain,layout);
    const Vector y=SelectValues(data.y,{layout.informative_rows.begin(),layout.informative_rows.end()});
    if(state) geometry=BuildTargetGeometry(domain,y,*state,context);
    else geometry.reason="missing-component-state";
    const bool usable=geometry.valid && !geometry.threshold_sensitive;
    bool targets_identified=true; std::size_t targets{};
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
    {
        const auto atom=layout.full_atoms[k];
        JointAtomIdentifiability evidence; evidence.atom=atom; evidence.reason=geometry.reason;
        if(usable)
        {
            const auto & svd=geometry.svd;
            const double leakage=svd.right_vectors.block(3*static_cast<Eigen::Index>(k),svd.rank,3,
                svd.right_vectors.cols()-svd.rank).norm();
            evidence.null_space_leakage=leakage;
            evidence.status=leakage<=1e-10 ? Status::Passed : Status::Failed;
            evidence.reason=evidence.status==Status::Passed ? "" : "nonunique-parameter-directions";
        }
        if(std::binary_search(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),atom))
        {++targets; targets_identified &= evidence.status==Status::Passed;}
        out.atoms.push_back(std::move(evidence));
    }
    for(const auto & group:layout.groups) for(auto atom:group.atoms)
        out.atoms.push_back({atom,Status::Unavailable,{},"observable-contribution-only"});
    std::sort(out.atoms.begin(),out.atoms.end(),[](const auto & a,const auto & b){return a.atom<b.atom;});
    if(geometry.valid)
    {
        out.column_scales=Values(geometry.scales);
        out.ranks.push_back({"column-scaled-observation-jacobian",scope,static_cast<std::size_t>(geometry.svd.rank),
            geometry.svd.threshold,Values(geometry.svd.singular_values)});
    }
    // A nuisance-only component has no target coordinates, but its state is still required.
    if(!targets)
    {
        check("target-identifiability",state ? Status::Passed : Status::Unavailable,state ? "no-target-coordinates" : geometry.reason);
        check("target-local-correction",state ? Status::Passed : Status::Unavailable,state ? "no-target-coordinates" : geometry.reason);
        return out;
    }
    check("target-identifiability",!usable ? Status::Unavailable : targets_identified ? Status::Passed : Status::Failed,
        !usable ? geometry.reason : targets_identified ? "" : "target-confounded-with-nuisance");
    if(!usable || !targets_identified)
    {check("target-local-correction",Status::Unavailable,!usable ? geometry.reason : "target-unidentified"); return out;}
    const auto & e=geometry.endpoint;
    const auto prepared=PrepareDerivative(e,context.scale,&context);
    const auto reduced=ReduceDerivative(prepared,e.residual,false);
    if(!reduced.valid) {check("target-local-correction",Status::Unavailable,reduced.reason); return out;}
    Matrix compact=SquareCompact(reduced.jacobian);
    Vector rhs=Vector::Zero(compact.rows()); rhs.head(reduced.response.size())=-reduced.response;
    const auto step=EvaluateRank(compact,{context.rank,compact.cols()},&rhs,CompactSvdVectors::Right);
    if(step.valid) out.ranks.push_back({"undamped-full-profile-step",scope,static_cast<std::size_t>(step.rank),step.threshold,Values(step.singular_values)});
    if(!step.valid || Sensitive(step))
    {check("target-local-correction",Status::Unavailable,step.valid ? "profile-step-rank-threshold-sensitive" : "profile-step-rank-unavailable"); return out;}
    double correction{},leakage{};
    for(std::size_t a=0;a<layout.full_atoms.size();++a)
    {
        if(!std::binary_search(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),layout.full_atoms[a])) continue;
        Matrix map=Matrix::Zero(3,compact.cols());
        for(int k=0;k<2;++k)
        {
            const auto col=2*static_cast<Eigen::Index>(a)+k;
            map.row(k)=-(prepared.coefficients.row(col)+prepared.correction.row(col))/
                (e.x.col(col).norm()*(1+std::abs(e.beta(col))));
        }
        map(2,static_cast<Eigen::Index>(a))=1;
        correction=std::max(correction,(map*step.solution).lpNorm<Eigen::Infinity>());
        const Matrix null=map*step.right_vectors.rightCols(step.right_vectors.cols()-step.rank);
        leakage=std::max(leakage,null.norm()/std::max(1.,map.norm()));
    }
    check("target-local-correction",leakage>1e-10 ? Status::Unavailable : correction<=1e-10 ? Status::Passed : Status::Failed,
        leakage>1e-10 ? "nonunique-target-profile-step" : "",correction,1e-10);
    check("target-step-null-space",leakage<=1e-10 ? Status::Passed : Status::Unavailable,
        leakage<=1e-10 ? "" : "nonunique-target-profile-step",leakage,1e-10);
    return out;
}
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> & evidence,bool state)
{
    if(!evidence) return Status::NotRun;
    if(!state) return Status::Unavailable;
    auto status=Status::Passed;
    for(const char * name:{"inner","kkt","width-stationarity","target-identifiability","target-local-correction"})
    {
        const auto found=std::find_if(evidence->checks.begin(),evidence->checks.end(),[&](const auto & c){return c.name==name;});
        status=MergeConvergenceStatus(status,found==evidence->checks.end() ? Status::Unavailable : found->status);
    }
    return status;
}
void AddTargetEvidence(JointFitResult & fit)
{
    if(!fit.problem || !fit.problem->Input().selection_domain) return;
    const auto started=std::chrono::steady_clock::now();
    const auto & problem=*fit.problem;
    for(auto & c:fit.components)
    {
        const JointParameterLayout layout=c.layout.value_or(JointParameterLayout{c.atoms,c.rows,{}});
        c.target_evidence=AssessTargetState(problem,layout,c.state,c.evidence,JointEvidenceScope::ComponentLocal,c.rows.size());
    }
    fit.target_evidence=AssessTargetState(problem,problem.ParameterLayout(),fit.assembled_state,fit.evidence,
        JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    fit.costs.assessment_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
JointFitResult FitWithSearchPolicy(const JointProblem & problem,const std::vector<double> & b,const SearchPolicy & policy)
{
    eigen_helper::ScopedEigenThreadCount eigen_thread_guard{1};
    auto out=FitWithSearchPolicyImpl(problem,b,policy); AddTargetEvidence(out); return out;
}
}
namespace rhbm_gem::core {
JointCheckStatus JointComponentResult::TargetRuntimeConvergence() const
{return joint_component::TargetConvergence(target_evidence,state.has_value());}
JointCheckStatus JointFitResult::TargetRuntimeConvergence() const
{
    if(!target_evidence) return JointCheckStatus::NotRun;
    if(!assembled_state || !prediction || !objective || components.empty() ||
        !std::all_of(available_row_mask.begin(),available_row_mask.end(),[](bool b){return b;})) return JointCheckStatus::Unavailable;
    auto status=joint_component::TargetConvergence(target_evidence,true);
    for(const auto & c:components) status=joint_component::MergeConvergenceStatus(status,c.TargetRuntimeConvergence());
    const auto found=std::find_if(evidence.begin(),evidence.end(),[](const auto & c){return c.name=="assembled-profile";});
    return joint_component::MergeConvergenceStatus(status,found==evidence.end() ? JointCheckStatus::Unavailable : found->status);
}
}
