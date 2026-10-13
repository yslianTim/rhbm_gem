#include "TargetEvidence.hpp"
#include "TiledQR.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <numeric>
#include <sstream>

namespace rhbm_gem::core::joint_component {
namespace {
using Status=JointCheckStatus;
std::vector<double> Values(const Vector & v) {return {v.data(),v.data()+v.size()};}
bool Sensitive(const CompactSvdResult & svd)
{
    return svd.threshold>0 && ((svd.singular_values.array()>=.5*svd.threshold) &&
        (svd.singular_values.array()<=2*svd.threshold)).any();
}
bool SameState(const JointState & a,const JointState & b)
{
    return a.ac==b.ac && a.b==b.b && a.log_b==b.log_b && a.width_gradient==b.width_gradient &&
        a.objective==b.objective && a.nuisance_amplitudes==b.nuisance_amplitudes;
}
bool FullRowIdentity(const std::vector<std::size_t> & rows,std::size_t count)
{
    if(rows.size()!=count) return false;
    for(std::size_t i=0;i<count;++i) if(rows[i]!=i) return false;
    return true;
}
void HashBytes(std::uint64_t & hash,const void * data,std::size_t size)
{
    const auto * bytes=static_cast<const unsigned char *>(data);
    for(std::size_t i=0;i<size;++i) {hash^=bytes[i]; hash*=1099511628211ULL;}
}
void HashWord(std::uint64_t & hash,std::uint64_t encoded)
{
    for(std::size_t i=0;i<sizeof(encoded);++i)
    {const auto byte=static_cast<unsigned char>(encoded&0xffU); HashBytes(hash,&byte,1); encoded>>=8;}
}
void HashSize(std::uint64_t & hash,std::size_t value)
{HashWord(hash,static_cast<std::uint64_t>(value));}
void HashText(std::uint64_t & hash,const std::string & value)
{HashSize(hash,value.size()); HashBytes(hash,value.data(),value.size());}
void HashDouble(std::uint64_t & hash,double value)
{HashWord(hash,std::bit_cast<std::uint64_t>(value));}
bool SameGroup(const JointNuisanceGroup & a,const JointNuisanceGroup & b)
{return a.row==b.row && a.atoms==b.atoms;}
bool SameLayout(const JointParameterLayout & a,const JointParameterLayout & b)
{
    if(a.full_atoms!=b.full_atoms || a.informative_rows!=b.informative_rows || a.groups.size()!=b.groups.size()) return false;
    for(std::size_t i=0;i<a.groups.size();++i) if(!SameGroup(a.groups[i],b.groups[i])) return false;
    return true;
}
bool LayoutMatchesScope(const JointProblem & problem,const JointParameterLayout & layout,
    JointEvidenceScope scope,std::size_t rows)
{
    const auto & data=JointProblemAccess::Get(problem);
    if(scope==JointEvidenceScope::AssembledGlobal)
        return rows==data.input->observations.size() && SameLayout(layout,data.layout);
    for(const auto & view:data.partition.components)
    {
        if(rows!=view.rows.size()) continue;
        JointParameterLayout expected;
        for(auto atom:data.layout.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(atom))>=0) expected.full_atoms.push_back(atom);
        for(auto row:data.layout.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(row))>=0) expected.informative_rows.push_back(row);
        for(const auto & group:data.layout.groups)
            if(view.LocalRow(static_cast<Eigen::Index>(group.row))>=0) expected.groups.push_back(group);
        if(SameLayout(layout,expected)) return true;
    }
    return false;
}
const JointCheck * FindCheck(const std::vector<JointCheck> & checks,const char * name,JointEvidenceScope scope)
{
    const auto found=std::find_if(checks.begin(),checks.end(),[&](const auto & check) {
        return check.name==name && check.scope==scope;
    });
    return found==checks.end() ? nullptr : &*found;
}
bool RankThresholdSensitive(const JointRankEvidence & evidence)
{
    if(evidence.singular_values.empty() || !(evidence.threshold>0)) return true;
    const double lower=std::nextafter(evidence.threshold,0.);
    const double upper=std::nextafter(evidence.threshold,std::numeric_limits<double>::infinity());
    return std::any_of(evidence.singular_values.begin(),evidence.singular_values.end(),[&](double value) {
        return value>=lower && value<=upper;
    });
}
bool ValidRankEvidence(const JointRankEvidence & evidence,std::size_t columns,std::size_t rows)
{
    if(evidence.rank!=columns || evidence.singular_values.size()!=columns || !columns || !rows ||
        !(evidence.threshold>0) || !std::isfinite(evidence.threshold)) return false;
    if(!std::all_of(evidence.singular_values.begin(),evidence.singular_values.end(),[](double value) {
        return std::isfinite(value) && value>0;
    })) return false;
    for(std::size_t i=1;i<evidence.singular_values.size();++i)
        if(evidence.singular_values[i]>evidence.singular_values[i-1]) return false;
    const double expected=std::numeric_limits<double>::epsilon()*static_cast<double>(std::max(rows,columns)) *
        evidence.singular_values.front();
    return std::abs(evidence.threshold-expected)<=8*std::numeric_limits<double>::epsilon()*
        std::max({std::numeric_limits<double>::min(),std::abs(expected),std::abs(evidence.threshold)});
}
Matrix SquareCompact(const Matrix & r)
{
    Matrix out=Matrix::Zero(r.cols(),r.cols()); out.topRows(r.rows())=r; return out;
}
MemoryEstimate TargetIdentifiabilityMemoryEstimateImpl(const JointProblemInput & input,const JointParameterLayout & layout,
    std::size_t rows)
{
    MemoryShape shape; shape.component_atoms=layout.full_atoms.size(); shape.component_observations=rows;
    const auto overflow=[&] {
        auto estimate=PlanMemory(shape,MemoryStage::TargetIdentifiability);
        estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
    };
    if(!CheckedMultiplyBytes(shape.component_atoms,2,shape.component_design_columns) ||
        !CheckedMultiplyBytes(shape.component_atoms,3,shape.component_jacobian_columns)) return overflow();
    shape.component_width_columns=shape.component_atoms;
    shape.maximum_tile_rows=derivative_tile_rows;
    shape.component_free_columns=shape.component_design_columns;
    for(const auto atom:layout.full_atoms)
        if(input.selection_domain && std::binary_search(input.selection_domain->target_indices.begin(),
            input.selection_domain->target_indices.end(),atom))
            if(!CheckedAddBytes(shape.component_target_atoms,1,shape.component_target_atoms)) return overflow();
    shape.component_halo_atoms=shape.component_atoms-shape.component_target_atoms;
    if(!CheckedMultiplyBytes(shape.component_target_atoms,3,shape.component_target_columns)) return overflow();
    for(const auto atom:layout.full_atoms)
        if(atom>=input.support.size() ||
            !CheckedAddBytes(shape.component_support_entries,input.support[atom].size(),shape.component_support_entries))
            return overflow();
    const auto snapshot=EstimateSnapshotBytes(input);
    std::size_t resident{};
    const auto add_product=[](std::size_t & total,std::size_t count,std::size_t width) {
        std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
    };
    if(!snapshot) return overflow();
    resident=*snapshot;
    if(!add_product(resident,rows,sizeof(double)) ||
        !add_product(resident,shape.component_atoms,4*sizeof(double)) ||
        !add_product(resident,rows,2*sizeof(Eigen::Index)) ||
        !add_product(resident,shape.component_atoms,2*sizeof(Eigen::Index)+sizeof(std::vector<LocalSupportIndex>)) ||
        !add_product(resident,shape.component_support_entries,sizeof(LocalSupportIndex))) return overflow();
    shape.resident_bytes=resident;
    return PlanMemory(shape,MemoryStage::TargetIdentifiability);
}
}
MemoryEstimate TargetIdentifiabilityMemoryEstimate(const JointProblemInput & input,const JointParameterLayout & layout,
    std::size_t rows)
{return TargetIdentifiabilityMemoryEstimateImpl(input,layout,rows);}
std::string TargetRankContextIdentity(const JointProblemInput & input,const JointParameterLayout & layout,
    const JointState & state,std::size_t rank_rows)
{
    const auto atoms=layout.full_atoms.size();
    if(atoms>std::numeric_limits<std::size_t>::max()/2 || state.ac.size()!=2*atoms ||
        state.b.size()!=atoms || state.log_b.size()!=atoms ||
        !std::all_of(state.ac.begin(),state.ac.end(),[](double value){return std::isfinite(value);}) ||
        !std::all_of(state.b.begin(),state.b.end(),[](double value){return std::isfinite(value) && value>0;}) ||
        !std::all_of(state.log_b.begin(),state.log_b.end(),[](double value){return std::isfinite(value);})) return {};
    std::uint64_t hash=14695981039346656037ULL;
    HashText(hash,"derived-target-rank-context-v1"); HashSize(hash,rank_rows);
    HashSize(hash,input.atom_ids.size()); HashSize(hash,input.row_ids.size());
    HashSize(hash,layout.full_atoms.size());
    for(const auto atom:layout.full_atoms)
    {
        if(atom>=input.atom_ids.size() || atom>=input.support.size()) return {};
        HashSize(hash,atom); HashText(hash,input.atom_ids[atom]);
        HashSize(hash,input.support[atom].size());
        for(const auto & support:input.support[atom])
        {HashSize(hash,support.row); HashDouble(hash,support.squared_distance);}
    }
    HashSize(hash,layout.informative_rows.size());
    for(const auto row:layout.informative_rows)
    {
        if(row>=input.row_ids.size() || row>=input.observations.size()) return {};
        HashSize(hash,row); HashText(hash,input.row_ids[row]); HashDouble(hash,input.observations[row]);
    }
    HashSize(hash,layout.groups.size());
    for(const auto & group:layout.groups)
    {
        if(group.row>=input.row_ids.size() || group.row>=input.observations.size()) return {};
        HashSize(hash,group.row); HashText(hash,input.row_ids[group.row]); HashDouble(hash,input.observations[group.row]);
        HashSize(hash,group.atoms.size());
        for(const auto atom:group.atoms)
        {
            if(atom>=input.atom_ids.size() || atom>=input.support.size()) return {};
            HashSize(hash,atom); HashText(hash,input.atom_ids[atom]); HashSize(hash,input.support[atom].size());
            for(const auto & support:input.support[atom])
            {HashSize(hash,support.row); HashDouble(hash,support.squared_distance);}
        }
    }
    HashSize(hash,state.ac.size()); for(const auto value:state.ac) HashDouble(hash,value);
    HashSize(hash,state.b.size()); for(const auto value:state.b) HashDouble(hash,value);
    HashSize(hash,state.log_b.size()); for(const auto value:state.log_b) HashDouble(hash,value);
    if(input.selection_domain)
    {
        HashText(hash,input.selection_domain->contract); HashText(hash,input.selection_domain->contributor_policy);
        HashDouble(hash,input.selection_domain->observation_radius); HashDouble(hash,input.selection_domain->support_radius);
        HashSize(hash,input.selection_domain->target_indices.size());
        for(const auto atom:input.selection_domain->target_indices) HashSize(hash,atom);
    }
    else HashText(hash,"no-selection-domain");
    std::ostringstream out; out<<"fnv1a64:"<<std::hex<<std::setw(16)<<std::setfill('0')<<hash; return out.str();
}
JointTargetEvidence DeriveTargetIdentifiability(const JointProblem & problem,const JointParameterLayout & layout,
    const std::optional<JointState> & state,const std::vector<JointCheck> & existing,
    const std::vector<JointRankEvidence> & ranks,JointEvidenceScope scope,std::size_t original_rows)
{
    const auto & input=problem.Input();
    JointTargetEvidence out; out.contract="target-estimability-v3"; out.original_rows=original_rows;
    for(const auto & check:existing)
        if(check.name=="conditional-ac-optimality" || check.name=="independent-replay" ||
            check.name=="width-stationarity" || check.name=="design-identifiability" ||
            check.name=="numerical-identifiability" || check.name=="assembly-mapping" ||
            check.name=="assembled-consistency" || check.name=="inner" || check.name=="kkt" ||
            check.name=="assembled-profile") out.checks.push_back(check);
    const auto targets=input.selection_domain ? &input.selection_domain->target_indices : nullptr;
    bool has_target=false;
    for(const auto atom:layout.full_atoms)
    {
        const bool target=targets && std::binary_search(targets->begin(),targets->end(),atom);
        has_target|=target;
        out.atoms.push_back({atom,JointCheckStatus::Unavailable,{},"missing-trusted-component-state"});
    }
    for(const auto & group:layout.groups) for(const auto atom:group.atoms)
        out.atoms.push_back({atom,JointCheckStatus::Unavailable,{},"observable-contribution-only"});
    std::sort(out.atoms.begin(),out.atoms.end(),[](const auto & a,const auto & b){return a.atom<b.atom;});
    auto target_check=[&](JointCheckStatus status,std::string reason) {
        out.checks.push_back({"target-identifiability",status,scope,{}, {},std::move(reason)});
    };
    auto mark_full=[&](JointCheckStatus status,const std::string & reason) {
        for(auto & atom:out.atoms) if(std::binary_search(layout.full_atoms.begin(),layout.full_atoms.end(),atom.atom))
        {atom.status=status; atom.reason=reason;}
    };
    if(!targets) {mark_full(JointCheckStatus::Unavailable,"missing-selection-metadata");
        target_check(JointCheckStatus::Unavailable,"missing-selection-metadata"); return out;}
    if(scope==JointEvidenceScope::AssembledGlobal && original_rows!=input.observations.size())
    {mark_full(JointCheckStatus::Unavailable,"observation-domain-mismatch");
        target_check(JointCheckStatus::Unavailable,"observation-domain-mismatch"); return out;}
    if(!LayoutMatchesScope(problem,layout,scope,original_rows))
    {mark_full(JointCheckStatus::Unavailable,"wrong-parameter-layout"); target_check(JointCheckStatus::Unavailable,"wrong-parameter-layout"); return out;}
    if(!state) {mark_full(JointCheckStatus::Unavailable,"missing-component-state");
        target_check(JointCheckStatus::Unavailable,"missing-component-state"); return out;}
    const auto columns=layout.full_atoms.size();
    if(!columns)
    {
        target_check(has_target ? JointCheckStatus::Unavailable : JointCheckStatus::Passed,
            has_target ? "target-not-in-full-abc-layout" : "no-target-coordinates");
        return out;
    }
    if(state->ac.size()!=2*columns || state->b.size()!=columns || state->log_b.size()!=columns ||
        state->width_gradient.size()!=columns || !std::isfinite(state->objective) ||
        !std::all_of(state->ac.begin(),state->ac.end(),[](double value){return std::isfinite(value);}) ||
        !std::all_of(state->b.begin(),state->b.end(),[](double value){return std::isfinite(value) && value>0;}) ||
        !std::all_of(state->log_b.begin(),state->log_b.end(),[](double value){return std::isfinite(value);}) ||
        !std::all_of(state->width_gradient.begin(),state->width_gradient.end(),[](double value){return std::isfinite(value);}))
    {mark_full(JointCheckStatus::Unavailable,"nonfinite-parameter-state");
        target_check(JointCheckStatus::Unavailable,"nonfinite-parameter-state"); return out;}
    for(std::size_t atom=0;atom<columns;++atom)
        if(std::abs(std::exp(state->log_b[atom])-state->b[atom])>
            8*std::numeric_limits<double>::epsilon()*state->b[atom])
        {mark_full(JointCheckStatus::Unavailable,"inconsistent-saved-widths");
            target_check(JointCheckStatus::Unavailable,"inconsistent-saved-widths"); return out;}
    for(std::size_t atom=0;atom<columns;++atom) if(state->ac[2*atom]==0)
    {mark_full(JointCheckStatus::Unavailable,"active-amplitude-boundary");
        target_check(JointCheckStatus::Unavailable,"active-amplitude-boundary"); return out;}
    for(std::size_t atom=0;atom<columns;++atom) if(state->ac[2*atom]<0)
    {mark_full(JointCheckStatus::Unavailable,"invalid-parameter-state");
        target_check(JointCheckStatus::Unavailable,"invalid-parameter-state"); return out;}
    const auto context_identity=TargetRankContextIdentity(input,layout,*state,original_rows);
    if(context_identity.empty()) {mark_full(JointCheckStatus::Unavailable,"missing-rank-context-identity");
        target_check(JointCheckStatus::Unavailable,"missing-rank-context-identity"); return out;}
    const auto find_rank=[&](const char * name)->const JointRankEvidence * {
        const JointRankEvidence * result{};
        for(const auto & rank:ranks) if(rank.name==name && rank.scope==scope)
        {if(result) return nullptr; result=&rank;}
        return result;
    };
    const auto * design=find_rank("design"); const auto * widths=find_rank("projected-width");
    if(!design || !widths) {mark_full(JointCheckStatus::Unavailable,"missing-rank-evidence");
        target_check(JointCheckStatus::Unavailable,"missing-rank-evidence"); return out;}
    out.ranks={*design,*widths};
    if(design->context_identity!=context_identity || widths->context_identity!=context_identity)
    {mark_full(JointCheckStatus::Unavailable,"rank-context-mismatch");
        target_check(JointCheckStatus::Unavailable,"rank-context-mismatch"); return out;}
    if(design->rank<2*columns || widths->rank<columns)
    {mark_full(JointCheckStatus::Unavailable,"component-rank-deficient");
        target_check(JointCheckStatus::Unavailable,"component-rank-deficient"); return out;}
    if(!ValidRankEvidence(*design,2*columns,original_rows) || !ValidRankEvidence(*widths,columns,original_rows))
    {mark_full(JointCheckStatus::Unavailable,"inconsistent-rank-evidence");
        target_check(JointCheckStatus::Unavailable,"inconsistent-rank-evidence"); return out;}
    if(RankThresholdSensitive(*design) || RankThresholdSensitive(*widths))
    {mark_full(JointCheckStatus::Unavailable,"rank-threshold-sensitive");
        target_check(JointCheckStatus::Unavailable,"rank-threshold-sensitive"); return out;}
    for(const char * name:{"conditional-ac-optimality","independent-replay","width-stationarity",
        "design-identifiability","numerical-identifiability"})
    {
        const auto * check=FindCheck(existing,name,scope);
        if(!check || check->status!=JointCheckStatus::Passed)
        {
            const std::string reason=check && check->reason.find("threshold-sensitive")!=std::string::npos ?
                "rank-threshold-sensitive" : !check ? "missing-component-qualification" : "component-not-qualified";
            mark_full(JointCheckStatus::Unavailable,reason); target_check(JointCheckStatus::Unavailable,reason); return out;
        }
    }
    if(scope==JointEvidenceScope::AssembledGlobal)
        for(const char * name:{"assembly-mapping","assembled-consistency"})
        {
            const auto * check=FindCheck(existing,name,scope);
            if(!check || check->status!=JointCheckStatus::Passed)
            {mark_full(JointCheckStatus::Unavailable,"assembly-not-qualified");
                target_check(JointCheckStatus::Unavailable,"assembly-not-qualified"); return out;}
        }
    for(auto & atom:out.atoms) if(std::binary_search(layout.full_atoms.begin(),layout.full_atoms.end(),atom.atom))
    {atom.status=JointCheckStatus::Passed; atom.reason="derived-from-component-identifiability";}
    target_check(JointCheckStatus::Passed,"derived-from-component-identifiability");
    return out;
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
    JointTargetEvidence out; out.contract="target-estimability-v2"; out.original_rows=original_rows;
    auto check=[&](const char * name,Status status,std::string reason,std::optional<double> value={},std::optional<double> threshold={}) {
        out.checks.push_back({name,status,scope,value,threshold,std::move(reason)});
    };
    for(const auto & c:existing)
        if(c.name=="conditional-ac-optimality" || c.name=="independent-replay" ||
            c.name=="width-stationarity" || c.name=="design-identifiability" ||
            c.name=="numerical-identifiability" ||
            c.name=="assembly-mapping" || c.name=="assembled-consistency" ||
            c.name=="inner" || c.name=="kkt" || c.name=="assembled-profile")
            out.checks.push_back(c);
    MemoryEstimate target_memory;
    if(state)
    {
        target_memory=TargetIdentifiabilityMemoryEstimate(input,layout,layout.informative_rows.size());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting preflight_stage("target-memory-preflight",
            static_cast<Eigen::Index>(layout.informative_rows.size()),
            static_cast<Eigen::Index>(target_memory.shape.component_jacobian_columns),&target_memory);
        preflight_stage.Finish();
#endif
        if(target_memory.status!=MemoryStatus::Pass || target_memory.estimated_peak_bytes>target_memory.budget_bytes/2)
            LogMemoryEstimate(target_memory);
        if(target_memory.status!=MemoryStatus::Pass)
        {
            const auto reason=std::string("resource-limit:target-identifiability");
            for(const auto atom:layout.full_atoms) out.atoms.push_back({atom,Status::Unavailable,{},reason});
            for(const auto & group:layout.groups) for(const auto atom:group.atoms)
                out.atoms.push_back({atom,Status::Unavailable,{},"observable-contribution-only"});
            std::sort(out.atoms.begin(),out.atoms.end(),[](const auto & a,const auto & b){return a.atom<b.atom;});
            check("target-identifiability",Status::Unavailable,reason);
            check("target-memory-preflight",Status::Unavailable,
                std::string("resource-limit:target-identifiability:")+target_memory.reason,
                static_cast<double>(target_memory.estimated_peak_bytes),static_cast<double>(target_memory.budget_bytes));
            return out;
        }
    }
    TargetGeometry geometry;
    auto context=ProfileContext(data.context,layout,static_cast<Eigen::Index>(original_rows));
    const auto domain=ProfileDomain(data.domain,layout);
    const Vector y=SelectValues(data.y,{layout.informative_rows.begin(),layout.informative_rows.end()});
    if(state)
    {
        const auto geometry_memory=PlanMemory(target_memory.shape,MemoryStage::TargetGeometry,target_memory.budget_bytes);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting geometry_stage("target-geometry",y.size(),
            static_cast<Eigen::Index>(target_memory.shape.component_jacobian_columns),&geometry_memory);
#endif
        geometry=BuildTargetGeometry(domain,y,*state,context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        geometry_stage.Finish();
#endif
    }
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
        return out;
    }
    check("target-identifiability",!usable ? Status::Unavailable : targets_identified ? Status::Passed : Status::Failed,
        !usable ? geometry.reason : targets_identified ? "" : "target-confounded-with-nuisance");
    return out;
}
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> & evidence,bool state)
{
    if(!evidence) return Status::NotRun;
    if(!state) return Status::Unavailable;
    if(evidence->contract!="target-estimability-v1" && evidence->contract!="target-estimability-v2" &&
        evidence->contract!="target-estimability-v3")
        return Status::Unavailable;
    auto status=Status::Passed;
    const bool reference_free=std::any_of(evidence->checks.begin(),evidence->checks.end(),[](const auto & check) {
        return check.name=="conditional-ac-optimality";
    });
    const std::vector<const char *> required=reference_free ?
        std::vector<const char *>{"conditional-ac-optimality","independent-replay","width-stationarity",
            "design-identifiability","numerical-identifiability",
            "target-identifiability"} :
        std::vector<const char *>{"inner","kkt","width-stationarity","target-identifiability"};
    for(const char * name:required)
    {
        const auto found=std::find_if(evidence->checks.begin(),evidence->checks.end(),[&](const auto & c){return c.name==name;});
        status=MergeConvergenceStatus(status,found==evidence->checks.end() ? Status::Unavailable : found->status);
    }
    if(evidence->contract=="target-estimability-v1")
        for(const char * name:{"local-correction","target-local-correction"})
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
    const auto & data=JointProblemAccess::Get(problem);
    for(auto & c:fit.components)
    {
        const JointParameterLayout layout=c.layout.value_or(JointParameterLayout{c.atoms,c.rows,{}});
        c.target_evidence=AssessTargetState(problem,layout,c.state,c.evidence,JointEvidenceScope::ComponentLocal,c.rows.size());
    }
    const auto & global_layout=problem.ParameterLayout();
    bool reused=false;
    if(fit.components.size()==1 && fit.components.front().state && fit.assembled_state &&
        fit.components.front().target_evidence && fit.components.front().layout &&
        fit.components.front().rows.size()==problem.Input().observations.size() &&
        FullRowIdentity(fit.components.front().rows,problem.Input().observations.size()) &&
        SameLayout(*fit.components.front().layout,global_layout) &&
        SameState(*fit.components.front().state,*fit.assembled_state) && problem.Input().selection_domain)
    {
        const auto component_context=ProfileContext(data.context,*fit.components.front().layout,
            static_cast<Eigen::Index>(fit.components.front().rows.size()));
        const auto global_context=ProfileContext(data.context,global_layout,
            static_cast<Eigen::Index>(problem.Input().observations.size()));
        if(SameAssessmentPolicy(component_context,global_context))
        {
            fit.target_evidence=*fit.components.front().target_evidence;
            fit.target_evidence->original_rows=problem.Input().observations.size();
            std::vector<JointCheck> global_checks;
            for(const auto & check:fit.evidence)
                if(check.name=="conditional-ac-optimality" || check.name=="independent-replay" ||
                    check.name=="width-stationarity" || check.name=="design-identifiability" ||
                    check.name=="numerical-identifiability" || check.name=="assembly-mapping" ||
                    check.name=="assembled-consistency" || check.name=="inner" || check.name=="kkt" ||
                    check.name=="assembled-profile")
                {
                    global_checks.push_back(check);
                    global_checks.back().scope=JointEvidenceScope::AssembledGlobal;
                }
            const auto target_check=std::find_if(fit.target_evidence->checks.begin(),fit.target_evidence->checks.end(),
                [](const auto & check){return check.name=="target-identifiability";});
            if(target_check!=fit.target_evidence->checks.end())
                global_checks.push_back({target_check->name,target_check->status,JointEvidenceScope::AssembledGlobal,
                    target_check->value,target_check->threshold,target_check->reason});
            fit.target_evidence->checks=std::move(global_checks);
            for(auto & rank:fit.target_evidence->ranks) rank.scope=JointEvidenceScope::AssembledGlobal;
            reused=true;
        }
    }
    if(!reused)
        fit.target_evidence=AssessTargetState(problem,global_layout,fit.assembled_state,fit.evidence,
            JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    fit.costs.assessment_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
}
JointFitResult FitFixedNeighborComponents(const JointProblem & problem,const std::vector<double> & b,
    const FixedNeighborSearchPolicy & policy,
    const JointProgressObserver & observer)
{
    auto out=FitFixedNeighborComponentsImpl(problem,b,policy,observer);
    eigen_helper::ScopedEigenThreadCount eigen_thread_guard{1};
    AddTargetEvidence(out);
    return out;
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
    for(const char * name:{"assembly-mapping","assembled-consistency"})
    {
        const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & c){return c.name==name;});
        status=joint_component::MergeConvergenceStatus(status,found==evidence.end() ? JointCheckStatus::Unavailable : found->status);
    }
    return status;
}
}
