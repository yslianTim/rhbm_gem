#include "TargetEvidence.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <sstream>

namespace rhbm_gem::core::joint_component {
namespace {
using Status=JointCheckStatus;
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
bool SameLayout(const JointParameterLayout & a,const JointParameterLayout & b)
{
    return a.full_atoms==b.full_atoms && a.informative_rows==b.informative_rows;
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
MemoryEstimate TargetIdentifiabilityMemoryEstimateImpl(const JointProblemInput & input,const JointParameterLayout & layout,
    std::size_t rows)
{
    MemoryShape shape; shape.component_atoms=layout.full_atoms.size(); shape.component_observations=rows;
    const auto overflow=[&] {
        auto estimate=PlanMemory(shape,MemoryStage::TargetIdentifiability);
        estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow"; return estimate;
    };
    if(!CheckedMultiplyBytes(shape.component_atoms,2,shape.component_design_columns)) return overflow();
    shape.component_width_columns=shape.component_atoms;
    shape.component_free_columns=shape.component_design_columns;
    for(const auto atom:layout.full_atoms)
        if(input.selection_domain && std::binary_search(input.selection_domain->target_indices.begin(),
            input.selection_domain->target_indices.end(),atom))
            if(!CheckedAddBytes(shape.component_target_atoms,1,shape.component_target_atoms)) return overflow();
    shape.component_halo_atoms=shape.component_atoms-shape.component_target_atoms;
    for(const auto atom:layout.full_atoms)
        if(atom>=input.support.size() ||
            !CheckedAddBytes(shape.component_support_entries,input.support[atom].size(),shape.component_support_entries))
            return overflow();
    const auto snapshot=EstimateSnapshotBytes(input);
    std::size_t resident{},rank_columns{};
    const auto add_product=[](std::size_t & total,std::size_t count,std::size_t width) {
        std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
    };
    if(!snapshot || !CheckedAddBytes(shape.component_free_columns,shape.component_width_columns,rank_columns))
        return overflow();
    resident=*snapshot;
    if(!add_product(resident,rows,sizeof(double)) ||
        !add_product(resident,shape.component_atoms,4*sizeof(double)) ||
        !add_product(resident,rows,2*sizeof(Eigen::Index)) ||
        !add_product(resident,shape.component_atoms,2*sizeof(Eigen::Index)+sizeof(std::vector<LocalSupportIndex>)) ||
        !add_product(resident,shape.component_support_entries,sizeof(LocalSupportIndex)) ||
        !add_product(resident,rank_columns,sizeof(double)) ||
        !add_product(resident,2,sizeof(JointRankEvidence)+64)) return overflow();
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
    JointTargetEvidence out; out.contract="target-estimability-v4"; out.original_rows=original_rows;
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
JointTargetEvidence AssessTargetState(const JointProblem & problem,const JointParameterLayout & layout,
    const std::optional<JointState> & state,const std::vector<JointCheck> & existing,
    const std::vector<JointRankEvidence> & ranks,JointEvidenceScope scope,std::size_t original_rows)
{
    const auto & input=problem.Input();
    const auto target_memory=state ? TargetIdentifiabilityMemoryEstimate(input,layout,layout.informative_rows.size()) :
        MemoryEstimate{};
    if(state)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting preflight_stage("derived-target-memory-preflight",
            static_cast<Eigen::Index>(layout.informative_rows.size()),
            static_cast<Eigen::Index>(layout.full_atoms.size()),&target_memory);
        preflight_stage.Finish();
#endif
        if(target_memory.status!=MemoryStatus::Pass || target_memory.estimated_peak_bytes>target_memory.budget_bytes/2)
            LogMemoryEstimate(target_memory);
        if(target_memory.status!=MemoryStatus::Pass)
        {
            JointTargetEvidence out; out.contract="target-estimability-v4"; out.original_rows=original_rows;
            for(const auto & check:existing)
                if(check.name=="conditional-ac-optimality" || check.name=="independent-replay" ||
                    check.name=="width-stationarity" || check.name=="design-identifiability" ||
                    check.name=="numerical-identifiability" || check.name=="assembly-mapping" ||
                    check.name=="assembled-consistency" || check.name=="inner" || check.name=="kkt" ||
                    check.name=="assembled-profile") out.checks.push_back(check);
            const auto reason=std::string("resource-limit:target-identifiability");
            for(const auto atom:layout.full_atoms) out.atoms.push_back({atom,Status::Unavailable,{},reason});
            std::sort(out.atoms.begin(),out.atoms.end(),[](const auto & a,const auto & b){return a.atom<b.atom;});
            out.checks.push_back({"target-identifiability",Status::Unavailable,scope,{}, {},reason});
            out.checks.push_back({"derived-target-memory-preflight",Status::Unavailable,scope,
                static_cast<double>(target_memory.estimated_peak_bytes),static_cast<double>(target_memory.budget_bytes),
                reason+":"+target_memory.reason});
            return out;
        }
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting derived_stage("derived-target-identifiability",
        static_cast<Eigen::Index>(layout.informative_rows.size()),static_cast<Eigen::Index>(layout.full_atoms.size()),
        state ? &target_memory : nullptr);
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    auto out=DeriveTargetIdentifiability(problem,layout,state,existing,ranks,scope,original_rows);
    derived_stage.Finish();
    return out;
#else
    return DeriveTargetIdentifiability(problem,layout,state,existing,ranks,scope,original_rows);
#endif
}
JointTargetEvidence AssessTargetPointState(const JointProblem & problem,const JointParameterLayout & layout,
    const std::optional<JointState> & state,JointCheckStatus convergence,
    const std::vector<JointBlockConditioning> & conditioning,
    const std::vector<JointBlockParameterStability> & stability,JointEvidenceScope scope,std::size_t original_rows)
{
    JointTargetEvidence out;
    out.contract=std::string(JointTargetPointValidationContractVersion);
    out.original_rows=original_rows;
    const auto & input=problem.Input();
    const bool layout_valid=LayoutMatchesScope(problem,layout,scope,original_rows);
    const bool dimensions=state && layout_valid && layout.full_atoms.size()<=std::numeric_limits<std::size_t>::max()/2 &&
        state->ac.size()==2*layout.full_atoms.size() && state->b.size()==layout.full_atoms.size() &&
        state->log_b.size()==layout.full_atoms.size();
    const auto parameter_validity=[&](std::size_t index) {
        if(!dimensions) return JointCheckStatus::Unavailable;
        const double amplitude=state->ac[2*index],width=state->b[index],log_width=state->log_b[index];
        if(!std::isfinite(amplitude) || amplitude<0 || !std::isfinite(state->ac[2*index+1]) ||
            !std::isfinite(width) || !(width>0) || !std::isfinite(log_width))
            return JointCheckStatus::Failed;
        const double exp_width=std::exp(log_width);
        return std::isfinite(exp_width) && exp_width>0 &&
            std::abs(exp_width-width)<=8*std::numeric_limits<double>::epsilon()*width ?
                JointCheckStatus::Passed : JointCheckStatus::Failed;
    };
    const auto local_status=[&](std::size_t atom,std::vector<std::string> & block_ids) {
        bool matched=false,unavailable=false,adequate=true;
        for(const auto & block:conditioning)
        {
            if(std::find(block.core_atoms.begin(),block.core_atoms.end(),atom)==block.core_atoms.end()) continue;
            matched=true; block_ids.push_back(block.block_id);
            if(block.status=="unavailable" || block.status=="not-assessed" ||
                block.free_ac_status=="unavailable" || block.projected_width_status=="unavailable" ||
                block.profile_status=="unavailable") unavailable=true;
            adequate &= block.status=="assessed" && block.active_amplitude_count==0 &&
                block.free_ac_status=="full-rank" && block.projected_width_status=="full-rank" &&
                block.profile_status=="full-rank";
        }
        if(!matched) return std::string("NotAssessed");
        if(unavailable) return std::string("Unavailable");
        return adequate ? std::string("Adequate") : std::string("Warning");
    };
    const auto stability_status=[&](std::size_t atom,std::vector<std::string> & block_ids) {
        bool matched=false,not_assessed=false;
        int worst{};
        std::string status="NotAssessed";
        for(const auto & block:stability)
        {
            if(std::find(block.core_atoms.begin(),block.core_atoms.end(),atom)==block.core_atoms.end()) continue;
            matched=true; block_ids.push_back(block.block_id);
            const int severity=block.status=="Unstable" ? 4 : block.status=="Sensitive" ? 3 :
                block.status=="Unavailable" ? 2 : block.status=="Stable" ? 1 : 0;
            if(!severity) not_assessed=true;
            if(severity>worst)
            {
                worst=severity; status=block.status;
            }
        }
        if(!matched || (not_assessed && worst<=1)) return std::string("NotAssessed");
        return status;
    };
    const auto targets=input.selection_domain ? &input.selection_domain->target_indices : nullptr;
    bool selection_targets_valid=targets!=nullptr;
    JointCheckStatus all_parameter_status=JointCheckStatus::Passed;
    if(!dimensions || layout.full_atoms.empty()) all_parameter_status=JointCheckStatus::Unavailable;
    for(std::size_t index=0;index<layout.full_atoms.size();++index)
    {
        JointAtomPointValidation point; point.atom=layout.full_atoms[index];
        point.component_convergence=convergence;
        point.parameter_validity=parameter_validity(index);
        point.local_conditioning=local_status(point.atom,point.local_conditioning_blocks);
        point.parameter_stability=stability_status(point.atom,point.parameter_stability_blocks);
        all_parameter_status=MergeConvergenceStatus(all_parameter_status,point.parameter_validity);
        if(point.component_convergence!=JointCheckStatus::Passed)
        {point.point_quality="Unavailable"; point.reason="component-not-converged";}
        else if(point.parameter_validity!=JointCheckStatus::Passed)
        {point.point_quality="Unavailable"; point.reason="invalid-or-unavailable-parameters";}
        else if(point.local_conditioning=="Warning" || point.local_conditioning=="Unavailable" ||
            point.parameter_stability=="Unstable" || point.parameter_stability=="Sensitive" ||
            point.parameter_stability=="Unavailable")
        {point.point_quality="Warning"; point.reason="conditioning-or-stability-warning";}
        else if(point.local_conditioning=="NotAssessed" || point.parameter_stability=="NotAssessed")
        {point.point_quality="NotAssessed"; point.reason="quality-diagnostics-not-assessed";}
        else
        {point.point_quality="Available"; point.reason="numerically-converged-point";}
        out.point_atoms.push_back(std::move(point));
    }
    JointCheckStatus target_status=targets ? convergence : JointCheckStatus::Unavailable;
    if(targets)
        for(const auto atom:*targets)
        {
            const auto found=std::lower_bound(out.point_atoms.begin(),out.point_atoms.end(),atom,
                [](const auto & point,std::size_t index){return point.atom<index;});
            if(found==out.point_atoms.end() || found->atom!=atom)
            {
                if(scope==JointEvidenceScope::ComponentLocal &&
                    !std::binary_search(layout.full_atoms.begin(),layout.full_atoms.end(),atom)) continue;
                selection_targets_valid=false;
                target_status=MergeConvergenceStatus(target_status,JointCheckStatus::Unavailable);
                continue;
            }
            target_status=MergeConvergenceStatus(target_status,found->component_convergence);
            target_status=MergeConvergenceStatus(target_status,found->parameter_validity);
        }
    if(!layout_valid || (scope==JointEvidenceScope::AssembledGlobal && original_rows!=input.observations.size()))
    {
        all_parameter_status=JointCheckStatus::Unavailable;
        target_status=JointCheckStatus::Unavailable;
    }
    if(!selection_targets_valid) target_status=JointCheckStatus::Unavailable;
    out.checks.push_back({"component-convergence",convergence,scope,{},{},
        convergence==JointCheckStatus::Passed ? "global-numerical-convergence" : "component-not-converged"});
    out.checks.push_back({"parameter-validity",all_parameter_status,scope,{},{},
        all_parameter_status==JointCheckStatus::Passed ? "all-full-abc-parameters-valid" : "invalid-or-unavailable-parameters"});
    out.checks.push_back({"target-point-validation",target_status,scope,{},{},
        target_status==JointCheckStatus::Passed ? "selected-target-points-valid" : "selected-target-point-unavailable"});
    return out;
}
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> & evidence,bool state)
{
    if(!evidence) return Status::NotRun;
    if(!state) return Status::Unavailable;
    if(evidence->contract==JointTargetPointValidationContractVersion)
    {
        auto status=Status::Passed;
        for(const char * name:{"component-convergence","parameter-validity","target-point-validation"})
        {
            const auto found=std::find_if(evidence->checks.begin(),evidence->checks.end(),[&](const auto & check) {
                return check.name==name;
            });
            status=MergeConvergenceStatus(status,found==evidence->checks.end() ? Status::Unavailable : found->status);
        }
        return status;
    }
    if(evidence->contract!="target-estimability-v1" && evidence->contract!="target-estimability-v2" &&
        evidence->contract!="target-estimability-v3" && evidence->contract!="target-estimability-v4")
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
    for(auto & c:fit.components)
    {
        const JointParameterLayout layout=c.layout.value_or(JointParameterLayout{c.atoms,c.rows,{}});
        static const std::vector<JointBlockParameterStability> no_stability;
        const auto & stability=c.parameter_stability ? c.parameter_stability->blocks : no_stability;
        c.target_evidence=AssessTargetPointState(problem,layout,c.state,c.RuntimeConvergence(),
            c.local_conditioning,stability,JointEvidenceScope::ComponentLocal,c.rows.size());
    }
    const auto & global_layout=problem.ParameterLayout();
    std::vector<JointBlockConditioning> conditioning;
    std::vector<JointBlockParameterStability> stability;
    for(const auto & component:fit.components)
    {
        conditioning.insert(conditioning.end(),component.local_conditioning.begin(),component.local_conditioning.end());
        if(component.parameter_stability)
            stability.insert(stability.end(),component.parameter_stability->blocks.begin(),component.parameter_stability->blocks.end());
    }
    fit.target_evidence=AssessTargetPointState(problem,global_layout,fit.assembled_state,fit.RuntimeConvergence(),
        conditioning,stability,JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
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
