#include "TargetEvidence.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <utility>

namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
TargetEvidenceWork & TargetEvidenceWorkForTesting() {static thread_local TargetEvidenceWork work; return work;}
#endif
namespace {
using Status=JointCheckStatus;
bool SameLayout(const JointParameterLayout & a,const JointParameterLayout & b)
{
    return a.full_atoms==b.full_atoms && a.informative_rows==b.informative_rows;
}
bool LayoutMatchesScope(const JointProblem & problem,const JointParameterLayout & layout,
    JointEvidenceScope scope,std::size_t rows)
{
    const auto & data=JointProblemAccess::Get(problem);
    const auto & optimization_layout=layout;
    if(scope==JointEvidenceScope::AssembledGlobal)
        return rows==data.input->observations.size() && SameLayout(optimization_layout,data.layout);
    for(const auto & view:data.partition.components)
    {
        if(rows!=view.rows.size()) continue;
        JointParameterLayout expected;
        for(auto atom:data.layout.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(atom))>=0) expected.full_atoms.push_back(atom);
        for(auto row:data.layout.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(row))>=0) expected.informative_rows.push_back(row);
        if(SameLayout(optimization_layout,expected)) return true;
    }
    return false;
}

template<class Blocks>
struct DiagnosticAtomIndex
{
    std::vector<std::size_t> offsets,block_indices;
    std::vector<bool> invalid_atoms;
    bool invalid{};
    std::string reason;

    std::size_t Count(std::size_t atom) const
    {return atom<invalid_atoms.size() ? offsets[atom+1]-offsets[atom] : 0;}
    bool Invalid(std::size_t atom) const
    {return atom>=invalid_atoms.size() || invalid_atoms[atom];}
    std::pair<std::size_t,std::size_t> Range(std::size_t atom) const
    {return atom<invalid_atoms.size() ? std::pair{offsets[atom],offsets[atom+1]} : std::pair<std::size_t,std::size_t>{0,0};}
};

void MappingError(bool & invalid,std::string & reason,const char * message)
{invalid=true; if(reason.empty()) reason=message;}

template<class Blocks>
DiagnosticAtomIndex<Blocks> BuildDiagnosticAtomIndex(std::size_t atom_count,const Blocks & blocks)
{
    DiagnosticAtomIndex<Blocks> index;
    index.offsets.assign(atom_count+1,0); index.invalid_atoms.assign(atom_count,false);
    std::unordered_map<std::string,std::size_t> identity_counts;
    identity_counts.reserve(blocks.size());
    for(const auto & block:blocks) if(!block.block_id.empty()) ++identity_counts[block.block_id];
    std::vector<bool> invalid_blocks(blocks.size());
    for(std::size_t block_index=0;block_index<blocks.size();++block_index)
    {
        const auto & block=blocks[block_index];
        if(block.block_id.empty())
        {
            invalid_blocks[block_index]=true;
            MappingError(index.invalid,index.reason,"empty-diagnostic-block-identity");
        }
        else if(identity_counts[block.block_id]>1)
        {
            invalid_blocks[block_index]=true;
            MappingError(index.invalid,index.reason,"duplicate-diagnostic-block-identity");
        }
        for(const auto atom:block.core_atoms)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++TargetEvidenceWorkForTesting().memberships_indexed;
#endif
            if(atom>=atom_count)
            {
                MappingError(index.invalid,index.reason,"diagnostic-core-atom-out-of-range");
                continue;
            }
            ++index.offsets[atom+1];
            if(invalid_blocks[block_index]) index.invalid_atoms[atom]=true;
        }
    }
    for(std::size_t atom=1;atom<index.offsets.size();++atom) index.offsets[atom]+=index.offsets[atom-1];
    index.block_indices.resize(index.offsets.back());
    auto next=index.offsets;
    for(std::size_t block_index=0;block_index<blocks.size();++block_index)
        for(const auto atom:blocks[block_index].core_atoms)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++TargetEvidenceWorkForTesting().memberships_indexed;
#endif
            if(atom<atom_count) index.block_indices[next[atom]++]=block_index;
        }
    for(std::size_t atom=0;atom<atom_count;++atom)
        if(index.Count(atom)>1)
        {
            index.invalid_atoms[atom]=true;
            MappingError(index.invalid,index.reason,"duplicate-diagnostic-core-membership");
        }
    return index;
}
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
    auto conditioning_index=BuildDiagnosticAtomIndex(input.atom_ids.size(),conditioning);
    auto stability_index=BuildDiagnosticAtomIndex(input.atom_ids.size(),stability);
    const auto mark_missing_membership=[&](auto & index,const auto & blocks) {
        if(blocks.empty()) return;
        for(const auto atom:layout.full_atoms)
            if(atom>=input.atom_ids.size())
                MappingError(index.invalid,index.reason,"diagnostic-layout-atom-out-of-range");
            else if(index.Count(atom)==0)
            {
                index.invalid_atoms[atom]=true;
                MappingError(index.invalid,index.reason,"missing-diagnostic-core-membership");
            }
    };
    mark_missing_membership(conditioning_index,conditioning);
    mark_missing_membership(stability_index,stability);
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
    const auto local_status=[&](std::size_t atom,std::vector<std::string> & block_ids,bool & mapping_invalid) {
        bool matched=false,unavailable=false,adequate=true;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++TargetEvidenceWorkForTesting().atom_lookups;
#endif
        const auto [first,last]=conditioning_index.Range(atom);
        mapping_invalid=conditioning_index.Invalid(atom) || (!conditioning.empty() && first==last);
        for(std::size_t position=first;position<last;++position)
        {
            const auto & block=conditioning[conditioning_index.block_indices[position]];
            matched=true; block_ids.push_back(block.block_id);
            if(block.status=="unavailable" || block.status=="not-assessed" ||
                block.conditioning_classification=="Unavailable" ||
                block.conditioning_classification=="NotAssessed") unavailable=true;
            adequate &= block.status=="assessed" &&
                block.conditioning_classification=="WellConditioned";
        }
        if(mapping_invalid) return std::string("Unavailable");
        if(!matched) return std::string("NotAssessed");
        if(unavailable) return std::string("Unavailable");
        return adequate ? std::string("Adequate") : std::string("Warning");
    };
    const auto stability_status=[&](std::size_t atom,std::vector<std::string> & block_ids,bool & mapping_invalid) {
        bool matched=false,not_assessed=false;
        int worst{};
        std::string status="NotAssessed";
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ++TargetEvidenceWorkForTesting().atom_lookups;
#endif
        const auto [first,last]=stability_index.Range(atom);
        mapping_invalid=stability_index.Invalid(atom) || (!stability.empty() && first==last);
        for(std::size_t position=first;position<last;++position)
        {
            const auto & block=stability[stability_index.block_indices[position]];
            matched=true; block_ids.push_back(block.block_id);
            const int severity=block.status=="Unstable" ? 4 : block.status=="Sensitive" ? 3 :
                block.status=="Unavailable" ? 2 : block.status=="Stable" ? 1 : 0;
            if(!severity) not_assessed=true;
            if(severity>worst)
            {
                worst=severity; status=block.status;
            }
        }
        if(mapping_invalid) return std::string("Unavailable");
        if(!matched || (not_assessed && worst<=1)) return std::string("NotAssessed");
        return status;
    };
    const auto targets=input.selection_domain ? &input.selection_domain->target_indices : nullptr;
    bool selection_targets_valid=targets!=nullptr;
    JointCheckStatus all_parameter_status=JointCheckStatus::Passed;
    if(!dimensions || layout.full_atoms.empty()) all_parameter_status=JointCheckStatus::Unavailable;
    bool point_mapping_invalid=false;
    for(std::size_t index=0;index<layout.full_atoms.size();++index)
    {
        JointAtomPointValidation point; point.atom=layout.full_atoms[index];
        point.component_convergence=convergence;
        point.parameter_validity=parameter_validity(index);
        bool invalid_conditioning_mapping{},invalid_stability_mapping{};
        point.local_conditioning=local_status(point.atom,point.local_conditioning_blocks,invalid_conditioning_mapping);
        point.parameter_stability=stability_status(point.atom,point.parameter_stability_blocks,invalid_stability_mapping);
        point_mapping_invalid|=invalid_conditioning_mapping || invalid_stability_mapping;
        all_parameter_status=MergeConvergenceStatus(all_parameter_status,point.parameter_validity);
        if(point.component_convergence!=JointCheckStatus::Passed)
        {point.point_quality="Unavailable"; point.reason="component-not-converged";}
        else if(point.parameter_validity!=JointCheckStatus::Passed)
        {point.point_quality="Unavailable"; point.reason="invalid-or-unavailable-parameters";}
        else if(invalid_conditioning_mapping || invalid_stability_mapping)
        {point.point_quality="Warning"; point.reason="invalid-diagnostic-block-mapping";}
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
    {
        std::size_t point_cursor{},layout_cursor{};
        if(!std::is_sorted(targets->begin(),targets->end())) selection_targets_valid=false;
        for(const auto atom:*targets)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            ++TargetEvidenceWorkForTesting().atom_lookups;
#endif
            while(point_cursor<out.point_atoms.size() && out.point_atoms[point_cursor].atom<atom) ++point_cursor;
            while(layout_cursor<layout.full_atoms.size() && layout.full_atoms[layout_cursor]<atom) ++layout_cursor;
            const bool target_in_scope=layout_cursor<layout.full_atoms.size() && layout.full_atoms[layout_cursor]==atom;
            if(point_cursor==out.point_atoms.size() || out.point_atoms[point_cursor].atom!=atom)
            {
                if(scope==JointEvidenceScope::ComponentLocal && !target_in_scope) continue;
                selection_targets_valid=false;
                target_status=MergeConvergenceStatus(target_status,JointCheckStatus::Unavailable);
                continue;
            }
            const auto & found=out.point_atoms[point_cursor];
            target_status=MergeConvergenceStatus(target_status,found.component_convergence);
            target_status=MergeConvergenceStatus(target_status,found.parameter_validity);
        }
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
        all_parameter_status==JointCheckStatus::Passed ? "all-optimized-parameters-valid" : "invalid-or-unavailable-parameters"});
    out.checks.push_back({"target-point-validation",target_status,scope,{},{},
        target_status==JointCheckStatus::Passed ? "selected-target-points-valid" : "selected-target-point-unavailable"});
    const bool diagnostic_mapping_valid=!conditioning_index.invalid && !stability_index.invalid && !point_mapping_invalid;
    const auto mapping_reason=!conditioning_index.reason.empty() ? conditioning_index.reason :
        !stability_index.reason.empty() ? stability_index.reason : "invalid-diagnostic-atom-mapping";
    out.checks.push_back({"diagnostic-block-mapping",diagnostic_mapping_valid ? JointCheckStatus::Passed :
        JointCheckStatus::Unavailable,scope,{},{},diagnostic_mapping_valid ? "unique-atom-to-block-mapping" : mapping_reason});
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
    return status;
}
void AddTargetEvidence(JointFitResult & fit)
{
    if(!fit.problem || !fit.problem->Input().selection_domain) return;
    const auto started=std::chrono::steady_clock::now();
    const auto & problem=*fit.problem;
    for(auto & c:fit.components)
    {
        const JointParameterLayout layout=c.layout.value_or(JointParameterLayout{c.atoms,c.rows});
        static const std::vector<JointBlockParameterStability> no_stability;
        const auto & stability=c.parameter_stability ? c.parameter_stability->blocks : no_stability;
        c.target_evidence=AssessTargetPointState(problem,layout,c.state,c.RuntimeConvergence(),
            c.local_conditioning,stability,JointEvidenceScope::ComponentLocal,c.rows.size());
    }
    auto global_layout=problem.ParameterLayout();
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
    const auto & input=problem.Input();
    if(input.selection_domain) AddTargetEvidence(out);
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
