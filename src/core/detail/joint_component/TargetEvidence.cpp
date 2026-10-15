#include "TargetEvidence.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <algorithm>
#include <chrono>

namespace rhbm_gem::core::joint_component {
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
                block.conditioning_classification=="Unavailable" ||
                block.conditioning_classification=="NotAssessed") unavailable=true;
            adequate &= block.status=="assessed" &&
                block.conditioning_classification=="WellConditioned";
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
