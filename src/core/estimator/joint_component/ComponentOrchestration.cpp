#include "Problem.hpp"
#include "GaussSeidelBcd.hpp"
#include "ResourceWork.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
std::vector<double> Values(const Vector & v) {return {v.data(),v.data()+v.size()};}
JointState ComponentState(const Endpoint & e,double scale)
{return {Values(e.beta),Values(e.eta.array().exp()),Values(e.eta),Values(e.gradient),e.certificate.objective/(scale*scale)};}
JointParameterLayout ComponentLayout(const JointParameterLayout & layout,const ComponentView & view)
{
    JointParameterLayout out;
    for(auto a:layout.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(a))>=0) out.full_atoms.push_back(a);
    for(auto r:layout.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(r))>=0) out.informative_rows.push_back(r);
    return out;
}
JointCandidateDiagnostic CandidateDiagnostic(const BcdSearchResult & search, std::size_t atoms)
{
    JointCandidateDiagnostic out;
    out.convergence_criteria_id="practical-v2";
    out.global_ac_kkt_threshold=PracticalJointConvergenceCriteria.global_ac_kkt;
    out.global_eta_gradient_threshold=PracticalJointConvergenceCriteria.global_eta_gradient;
    const auto & state=search.state;
    bool finite_widths=state.eta.size()==static_cast<Eigen::Index>(atoms) && state.eta.allFinite();
    if(finite_widths) for(Eigen::Index atom=0;atom<state.eta.size();++atom)
    {
        const double tau=std::exp(state.eta(atom));
        finite_widths &= std::isfinite(tau) && tau>0;
    }
    const bool finite_ac=atoms<=std::numeric_limits<std::size_t>::max()/2 &&
        state.beta.size()==static_cast<Eigen::Index>(2*atoms) && state.beta.allFinite();
    out.candidate_available=finite_widths && finite_ac && std::isfinite(state.objective);
    if(out.candidate_available)
    {
        out.candidate_eta=Values(state.eta);
        out.candidate_ac=Values(state.beta);
        out.candidate_objective=state.objective;
    }
    if(search.final_sweep)
    {
        const auto & sweep=*search.final_sweep;
        if(std::isfinite(sweep.global_ac_kkt)) out.candidate_ac_kkt=sweep.global_ac_kkt;
        if(std::isfinite(sweep.global_eta_gradient_inf_norm))
            out.candidate_eta_gradient_inf_norm=sweep.global_eta_gradient_inf_norm;
        if(std::isfinite(sweep.global_a_feasibility)) out.candidate_feasible=sweep.global_a_feasibility==0;
    }
    if(search.certification)
    {
        const auto & certification=*search.certification;
        out.endpoint_verification_status=certification.verified() ? JointCheckStatus::Passed :
            certification.verification_memory_estimate.status!=MemoryStatus::Pass ||
                certification.assessment.failure.starts_with("resource-limit:") ?
                    JointCheckStatus::Unavailable : JointCheckStatus::Failed;
        out.certification_status=certification.qualified() ? JointCheckStatus::Passed :
            certification.resource_limited() ? JointCheckStatus::Unavailable : JointCheckStatus::Failed;
        out.numerical_qualification_status=out.certification_status;
        out.verification_reason=!certification.trust.reason.empty() ? certification.trust.reason :
            certification.assessment.failure;
        if(certification.verification_memory_estimate.status!=MemoryStatus::Pass)
            out.resource_limit_reason="resource-limit:endpoint-verification:"+certification.verification_memory_estimate.reason;
    }
    else if(out.candidate_available && search.reason.starts_with("resource-limit:"))
    {
        out.endpoint_verification_status=JointCheckStatus::Unavailable;
        out.numerical_qualification_status=JointCheckStatus::Unavailable;
        out.verification_reason=search.reason;
    }
    if(out.resource_limit_reason.empty() && search.memory_estimate.status!=MemoryStatus::Pass)
        out.resource_limit_reason=std::string("resource-limit:")+MemoryStageName(search.memory_estimate.stage)+":"+
            search.memory_estimate.reason;
    if(out.resource_limit_reason.empty() && search.reason.starts_with("resource-limit:"))
        out.resource_limit_reason=search.reason;
    return out;
}
}
BcdComponentBatch OrchestrateBcdComponents(
    const JointProblem & problem,const std::vector<double> & initial_tau,
    const BcdExecutionPolicy & production_policy,const JointProgressObserver & observer)
{
    const auto & data=JointProblemAccess::Get(problem);
    BcdComponentBatch out;
    const auto component_count=data.partition.components.size();
    out.components.reserve(component_count);
    const VectorMap tau_values(initial_tau.data(),static_cast<Eigen::Index>(initial_tau.size()));
    for(std::size_t component_index=0;component_index<component_count;++component_index)
    {
        const auto & view=data.partition.components[component_index];
        const auto component_start=Clock::now();
        std::optional<JointProgressComponent> progress_component;
        if(observer)
        {
            progress_component.emplace(JointProgressComponent{component_index+1,component_count,view.id,
                view.atoms.size(),view.rows.size()});
            NotifyJointProgress(observer,JointProgressPhase::ComponentStarted,*progress_component,
                0,data.context.profile_budget,0,data.context.update_budget);
        }
        const auto layout=ComponentLayout(data.layout,view);
        JointComponentResult component; component.id=view.id;
        component.atoms.assign(view.atoms.begin(),view.atoms.end());
        component.rows.assign(view.rows.begin(),view.rows.end());
        BcdComponentResult fit;
        std::optional<double> accepted_objective,accepted_gradient_inf_norm;
        const bool valid=std::all_of(layout.full_atoms.begin(),layout.full_atoms.end(),[&](auto a) {
            return a<initial_tau.size() && std::isfinite(initial_tau[a]) && initial_tau[a]>0;
        });
        if(valid && !layout.informative_rows.empty())
        {
            fit=SolveBcdComponentView(view,layout,tau_values,data.context,production_policy,observer,
                progress_component ? &*progress_component : nullptr);
            accepted_objective=LastBcdObjective(fit.search);
            accepted_gradient_inf_norm=LastBcdGradient(fit.search);
            component.search_completed=fit.succeeded; component.stop_reason=fit.search.reason;
            component.search_sweeps=fit.search.sweep_count;
            component.first_order_stationarity_sweep=fit.search.first_order_stationarity_sweep;
            component.confirmed_stationarity_sweep=fit.search.confirmed_stationarity_sweep;
            component.endpoint_certification_sweep=fit.search.certification_sweep;
            component.profile_evaluations=static_cast<int>(fit.search.total_profile_evaluations);
            component.accepted_updates=static_cast<int>(fit.search.total_accepted_local_updates);
            component.candidate_diagnostic=CandidateDiagnostic(fit.search, layout.full_atoms.size());
            fit.search.state=joint_component::BlockCoordinateState{};
            fit.search.blocks={}; fit.search.final_sweep.reset();
            component.local_conditioning=std::move(fit.search.local_conditioning);
            const auto remap_core_atoms=[&](auto & blocks) {
                for(auto & block:blocks)
                    for(auto & atom:block.core_atoms)
                    {
                        if(atom>=component.atoms.size())
                        {
                            block.status="unavailable";
                            block.reason="invalid-core-atom-mapping";
                            continue;
                        }
                        atom=component.atoms[atom];
                    }
            };
            remap_core_atoms(component.local_conditioning);
            component.local_conditioning_seconds=fit.search.work.local_conditioning_seconds;
            if(fit.search.work.local_conditioning_peak_memory_estimate_bytes)
                component.local_conditioning_peak_memory_estimate_bytes=
                    fit.search.work.local_conditioning_peak_memory_estimate_bytes;
            if(fit.search.certification)
            {
                const auto & certification=*fit.search.certification;
                component.evidence=AssessmentEvidence(certification.assessment,JointEvidenceScope::ComponentLocal);
                if(fit.search.confirmed_stationarity_sweep)
                    component.evidence.push_back({"complete-sweep-stationarity",JointCheckStatus::Passed,
                        JointEvidenceScope::ComponentLocal,
                        static_cast<double>(fit.search.confirmed_stationarity_sweep),{},"confirmed"});
                if(certification.verification_memory_estimate.status!=MemoryStatus::Pass)
                    component.evidence.push_back({"endpoint-verification-memory",JointCheckStatus::Unavailable,
                        JointEvidenceScope::ComponentLocal,
                        static_cast<double>(certification.verification_memory_estimate.estimated_peak_bytes),
                        static_cast<double>(certification.verification_memory_estimate.budget_bytes),
                        std::string("resource-limit:endpoint-verification:")+certification.verification_memory_estimate.reason});
                if(certification.qualified() && fit.search.accepted_state_integrity)
                    component.state=ComponentState(certification.assessment.primary,data.context.scale);
            }
        }
        else component.stop_reason=valid ? "unobserved-full-parameters" : "invalid-initial-widths";

        if(!component.state && !fit.search.certification)
        {
            Assessment missing; missing.failure=component.stop_reason;
            component.evidence=AssessmentEvidence(missing,JointEvidenceScope::ComponentLocal);
        }
        component.evidence.push_back({"accepted-state-integrity",fit.search.accepted_state_integrity ?
            JointCheckStatus::Passed : JointCheckStatus::Failed,JointEvidenceScope::ComponentLocal,
            static_cast<double>(fit.search.unresolved_failures),0.0,
            "recovery-attempts="+std::to_string(fit.search.recovery_attempts)+
                "; recovery-successes="+std::to_string(fit.search.recovery_successes)+
                "; unresolved="+std::to_string(fit.search.unresolved_failures)});
        const auto & work=fit.search.work.bounded_lm;
        const bool certified_without_steps=work.subproblem_solves==0 && fit.search.certification &&
            fit.search.certification->qualified();
        const auto status=work.subproblem_solves==0 ?
            (certified_without_steps ? JointCheckStatus::Passed : JointCheckStatus::Unavailable) :
            work.failed_subproblems ? JointCheckStatus::Failed : JointCheckStatus::Passed;
        component.evidence.push_back({"active-set-lm-subproblems",status,
            JointEvidenceScope::ComponentLocal,
            std::isfinite(work.maximum_subproblem_kkt) ?
                std::optional<double>(work.maximum_subproblem_kkt) : std::nullopt,{},
            "solves="+std::to_string(work.subproblem_solves)+
                "; failed="+std::to_string(work.failed_subproblems)+
                "; damping-solves="+std::to_string(work.damping_solves)+
                "; bracket-solves="+std::to_string(work.bracket_solves)+
                "; qr-factorizations="+std::to_string(work.damping_qr_factorizations)+
                "; active-set-iterations="+std::to_string(work.active_set_iterations)+
                "; releases="+std::to_string(work.working_set_releases)+
                "; active-bound-occurrences="+std::to_string(work.active_bound_occurrences)+
                (certified_without_steps ? "; no-step-endpoint-certified" : "")});
        if(fit.search.memory_estimate.status!=MemoryStatus::Pass)
            component.evidence.push_back({"memory-preflight",JointCheckStatus::Unavailable,JointEvidenceScope::ComponentLocal,
                static_cast<double>(fit.search.memory_estimate.estimated_peak_bytes),
                static_cast<double>(fit.search.memory_estimate.budget_bytes),
                std::string("resource-limit:")+MemoryStageName(fit.search.memory_estimate.stage)+":"+
                    fit.search.memory_estimate.reason});
        fit.search.certification.reset();
        const double elapsed=Seconds(component_start);
        out.costs.search_seconds+=fit.search.search_seconds;
        out.costs.assessment_seconds+=elapsed-fit.search.search_seconds;
        if(progress_component)
            NotifyJointProgress(observer,JointProgressPhase::ComponentCompleted,*progress_component,
                component.profile_evaluations,data.context.profile_budget,component.accepted_updates,
                data.context.update_budget,elapsed,component.stop_reason,component.state.has_value(),
                accepted_objective,accepted_gradient_inf_norm);
        out.components.push_back(std::move(component));
    }
    return out;
}
BcdComponentResult SolveBcdComponentView(
    const ComponentView & view,const JointParameterLayout & component_layout,const Vector & initial_tau,
    const EvaluationContext & parent_context,const BcdExecutionPolicy & production_policy,
    const JointProgressObserver & observer,const JointProgressComponent * progress_component)
{
    BcdComponentResult out;
    MemoryShape shape; shape.component_atoms=component_layout.full_atoms.size();
    shape.component_observations=component_layout.informative_rows.size();
    shape.component_blocks=shape.component_atoms;
    std::size_t support_entries{};
    for(const auto atom:component_layout.full_atoms)
    {
        const auto local=view.LocalAtom(static_cast<Eigen::Index>(atom));
        if(local<0 || !CheckedAddBytes(support_entries,view.domain.atoms[static_cast<std::size_t>(local)].size(),support_entries))
        {out.search.memory_estimate.stage=MemoryStage::BcdPreparation;
            out.search.memory_estimate.status=MemoryStatus::Overflow; out.search.memory_estimate.reason="size-overflow";
            out.search.reason="resource-limit:bcd-preparation:size-overflow";
            return out;}
    }
    for(const auto row:component_layout.informative_rows)
        if(view.LocalRow(static_cast<Eigen::Index>(row))<0)
        {out.search.reason="block-invalid-partition"; return out;}
    shape.component_support_entries=support_entries; shape.affected_row_entries=support_entries;
    shape.maximum_block_rows=shape.component_observations;
    shape.maximum_block_support_entries=support_entries;
    if(!CheckedMultiplyBytes(shape.component_atoms,2,shape.component_free_columns))
    {out.search.memory_estimate.stage=MemoryStage::BcdPreparation;
        out.search.memory_estimate.status=MemoryStatus::Overflow; out.search.memory_estimate.reason="size-overflow";
        out.search.reason="resource-limit:bcd-preparation:size-overflow";
        return out;}
    const auto * input=view.domain.atoms.Storage();
    const auto snapshot=input ? EstimateSnapshotBytes(*input) : std::optional<std::size_t>{0};
    std::size_t resident=snapshot.value_or(0),root_entries{},root_mappings{},view_entries{};
    const auto support_index_bytes=view.domain.atoms.LocalSupportIndexBytes();
    if(!snapshot || support_index_bytes==std::numeric_limits<std::size_t>::max() ||
        (input && (!CheckedAddBytes(input->atom_ids.size(),input->observations.size(),root_entries) ||
            !CheckedMultiplyBytes(root_entries,4*sizeof(Eigen::Index),root_mappings) ||
            !CheckedAddBytes(resident,root_mappings,resident))) ||
        !CheckedAddBytes(resident,support_index_bytes,resident) ||
        !CheckedAddBytes(view.atoms.capacity(),view.rows.capacity(),view_entries) ||
        !CheckedMultiplyBytes(view_entries,sizeof(Eigen::Index),root_mappings) ||
        !CheckedAddBytes(resident,root_mappings,resident))
    {out.search.memory_estimate.shape=shape; out.search.memory_estimate.stage=MemoryStage::BcdPreparation;
        out.search.memory_estimate.status=MemoryStatus::Overflow; out.search.memory_estimate.reason="size-overflow";
        out.search.reason="resource-limit:bcd-preparation:size-overflow";
        return out;}
    shape.resident_bytes=resident;
    out.search.memory_estimate=PlanMemory(shape,MemoryStage::BcdPreparation);
    if(out.search.memory_estimate.status!=MemoryStatus::Pass)
    {
        out.search.reason="resource-limit:bcd-preparation:"+out.search.memory_estimate.reason;
        return out;
    }
    const auto prepared=MakeComponentProblemView(component_layout,view);
    auto diagnostic_policy=production_policy;
    diagnostic_policy.component_index=static_cast<std::size_t>(view.component_index);
    diagnostic_policy.component_id=view.id;
    out=SolveBcdComponent(prepared,initial_tau,parent_context,
        diagnostic_policy,observer,progress_component);
    return out;
}
}
