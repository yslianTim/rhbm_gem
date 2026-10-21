#include "Numerics.hpp"
#include "FixedNeighborBlockCoordinate.hpp"
#include "Problem.hpp"
#include "ResourceWork.hpp"
#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <cmath>
#include <chrono>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
std::size_t Index(Eigen::Index i) {return static_cast<std::size_t>(i);}
}
ComponentPartition Partition(const Domain & domain,const Identities & ids)
{
    const std::vector<std::string> identities=ids;
    const auto count=static_cast<Eigen::Index>(domain.atoms.size());
    if(domain.rows<0 || ids.size()!=domain.atoms.size() || std::set<std::string>(identities.begin(),identities.end()).size()!=ids.size())
        throw std::invalid_argument("Invalid structural partition identities.");
    std::vector<Eigen::Index> parent(static_cast<std::size_t>(count)); std::iota(parent.begin(),parent.end(),0);
    auto root=[&](Eigen::Index a) {while(parent[Index(a)]!=a) {parent[Index(a)]=parent[Index(parent[Index(a)])]; a=parent[Index(a)];} return a;};
    std::vector<std::vector<Eigen::Index>> contributors(static_cast<std::size_t>(domain.rows));
    ComponentPartition out; auto mappings=std::make_shared<PartitionMappings>(); out.mappings=mappings;
    mappings->atom_component.assign(ids.size(),-1); mappings->row_component.assign(static_cast<std::size_t>(domain.rows),-1);
    mappings->atom_to_local.assign(ids.size(),-1); mappings->row_to_local.assign(static_cast<std::size_t>(domain.rows),-1);
    for(Eigen::Index a=0;a<count;++a)
    {
        std::set<Eigen::Index> seen;
        for(const auto & s:domain.atoms[Index(a)])
        {
            if(s.row<0 || s.row>=domain.rows || !std::isfinite(s.square) || s.square<0 || s.square>6.25 || !seen.insert(s.row).second)
                throw std::invalid_argument("Invalid structural membership.");
            contributors[Index(s.row)].push_back(a);
        }
        if(domain.atoms[Index(a)].empty()) out.unobserved_atoms.push_back(a);
    }
    for(const auto & row:contributors) for(std::size_t k=1;k<row.size();++k) parent[Index(root(row[k]))]=root(row[0]);
    std::map<Eigen::Index,std::vector<Eigen::Index>> groups;
    for(Eigen::Index a=0;a<count;++a) groups[root(a)].push_back(a);
    for(const auto & [key,atoms]:groups)
    {
        (void)key; ComponentView v; v.atoms=atoms;
        v.id=ids[Index(*std::min_element(atoms.begin(),atoms.end(),[&](auto a,auto b){return ids[Index(a)]<ids[Index(b)];}))];
        out.components.push_back(std::move(v));
    }
    std::sort(out.components.begin(),out.components.end(),[](const auto & a,const auto & b){return a.id<b.id;});
    for(std::size_t c=0;c<out.components.size();++c) for(auto a:out.components[c].atoms) mappings->atom_component[Index(a)]=static_cast<Eigen::Index>(c);
    for(Eigen::Index r=0;r<domain.rows;++r)
    {
        if(contributors[Index(r)].empty()) {out.constant_rows.push_back(r); continue;}
        const auto c=mappings->atom_component[Index(contributors[Index(r)][0])]; mappings->row_component[Index(r)]=c; out.components[Index(c)].rows.push_back(r);
    }
    for(std::size_t c=0;c<out.components.size();++c)
    {
        auto & v=out.components[c]; v.mappings=mappings; v.component_index=static_cast<Eigen::Index>(c);
        for(std::size_t r=0;r<v.rows.size();++r) mappings->row_to_local[Index(v.rows[r])]=static_cast<Eigen::Index>(r);
        for(std::size_t a=0;a<v.atoms.size();++a) mappings->atom_to_local[Index(v.atoms[a])]=static_cast<Eigen::Index>(a);
        const auto rows=CompactRowMapping::FromParentRows(v.rows);
        v.domain=domain.Select(v.atoms,static_cast<Eigen::Index>(v.rows.size()),rows);
    }
    return out;
}
ComponentProblemView MakeComponentProblemView(const JointParameterLayout & parent_layout,const ComponentView & view)
{
    const auto parent=view.domain.atoms.Snapshot();
    ComponentProblemView out; out.parent_row_count=static_cast<Eigen::Index>(view.rows.size());
    Indices view_to_local_rows;
    Indices local_to_view_atoms;
    Indices identity_parent_atoms;
    out.local_to_parent_atoms.reserve(parent_layout.full_atoms.size());
    out.local_to_parent_rows.reserve(parent_layout.informative_rows.size());
    local_to_view_atoms.reserve(parent_layout.full_atoms.size());
    identity_parent_atoms.reserve(parent_layout.full_atoms.size());
    for(std::size_t k=0;k<parent_layout.informative_rows.size();++k)
    {
        const auto row=parent_layout.informative_rows[k];
        const auto component_row=view.LocalRow(static_cast<Eigen::Index>(row));
        if(component_row<0) throw std::invalid_argument("Component row is outside its view.");
        view_to_local_rows.push_back(component_row);
        out.local_to_parent_rows.push_back(static_cast<Eigen::Index>(row));
        out.local_layout.informative_rows.push_back(k);
    }
    for(std::size_t k=0;k<parent_layout.full_atoms.size();++k)
    {
        const auto atom=parent_layout.full_atoms[k];
        const auto component_atom=view.LocalAtom(static_cast<Eigen::Index>(atom));
        if(component_atom<0) throw std::invalid_argument("Component atom is outside its view.");
        out.local_to_parent_atoms.push_back(static_cast<Eigen::Index>(atom));
        local_to_view_atoms.push_back(component_atom);
        identity_parent_atoms.push_back(static_cast<Eigen::Index>(view.domain.atoms.ParentAtom(static_cast<std::size_t>(component_atom))));
        out.local_layout.full_atoms.push_back(k);
    }
    const auto row_mapping=CompactRowMapping::FromParentRows(view_to_local_rows);
    out.domain=view.domain.Select(local_to_view_atoms,static_cast<Eigen::Index>(out.local_layout.informative_rows.size()),row_mapping);
    const Identities parent_atoms(std::shared_ptr<const std::vector<std::string>>(parent,&parent->atom_ids));
    const Identities parent_rows(std::shared_ptr<const std::vector<std::string>>(parent,&parent->row_ids));
    out.atom_ids=parent_atoms.Select(identity_parent_atoms);
    out.row_ids=parent_rows.Select(out.local_to_parent_rows);
    return out;
}
Eigen::VectorXd SelectValues(VectorRef v,const std::vector<Eigen::Index> & indices)
{Eigen::VectorXd out(static_cast<Eigen::Index>(indices.size())); for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=v(indices[k]); return out;}
EvaluationContext ChildContext(const EvaluationContext & parent,const ComponentView & view,bool independent)
{
    auto c=parent; c.atom_ids=parent.atom_ids.Select(view.atoms); c.row_ids=parent.row_ids.Select(view.rows);
    if(parent.audit.directions.size())
    {
        c.audit.directions.resize(static_cast<Eigen::Index>(view.atoms.size()),parent.audit.directions.cols());
        for(std::size_t a=0;a<view.atoms.size();++a) c.audit.directions.row(static_cast<Eigen::Index>(a))=parent.audit.directions.row(view.atoms[a]);
    }
    c.audit.boundary_atoms.clear();
    for(auto a:parent.audit.boundary_atoms) if(view.LocalAtom(a)>=0) c.audit.boundary_atoms.push_back(view.LocalAtom(a));
    c.audit.boundary=!c.audit.boundary_atoms.empty();
    c.independent_search=independent;
    if(independent)
    {
        const auto n=static_cast<Eigen::Index>(view.atoms.size()); c.rank={view.domain.rows,2*n,n};
        c.linear.rank_relative=c.rank.Relative(2*n);
    }
    return c;
}
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
JointCandidateDiagnostic CandidateDiagnostic(const FixedNeighborSearchResult & search,std::size_t atoms,
    const JointConvergenceCriteria & criteria)
{
    JointCandidateDiagnostic out;
    out.convergence_criteria_id=criteria.global_ac_kkt==PracticalJointConvergenceCriteria.global_ac_kkt &&
        criteria.global_width_gradient==PracticalJointConvergenceCriteria.global_width_gradient &&
        criteria.confirmed_log_width_change==PracticalJointConvergenceCriteria.confirmed_log_width_change ?
            "practical-v1" : criteria.global_ac_kkt==StrictJointConvergenceCriteria.global_ac_kkt &&
                criteria.global_width_gradient==StrictJointConvergenceCriteria.global_width_gradient &&
                criteria.confirmed_log_width_change==StrictJointConvergenceCriteria.confirmed_log_width_change ?
                    "strict-v1" : "custom-v1";
    out.global_ac_kkt_threshold=criteria.global_ac_kkt;
    out.global_width_gradient_threshold=criteria.global_width_gradient;
    out.confirmed_log_width_change_threshold=criteria.confirmed_log_width_change;
    const auto & state=search.state;
    bool finite_widths=state.eta.size()==static_cast<Eigen::Index>(atoms) && state.eta.allFinite();
    if(finite_widths) for(Eigen::Index atom=0;atom<state.eta.size();++atom)
    {
        const double width=std::exp(state.eta(atom));
        finite_widths &= std::isfinite(width) && width>0;
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
        if(std::isfinite(sweep.global_width_gradient_inf_norm))
            out.candidate_width_gradient_inf_norm=sweep.global_width_gradient_inf_norm;
        if(std::isfinite(sweep.eta_change_inf)) out.candidate_log_width_change_inf_norm=sweep.eta_change_inf;
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
FixedNeighborComponentBatch OrchestrateFixedNeighborComponents(
    const JointProblem & problem,const std::vector<double> & initial_b,
    const FixedNeighborSearchPolicy & production_policy,const JointProgressObserver & observer)
{
    const auto & data=JointProblemAccess::Get(problem);
    FixedNeighborComponentBatch out;
    const auto component_count=data.partition.components.size();
    out.components.reserve(component_count);
    const VectorMap widths(initial_b.data(),static_cast<Eigen::Index>(initial_b.size()));
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
        FixedNeighborComponentResult fit;
        std::optional<double> accepted_objective,accepted_gradient_inf_norm;
        const bool valid=std::all_of(layout.full_atoms.begin(),layout.full_atoms.end(),[&](auto a) {
            return a<initial_b.size() && std::isfinite(initial_b[a]) && initial_b[a]>0;
        });
        if(valid && !layout.informative_rows.empty())
        {
            fit=SolveFixedNeighborComponentView(view,layout,widths,data.context,production_policy,observer,
                progress_component ? &*progress_component : nullptr);
            if(production_policy.search_result_observer)
                production_policy.search_result_observer(view.id,fit.search);
            accepted_objective=LastFixedNeighborObjective(fit.search);
            accepted_gradient_inf_norm=LastFixedNeighborGradient(fit.search);
            component.search_completed=fit.succeeded; component.stop_reason=fit.stop_reason;
            component.search_sweeps=fit.search.sweep_count;
            component.first_order_stationarity_sweep=fit.search.first_order_stationarity_sweep;
            component.confirmed_stationarity_sweep=fit.search.confirmed_stationarity_sweep;
            component.endpoint_certification_sweep=fit.search.certification_sweep;
            component.profile_evaluations=static_cast<int>(fit.search.total_profile_evaluations);
            component.accepted_updates=static_cast<int>(fit.search.total_accepted_local_updates);
            component.candidate_diagnostic=CandidateDiagnostic(fit.search,layout.full_atoms.size(),
                production_policy.convergence_criteria);
            component.local_conditioning=fit.search.local_conditioning;
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
            component.parameter_stability=fit.search.parameter_stability;
            if(component.parameter_stability) remap_core_atoms(component.parameter_stability->blocks);
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
                if(certification.qualified())
                    component.state=ComponentState(certification.assessment.primary,data.context.scale);
            }
        }
        else component.stop_reason=valid ? "unobserved-full-parameters" : "invalid-initial-widths";

        if(!component.state && !fit.search.certification)
        {
            Assessment missing; missing.failure=component.stop_reason;
            component.evidence=AssessmentEvidence(missing,JointEvidenceScope::ComponentLocal);
        }
        if(production_policy.eta_upper_bound)
        {
            const auto & work=fit.search.work.local_profile_work.bounded_lm;
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
                    "; active-set-iterations="+std::to_string(work.active_set_iterations)+
                    "; releases="+std::to_string(work.releases)+
                    "; active-bound-occurrences="+std::to_string(work.active_bound_occurrences)+
                    (certified_without_steps ? "; no-step-endpoint-certified" : "")});
        }
        if(fit.search.memory_estimate.status!=MemoryStatus::Pass)
            component.evidence.push_back({"memory-preflight",JointCheckStatus::Unavailable,JointEvidenceScope::ComponentLocal,
                static_cast<double>(fit.search.memory_estimate.estimated_peak_bytes),
                static_cast<double>(fit.search.memory_estimate.budget_bytes),
                std::string("resource-limit:")+MemoryStageName(fit.search.memory_estimate.stage)+":"+
                    fit.search.memory_estimate.reason});
        const double elapsed=Seconds(component_start);
        out.costs.search_seconds+=fit.search.search_seconds;
        out.costs.assessment_seconds+=elapsed-fit.search.search_seconds;
        if(progress_component)
            NotifyJointProgress(observer,JointProgressPhase::ComponentCompleted,*progress_component,
                component.profile_evaluations,data.context.profile_budget,component.accepted_updates,
                data.context.update_budget,elapsed,component.stop_reason,component.state.has_value(),
                accepted_objective,accepted_gradient_inf_norm);
        if(component_count==1 && data.partition.constant_rows.empty() &&
            fit.search.certification && fit.search.certification->qualified())
            out.single_component_assessment=std::move(fit.search.certification->assessment);
        out.components.push_back(std::move(component));
    }
    return out;
}
FixedNeighborComponentResult SolveFixedNeighborComponentView(
    const ComponentView & view,const JointParameterLayout & component_layout,const Vector & initial_b,
    const EvaluationContext & parent_context,const FixedNeighborSearchPolicy & production_policy,
    const JointProgressObserver & observer,const JointProgressComponent * progress_component)
{
    FixedNeighborComponentResult out;
    if(production_policy.core_atoms==0)
    {out.stop_reason="block-invalid-partition"; out.search.reason=out.stop_reason; return out;}
    MemoryShape shape; shape.component_atoms=component_layout.full_atoms.size();
    shape.component_observations=component_layout.informative_rows.size();
    shape.component_blocks=shape.component_atoms;
    std::size_t support_entries{};
    for(const auto atom:component_layout.full_atoms)
    {
        const auto local=view.LocalAtom(static_cast<Eigen::Index>(atom));
        if(local<0 || !CheckedAddBytes(support_entries,view.domain.atoms[static_cast<std::size_t>(local)].size(),support_entries))
        {out.search.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
            out.search.memory_estimate.status=MemoryStatus::Overflow; out.search.memory_estimate.reason="size-overflow";
            out.search.reason="resource-limit:fixed-neighbor-preparation:size-overflow";
            out.stop_reason=out.search.reason; return out;}
    }
    for(const auto row:component_layout.informative_rows)
        if(view.LocalRow(static_cast<Eigen::Index>(row))<0)
        {out.stop_reason="block-invalid-partition"; out.search.reason=out.stop_reason; return out;}
    shape.component_support_entries=support_entries; shape.affected_row_entries=support_entries;
    shape.maximum_block_rows=shape.component_observations;
    shape.maximum_block_support_entries=support_entries;
    if(!CheckedMultiplyBytes(shape.component_atoms,2,shape.component_free_columns))
    {out.search.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
        out.search.memory_estimate.status=MemoryStatus::Overflow; out.search.memory_estimate.reason="size-overflow";
        out.search.reason="resource-limit:fixed-neighbor-preparation:size-overflow";
        out.stop_reason=out.search.reason; return out;}
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
    {out.search.memory_estimate.shape=shape; out.search.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
        out.search.memory_estimate.status=MemoryStatus::Overflow; out.search.memory_estimate.reason="size-overflow";
        out.search.reason="resource-limit:fixed-neighbor-preparation:size-overflow";
        out.stop_reason=out.search.reason; return out;}
    shape.resident_bytes=resident;
    out.search.memory_estimate=PlanMemory(shape,MemoryStage::FixedNeighborPreparation);
    if(out.search.memory_estimate.status!=MemoryStatus::Pass)
    {
        out.search.reason="resource-limit:fixed-neighbor-preparation:"+out.search.memory_estimate.reason;
        out.stop_reason=out.search.reason; return out;
    }
    const auto prepared=MakeComponentProblemView(component_layout,view);
    auto diagnostic_policy=production_policy;
    diagnostic_policy.component_index=static_cast<std::size_t>(view.component_index);
    diagnostic_policy.component_id=view.id;
    const auto input_owner=view.domain.atoms.Snapshot();
    diagnostic_policy.halo_census=[input_owner,&component_layout,&view](
        std::size_t & halo_count,std::size_t & minimum_support) {
        const auto census=DiagnoseHaloSupport(*input_owner,component_layout,&view.atoms);
        halo_count=census.size(); minimum_support=census.empty() ? 0 : std::numeric_limits<std::size_t>::max();
        for(const auto & halo:census)
            minimum_support=std::min(minimum_support,halo.original_observation_support_count);
    };
    out=SolveFixedNeighborComponent(prepared,initial_b,parent_context,
        diagnostic_policy,observer,progress_component);
    return out;
}
AssemblyResult AssembleComponents(const Domain & domain,VectorRef y,const ComponentPartition & partition,
    const EvaluationContext & context,const std::vector<JointComponentResult> & components,const AssessmentReuse * reuse,
    const JointConvergenceCriteria & convergence_criteria)
{
    ResourcePhase phase("assembly");
    AssemblyResult out; out.assessment.convergence_criteria=convergence_criteria;
    if(components.size()!=partition.components.size())
    {out.mapping_reason="missing-component-result"; return out;}
    out.available=true;
    out.row_mask.assign(static_cast<std::size_t>(domain.rows),false);
    out.eta=Vector::Constant(static_cast<Eigen::Index>(domain.atoms.size()),unavailable);
    out.beta=Vector::Constant(2*out.eta.size(),unavailable);
    if(domain.rows<0 || y.size()!=domain.rows || !partition.mappings ||
        partition.mappings->atom_component.size()!=domain.atoms.size() ||
        partition.mappings->atom_to_local.size()!=domain.atoms.size() ||
        partition.mappings->row_component.size()!=static_cast<std::size_t>(domain.rows) ||
        partition.mappings->row_to_local.size()!=static_cast<std::size_t>(domain.rows))
    {out.available=false; out.mapping_reason="invalid-partition-mapping"; return out;}
    std::vector<bool> mapped_atoms(domain.atoms.size(),false),mapped_rows(static_cast<std::size_t>(domain.rows),false);
    for(auto r:partition.constant_rows)
    {
        if(r<0 || r>=domain.rows || mapped_rows[static_cast<std::size_t>(r)] ||
            partition.mappings->row_component[static_cast<std::size_t>(r)]!=-1)
        {out.available=false; out.mapping_reason="invalid-constant-row-mapping"; return out;}
        mapped_rows[static_cast<std::size_t>(r)]=true;
        out.row_mask[static_cast<std::size_t>(r)]=true;
    }
    for(std::size_t k=0;k<components.size();++k)
    {
        const auto & component=components[k]; const auto & view=partition.components[k];
        out.completed &= component.search_completed;
        const bool same_atoms=component.atoms.size()==view.atoms.size() &&
            std::equal(component.atoms.begin(),component.atoms.end(),view.atoms.begin(),[](std::size_t a,Eigen::Index b) {
                return b>=0 && a==static_cast<std::size_t>(b);
            });
        const bool same_rows=component.rows.size()==view.rows.size() &&
            std::equal(component.rows.begin(),component.rows.end(),view.rows.begin(),[](std::size_t a,Eigen::Index b) {
                return b>=0 && a==static_cast<std::size_t>(b);
            });
        if(component.id!=view.id || view.component_index!=static_cast<Eigen::Index>(k) || !same_atoms || !same_rows)
        {
            out.available=false; out.mapping_reason="component-identity-mapping-mismatch"; return out;
        }
        for(std::size_t atom_index=0;atom_index<view.atoms.size();++atom_index)
        {
            const auto atom=view.atoms[atom_index];
            if(atom<0 || atom>=static_cast<Eigen::Index>(domain.atoms.size()) ||
                partition.mappings->atom_component[static_cast<std::size_t>(atom)]!=static_cast<Eigen::Index>(k) ||
                partition.mappings->atom_to_local[static_cast<std::size_t>(atom)]!=static_cast<Eigen::Index>(atom_index) ||
                mapped_atoms[static_cast<std::size_t>(atom)])
            {out.available=false; out.mapping_reason="component-atom-mapping-mismatch"; return out;}
            mapped_atoms[static_cast<std::size_t>(atom)]=true;
            for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
                if(support.row<0 || support.row>=domain.rows ||
                    partition.mappings->row_component[static_cast<std::size_t>(support.row)]!=static_cast<Eigen::Index>(k))
                {out.available=false; out.mapping_reason="component-support-row-mapping-mismatch"; return out;}
        }
        for(std::size_t row_index=0;row_index<view.rows.size();++row_index)
        {
            const auto r=view.rows[row_index];
            if(r<0 || r>=domain.rows || mapped_rows[static_cast<std::size_t>(r)] ||
                partition.mappings->row_component[static_cast<std::size_t>(r)]!=static_cast<Eigen::Index>(k) ||
                partition.mappings->row_to_local[static_cast<std::size_t>(r)]!=static_cast<Eigen::Index>(row_index))
            {out.available=false; out.mapping_reason="component-row-mapping-mismatch"; return out;}
            mapped_rows[static_cast<std::size_t>(r)]=true;
        }
    }
    for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        if(!mapped_atoms[atom] || partition.mappings->atom_component[atom]<0)
        {out.available=false; out.mapping_reason="unmapped-component-atom"; return out;}
    if(!std::all_of(mapped_rows.begin(),mapped_rows.end(),[](bool mapped){return mapped;}))
    {out.available=false; out.mapping_reason="unmapped-observation-row"; return out;}
    for(const auto & support:domain.atoms)
        for(const auto & entry:support)
            if(entry.row<0 || entry.row>=domain.rows || !mapped_rows[static_cast<std::size_t>(entry.row)] ||
                std::binary_search(partition.constant_rows.begin(),partition.constant_rows.end(),entry.row))
            {out.available=false; out.mapping_reason="constant-row-mapping-mismatch"; return out;}
    out.mapping_valid=true;
    for(std::size_t k=0;k<components.size();++k)
    {
        const auto & component=components[k]; const auto & view=partition.components[k];
        if(!component.state)
        {
            out.available=false;
            if(out.mapping_reason.empty()) out.mapping_reason="missing-component-state";
            continue;
        }
        if(component.RuntimeConvergence()!=JointCheckStatus::Passed)
        {
            out.available=false;
            if(out.mapping_reason.empty()) out.mapping_reason="component-not-converged";
            continue;
        }
        const auto & state=*component.state;
        if(state.ac.size()!=2*view.atoms.size() || state.b.size()!=view.atoms.size() ||
            state.log_b.size()!=view.atoms.size() || state.width_gradient.size()!=view.atoms.size())
        {
            out.available=false;
            if(out.mapping_reason.empty()) out.mapping_reason="invalid-component-state-dimensions";
            continue;
        }
        for(std::size_t a=0;a<view.atoms.size();++a)
        {
            out.eta(view.atoms[a])=state.log_b[a];
            out.beta(2*view.atoms[a])=state.ac[2*a];
            out.beta(2*view.atoms[a]+1)=state.ac[2*a+1];
        }
        for(const auto row:view.rows) out.row_mask[static_cast<std::size_t>(row)]=true;
    }
    if(!out.available) return out;
    out.row_mask=std::move(mapped_rows);
    if(!out.eta.allFinite() || !out.beta.allFinite())
    {out.available=false; out.mapping_reason="invalid-assembled-parameter-state"; return out;}
    out.verification_memory_estimate=EndpointVerificationMemoryEstimate(domain,y,context);
    if(out.verification_memory_estimate.status!=MemoryStatus::Pass)
    {
        LogMemoryEstimate(out.verification_memory_estimate);
        out.assessment.failure="resource-limit:endpoint-verification:"+out.verification_memory_estimate.reason;
        out.assessment.verification_memory_estimate=out.verification_memory_estimate;
        out.assessment.primary.eta=out.eta; out.assessment.primary.reason=out.assessment.failure;
        out.replay.reason=out.assessment.failure;
        return out;
    }
    if(out.verification_memory_estimate.estimated_peak_bytes>out.verification_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.verification_memory_estimate);
    const auto primary_started=Clock::now();
    const auto raw=EvaluateStateStreaming(domain,y,out.eta,out.beta,context);
    const double primary_seconds=Seconds(primary_started);
    out.raw=raw;
    const auto replay_started=Clock::now();
    out.replay=CheckIndependentReplay(domain,y,raw,context,nullptr,nullptr,convergence_criteria);
    const double replay_seconds=Seconds(replay_started);
    out.consistency_passed=out.replay.passed;
    bool same=reuse && SameAssessmentPolicy(reuse->context,context) && reuse->domain.rows==domain.rows &&
        reuse->domain.atoms.size()==domain.atoms.size() && reuse->observations.size()==y.size() &&
        (reuse->observations.array()==y.array()).all() && reuse->assessment.primary.eta.size()==out.eta.size() &&
        reuse->assessment.primary.beta.size()==out.beta.size() &&
        (reuse->assessment.primary.eta.array()==out.eta.array()).all() && (reuse->assessment.primary.beta.array()==out.beta.array()).all();
    if(same) for(std::size_t a=0;a<domain.atoms.size();++a)
    {
        const auto & left=reuse->domain.atoms[a], & right=domain.atoms[a];
        same &= left.size()==right.size();
        if(same) for(std::size_t k=0;k<left.size();++k) same &= left[k].row==right[k].row && left[k].square==right[k].square;
    }
    if(same)
    {
        out.assessment=VerifyJointEndpoint(domain,y,raw,out.replay,context,convergence_criteria);
    }
    else
        out.assessment=VerifyJointEndpoint(domain,y,raw,out.replay,context,convergence_criteria);
    out.assessment.verification_memory_estimate=out.verification_memory_estimate;
    out.assessment.primary_endpoint_seconds=primary_seconds;
    out.assessment.independent_replay_seconds=replay_seconds;
    if(raw.valid) {out.prediction=raw.residual+y; out.objective=raw.certificate.objective;}
    return out;
}

}
