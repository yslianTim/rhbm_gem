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
    out.local_to_parent_atoms.reserve(parent_layout.full_atoms.size());
    out.local_to_parent_rows.reserve(parent_layout.informative_rows.size());
    local_to_view_atoms.reserve(parent_layout.full_atoms.size());
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
        out.local_layout.full_atoms.push_back(k);
    }
    const auto row_mapping=CompactRowMapping::FromParentRows(view_to_local_rows);
    out.domain=view.domain.Select(local_to_view_atoms,static_cast<Eigen::Index>(out.local_layout.informative_rows.size()),row_mapping);
    const Identities parent_atoms(std::shared_ptr<const std::vector<std::string>>(parent,&parent->atom_ids));
    const Identities parent_rows(std::shared_ptr<const std::vector<std::string>>(parent,&parent->row_ids));
    out.atom_ids=parent_atoms.Select(out.local_to_parent_atoms);
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
{return {Values(e.beta),Values(e.eta.array().exp()),Values(e.eta),Values(e.gradient),e.certificate.objective/(scale*scale),{}};}
JointParameterLayout ComponentLayout(const JointParameterLayout & layout,const ComponentView & view)
{
    JointParameterLayout out;
    for(auto a:layout.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(a))>=0) out.full_atoms.push_back(a);
    for(auto r:layout.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(r))>=0) out.informative_rows.push_back(r);
    for(const auto & group:layout.groups) if(view.LocalRow(static_cast<Eigen::Index>(group.row))>=0) out.groups.push_back(group);
    return out;
}
JointCandidateDiagnostic CandidateDiagnostic(const FixedNeighborSearchResult & search,std::size_t atoms)
{
    JointCandidateDiagnostic out;
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
        if(std::isfinite(sweep.global_a_feasibility)) out.candidate_feasible=sweep.global_a_feasibility==0;
    }
    if(search.certification)
    {
        const auto & certification=*search.certification;
        out.reference_status=JointCheckStatus::NotRun;
        out.endpoint_verification_status=certification.verified() ? JointCheckStatus::Passed :
            certification.verification_memory_estimate.status!=MemoryStatus::Pass ||
                certification.assessment.failure.starts_with("resource-limit:") ?
                    JointCheckStatus::Unavailable : JointCheckStatus::Failed;
        const bool rank_unavailable=certification.assessment.rank_memory_estimate.status!=MemoryStatus::Pass ||
            !certification.assessment.design || !certification.assessment.widths ||
            !certification.assessment.normalized_widths || !certification.assessment.jacobian;
        const bool rank_identified=!rank_unavailable &&
            certification.assessment.design->rank==certification.assessment.design->columns &&
            certification.assessment.identified;
        out.rank_assessment_status=rank_unavailable ? JointCheckStatus::Unavailable :
            rank_identified ? JointCheckStatus::Passed : JointCheckStatus::Failed;
        out.certification_status=certification.qualified() ? JointCheckStatus::Passed :
            certification.verification_memory_estimate.status!=MemoryStatus::Pass ||
                certification.assessment.rank_memory_estimate.status!=MemoryStatus::Pass ?
                JointCheckStatus::Unavailable : JointCheckStatus::Failed;
        out.numerical_qualification_status=out.certification_status;
        out.verification_reason=!certification.trust.reason.empty() ? certification.trust.reason :
            certification.assessment.failure;
        if(certification.verification_memory_estimate.status!=MemoryStatus::Pass)
            out.resource_limit_reason="resource-limit:endpoint-verification:"+certification.verification_memory_estimate.reason;
        else if(certification.assessment.rank_memory_estimate.status!=MemoryStatus::Pass)
            out.resource_limit_reason="resource-limit:rank-assessment:"+certification.assessment.rank_memory_estimate.reason;
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
std::vector<JointCheck> AnalyticNuisanceEvidence(JointEvidenceScope scope)
{
    constexpr auto passed=JointCheckStatus::Passed;
    constexpr auto reason="analytic-nuisance-only";
    return {
        {"ac-feasibility",passed,scope,0.0,0.0,reason},
        {"conditional-ac-optimality",passed,scope,0.0,1e-10,reason},
        {"independent-ac-optimality",passed,scope,0.0,1e-10,reason},
        {"independent-kkt-consistency",passed,scope,0.0,1e-13,reason},
        {"independent-replay",passed,scope,0.0,{},reason},
        {"width-gradient-consistency",passed,scope,0.0,{},reason},
        {"width-stationarity",passed,scope,0.0,1e-12,reason},
        {"design-identifiability",passed,scope,0.0,0.0,reason},
        {"local-correction",passed,scope,0.0,1e-10,reason},
        {"numerical-identifiability",passed,scope,0.0,0.0,reason}};
}
FixedNeighborComponentBatch OrchestrateFixedNeighborComponents(
    const JointProblem & problem,const std::vector<double> & initial_b,
    const FixedNeighborSearchPolicy & production_policy,const JointProgressObserver & observer)
{
    const auto & data=JointProblemAccess::Get(problem);
    const bool has_observable_groups=!data.layout.groups.empty();
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
        if(has_observable_groups) component.layout=layout;
        FixedNeighborComponentResult fit;
        std::optional<double> accepted_objective,accepted_gradient_inf_norm;
        const bool valid=std::all_of(layout.full_atoms.begin(),layout.full_atoms.end(),[&](auto a) {
            return a<initial_b.size() && std::isfinite(initial_b[a]) && initial_b[a]>0;
        });
        if(has_observable_groups && layout.full_atoms.empty())
        {
            component.state.emplace(); component.search_completed=true; component.stop_reason="analytic-nuisance-only";
            component.evidence=AnalyticNuisanceEvidence(JointEvidenceScope::ComponentLocal);
        }
        else if(valid && !layout.informative_rows.empty())
        {
            fit=SolveFixedNeighborComponentView(view,layout,widths,data.context,production_policy,observer,
                progress_component ? &*progress_component : nullptr);
            accepted_objective=LastFixedNeighborObjective(fit.search);
            accepted_gradient_inf_norm=LastFixedNeighborGradient(fit.search);
            component.search_completed=fit.succeeded; component.stop_reason=fit.stop_reason;
            component.search_sweeps=fit.search.sweep_count;
            component.first_order_stationarity_sweep=fit.search.first_order_stationarity_sweep;
            component.confirmed_stationarity_sweep=fit.search.confirmed_stationarity_sweep;
            component.endpoint_certification_sweep=fit.search.certification_sweep;
            component.profile_evaluations=static_cast<int>(fit.search.total_profile_evaluations);
            component.accepted_updates=static_cast<int>(fit.search.total_accepted_local_updates);
            component.candidate_diagnostic=CandidateDiagnostic(fit.search,layout.full_atoms.size());
            if(fit.search.certification)
            {
                const auto & certification=*fit.search.certification;
                component.evidence=AssessmentEvidence(certification.assessment,JointEvidenceScope::ComponentLocal);
                component.ranks=AssessmentRanks(certification.assessment,JointEvidenceScope::ComponentLocal);
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
        if(!has_observable_groups && component_count==1 && data.partition.constant_rows.empty() &&
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
    shape.component_target_columns=parent_context.target_columns;
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
    out=SolveFixedNeighborComponent(prepared,initial_b,parent_context,
        production_policy,observer,progress_component);
    return out;
}
AssemblyResult AssembleComponents(const Domain & domain,VectorRef y,const ComponentPartition & partition,
    const EvaluationContext & context,const std::vector<JointComponentResult> & components,const AssessmentReuse * reuse)
{
    ResourcePhase phase("assembly");
    if(components.size()!=partition.components.size()) throw std::invalid_argument("Missing component result.");
    AssemblyResult out; out.available=true; out.row_mask.assign(static_cast<std::size_t>(domain.rows),false);
    out.eta=Vector::Constant(static_cast<Eigen::Index>(domain.atoms.size()),unavailable);
    out.beta=Vector::Constant(2*out.eta.size(),unavailable);
    for(auto r:partition.constant_rows) out.row_mask[static_cast<std::size_t>(r)]=true;
    for(std::size_t k=0;k<components.size();++k)
    {
        const auto & component=components[k]; const auto & view=partition.components[k];
        out.completed &= component.search_completed;
        if(!component.state) {out.available=false; continue;}
        const auto & state=*component.state;
        if(state.ac.size()!=2*view.atoms.size() || state.log_b.size()!=view.atoms.size())
            throw std::invalid_argument("Invalid certified component state.");
        for(std::size_t a=0;a<view.atoms.size();++a)
        {
            out.eta(view.atoms[a])=state.log_b[a];
            out.beta(2*view.atoms[a])=state.ac[2*a];
            out.beta(2*view.atoms[a]+1)=state.ac[2*a+1];
        }
        for(auto r:view.rows) out.row_mask[static_cast<std::size_t>(r)]=true;
    }
    out.mapping_valid=out.available && out.eta.allFinite() && out.beta.allFinite() &&
        std::all_of(out.row_mask.begin(),out.row_mask.end(),[](bool mapped){return mapped;});
    if(!out.mapping_valid) {out.available=false; return out;}
    out.verification_memory_estimate=EndpointVerificationMemoryEstimate(domain,y,context);
    if(out.verification_memory_estimate.status!=MemoryStatus::Pass)
    {
        LogMemoryEstimate(out.verification_memory_estimate);
        out.assessment.reference_free=true;
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
    out.replay=CheckIndependentReplay(domain,y,raw,context);
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
        out.assessment=reuse->assessment;
        out.assessment.independent_replay=out.replay;
        out.assessment.endpoint_verified=out.replay.passed && out.replay.primary_width_gradient<=1e-12;
        out.assessment.primary=raw;
    }
    else
        out.assessment=AssessEvaluated(domain,y,raw,out.replay,context,true);
    out.assessment.verification_memory_estimate=out.verification_memory_estimate;
    out.assessment.primary_endpoint_seconds=primary_seconds;
    out.assessment.independent_replay_seconds=replay_seconds;
    if(raw.valid) {out.prediction=raw.residual+y; out.objective=raw.certificate.objective;}
    return out;
}

}
