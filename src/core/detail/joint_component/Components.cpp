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
}
std::vector<JointCheck> AnalyticNuisanceEvidence(JointEvidenceScope scope)
{
    auto out=AssessmentEvidence(Assessment{},scope);
    for(auto & check:out) if(check.status==JointCheckStatus::Unavailable)
    {check.status=JointCheckStatus::Passed; check.reason="analytic-nuisance-only"; check.value=0;}
    return out;
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
            if(fit.search.certification)
            {
                const auto & certification=*fit.search.certification;
                component.evidence=AssessmentEvidence(certification.assessment,JointEvidenceScope::ComponentLocal);
                component.ranks=AssessmentRanks(certification.assessment,JointEvidenceScope::ComponentLocal);
                component.evidence.push_back({"endpoint-trust",certification.passed() ? JointCheckStatus::Passed :
                    JointCheckStatus::Failed,JointEvidenceScope::ComponentLocal,{},{},certification.trust.reason});
                const auto add_width_gradient=[&](const char * name,const Endpoint & evaluation) {
                    const bool available=evaluation.valid && evaluation.gradient.size()>0 && evaluation.gradient.allFinite();
                    const double value=available ? evaluation.gradient.lpNorm<Eigen::Infinity>() : 0.0;
                    component.evidence.push_back({name,!available ? JointCheckStatus::Unavailable :
                        value<=1e-12 ? JointCheckStatus::Passed : JointCheckStatus::Failed,
                        JointEvidenceScope::ComponentLocal,available ? std::optional<double>(value) : std::nullopt,
                        1e-12,available ? "" : evaluation.reason});
                };
                add_width_gradient("endpoint-primary-width-gradient",certification.assessment.primary);
                add_width_gradient("endpoint-reference-width-gradient",certification.assessment.reference);
                if(certification.passed())
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
            fit.search.certification && fit.search.certification->passed())
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
    const auto prepared=MakeComponentProblemView(component_layout,view);
    return SolveFixedNeighborComponent(prepared,initial_b,parent_context,
        production_policy,observer,progress_component);
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
    if(!out.available) return out;
    const auto raw=EvaluateState(domain,y,out.eta,out.beta,context);
    out.raw=raw;
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
        // Preserve raw-state metadata; the shared numerical evidence is unchanged.
        out.assessment.primary=raw;
    }
    else
    {
        const auto reference=EvaluateProfile(domain,y,out.eta,true,&context);
        out.assessment=AssessEvaluated(domain,y,raw,reference,context,true);
    }
    if(raw.valid) {out.prediction=raw.x*out.beta; out.objective=.5*raw.residual.squaredNorm();}
    if(raw.valid && out.assessment.primary.valid)
    {
        const auto profile=EvaluateProfile(domain,y,out.eta,false,&context);
        out.profile_control=profile; out.profile_evaluated=true;
        if(!profile.valid) return out;
        out.profile_difference=((out.beta-profile.beta).array().abs()/(1+out.beta.array().abs().max(profile.beta.array().abs()))).maxCoeff();
        out.profile_agrees=out.profile_difference<=1e-10;
    }
    return out;
}

}
