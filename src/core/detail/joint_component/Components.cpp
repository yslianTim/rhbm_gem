#include "Numerics.hpp"
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
        const std::shared_ptr<const Indices> rows(mappings,&mappings->row_to_local);
        v.domain=domain.Select(v.atoms,static_cast<Eigen::Index>(v.rows.size()),rows);
    }
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
ComponentResult SolveComponent(const ComponentView & view,VectorRef y,const Vector & initial_b,const EvaluationContext & parent)
{
    const auto context=ChildContext(parent,view,true);
    const Vector start=SelectValues(initial_b,view.atoms);
    const bool contiguous=!view.rows.empty() && view.rows.back()-view.rows.front()+1==static_cast<Eigen::Index>(view.rows.size());
    if(contiguous)
    {
        const auto local_y=y.segment(view.rows.front(),static_cast<Eigen::Index>(view.rows.size()));
        return AssessComponentSearch(view.domain,local_y,context,SearchProfile(view.domain,local_y,start,context));
    }
    const Vector local_y=SelectValues(y,view.rows);
    return AssessComponentSearch(view.domain,local_y,context,SearchProfile(view.domain,local_y,start,context));
}
ComponentResult AssessComponentSearch(const Domain & domain,VectorRef y,const EvaluationContext & context,SearchResult search)
{
    ComponentResult out; out.search=std::move(search);
    const auto audit_start=std::chrono::steady_clock::now();
    const auto endpoint=EvaluateProfile(domain,y,out.search.eta,false,&context);
    const auto reference=EvaluateProfile(domain,y,out.search.eta,true,&context);
    out.assessment=AssessEvaluated(domain,y,endpoint,reference,context);
    const bool accepted=std::any_of(out.search.trials.begin(),out.search.trials.end(),[](const Trial & trial) {
        return trial.accepted && trial.trust && trial.trust->passed;
    });
    if(accepted)
    {
        out.endpoint_trust=CheckTrust(domain,y,endpoint,context,reference);
        if(out.endpoint_trust->passed)
        {out.trusted_state=endpoint; out.trusted_assessment=out.assessment;}
        else
        {
            if(!out.search.stopped) out.search.stop_reason="endpoint-certification-failed";
            out.search.stopped=true;
            auto fallback_reference=reference;
            for(std::size_t k=out.search.trials.size();k>0;--k)
            {
                const auto & trial=out.search.trials[k-1];
                if(!trial.accepted || !trial.trust || !trial.trust->passed) continue;
                const auto & state=trial.endpoint;
                const auto candidate=EvaluateState(domain,y,state.eta,state.beta,context);
                if(!CheckReplay(domain,y,candidate,context).passed) continue;
                if(state.eta.size()!=fallback_reference.eta.size() ||
                    !(state.eta.array()==fallback_reference.eta.array()).all())
                    fallback_reference=EvaluateProfile(domain,y,state.eta,true,&context);
                if(!CheckTrust(domain,y,candidate,context,fallback_reference).passed) continue;
                out.trusted_state=state; out.trusted_trial=k-1;
                out.trusted_assessment=AssessEvaluated(domain,y,candidate,fallback_reference,context,true);
                break;
            }
        }
    }
    out.assessment_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-audit_start).count();
    out.search_success=out.trusted_state.has_value() && !out.search.stopped; return out;
}
AssemblyResult AssembleComponents(const Domain & domain,VectorRef y,const ComponentPartition & partition,
    const EvaluationContext & context,const std::vector<ComponentResult> & fits,const AssessmentReuse * reuse)
{
    ResourcePhase phase("assembly");
    if(fits.size()!=partition.components.size()) throw std::invalid_argument("Missing component result.");
    AssemblyResult out; out.available=true; out.row_mask.assign(static_cast<std::size_t>(domain.rows),false);
    out.eta=Vector::Constant(static_cast<Eigen::Index>(domain.atoms.size()),unavailable);
    out.beta=Vector::Constant(2*out.eta.size(),unavailable);
    for(auto r:partition.constant_rows) out.row_mask[static_cast<std::size_t>(r)]=true;
    for(std::size_t k=0;k<fits.size();++k)
    {
        const auto & fit=fits[k]; const auto & view=partition.components[k]; out.completed &= fit.search_success;
        if(!fit.trusted_state) {out.available=false; continue;}
        const auto & state=*fit.trusted_state;
        for(std::size_t a=0;a<view.atoms.size();++a)
        {
            out.eta(view.atoms[a])=state.eta(static_cast<Eigen::Index>(a));
            out.beta(2*view.atoms[a])=state.beta(static_cast<Eigen::Index>(2*a));
            out.beta(2*view.atoms[a]+1)=state.beta(static_cast<Eigen::Index>(2*a+1));
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
