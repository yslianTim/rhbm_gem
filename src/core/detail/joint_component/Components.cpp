#include "Numerics.hpp"
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
ComponentPartition Partition(const Domain & domain,const std::vector<std::string> & ids)
{
    const auto count=static_cast<Eigen::Index>(domain.atoms.size());
    if(domain.rows<0 || ids.size()!=domain.atoms.size() || std::set<std::string>(ids.begin(),ids.end()).size()!=ids.size())
        throw std::invalid_argument("Invalid structural partition identities.");
    std::vector<Eigen::Index> parent(static_cast<std::size_t>(count)); std::iota(parent.begin(),parent.end(),0);
    auto root=[&](Eigen::Index a) {while(parent[Index(a)]!=a) {parent[Index(a)]=parent[Index(parent[Index(a)])]; a=parent[Index(a)];} return a;};
    std::vector<std::vector<Eigen::Index>> contributors(static_cast<std::size_t>(domain.rows));
    ComponentPartition out; out.atom_component.assign(ids.size(),-1); out.row_component.assign(static_cast<std::size_t>(domain.rows),-1);
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
    for(std::size_t c=0;c<out.components.size();++c) for(auto a:out.components[c].atoms) out.atom_component[Index(a)]=static_cast<Eigen::Index>(c);
    for(Eigen::Index r=0;r<domain.rows;++r)
    {
        if(contributors[Index(r)].empty()) {out.constant_rows.push_back(r); continue;}
        const auto c=out.atom_component[Index(contributors[Index(r)][0])]; out.row_component[Index(r)]=c; out.components[Index(c)].rows.push_back(r);
    }
    for(auto & v:out.components)
    {
        v.atom_to_local.assign(ids.size(),-1); v.row_to_local.assign(static_cast<std::size_t>(domain.rows),-1);
        for(std::size_t r=0;r<v.rows.size();++r) v.row_to_local[Index(v.rows[r])]=static_cast<Eigen::Index>(r);
        v.domain.rows=static_cast<Eigen::Index>(v.rows.size()); v.domain.atoms.resize(v.atoms.size());
        for(std::size_t a=0;a<v.atoms.size();++a)
        {
            v.atom_to_local[Index(v.atoms[a])]=static_cast<Eigen::Index>(a);
            for(const auto & s:domain.atoms[Index(v.atoms[a])]) v.domain.atoms[a].push_back({v.row_to_local[Index(s.row)],s.square});
        }
    }
    return out;
}
Eigen::VectorXd SelectValues(const Eigen::VectorXd & v,const std::vector<Eigen::Index> & indices)
{Eigen::VectorXd out(static_cast<Eigen::Index>(indices.size())); for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=v(indices[k]); return out;}
EvaluationContext ChildContext(const EvaluationContext & parent,const ComponentView & view,bool independent)
{
    auto c=parent; c.atom_ids.clear(); c.row_ids.clear();
    for(auto a:view.atoms) c.atom_ids.push_back(parent.atom_ids.at(static_cast<std::size_t>(a)));
    for(auto r:view.rows) c.row_ids.push_back(parent.row_ids.at(static_cast<std::size_t>(r)));
    if(parent.audit.directions.size())
    {
        c.audit.directions.resize(static_cast<Eigen::Index>(view.atoms.size()),parent.audit.directions.cols());
        for(std::size_t a=0;a<view.atoms.size();++a) c.audit.directions.row(static_cast<Eigen::Index>(a))=parent.audit.directions.row(view.atoms[a]);
    }
    c.audit.boundary_atoms.clear();
    for(auto a:parent.audit.boundary_atoms) if(view.atom_to_local.at(static_cast<std::size_t>(a))>=0) c.audit.boundary_atoms.push_back(view.atom_to_local[Index(a)]);
    c.audit.boundary=!c.audit.boundary_atoms.empty();
    c.independent_search=independent;
    if(independent)
    {
        const auto n=static_cast<Eigen::Index>(view.atoms.size()); c.rank={view.domain.rows,2*n,n};
        c.linear.rank_relative=c.rank.Relative(2*n);
    }
    return c;
}
ComponentResult SolveComponent(const ComponentView & view,const Vector & y,const Vector & initial_b,const EvaluationContext & parent)
{
    const auto context=ChildContext(parent,view,true);
    const Vector local_y=SelectValues(y,view.rows),start=SelectValues(initial_b,view.atoms);
    return AssessComponentSearch(view.domain,local_y,context,SearchProfile(view.domain,local_y,start,context));
}
ComponentResult AssessComponentSearch(const Domain & domain,const Vector & y,const EvaluationContext & context,SearchResult search)
{
    ComponentResult out; out.search=std::move(search);
    const auto audit_start=std::chrono::steady_clock::now();
    const auto endpoint=EvaluateProfile(domain,y,out.search.eta,false,&context);
    const auto reference=EvaluateProfile(domain,y,out.search.eta,true,&context);
    out.assessment=AssessEvaluated(domain,y,endpoint,reference,context);
    for(std::size_t k=0;k<out.search.trials.size();++k)
    {
        const auto & trial=out.search.trials[k];
        if(trial.accepted && trial.trust && trial.trust->passed) {out.trusted_state=trial.endpoint; out.trusted_trial=k;}
    }
    if(out.trusted_state && out.assessment.primary.valid)
    {
        out.endpoint_trust=CheckTrust(domain,y,endpoint,context,reference);
        if(out.endpoint_trust->passed) {out.trusted_state=out.assessment.primary; out.trusted_trial.reset();}
    }
    if(out.trusted_state)
    {
        const auto & state=*out.trusted_state;
        const bool same=state.eta.size()==endpoint.eta.size() && state.beta.size()==endpoint.beta.size() &&
            (state.eta.array()==endpoint.eta.array()).all() && (state.beta.array()==endpoint.beta.array()).all();
        if(same) out.trusted_assessment=out.assessment;
        else out.trusted_assessment=AssessProfile(domain,y,state.eta,context,&state.beta);
    }
    out.assessment_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-audit_start).count();
    out.search_success=out.trusted_state.has_value() && !out.search.stopped; return out;
}
AssemblyResult AssembleComponents(const Domain & domain,const Vector & y,const ComponentPartition & partition,
    const EvaluationContext & context,const std::vector<ComponentResult> & fits,const AssessmentReuse * reuse)
{
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
