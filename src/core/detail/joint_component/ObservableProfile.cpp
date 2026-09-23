#include "Problem.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
Indices IndicesOf(const std::vector<std::size_t> & values) {return {values.begin(),values.end()};}
std::vector<double> Values(const Vector & v) {return {v.data(),v.data()+v.size()};}
JointState State(const Endpoint & e,double scale)
{return {Values(e.beta),Values(e.eta.array().exp()),Values(e.eta),Values(e.gradient),e.certificate.objective/(scale*scale),{}};}
std::vector<JointCheck> AnalyticEvidence(JointEvidenceScope scope)
{
    auto out=AssessmentEvidence(Assessment{},scope);
    for(auto & check:out) if(check.status==JointCheckStatus::Unavailable)
    {check.status=JointCheckStatus::Passed; check.reason="analytic-nuisance-only"; check.value=0;}
    return out;
}
}
static JointParameterLayout ComponentLayout(const JointParameterLayout & layout,const ComponentView & view)
{
    JointParameterLayout out;
    for(auto a:layout.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(a))>=0) out.full_atoms.push_back(a);
    for(auto r:layout.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(r))>=0) out.informative_rows.push_back(r);
    for(const auto & group:layout.groups) if(view.LocalRow(static_cast<Eigen::Index>(group.row))>=0) out.groups.push_back(group);
    return out;
}
namespace {
void Reconstruct(const JointProblemInput & input,const JointParameterLayout & layout,JointState & state,Vector & prediction)
{
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
        for(const auto & support:input.support[layout.full_atoms[k]])
        {
            const auto basis=EvaluateKernel(support.squared_distance,state.b[k],2.5);
            prediction(static_cast<Eigen::Index>(support.row))+=state.ac[2*k]*basis.gaussian+state.ac[2*k+1]*basis.charge;
        }
    for(const auto & group:layout.groups)
    {
        const auto row=static_cast<Eigen::Index>(group.row);
        const double lambda=input.observations[group.row]-prediction(row);
        state.nuisance_amplitudes.push_back(lambda); prediction(row)+=lambda;
    }
}
}
JointParameterLayout BuildParameterLayout(const JointProblemInput & input)
{
    JointParameterLayout out; std::map<std::size_t,std::vector<std::size_t>> groups;
    for(std::size_t a=0;a<input.atom_ids.size();++a)
    {
        const bool halo=input.selection_domain && !std::binary_search(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),a);
        if(halo && input.support[a].size()==1) groups[input.support[a][0].row].push_back(a);
        else out.full_atoms.push_back(a);
    }
    for(auto & [row,atoms]:groups) out.groups.push_back({row,std::move(atoms)});
    for(std::size_t r=0;r<input.observations.size();++r) if(!groups.contains(r)) out.informative_rows.push_back(r);
    return out;
}
Domain ProfileDomain(const Domain & domain,const JointParameterLayout & layout)
{
    auto mapping=std::make_shared<Indices>(static_cast<std::size_t>(domain.rows),-1);
    for(std::size_t r=0;r<layout.informative_rows.size();++r) mapping->at(layout.informative_rows[r])=static_cast<Eigen::Index>(r);
    return domain.Select(IndicesOf(layout.full_atoms),static_cast<Eigen::Index>(layout.informative_rows.size()),mapping);
}
EvaluationContext ProfileContext(const EvaluationContext & parent,const JointParameterLayout & layout,Eigen::Index original_rows)
{
    auto context=parent;
    context.atom_ids=parent.atom_ids.Select(IndicesOf(layout.full_atoms));
    context.row_ids=parent.row_ids.Select(IndicesOf(layout.informative_rows));
    const auto atoms=static_cast<Eigen::Index>(layout.full_atoms.size());
    context.rank={original_rows,2*atoms,atoms}; context.linear.rank_relative=context.rank.Relative(2*atoms);
    return context;
}
JointFitResult FitObservableComponents(const JointProblem & problem,const std::vector<double> & initial_b)
{
    const auto & data=JointProblemAccess::Get(problem); const auto & input=*data.input;
    JointFitResult out; out.problem=problem; out.layout=data.layout; out.observation_scale=data.context.scale;
    out.initialization.b=initial_b;
    out.initialization.valid=initial_b.size()==input.atom_ids.size() && std::all_of(data.layout.full_atoms.begin(),data.layout.full_atoms.end(),[&](auto a){return std::isfinite(initial_b[a]) && initial_b[a]>0;});
    out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
    out.available_row_mask.assign(input.observations.size(),false);
    for(auto r:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(r)]=true;
    if(initial_b.size()!=input.atom_ids.size()) return out;
    const VectorMap widths(initial_b.data(),static_cast<Eigen::Index>(initial_b.size()));
    Vector prediction=Vector::Zero(data.y.size());
    for(const auto & view:data.partition.components)
    {
        const auto started=Clock::now();
        const double assembly_before=out.costs.assembly_seconds;
        JointComponentResult component; component.id=view.id;
        component.atoms.assign(view.atoms.begin(),view.atoms.end()); component.rows.assign(view.rows.begin(),view.rows.end());
        component.layout=ComponentLayout(data.layout,view); const auto & layout=*component.layout;
        const bool valid=std::all_of(layout.full_atoms.begin(),layout.full_atoms.end(),[&](auto a){return std::isfinite(initial_b[a]) && initial_b[a]>0;});
        if(layout.full_atoms.empty())
        {
            component.state.emplace(); component.search_completed=true; component.stop_reason="analytic-nuisance-only";
            component.evidence=AnalyticEvidence(JointEvidenceScope::ComponentLocal);
        }
        else if(valid && !layout.informative_rows.empty())
        {
            const auto domain=ProfileDomain(data.domain,layout);
            auto context=ProfileContext(data.context,layout,static_cast<Eigen::Index>(view.rows.size())); context.independent_search=true;
            const Vector y=SelectValues(data.y,IndicesOf(layout.informative_rows));
            const Vector start=SelectValues(widths,IndicesOf(layout.full_atoms));
            const auto result=AssessComponentSearch(domain,y,context,SearchProfile(domain,y,start,context));
            component.search_completed=result.search_success; component.stop_reason=result.search.stop_reason;
            component.profile_evaluations=result.search.evaluations; component.reference_evaluations=result.search.references;
            component.accepted_updates=result.search.accepted; component.native_status=result.search.lm_status;
            out.costs.search_seconds+=result.search.seconds; out.costs.search_reference_seconds+=result.search.reference_seconds;
            if(result.trusted_state)
            {
                component.state=State(*result.trusted_state,data.context.scale);
                component.evidence=AssessmentEvidence(*result.trusted_assessment,JointEvidenceScope::ComponentLocal);
                component.ranks=AssessmentRanks(*result.trusted_assessment,JointEvidenceScope::ComponentLocal);
            }
        }
        else component.stop_reason=valid ? "unobserved-full-parameters" : "invalid-initial-widths";
        if(component.state)
        {
            const auto reconstruction_start=Clock::now();
            Reconstruct(input,layout,*component.state,prediction);
            for(auto row:component.rows) out.available_row_mask[row]=true;
            out.costs.assembly_seconds+=Seconds(reconstruction_start);
        }
        else
        {
            Assessment missing; missing.failure=component.stop_reason;
            component.evidence=AssessmentEvidence(missing,JointEvidenceScope::ComponentLocal);
        }
        out.costs.assessment_seconds+=Seconds(started)-(out.costs.assembly_seconds-assembly_before);
        out.components.push_back(std::move(component));
    }
    // Search is included in the component elapsed time above.
    out.costs.assessment_seconds-=out.costs.search_seconds;
    out.search_completed=std::all_of(out.components.begin(),out.components.end(),[](const auto & c){return c.search_completed;});
    if(std::any_of(out.components.begin(),out.components.end(),[](const auto & c){return !c.state;}))
    {
        out.evidence=AssessmentEvidence(Assessment{},JointEvidenceScope::AssembledGlobal); return out;
    }
    const auto started=Clock::now(); const auto & layout=data.layout;
    JointState state; const auto count=layout.full_atoms.size();
    state.ac.resize(2*count); state.b.resize(count); state.log_b.resize(count); state.width_gradient.resize(count);
    std::vector<std::size_t> position(input.atom_ids.size());
    for(std::size_t k=0;k<count;++k) position[layout.full_atoms[k]]=k;
    std::map<std::size_t,double> amplitudes;
    for(const auto & component:out.components)
    {
        const auto & local=*component.state;
        for(std::size_t a=0;a<component.layout->full_atoms.size();++a)
        {
            const auto k=position[component.layout->full_atoms[a]];
            state.ac[2*k]=local.ac[2*a]; state.ac[2*k+1]=local.ac[2*a+1];
            state.b[k]=local.b[a]; state.log_b[k]=local.log_b[a]; state.width_gradient[k]=local.width_gradient[a];
        }
        for(std::size_t k=0;k<component.layout->groups.size();++k) amplitudes[component.layout->groups[k].row]=local.nuisance_amplitudes[k];
    }
    for(const auto & group:layout.groups) state.nuisance_amplitudes.push_back(amplitudes.at(group.row));
    state.objective=.5*(prediction-data.y).squaredNorm()/(data.context.scale*data.context.scale);
    bool agrees=true; double difference=0;
    Vector expected=Vector::Zero(data.y.size());
    for(const auto & group:layout.groups) expected(static_cast<Eigen::Index>(group.row))=data.y(static_cast<Eigen::Index>(group.row));
    if(count)
    {
        const auto domain=ProfileDomain(data.domain,layout); const auto context=ProfileContext(data.context,layout,data.y.size());
        const Vector y=SelectValues(data.y,IndicesOf(layout.informative_rows));
        const Vector eta=Eigen::Map<const Vector>(state.log_b.data(),static_cast<Eigen::Index>(count));
        const Vector beta=Eigen::Map<const Vector>(state.ac.data(),static_cast<Eigen::Index>(2*count));
        const auto assessment=AssessProfile(domain,y,eta,context,&beta);
        out.evidence=AssessmentEvidence(assessment,JointEvidenceScope::AssembledGlobal);
        out.ranks=AssessmentRanks(assessment,JointEvidenceScope::AssembledGlobal);
        const auto control=EvaluateProfile(domain,y,eta,false,&context);
        agrees=control.valid;
        if(agrees)
        {
            const Vector full_prediction=control.x*beta;
            for(std::size_t r=0;r<layout.informative_rows.size();++r)
                expected(static_cast<Eigen::Index>(layout.informative_rows[r]))=full_prediction(static_cast<Eigen::Index>(r));
            difference=((beta-control.beta).array().abs()/(1+beta.array().abs().max(control.beta.array().abs()))).maxCoeff();
            agrees=difference<=1e-10;
        }
    }
    else out.evidence=AnalyticEvidence(JointEvidenceScope::AssembledGlobal);
    const double prediction_difference=((prediction-expected).array().abs()/
        (1+prediction.array().abs().max(expected.array().abs()))).maxCoeff();
    const double profile_objective=.5*(expected-data.y).squaredNorm()/(data.context.scale*data.context.scale);
    const bool reconstruction=prediction.allFinite() && expected.allFinite() && prediction_difference<=1e-10 &&
        std::abs(state.objective-profile_objective)<=1e-12;
    difference=std::max(difference,prediction_difference); agrees &= reconstruction;
    out.evidence.push_back({"assembled-profile",agrees ? JointCheckStatus::Passed : JointCheckStatus::Failed,
        JointEvidenceScope::AssembledGlobal,difference,1e-10,reconstruction ? "" : "full-domain-reconstruction-failed"});
    out.assembled_state=std::move(state); out.prediction=Values(prediction); out.objective=out.assembled_state->objective;
    out.costs.assembly_seconds+=Seconds(started); return out;
}
}
