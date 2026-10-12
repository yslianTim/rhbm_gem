#include "Problem.hpp"
#include "FixedNeighborBlockCoordinate.hpp"
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
    Indices rows; rows.reserve(layout.informative_rows.size());
    for(const auto row:layout.informative_rows) rows.push_back(static_cast<Eigen::Index>(row));
    const auto mapping=CompactRowMapping::FromParentRows(rows);
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
JointFitResult FitObservableComponents(const JointProblem & problem,const std::vector<double> & initial_b,
    const FixedNeighborSearchPolicy & production_policy,const JointProgressObserver & observer)
{
    const auto & data=JointProblemAccess::Get(problem); const auto & input=*data.input;
    JointFitResult out; out.problem=problem; out.layout=data.layout; out.observation_scale=data.context.scale;
    out.initialization.b=initial_b;
    out.initialization.valid=initial_b.size()==input.atom_ids.size() && std::all_of(data.layout.full_atoms.begin(),data.layout.full_atoms.end(),[&](auto a){return std::isfinite(initial_b[a]) && initial_b[a]>0;});
    out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
    out.available_row_mask.assign(input.observations.size(),false);
    for(auto r:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(r)]=true;
    if(initial_b.size()!=input.atom_ids.size()) return out;
    Vector prediction=Vector::Zero(data.y.size());
    const auto component_count=data.partition.components.size();
    auto component_batch=OrchestrateFixedNeighborComponents(problem,initial_b,production_policy,observer);
    out.components=std::move(component_batch.components);
    out.costs.search_seconds=component_batch.costs.search_seconds;
    out.costs.assessment_seconds=component_batch.costs.assessment_seconds;
    out.search_completed=std::all_of(out.components.begin(),out.components.end(),[](const auto & c){return c.search_completed;});
    const auto assembly_started=Clock::now();
    NotifyJointAssemblyProgress(observer,JointProgressPhase::AssemblyStarted,component_count);
    for(auto & component:out.components) if(component.state)
    {
        const auto started=Clock::now();
        Reconstruct(input,*component.layout,*component.state,prediction);
        for(auto row:component.rows) out.available_row_mask[row]=true;
        out.costs.assembly_seconds+=Seconds(started);
    }
    if(std::any_of(out.components.begin(),out.components.end(),[](const auto & c){return !c.state;}))
    {
        out.evidence=AssessmentEvidence(Assessment{},JointEvidenceScope::AssembledGlobal);
        NotifyJointAssemblyProgress(observer,JointProgressPhase::AssemblyCompleted,component_count,
            Seconds(assembly_started),false);
        return out;
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
    Vector expected=Vector::Zero(data.y.size());
    for(const auto & group:layout.groups) expected(static_cast<Eigen::Index>(group.row))=data.y(static_cast<Eigen::Index>(group.row));
    Assessment assessment;
    bool profile_replay_passed=true;
    if(count)
    {
        const auto domain=ProfileDomain(data.domain,layout); const auto context=ProfileContext(data.context,layout,data.y.size());
        const Vector y=SelectValues(data.y,IndicesOf(layout.informative_rows));
        const Vector eta=Eigen::Map<const Vector>(state.log_b.data(),static_cast<Eigen::Index>(count));
        const Vector beta=Eigen::Map<const Vector>(state.ac.data(),static_cast<Eigen::Index>(2*count));
        assessment=AssessProfile(domain,y,eta,context,&beta);
        out.evidence=AssessmentEvidence(assessment,JointEvidenceScope::AssembledGlobal);
        out.ranks=AssessmentRanks(assessment,JointEvidenceScope::AssembledGlobal);
        profile_replay_passed=assessment.independent_replay.passed;
        if(assessment.primary.valid)
        {
            for(std::size_t r=0;r<layout.informative_rows.size();++r)
                expected(static_cast<Eigen::Index>(layout.informative_rows[r]))=
                    assessment.primary.residual(static_cast<Eigen::Index>(r))+y(static_cast<Eigen::Index>(r));
            for(std::size_t k=0;k<count;++k)
                state.width_gradient[k]=assessment.primary.gradient(static_cast<Eigen::Index>(k));
        }
    }
    else out.evidence=AnalyticNuisanceEvidence(JointEvidenceScope::AssembledGlobal);
    double prediction_difference{},prediction_error_ratio{}; bool prediction_passed=true;
    for(Eigen::Index row=0;row<data.y.size();++row)
    {
        const double error=std::abs(prediction(row)-expected(row));
        prediction_difference=std::max(prediction_difference,error);
        const double tolerance=2e-12+2e-13*std::abs(expected(row));
        prediction_error_ratio=std::max(prediction_error_ratio,error/tolerance);
        prediction_passed &= error<=tolerance;
    }
    const double raw_objective=.5*(prediction-data.y).squaredNorm();
    const double expected_raw_objective=.5*(expected-data.y).squaredNorm();
    const double objective_difference=std::abs(raw_objective-expected_raw_objective);
    const bool reconstruction=prediction.allFinite() && expected.allFinite() && prediction_passed &&
        objective_difference<=1e-12+2e-12*std::abs(expected_raw_objective);
    const bool mapping=std::all_of(out.available_row_mask.begin(),out.available_row_mask.end(),[](bool value){return value;}) &&
        std::all_of(state.ac.begin(),state.ac.end(),[](double value){return std::isfinite(value);}) &&
        std::all_of(state.log_b.begin(),state.log_b.end(),[](double value){return std::isfinite(value);});
    const bool consistency=mapping && reconstruction && profile_replay_passed;
    out.evidence.push_back({"assembly-mapping",mapping ? JointCheckStatus::Passed : JointCheckStatus::Failed,
        JointEvidenceScope::AssembledGlobal,{},{},mapping ? "exact-component-identity-mapping" : "incomplete-component-mapping"});
    out.evidence.push_back({"assembly-prediction-consistency",prediction_passed ? JointCheckStatus::Passed :
        JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,prediction_error_ratio,1.0,{}});
    out.evidence.push_back({"assembly-objective-consistency",
        objective_difference<=1e-12+2e-12*std::abs(expected_raw_objective) ? JointCheckStatus::Passed : JointCheckStatus::Failed,
        JointEvidenceScope::AssembledGlobal,objective_difference,1e-12+2e-12*std::abs(expected_raw_objective),{}});
    out.evidence.push_back({"assembled-consistency",consistency ? JointCheckStatus::Passed :
        !prediction.allFinite() || !expected.allFinite() ? JointCheckStatus::Unavailable : JointCheckStatus::Failed,
        JointEvidenceScope::AssembledGlobal,prediction_difference,2e-12,
        !profile_replay_passed ? assessment.independent_replay.reason :
            reconstruction ? "" : "full-domain-prediction-objective-mismatch"});
    if(!assessment.verification_memory_estimate.reason.empty())
        out.evidence.push_back({"assembly-verification-memory",JointCheckStatus::Unavailable,JointEvidenceScope::AssembledGlobal,
            static_cast<double>(assessment.verification_memory_estimate.estimated_peak_bytes),
            static_cast<double>(assessment.verification_memory_estimate.budget_bytes),
            "resource-limit:endpoint-verification:"+assessment.verification_memory_estimate.reason});
    if(!count || (assessment.qualified && consistency))
    {
        out.assembled_state=std::move(state); out.prediction=Values(prediction); out.objective=out.assembled_state->objective;
    }
    out.costs.assembly_seconds+=Seconds(started);
    NotifyJointAssemblyProgress(observer,JointProgressPhase::AssemblyCompleted,component_count,
        Seconds(assembly_started),out.assembled_state.has_value());
    return out;
}
}
