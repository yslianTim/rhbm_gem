#include "support/JointFixedNeighborDiagnostic.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <algorithm>
#include <cmath>

namespace rhbm_gem::core::joint_component {
namespace {
Vector Select(VectorRef values,const Indices & indices)
{
    Vector out(static_cast<Eigen::Index>(indices.size()));
    for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=values(indices[k]);
    return out;
}
Indices IndicesOf(const std::vector<std::size_t> & values)
{return {values.begin(),values.end()};}
std::vector<double> Values(const Vector & values)
{return {values.data(),values.data()+values.size()};}
Vector SelectBeta(VectorRef beta,const std::vector<std::size_t> & atoms)
{
    Vector out(2*static_cast<Eigen::Index>(atoms.size()));
    for(std::size_t k=0;k<atoms.size();++k)
        out.segment<2>(2*static_cast<Eigen::Index>(k))=beta.segment<2>(2*static_cast<Eigen::Index>(atoms[k]));
    return out;
}
JointState State(const Endpoint & endpoint,double scale)
{
    return {Values(endpoint.beta),Values(endpoint.eta.array().exp()),Values(endpoint.eta),
        Values(endpoint.gradient),endpoint.certificate.objective/(scale*scale)};
}
void AssembleResult(const JointProblem & problem,const ProblemData & data,const JointParameterLayout & layout,
    const Domain & domain,VectorRef y,const EvaluationContext & context,
    const FixedNeighborEndpointCertification & certification,const JointConvergenceCriteria & convergence_criteria,
    FixedNeighborDiagnosticResult & out,bool accepted_state_integrity)
{
    out.assessment=certification.assessment;
    out.endpoint_trust=certification.trust;
    out.endpoint_certified=certification.verified();
    const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
    const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);
    JointCheck assembly_mapping,assembled_consistency;
    bool exact_identity_mapping=layout.full_atoms.size()==domain.atoms.size() &&
        layout.informative_rows.size()==static_cast<std::size_t>(domain.rows);
    for(std::size_t k=0;k<layout.full_atoms.size();++k) exact_identity_mapping &= layout.full_atoms[k]==k;
    for(std::size_t k=0;k<layout.informative_rows.size();++k) exact_identity_mapping &= layout.informative_rows[k]==k;

    const auto & component_view=data.partition.components.front();
    JointComponentResult component; component.id=component_view.id;
    component.atoms.assign(component_view.atoms.begin(),component_view.atoms.end());
    component.rows.assign(component_view.rows.begin(),component_view.rows.end());
    component.layout=layout; component.search_completed=out.search_converged; component.stop_reason=out.reason;
    component.search_sweeps=out.sweeps.size();
    component.first_order_stationarity_sweep=out.first_order_stationarity_sweep;
    component.confirmed_stationarity_sweep=out.confirmed_stationarity_sweep;
    component.endpoint_certification_sweep=out.sweeps.size();
    const bool qualified=certification.qualified() && accepted_state_integrity;
    if(qualified) component.state=State(out.assessment.primary,context.scale);
    component.evidence=AssessmentEvidence(out.assessment,JointEvidenceScope::ComponentLocal);
    component.evidence.push_back({"accepted-state-integrity",accepted_state_integrity ? JointCheckStatus::Passed :
        JointCheckStatus::Failed,JointEvidenceScope::ComponentLocal,{},{},{}});
    if(out.confirmed_stationarity_sweep)
        component.evidence.push_back({"complete-sweep-stationarity",JointCheckStatus::Passed,
            JointEvidenceScope::ComponentLocal,static_cast<double>(out.confirmed_stationarity_sweep),{},"confirmed"});
    component.evidence.push_back({"endpoint-trust",certification.verified() ? JointCheckStatus::Passed :
        JointCheckStatus::Failed,JointEvidenceScope::ComponentLocal,{},{},certification.trust.reason});

    out.fit.problem=problem; out.fit.initialization.valid=true;
    out.fit.initialization.reason="valid-widths"; out.fit.initialization.b=Values(out.state.eta.array().exp());
    out.fit.components.push_back(std::move(component)); out.fit.search_completed=out.search_converged;
    out.fit.observation_scale=context.scale;
    out.fit.available_row_mask.assign(data.input->observations.size(),false);
    if(qualified && exact_identity_mapping)
    {
        std::fill(out.fit.available_row_mask.begin(),out.fit.available_row_mask.end(),true);
        out.fit.assembled_state=out.fit.components.front().state;
        const auto assembled=EvaluateStateStreaming(domain,y,eta,beta,context);
        const auto replay=CheckIndependentReplay(domain,y,assembled,context,nullptr,nullptr,convergence_criteria);
        bool prediction_agrees=assembled.valid && assembled.residual.size()==out.state.prediction.size();
        if(prediction_agrees) for(Eigen::Index row=0;row<assembled.residual.size();++row)
        {
            const double prediction=assembled.residual(row)+y(row);
            const double error=std::abs(prediction-out.state.prediction(row));
            prediction_agrees &= error<=2e-12+2e-13*std::abs(prediction);
        }
        const double objective=assembled.certificate.objective/(context.scale*context.scale);
        const bool objective_agrees=std::isfinite(objective) &&
            std::abs(objective-out.state.objective)<=1e-12+2e-12*std::abs(objective);
        if(assembled.valid)
        {
            out.fit.prediction=Values(assembled.residual+y);
            out.fit.objective=objective;
        }
        assembly_mapping={"assembly-mapping",JointCheckStatus::Passed,
            JointEvidenceScope::AssembledGlobal,{}, {},"exact-component-identity-mapping"};
        const bool consistency=prediction_agrees && objective_agrees && replay.valid && replay.passed;
        assembled_consistency={"assembled-consistency",consistency ? JointCheckStatus::Passed :
            JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
            replay.available ? std::optional<double>(replay.prediction_difference) : std::nullopt,2e-12,
            consistency ? "verified" : "assembled-streaming-replay-failed"};
    }
    else if(!qualified)
    {
        assembly_mapping={"assembly-mapping",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},"missing-qualified-component-state"};
        assembled_consistency={"assembled-consistency",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},"missing-qualified-component-state"};
    }
    else
    {
        assembly_mapping={"assembly-mapping",JointCheckStatus::Failed,
            JointEvidenceScope::AssembledGlobal,{}, {},"component-layout-is-not-identity"};
        assembled_consistency={"assembled-consistency",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},"component-layout-is-not-identity"};
    }
    out.fit.evidence=AssessmentEvidence(out.assessment,JointEvidenceScope::AssembledGlobal);
    out.fit.evidence.push_back(std::move(assembly_mapping));
    out.fit.evidence.push_back(std::move(assembled_consistency));

}
}

FixedNeighborDiagnosticResult RunFixedNeighborDiagnostic(const JointProblem & problem,VectorRef initial_eta,
    const FixedNeighborPolicy & policy)
{
    FixedNeighborDiagnosticResult out;
    const auto & data=JointProblemAccess::Get(problem); const auto & input=*data.input; const auto & layout=data.layout;
    if(initial_eta.size()!=static_cast<Eigen::Index>(input.atom_ids.size()) || !initial_eta.allFinite() ||
        !(data.context.scale>0) || !std::isfinite(data.context.scale) || policy.core_atoms==0 ||
        layout.full_atoms.size()!=input.atom_ids.size() ||
        layout.informative_rows.size()!=input.observations.size() || data.partition.components.size()!=1)
    {out.reason="block-invalid-partition"; return out;}
    for(Eigen::Index a=0;a<initial_eta.size();++a)
        if(!(std::exp(initial_eta(a))>0) || !std::isfinite(std::exp(initial_eta(a))))
        {out.reason="block-invalid-partition"; return out;}
    if(!policy.eta_upper_bound || !std::isfinite(*policy.eta_upper_bound))
    {out.reason="invalid-width-upper-bound"; return out;}

    const Domain domain=data.domain;
    auto context=data.context;
    context.eta_upper_bound=policy.eta_upper_bound;
    const Vector y=Select(data.y,IndicesOf(layout.informative_rows));
    auto run_policy=policy;
    run_policy.collect_diagnostics=true;
    auto search=SearchFixedNeighborComponent(domain,layout,y,initial_eta,context,run_policy);
    out.state=std::move(search.state); out.blocks=std::move(search.blocks); out.sweeps=std::move(search.sweeps);
    out.search_converged=search.search_converged; out.reason=search.reason;
    out.first_order_stationarity_sweep=search.first_order_stationarity_sweep;
    out.confirmed_stationarity_sweep=search.confirmed_stationarity_sweep;
    out.sweep_count=search.sweep_count;
    out.total_block_solves=search.total_block_solves;
    out.total_profile_evaluations=search.total_profile_evaluations;
    out.total_accepted_local_updates=search.total_accepted_local_updates;
    out.block_preparations=search.block_preparations;
    out.domain_preparations=search.domain_preparations;
    out.mapping_preparations=search.mapping_preparations;
    out.prepared_mapping_entries=search.prepared_mapping_entries;
    out.prepared_mapping_bytes=search.prepared_mapping_bytes;
    out.memory_estimate=search.memory_estimate;
    out.global_state_memory_estimate=search.global_state_memory_estimate;
    out.work=search.work;
    out.search_seconds=search.search_seconds;
    out.final_sweep=search.final_sweep;
    out.local_conditioning=std::move(search.local_conditioning);
    out.certification_attempted=search.certification.has_value();
    if(policy.certify_endpoint)
    {
        if(search.certification)
        {
            out.work=search.work;
            AssembleResult(problem,data,layout,domain,y,context,*search.certification,
                policy.convergence_criteria,out,search.accepted_state_integrity);
        }
        else if(!policy.stop_after_stationarity && search.search_converged &&
            search.state.eta.size()==static_cast<Eigen::Index>(domain.atoms.size()))
        {
            const auto certification=CertifyFixedNeighborEndpoint(domain,y,search.state,context,
                0,
                policy.convergence_criteria);
            out.work=search.work;
            AssembleResult(problem,data,layout,domain,y,context,certification,
                policy.convergence_criteria,out,search.accepted_state_integrity);
        }
    }
    return out;
}
}
