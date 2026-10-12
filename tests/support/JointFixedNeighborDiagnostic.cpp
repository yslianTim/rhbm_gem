#include "support/JointFixedNeighborDiagnostic.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

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
        Values(endpoint.gradient),endpoint.certificate.objective/(scale*scale),{}};
}
double ScaledCoefficientDifference(VectorRef lhs,VectorRef rhs)
{
    if(lhs.size()!=rhs.size() || lhs.size()==0) return std::numeric_limits<double>::infinity();
    return ((lhs-rhs).array().abs()/(1.0+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
}
void AssembleResult(const JointProblem & problem,const ProblemData & data,const JointParameterLayout & layout,
    const Domain & domain,VectorRef y,const EvaluationContext & context,
    const FixedNeighborEndpointCertification & certification,FixedNeighborDiagnosticResult & out)
{
    out.assessment=certification.assessment;
    out.endpoint_trust=certification.trust;
    out.endpoint_certified=certification.passed();
    const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
    const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);

    const auto & component_view=data.partition.components.front();
    JointComponentResult component; component.id=component_view.id;
    component.atoms.assign(component_view.atoms.begin(),component_view.atoms.end());
    component.rows.assign(component_view.rows.begin(),component_view.rows.end());
    component.layout=layout; component.search_completed=out.search_converged; component.stop_reason=out.reason;
    component.search_sweeps=out.sweeps.size();
    component.first_order_stationarity_sweep=out.first_order_stationarity_sweep;
    component.confirmed_stationarity_sweep=out.confirmed_stationarity_sweep;
    component.endpoint_certification_sweep=out.sweeps.size();
    if(certification.passed()) component.state=State(out.assessment.primary,context.scale);
    component.evidence=AssessmentEvidence(out.assessment,JointEvidenceScope::ComponentLocal);
    component.ranks=AssessmentRanks(out.assessment,JointEvidenceScope::ComponentLocal);
    component.evidence.push_back({"endpoint-trust",certification.passed() ? JointCheckStatus::Passed :
        JointCheckStatus::Failed,JointEvidenceScope::ComponentLocal,{},{},certification.trust.reason});

    out.fit.problem=problem; out.fit.layout=layout; out.fit.initialization.valid=true;
    out.fit.initialization.reason="valid-widths"; out.fit.initialization.b=Values(out.state.eta.array().exp());
    out.fit.components.push_back(std::move(component)); out.fit.search_completed=out.search_converged;
    out.fit.observation_scale=context.scale; out.fit.available_row_mask.assign(data.input->observations.size(),true);
    if(certification.passed())
    {
        out.fit.assembled_state=out.fit.components.front().state;
        out.fit.prediction=Values(out.state.prediction); out.fit.objective=out.state.objective;
    }
    out.fit.evidence=AssessmentEvidence(out.assessment,JointEvidenceScope::AssembledGlobal);
    out.fit.ranks=AssessmentRanks(out.assessment,JointEvidenceScope::AssembledGlobal);

    const auto control=EvaluateProfile(domain,y,eta,false,&context);
    bool agrees=control.valid; double coefficient_difference=std::numeric_limits<double>::infinity();
    Vector expected=Vector::Zero(y.size());
    if(agrees)
    {
        expected=control.x*beta;
        coefficient_difference=ScaledCoefficientDifference(beta,control.beta);
        agrees=coefficient_difference<=1e-10;
    }
    const Vector actual=Select(out.state.prediction,IndicesOf(layout.informative_rows));
    const double prediction_difference=((actual-expected).array().abs()/(1+actual.array().abs().max(expected.array().abs()))).maxCoeff();
    const double profile_objective=.5*(expected-y).squaredNorm()/(context.scale*context.scale);
    const bool reconstruction=actual.allFinite() && expected.allFinite() && prediction_difference<=1e-10 &&
        std::abs(out.state.objective-profile_objective)<=1e-12;
    agrees &= reconstruction;
    out.fit.evidence.push_back({"assembled-profile",agrees ? JointCheckStatus::Passed : JointCheckStatus::Failed,
        JointEvidenceScope::AssembledGlobal,std::max(coefficient_difference,prediction_difference),1e-10,
        reconstruction ? "" : "full-domain-reconstruction-failed"});
}
}

FixedNeighborDiagnosticResult RunFixedNeighborDiagnostic(const JointProblem & problem,VectorRef initial_eta,
    const FixedNeighborPolicy & policy,Vector initial_beta)
{
    FixedNeighborDiagnosticResult out;
    const auto & data=JointProblemAccess::Get(problem); const auto & input=*data.input; const auto & layout=data.layout;
    if(initial_eta.size()!=static_cast<Eigen::Index>(input.atom_ids.size()) || !initial_eta.allFinite() ||
        !(data.context.scale>0) || !std::isfinite(data.context.scale) || policy.core_atoms==0 ||
        layout.groups.size()!=0 || layout.full_atoms.size()!=input.atom_ids.size() ||
        layout.informative_rows.size()!=input.observations.size() || data.partition.components.size()!=1)
    {out.reason="block-invalid-partition"; return out;}
    for(Eigen::Index a=0;a<initial_eta.size();++a)
        if(!(std::exp(initial_eta(a))>0) || !std::isfinite(std::exp(initial_eta(a))))
        {out.reason="block-invalid-partition"; return out;}

    const Domain domain=ProfileDomain(data.domain,layout);
    const EvaluationContext context=ProfileContext(data.context,layout,data.domain.rows);
    const Vector y=Select(data.y,IndicesOf(layout.informative_rows));
    auto run_policy=policy; run_policy.collect_diagnostics=true;
    auto search=SearchFixedNeighborComponent(domain,layout,y,initial_eta,context,
        run_policy,std::move(initial_beta));
    out.state=search.state; out.blocks=search.blocks; out.sweeps=search.sweeps;
    out.search_converged=search.search_converged; out.reason=search.reason;
    out.first_order_stationarity_sweep=search.first_order_stationarity_sweep;
    out.confirmed_stationarity_sweep=search.confirmed_stationarity_sweep;
    out.block_preparations=search.block_preparations;
    out.domain_preparations=search.domain_preparations;
    out.mapping_preparations=search.mapping_preparations;
    out.prepared_mapping_entries=search.prepared_mapping_entries;
    out.prepared_mapping_bytes=search.prepared_mapping_bytes;
    out.memory_estimate=search.memory_estimate;
    out.work=search.work;
    if(policy.certify_endpoint)
    {
        if(search.certification)
        {
            out.work=search.work;
            AssembleResult(problem,data,layout,domain,y,context,*search.certification,out);
        }
        else if(!policy.stop_after_stationarity && search.search_converged &&
            search.state.eta.size()==static_cast<Eigen::Index>(domain.atoms.size()))
        {
            const auto certification=CertifyFixedNeighborEndpoint(domain,y,search.state,context,
                policy.collect_telemetry ? &search.work.local_profile_work : nullptr);
            out.work=search.work;
            AssembleResult(problem,data,layout,domain,y,context,certification,out);
        }
    }
    return out;
}
}
