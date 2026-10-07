#include "FixedNeighborBlockCoordinate.hpp"
#include "Problem.hpp"
#include "SparseFactor.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
enum class StationarityState {NotStationary,CandidateStationary,ConfirmedStationary};
constexpr double EtaChangeConfirmationThreshold=1e-10;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
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
struct ReplayState {Vector prediction,residual; double objective{};};
ReplayState Replay(const JointProblemInput & input,const JointParameterLayout & layout,VectorRef observations,
    VectorRef eta,VectorRef beta,double scale)
{
    ReplayState out; out.prediction=Vector::Zero(observations.size()); out.residual=Vector::Zero(observations.size());
    std::vector<std::uint8_t> informative(input.observations.size());
    for(auto row:layout.informative_rows) informative.at(row)=1;
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); const double width=std::exp(eta(a));
        for(const auto & support:input.support.at(atom)) if(informative.at(support.row))
        {
            const auto basis=EvaluateKernel(support.squared_distance,width,2.5);
            out.prediction(static_cast<Eigen::Index>(support.row))+=beta(2*a)*basis.gaussian+beta(2*a+1)*basis.charge;
        }
    }
    double squared{};
    for(auto row:layout.informative_rows)
    {
        const auto r=static_cast<Eigen::Index>(row); out.residual(r)=out.prediction(r)-observations(r);
        squared+=out.residual(r)*out.residual(r);
    }
    out.objective=.5*squared/(scale*scale); return out;
}
Sparse Design(const Domain & domain,VectorRef eta)
{
    Sparse out(domain.rows,2*eta.size()); std::vector<Eigen::Triplet<double>> entries;
    for(Eigen::Index a=0;a<eta.size();++a)
    {
        const double width=std::exp(eta(a));
        for(const auto & support:domain.atoms[static_cast<std::size_t>(a)])
        {
            const auto basis=EvaluateKernel(support.square,width,2.5);
            if(basis.gaussian!=0) entries.emplace_back(support.row,2*a,basis.gaussian);
            if(basis.charge!=0) entries.emplace_back(support.row,2*a+1,basis.charge);
        }
    }
    out.setFromTriplets(entries.begin(),entries.end()); return out;
}
double PredictionDifference(VectorRef a,VectorRef b)
{return (a-b).lpNorm<Eigen::Infinity>();}
double ScaledCoefficientDifference(VectorRef lhs,VectorRef rhs)
{
    if(lhs.size()!=rhs.size() || lhs.size()==0) return std::numeric_limits<double>::infinity();
    return ((lhs-rhs).array().abs()/(1.0+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
}
JointState State(const Endpoint & endpoint,double scale)
{
    return {Values(endpoint.beta),Values(endpoint.eta.array().exp()),Values(endpoint.eta),
        Values(endpoint.gradient),endpoint.certificate.objective/(scale*scale),{}};
}
void BuildEndpointResult(const JointProblem & problem,const ProblemData & data,const JointParameterLayout & layout,
    const Domain & domain,VectorRef y,const EvaluationContext & context,
    const FixedNeighborSearchResult & search,FixedNeighborResult & out)
{
    out.assessment=search.assessment;
    out.endpoint_trust=search.endpoint_trust;
    out.endpoint_certified=search.endpoint_certified;
    const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
    const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);

    const auto & component_view=data.partition.components.front();
    JointComponentResult component; component.id=component_view.id;
    component.atoms.assign(component_view.atoms.begin(),component_view.atoms.end());
    component.rows.assign(component_view.rows.begin(),component_view.rows.end());
    component.layout=layout; component.search_completed=out.search_converged; component.stop_reason=out.reason;
    component.state=State(out.assessment.primary,context.scale);
    component.evidence=AssessmentEvidence(out.assessment,JointEvidenceScope::ComponentLocal);
    component.ranks=AssessmentRanks(out.assessment,JointEvidenceScope::ComponentLocal);

    out.fit.problem=problem; out.fit.layout=layout; out.fit.initialization.valid=true;
    out.fit.initialization.reason="valid-widths"; out.fit.initialization.b=Values(out.state.eta.array().exp());
    out.fit.components.push_back(std::move(component)); out.fit.search_completed=out.search_converged;
    out.fit.observation_scale=context.scale; out.fit.available_row_mask.assign(data.input->observations.size(),true);
    out.fit.assembled_state=out.fit.components.front().state;
    out.fit.prediction=Values(out.state.prediction); out.fit.objective=out.state.objective;
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
bool IsCertifiedLocalEndpoint(const Assessment & assessment,const TrustEvidence & trust)
{return assessment.inner && assessment.gradient && assessment.local && assessment.identified && trust.passed;}
bool IsFixedNeighborEtaChangeConfirmed(double eta_change_inf,bool has_previous_complete_sweep)
{
    return has_previous_complete_sweep && std::isfinite(eta_change_inf) && eta_change_inf>=0.0 &&
        eta_change_inf<=EtaChangeConfirmationThreshold;
}
int LocalUpdateBudget(FixedNeighborLocalWork work)
{
    switch(work)
    {
    case FixedNeighborLocalWork::Full: return 0;
    case FixedNeighborLocalWork::OneAcceptedUpdate: return 1;
    case FixedNeighborLocalWork::TwoAcceptedUpdates: return 2;
    }
    return 0;
}
FixedNeighborSearchResult SearchFixedNeighborComponent(
    const JointProblemInput & input,const JointParameterLayout & layout,const Domain & domain,
    VectorRef observations,VectorRef y,VectorRef initial_eta,const EvaluationContext & context,
    const FixedNeighborPolicy & policy)
{
    FixedNeighborSearchResult out;
    if(initial_eta.size()!=static_cast<Eigen::Index>(input.atom_ids.size()) || !initial_eta.allFinite() ||
        observations.size()!=static_cast<Eigen::Index>(input.observations.size()) || !observations.allFinite() ||
        y.size()!=domain.rows || !y.allFinite() || !(context.scale>0) || !std::isfinite(context.scale) ||
        policy.core_atoms==0 || layout.groups.size()!=0 || layout.full_atoms.size()!=input.atom_ids.size() ||
        layout.informative_rows.size()!=input.observations.size() || domain.atoms.size()!=input.atom_ids.size())
    {out.reason="block-invalid-partition"; return out;}
    for(Eigen::Index a=0;a<initial_eta.size();++a)
        if(!(std::exp(initial_eta(a))>0) || !std::isfinite(std::exp(initial_eta(a))))
        {out.reason="block-invalid-partition"; return out;}

    StructuralBlockPartition partition;
    try {partition=BuildStructuralBlockPartition(input,layout,policy.core_atoms);}
    catch(const std::exception &) {out.reason="block-invalid-partition"; return out;}
    if(partition.cores.empty()) {out.reason="block-invalid-partition"; return out;}

    out.state.eta=initial_eta; out.state.beta=Vector::Zero(2*static_cast<Eigen::Index>(input.atom_ids.size()));
    const auto initial=Replay(input,layout,observations,out.state.eta,out.state.beta,context.scale);
    out.state.prediction=initial.prediction; out.state.residual=initial.residual; out.state.objective=initial.objective;
    Indices atom_position(input.atom_ids.size(),-1),row_position(input.observations.size(),-1);
    for(std::size_t k=0;k<layout.full_atoms.size();++k) atom_position.at(layout.full_atoms[k])=static_cast<Eigen::Index>(k);
    for(std::size_t k=0;k<layout.informative_rows.size();++k) row_position.at(layout.informative_rows[k])=static_cast<Eigen::Index>(k);
    std::vector<PreparedFixedNeighborBlock> prepared_blocks;
    prepared_blocks.reserve(partition.cores.size());
    const int local_update_budget=LocalUpdateBudget(policy.local_work);
    for(std::size_t block_index=0;block_index<partition.cores.size();++block_index)
    {
        const auto & core=partition.cores[block_index]; PreparedFixedNeighborBlock prepared;
        prepared.block=block_index; prepared.core_atoms=core.atoms; prepared.affected_rows=core.affected_rows;
        prepared.profile_atoms.reserve(core.atoms.size()); prepared.profile_rows.reserve(core.affected_rows.size());
        for(auto atom:core.atoms)
        {
            const auto position=atom_position.at(static_cast<std::size_t>(atom));
            if(position<0) {out.reason="block-invalid-partition"; return out;}
            prepared.profile_atoms.push_back(position);
        }
        for(auto row:core.affected_rows)
        {
            const auto position=row_position.at(static_cast<std::size_t>(row));
            if(position<0) {out.reason="block-invalid-partition"; return out;}
            prepared.profile_rows.push_back(position);
        }
        std::vector<Eigen::Index> row_mapping(static_cast<std::size_t>(domain.rows),-1);
        for(std::size_t k=0;k<prepared.profile_rows.size();++k)
            row_mapping.at(static_cast<std::size_t>(prepared.profile_rows[k]))=static_cast<Eigen::Index>(k);
        prepared.row_mapping=std::make_shared<Indices>(std::move(row_mapping));
        prepared.domain=domain.Select(prepared.profile_atoms,
            static_cast<Eigen::Index>(prepared.profile_rows.size()),prepared.row_mapping);
        prepared.workspace=std::make_shared<LinearWorkspace>();
        prepared.context=context;
        prepared.context.atom_ids=context.atom_ids.Select(prepared.profile_atoms);
        prepared.context.row_ids=context.row_ids.Select(prepared.profile_rows);
        prepared.context.independent_search=true;
        prepared.local_atoms=static_cast<Eigen::Index>(prepared.profile_atoms.size());
        prepared.context.rank={prepared.domain.rows,2*prepared.local_atoms,prepared.local_atoms};
        prepared.context.linear.rank_relative=prepared.context.rank.Relative(2*prepared.local_atoms);
        prepared.context.search=SearchPolicy{};
        prepared.context.search.method=SearchMethod::LegacyCompact;
        if(local_update_budget>0) prepared.context.update_budget=local_update_budget;
        prepared_blocks.push_back(std::move(prepared));
    }
    out.block_preparations=prepared_blocks.size();
    out.domain_preparations=prepared_blocks.size();
    out.mapping_preparations=prepared_blocks.size();
    std::vector<std::size_t> order(partition.cores.size()); std::iota(order.begin(),order.end(),0);
    if(policy.order==FixedNeighborBlockOrder::Reverse) std::reverse(order.begin(),order.end());
    const bool retain_diagnostics=policy.collect_diagnostics || policy.capture_local_trajectory ||
        policy.certify_local_candidates || static_cast<bool>(policy.state_observer);
    bool stop=false;
    for(std::size_t sweep_index=0;sweep_index<policy.maximum_sweeps && !stop;++sweep_index)
    {
        const auto sweep_started=Clock::now(); FixedNeighborBlockSweep sweep; sweep.sweep=sweep_index+1;
        const Vector eta_before=out.state.eta,beta_before=out.state.beta;
        sweep.objective_before=out.state.objective;
        auto finish_record=[&](FixedNeighborBlockRecord record) {
            if(record.accepted)
            {
                ++sweep.accepted_blocks;
                sweep.accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
                out.total_accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
            }
            else if(record.status=="unchanged") ++sweep.unchanged_blocks;
            if(policy.diagnostic_sink && policy.diagnostic_sink->block) policy.diagnostic_sink->block(record);
            if(retain_diagnostics) out.blocks.push_back(std::move(record));
        };
        for(auto block_index:order)
        {
            const auto & prepared=prepared_blocks.at(block_index);
            const auto & core_atoms=prepared.core_atoms; const auto & affected_rows=prepared.affected_rows;
            const auto & local_domain=prepared.domain; const auto & local_context=prepared.context;
            const auto local_atoms=prepared.local_atoms; FixedNeighborBlockRecord record;
            record.sweep=sweep_index+1; record.block=block_index+1; record.core_atoms=core_atoms;
            record.affected_rows=affected_rows.size(); record.objective_before=out.state.objective;

            Vector old_eta(local_atoms),old_beta(2*local_atoms);
            for(std::size_t k=0;k<core_atoms.size();++k)
            {
                const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                old_eta(static_cast<Eigen::Index>(k))=out.state.eta(atom);
                old_beta.segment<2>(2*static_cast<Eigen::Index>(k))=out.state.beta.segment<2>(2*atom);
            }
            const Sparse old_design=Design(local_domain,old_eta);
            const Vector old_core=(old_design*old_beta).eval();
            Vector local_y(static_cast<Eigen::Index>(affected_rows.size()));
            double local_before_squared{};
            for(std::size_t k=0;k<affected_rows.size();++k)
            {
                const auto row=static_cast<Eigen::Index>(affected_rows[k]);
                local_y(static_cast<Eigen::Index>(k))=old_core(static_cast<Eigen::Index>(k))-out.state.residual(row);
                local_before_squared+=out.state.residual(row)*out.state.residual(row);
            }
            record.local_objective_before=.5*local_before_squared/(context.scale*context.scale);
            const auto search_started=Clock::now();
            const auto old_widths=old_eta.array().exp().eval();
            const auto local_search=SearchProfile(local_domain,local_y,old_widths,local_context,{},nullptr,
                prepared.workspace.get(),&prepared);
            record.search_seconds=Seconds(search_started); record.profile_evaluations=local_search.evaluations;
            record.accepted_updates=local_search.accepted; record.local_search_stop_reason=local_search.stop_reason;
            if(policy.capture_local_trajectory)
            {
                record.profile_trials.reserve(local_search.trials.size());
                double accepted_objective=record.local_objective_before;
                Vector accepted_eta=old_eta; double cumulative_factor_seconds{};
                for(std::size_t trial_index=0;trial_index<local_search.trials.size();++trial_index)
                {
                    const auto & trial=local_search.trials[trial_index];
                    const double objective=trial.endpoint.certificate.evaluated && trial.endpoint.certificate.available ?
                        trial.endpoint.certificate.objective/(context.scale*context.scale) : unavailable;
                    FixedNeighborProfileTrial telemetry;
                    telemetry.trial_index=trial_index; telemetry.profile_evaluation=trial.evaluation;
                    telemetry.accepted=trial.accepted; telemetry.accepted_update=trial.accepted_update;
                    telemetry.local_objective_before=accepted_objective;
                    telemetry.local_objective_after=objective;
                    telemetry.objective_reduction=std::isfinite(accepted_objective) && std::isfinite(objective) ?
                        accepted_objective-objective : unavailable;
                    telemetry.eta_change_inf=trial.endpoint.eta.size()==accepted_eta.size() ?
                        (trial.endpoint.eta-accepted_eta).lpNorm<Eigen::Infinity>() : unavailable;
                    telemetry.gradient_inf_norm=trial.endpoint.gradient.size() ?
                        trial.endpoint.gradient.lpNorm<Eigen::Infinity>() : unavailable;
                    telemetry.profile_seconds=trial.seconds; telemetry.factor_seconds=trial.factor_seconds;
                    cumulative_factor_seconds+=trial.factor_seconds;
                    telemetry.cumulative_factor_seconds=cumulative_factor_seconds;
                    record.profile_trials.push_back(std::move(telemetry));
                    if(trial.accepted && std::isfinite(objective))
                    {
                        accepted_objective=objective;
                        if(trial.endpoint.eta.size()==accepted_eta.size()) accepted_eta=trial.endpoint.eta;
                    }
                }
                record.profile_factor_seconds=cumulative_factor_seconds;
            }
            ++sweep.block_solves; ++out.total_block_solves;
            sweep.profile_evaluations+=static_cast<std::size_t>(local_search.evaluations);
            out.total_profile_evaluations+=static_cast<std::size_t>(local_search.evaluations);
            sweep.maximum_block_rows=std::max(sweep.maximum_block_rows,static_cast<std::size_t>(local_domain.rows));
            sweep.maximum_block_columns=std::max(sweep.maximum_block_columns,static_cast<std::size_t>(2*local_atoms));
            const Trial * accepted=nullptr;
            for(auto it=local_search.trials.rbegin();it!=local_search.trials.rend();++it)
                if(it->accepted && it->endpoint.valid && it->trust && it->trust->passed) {accepted=&*it; break;}
            if(!accepted)
            {
                record.status="failed"; record.reason="block-search-failed"; out.reason=record.reason;
                finish_record(std::move(record)); stop=true; break;
            }
            const auto local_state=EvaluateState(local_domain,local_y,accepted->endpoint.eta,accepted->endpoint.beta,local_context);
            if(!local_state.valid || !local_state.certificate.available || !local_state.certificate.feasible ||
                !local_state.certificate.kkt_passed)
            {
                record.status="failed"; record.reason="block-inner-invalid"; out.reason=record.reason;
                finish_record(std::move(record)); stop=true; break;
            }
            record.local_final_gradient_inf_norm=local_state.gradient.lpNorm<Eigen::Infinity>();
            record.local_final_ac_kkt=local_state.certificate.projected_kkt;
            if(policy.certify_local_candidates)
            {
                record.local_assessment_attempted=true;
                record.local_assessment_rows=local_state.x.rows();
                record.local_assessment_columns=local_state.x.cols();
                const auto assessment_started=Clock::now();
                const auto local_reference=EvaluateProfile(local_domain,local_y,local_state.eta,true,&local_context);
                const auto local_assessment=AssessEvaluated(local_domain,local_y,local_state,local_reference,local_context,true);
                const auto local_trust=CheckTrust(local_domain,local_y,local_state,local_context,local_reference);
                record.local_assessment_seconds=Seconds(assessment_started);
                sweep.local_assessment_seconds+=record.local_assessment_seconds;
                ++sweep.local_assessments; ++out.total_local_assessments;
                sweep.maximum_local_assessment_rows=std::max(sweep.maximum_local_assessment_rows,
                    static_cast<std::size_t>(record.local_assessment_rows));
                sweep.maximum_local_assessment_columns=std::max(sweep.maximum_local_assessment_columns,
                    static_cast<std::size_t>(record.local_assessment_columns));
                record.local_inner_passed=local_assessment.inner;
                record.local_gradient_passed=local_assessment.gradient;
                record.local_correction_passed=local_assessment.local;
                record.local_identified=local_assessment.identified;
                record.local_trust_passed=local_trust.passed;
                record.local_profile_gradient_inf_norm=local_assessment.primary.gradient.size() ?
                    local_assessment.primary.gradient.lpNorm<Eigen::Infinity>() : std::numeric_limits<double>::infinity();
                record.local_reference_gradient_inf_norm=local_assessment.reference.gradient.size() ?
                    local_assessment.reference.gradient.lpNorm<Eigen::Infinity>() : std::numeric_limits<double>::infinity();
                record.local_correction_inf_norm=local_assessment.correction.size() ?
                    local_assessment.correction.lpNorm<Eigen::Infinity>() : std::numeric_limits<double>::infinity();
                if(local_assessment.widths)
                {
                    record.local_projected_width_rank=local_assessment.widths->rank;
                    record.local_projected_width_minimum=local_assessment.widths->minimum;
                }
                if(local_assessment.jacobian)
                {
                    record.local_corrected_jacobian_rank=local_assessment.jacobian->rank;
                    record.local_corrected_jacobian_minimum=local_assessment.jacobian->minimum;
                }
                if(local_assessment.normalized_widths)
                {
                    record.local_normalized_width_rank=local_assessment.normalized_widths->rank;
                    record.local_normalized_width_minimum=local_assessment.normalized_widths->minimum;
                }
                record.local_assessment_failure=local_assessment.failure;
                record.local_trust_reason=local_trust.reason;
                record.local_assessment_passed=IsCertifiedLocalEndpoint(local_assessment,local_trust);
                if(!record.local_assessment_passed)
                {
                    record.status="unchanged"; record.reason="local-endpoint-uncertified";
                    record.local_objective_after=record.local_objective_before;
                    record.objective_after=record.objective_before;
                    finish_record(std::move(record)); continue;
                }
                ++sweep.certified_local_candidates;
            }
            record.local_objective_after=local_state.certificate.objective/(context.scale*context.scale);
            if(record.local_objective_after>record.local_objective_before)
            {
                record.status="unchanged"; record.reason="block-objective-increase";
                record.objective_after=record.objective_before; finish_record(std::move(record)); continue;
            }
            Vector candidate_eta=out.state.eta,candidate_beta=out.state.beta;
            for(std::size_t k=0;k<core_atoms.size();++k)
            {
                const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                candidate_eta(atom)=accepted->endpoint.eta(static_cast<Eigen::Index>(k));
                candidate_beta.segment<2>(2*atom)=accepted->endpoint.beta.segment<2>(2*static_cast<Eigen::Index>(k));
            }
            const auto replay=Replay(input,layout,observations,candidate_eta,candidate_beta,context.scale);
            record.objective_after=replay.objective; record.global_replay_delta=replay.objective-out.state.objective;
            const auto replay_reference=std::max(std::abs(replay.objective),std::abs(out.state.objective));
            record.objective_replay_enclosure=BlockObjectiveReplayEnclosure(replay_reference);
            record.local_global_delta_error=std::abs(record.global_replay_delta+
                record.local_objective_before-record.local_objective_after);
            const bool delta_enclosed=WithinBlockObjectiveReplay(record.local_global_delta_error,replay_reference);
            if(!delta_enclosed)
            {
                record.status="failed"; record.reason="block-objective-replay-failed"; out.reason=record.reason;
                finish_record(std::move(record)); stop=true; break;
            }
            if(record.global_replay_delta>0 && !(record.local_objective_after<=record.local_objective_before &&
                record.global_replay_delta<=record.objective_replay_enclosure))
            {
                record.status="unchanged"; record.reason="block-objective-increase";
                finish_record(std::move(record)); continue;
            }
            const bool unchanged=(out.state.eta.array()==candidate_eta.array()).all() &&
                (out.state.beta.array()==candidate_beta.array()).all();
            const Vector new_core=(local_state.x*local_state.beta).eval(); double old_squared{},new_squared{};
            for(std::size_t k=0;k<affected_rows.size();++k)
            {
                const auto row=static_cast<Eigen::Index>(affected_rows[k]); const double before=out.state.residual(row);
                const double change=new_core(static_cast<Eigen::Index>(k))-old_core(static_cast<Eigen::Index>(k));
                out.state.prediction(row)+=change; out.state.residual(row)+=change;
                old_squared+=before*before; new_squared+=out.state.residual(row)*out.state.residual(row);
            }
            out.state.objective+=(new_squared-old_squared)/(2*context.scale*context.scale);
            out.state.eta=std::move(candidate_eta); out.state.beta=std::move(candidate_beta);
            record.accepted=!unchanged; record.status=unchanged ? "unchanged" : "accepted";
            record.reason=unchanged ? "block-unchanged" : "";
            finish_record(std::move(record));
        }
        if(stop) break;
        sweep.objective_after=out.state.objective;
        const auto replay=Replay(input,layout,observations,out.state.eta,out.state.beta,context.scale);
        sweep.cache_replay_error=PredictionDifference(out.state.prediction,replay.prediction);
        sweep.objective_replay_error=std::abs(out.state.objective-replay.objective);
        const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
        const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);
        const auto global=EvaluateState(domain,y,eta,beta,context);
        if(!global.valid || !global.certificate.available)
        {out.reason="block-inner-invalid"; stop=true; break;}
        sweep.global_a_feasibility=global.certificate.feasible ? 0.0 : 1.0;
        sweep.global_ac_kkt=global.certificate.projected_kkt;
        sweep.global_width_gradient_inf_norm=global.gradient.lpNorm<Eigen::Infinity>();
        sweep.eta_change_inf=(out.state.eta-eta_before).lpNorm<Eigen::Infinity>();
        sweep.beta_scaled_change=ScaledCoefficientDifference(out.state.beta,beta_before);
        sweep.coordinate_confirmation_available=sweep_index>0;
        sweep.wall_seconds=Seconds(sweep_started); out.search_seconds+=sweep.wall_seconds;
        ++out.sweep_count; out.final_sweep=sweep;
        if(policy.diagnostic_sink && policy.diagnostic_sink->sweep) policy.diagnostic_sink->sweep(sweep);
        if(policy.production_progress) policy.production_progress(sweep);
        if(retain_diagnostics) out.sweeps.push_back(sweep);
        if(policy.sweep_observer)
        {
            const auto & observed=retain_diagnostics ? out.sweeps.back() : sweep;
            policy.sweep_observer(observed);
        }
        if(policy.state_observer) policy.state_observer(sweep_index+1,out.state,out.sweeps.back(),out.blocks);
        if(sweep.cache_replay_error>2e-12+2e-13*std::max(1.0,replay.prediction.cwiseAbs().maxCoeff()) ||
            !WithinBlockObjectiveReplay(sweep.objective_replay_error,replay.objective))
        {out.reason="block-cache-replay-failed"; stop=true; break;}
        if(policy.stop_after_no_certified_update && sweep.local_assessments>0 && sweep.certified_local_candidates==0)
        {out.reason="local-endpoint-uncertified"; break;}
        StationarityState stationarity=StationarityState::NotStationary;
        if(sweep.global_ac_kkt<=1e-10 && sweep.global_width_gradient_inf_norm<=1e-12)
        {
            if(out.first_order_stationarity_sweep==0) out.first_order_stationarity_sweep=sweep_index+1;
            stationarity=StationarityState::CandidateStationary;
            if(IsFixedNeighborEtaChangeConfirmed(sweep.eta_change_inf,sweep.coordinate_confirmation_available))
                stationarity=StationarityState::ConfirmedStationary;
        }
        if(stationarity==StationarityState::ConfirmedStationary)
        {
            if(out.confirmed_stationarity_sweep==0) out.confirmed_stationarity_sweep=sweep_index+1;
            out.search_converged=true; out.reason="block-stationary";
            if(policy.stop_after_stationarity) break;
        }
    }
    if(!out.search_converged && out.reason.empty()) out.reason="block-sweep-budget";
    if(policy.assess_final_endpoint)
    {
        const auto endpoint=EvaluateState(domain,y,out.state.eta,out.state.beta,context);
        const auto reference=EvaluateProfile(domain,y,out.state.eta,true,&context);
        out.assessment=AssessEvaluated(domain,y,endpoint,reference,context,true);
        out.endpoint_trust=CheckTrust(domain,y,endpoint,context,reference);
        out.endpoint_certified=out.endpoint_trust.passed;
    }
    return out;
}

namespace {
std::optional<double> LastFixedNeighborObjective(const FixedNeighborSearchResult & result)
{
    if(result.final_sweep) return result.final_sweep->objective_after;
    return std::nullopt;
}

std::optional<double> LastFixedNeighborGradient(const FixedNeighborSearchResult & result)
{
    if(result.final_sweep) return result.final_sweep->global_width_gradient_inf_norm;
    return std::nullopt;
}
}

ComponentResult SolveFixedNeighborComponent(
    const PreparedComponent & prepared,VectorRef initial_b,const EvaluationContext & parent_context,
    const FixedNeighborSearchPolicy & production_policy,const JointProgressObserver & observer,
    const JointProgressComponent * progress_component)
{
    ComponentResult out;
    const auto & component_layout=prepared.parent_layout; const auto & view=prepared.view;
    if(component_layout.full_atoms.empty() || component_layout.informative_rows.empty())
    {
        out.search.stop_reason=component_layout.full_atoms.empty() ? "analytic-nuisance-only" : "unobserved-full-parameters";
        out.search.stopped=true;
        return out;
    }
    for(const auto atom:component_layout.full_atoms)
        if(atom>=static_cast<std::size_t>(initial_b.size()) || !std::isfinite(initial_b(static_cast<Eigen::Index>(atom))) ||
            !(initial_b(static_cast<Eigen::Index>(atom))>0))
        {
            out.search.stop_reason="invalid-initial-widths"; out.search.stopped=true; return out;
        }

    const auto & local_input=prepared.input; const auto & local_layout=prepared.local_layout;
    const auto & local_domain=prepared.domain;
    const VectorMap local_observations(local_input->observations.data(),
        static_cast<Eigen::Index>(local_input->observations.size()));
    const Vector local_y=local_observations;
    auto context=parent_context;
    context.atom_ids=prepared.atom_ids; context.row_ids=prepared.row_ids;
    context.independent_search=true;
    const auto atom_count=static_cast<Eigen::Index>(component_layout.full_atoms.size());
    context.rank={static_cast<Eigen::Index>(view.rows.size()),2*atom_count,atom_count};
    context.linear.rank_relative=context.rank.Relative(2*atom_count);

    Vector initial_eta(atom_count);
    for(std::size_t k=0;k<prepared.local_to_parent_atoms.size();++k)
        initial_eta(static_cast<Eigen::Index>(k))=std::log(initial_b(prepared.local_to_parent_atoms[k]));

    FixedNeighborPolicy policy;
    policy.core_atoms=production_policy.core_atoms;
    policy.maximum_sweeps=production_policy.maximum_sweeps;
    policy.order=production_policy.order;
    policy.local_work=production_policy.local_work;
    policy.assess_final_endpoint=true;
    if(progress_component && observer)
    {
        policy.production_progress=[&observer,progress_component,&production_policy,context](
            const FixedNeighborBlockSweep & sweep) {
            auto event=MakeJointProgressEvent(JointProgressPhase::SearchProgress,*progress_component);
            event.profile_evaluations=static_cast<int>(sweep.profile_evaluations);
            event.profile_budget=context.profile_budget;
            event.accepted_updates=static_cast<int>(sweep.accepted_blocks);
            event.update_budget=static_cast<int>(sweep.block_solves);
            event.elapsed_seconds=sweep.wall_seconds;
            event.accepted_objective=sweep.objective_after;
            event.accepted_gradient_inf_norm=sweep.global_width_gradient_inf_norm;
            event.fixed_neighbor=JointFixedNeighborProgress{
                sweep.sweep,production_policy.maximum_sweeps,sweep.block_solves,sweep.accepted_blocks,
                sweep.accepted_local_updates,
                sweep.objective_after,sweep.global_ac_kkt,sweep.global_width_gradient_inf_norm,sweep.eta_change_inf};
            observer(event);
        };
    }
    const auto result=SearchFixedNeighborComponent(*local_input,local_layout,local_domain,
        local_observations,local_y,initial_eta,context,policy);
    out.search.initial.eta=initial_eta; out.search.initial.beta=Vector::Zero(2*atom_count);
    out.search.initial.valid=true; out.search.initial.certificate.evaluated=true;
    out.search.eta=result.state.eta; out.search.stopped=!result.search_converged;
    out.search.stop_reason=result.reason; out.search.lm_status=0;
    out.search.accepted_objective=LastFixedNeighborObjective(result);
    out.search.accepted_gradient_inf_norm=LastFixedNeighborGradient(result);
    out.search.evaluations=static_cast<int>(result.total_profile_evaluations);
    out.search.seconds=result.search_seconds;
    out.search.accepted=static_cast<int>(result.total_accepted_local_updates);
    out.search.references=static_cast<int>(result.total_local_assessments);
    out.assessment=result.assessment;
    if(result.search_converged && !result.endpoint_certified)
    {
        out.search.stop_reason="endpoint-certification-failed";
        out.search.stopped=true;
    }
    if(result.endpoint_certified)
    {
        out.endpoint_trust=result.endpoint_trust;
        out.trusted_state=result.assessment.primary;
        out.trusted_assessment=result.assessment;
    }
    out.search_success=result.search_converged && result.endpoint_certified;
    if(progress_component)
        NotifyJointProgress(observer,JointProgressPhase::CertificationStarted,*progress_component,
            out.search.evaluations,parent_context.profile_budget,out.search.accepted,parent_context.update_budget,
            out.search.seconds,out.search.stop_reason,out.trusted_state.has_value(),
            out.search.accepted_objective,out.search.accepted_gradient_inf_norm);
    return out;
}

FixedNeighborResult SearchFixedNeighbor(const JointProblem & problem,VectorRef initial_eta,const FixedNeighborPolicy & policy)
{
    FixedNeighborResult out;
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
    const auto search=SearchFixedNeighborComponent(*data.input,layout,domain,data.y,y,initial_eta,context,run_policy);
    out.state=search.state; out.blocks=search.blocks; out.sweeps=search.sweeps;
    out.search_converged=search.search_converged; out.reason=search.reason;
    out.first_order_stationarity_sweep=search.first_order_stationarity_sweep;
    out.confirmed_stationarity_sweep=search.confirmed_stationarity_sweep;
    out.block_preparations=search.block_preparations;
    out.domain_preparations=search.domain_preparations;
    out.mapping_preparations=search.mapping_preparations;
    if(policy.assess_final_endpoint) BuildEndpointResult(problem,data,layout,domain,y,context,search,out);
    return out;
}
}
