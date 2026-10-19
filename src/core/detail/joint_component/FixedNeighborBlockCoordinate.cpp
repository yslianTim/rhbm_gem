#include "FixedNeighborBlockCoordinate.hpp"
#include "LocalConditioning.hpp"
#include "ParameterStability.hpp"
#include "Problem.hpp"
#include "ResourceWork.hpp"
#include "StructuralPartition.hpp"
#include "TiledDerivative.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>

namespace rhbm_gem::core::joint_component {
double BlockObjectiveReplayEnclosure(double reference)
{return 1e-12+2e-12*std::abs(reference);}
bool WithinBlockObjectiveReplay(double error,double reference)
{return error<=BlockObjectiveReplayEnclosure(reference);}
IncrementalCandidateReplay EvaluateCandidateDelta(VectorRef current_residual,
    const Indices & affected_rows,VectorRef old_core,VectorRef new_core,
    double current_objective,double scale)
{
    if(affected_rows.size()!=static_cast<std::size_t>(old_core.size()) || old_core.size()!=new_core.size() ||
        !(scale>0) || !std::isfinite(scale)) throw std::invalid_argument("Invalid incremental candidate dimensions.");
    IncrementalCandidateReplay out; out.residual_delta.resize(new_core.size());
    long double delta{},correction{},absolute_delta{};
    const long double normalization=2.0L*static_cast<long double>(scale)*static_cast<long double>(scale);
    for(Eigen::Index k=0;k<new_core.size();++k)
    {
        const auto row=affected_rows[static_cast<std::size_t>(k)];
        const double change=new_core(k)-old_core(k); out.residual_delta(k)=change;
        const long double residual=current_residual(row),difference=change;
        const long double term=(2.0L*residual*difference+difference*difference)/normalization;
        const long double adjusted=term-correction,next=delta+adjusted;
        correction=(next-delta)-adjusted; delta=next; absolute_delta+=std::abs(term);
    }
    out.objective_delta=static_cast<double>(delta);
    out.objective=current_objective+out.objective_delta;
    out.delta_condition=absolute_delta==0 ? 0.0 : delta==0 ? std::numeric_limits<double>::infinity() :
        static_cast<double>(absolute_delta/std::abs(delta));
    return out;
}

namespace {
using Clock=std::chrono::steady_clock;
constexpr std::size_t EndpointCertificationRetryInterval=10;
constexpr std::size_t MaximumBlockFailureVisits=2;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
struct OptionalPhaseTimer
{
    double * target{};
    Clock::time_point started{Clock::now()};
    explicit OptionalPhaseTimer(double * target):target(target) {}
    ~OptionalPhaseTimer() {if(target) *target+=Seconds(started);}
};
OptionalPhaseTimer Measure(bool enabled,double & seconds)
{return OptionalPhaseTimer{enabled ? &seconds : nullptr};}
Vector Select(VectorRef values,const Indices & indices)
{
    Vector out(static_cast<Eigen::Index>(indices.size()));
    for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=values(indices[k]);
    return out;
}
Indices IndicesOf(const std::vector<std::size_t> & values)
{return {values.begin(),values.end()};}
bool AddMemoryBytes(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
}
bool PartitionStorageBytes(const StructuralBlockPartition & partition,bool include_metadata,std::size_t & bytes)
{
    bytes=0;
    if(!AddMemoryBytes(bytes,partition.cores.capacity(),sizeof(StructuralCore))) return false;
    for(const auto & core:partition.cores)
    {
        std::size_t id_bytes{};
        if(!CheckedAddBytes(core.id.capacity(),1,id_bytes) || !CheckedAddBytes(id_bytes,32,id_bytes) ||
            !CheckedAddBytes(bytes,id_bytes,bytes) || !AddMemoryBytes(bytes,core.atoms.capacity(),sizeof(Eigen::Index)) ||
            !AddMemoryBytes(bytes,1,4*32)) return false;
        if(include_metadata && (!AddMemoryBytes(bytes,core.affected_rows.capacity(),sizeof(Eigen::Index)) ||
            !AddMemoryBytes(bytes,core.neighbor_atoms.capacity(),sizeof(Eigen::Index)) ||
            !AddMemoryBytes(bytes,core.context_atoms.capacity(),sizeof(Eigen::Index)))) return false;
    }
    return true;
}
Vector SelectBeta(VectorRef beta,const std::vector<std::size_t> & atoms)
{
    Vector out(2*static_cast<Eigen::Index>(atoms.size()));
    for(std::size_t k=0;k<atoms.size();++k)
        out.segment<2>(2*static_cast<Eigen::Index>(k))=beta.segment<2>(2*static_cast<Eigen::Index>(atoms[k]));
    return out;
}
struct ReplayState {Vector prediction,residual; double objective{};};
ReplayState Replay(const Domain & domain,const JointParameterLayout & layout,VectorRef observations,
    VectorRef eta,VectorRef beta,double scale)
{
    ReplayState out; out.prediction=Vector::Zero(observations.size()); out.residual=Vector::Zero(observations.size());
    std::vector<std::uint8_t> informative(static_cast<std::size_t>(domain.rows));
    for(auto row:layout.informative_rows) informative.at(row)=1;
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); const double width=std::exp(eta(a));
        for(const auto & support:domain.atoms[atom]) if(informative.at(static_cast<std::size_t>(support.row)))
        {
            const auto basis=EvaluateKernel(support.square,width,2.5);
            out.prediction(support.row)+=beta(2*a)*basis.gaussian+beta(2*a+1)*basis.charge;
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
double SparseColumnDot(const Sparse & x,Eigen::Index left,Eigen::Index right)
{
    double out{};
    Sparse::InnerIterator a(x,left),b(x,right);
    while(a && b)
    {
        if(a.row()==b.row()) {out+=a.value()*b.value(); ++a; ++b;}
        else if(a.row()<b.row()) ++a;
        else ++b;
    }
    return out;
}
void CaptureWidthTrajectory(FixedNeighborBlockRecord & record,const SearchResult & search,
    const Domain & domain,VectorRef response,VectorRef block_entry_eta,const EvaluationContext & context,
    FixedNeighborWidthTrajectoryCapture & capture)
{
    if(capture.maximum_events==0 || capture.maximum_atoms==0 || capture.atom_initial_widths.empty()) return;
    std::vector<std::pair<Eigen::Index,double>> tracked;
    const auto count=std::min(capture.maximum_atoms,capture.atom_initial_widths.size());
    for(std::size_t atom_index=0;atom_index<context.atom_ids.size();++atom_index)
    {
        const auto atom=static_cast<Eigen::Index>(atom_index);
        const auto & id=context.atom_ids[atom_index];
        const auto initial=std::find_if(capture.atom_initial_widths.begin(),
            capture.atom_initial_widths.begin()+static_cast<std::ptrdiff_t>(count),
            [&](const auto & value){return value.first==id;});
        if(initial!=capture.atom_initial_widths.begin()+static_cast<std::ptrdiff_t>(count))
            tracked.emplace_back(atom,initial->second);
    }
    if(tracked.empty()) return;
    double previous_accepted_objective=record.local_objective_before;
    for(std::size_t trial_index=0;trial_index<search.trials.size();++trial_index)
    {
        if(capture.events_emitted>=capture.maximum_events) break;
        const auto & trial=search.trials[trial_index];
        if(capture.accepted_updates_only &&
            (!trial.accepted || trial.accepted_update.value_or(0)==0)) continue;
        FixedNeighborProfileTrial * destination=nullptr;
        const auto existing=std::find_if(record.profile_trials.begin(),record.profile_trials.end(),
            [&](const auto & value){return value.profile_evaluation==trial.evaluation;});
        if(existing!=record.profile_trials.end()) destination=&*existing;
        else
        {
            record.profile_trials.emplace_back(); destination=&record.profile_trials.back();
            destination->trial_index=trial_index; destination->profile_evaluation=trial.evaluation;
            destination->accepted=trial.accepted; destination->accepted_update=trial.accepted_update;
            destination->profile_seconds=trial.seconds; destination->factor_seconds=trial.factor_seconds;
            destination->local_objective_before=previous_accepted_objective;
            destination->local_objective_after=trial.endpoint.certificate.available ?
                trial.endpoint.certificate.objective/(context.scale*context.scale) : unavailable;
            destination->objective_reduction=std::isfinite(destination->local_objective_before) &&
                std::isfinite(destination->local_objective_after) ?
                destination->local_objective_before-destination->local_objective_after : unavailable;
            destination->eta_change_inf=trial.endpoint.eta.size()==block_entry_eta.size() ?
                (trial.endpoint.eta-block_entry_eta).lpNorm<Eigen::Infinity>() : unavailable;
            destination->gradient_inf_norm=trial.endpoint.gradient.size() ?
                trial.endpoint.gradient.lpNorm<Eigen::Infinity>() : unavailable;
        }
        destination->tracked_atoms.clear();
        const auto objective=trial.endpoint.certificate.available ?
            trial.endpoint.certificate.objective/(context.scale*context.scale) : unavailable;
        if(trial.lm)
        {
            destination->lm_actual_decrease=trial.lm->actual_decrease;
            destination->lm_predicted_decrease=trial.lm->predicted_decrease;
            destination->lm_ratio=trial.lm->ratio;
            destination->lm_proposed_acceptance=trial.lm->proposed_acceptance;
        }
        destination->weak_width_step_fraction=trial.weak_width_step_fraction;
        destination->weak_width_guard_rejected=trial.rejection_reason=="weak-width-identifiability-limited";
        destination->acceptance_reason=!trial.endpoint.valid ? trial.endpoint.reason :
            !trial.rejection_reason.empty() ? trial.rejection_reason :
            trial.trust && !trial.trust->passed ? "untrusted:"+trial.trust->reason :
            trial.lm && !trial.lm->proposed_acceptance ? "lm-rejected" :
            trial.accepted && trial.accepted_update.value_or(0)==0 ? "initial-profile" :
            trial.accepted ? "accepted" : "not-selected";
        if(trial.accepted)
        {
            destination->accepted_local_objective_improvement=trial.accepted_update.value_or(0)>0 &&
                std::isfinite(previous_accepted_objective) && std::isfinite(objective) ?
                previous_accepted_objective-objective : unavailable;
            if(std::isfinite(objective)) previous_accepted_objective=objective;
        }

        Evaluation state;
        Spectrum spectrum;
        const bool within_row_limit=static_cast<std::size_t>(domain.rows)<=capture.maximum_rows;
        if(within_row_limit && trial.endpoint.eta.size()==block_entry_eta.size() &&
            trial.endpoint.beta.size()==2*block_entry_eta.size())
        {
            state=EvaluateState(domain,response,trial.endpoint.eta,trial.endpoint.beta,context);
            if(state.valid) spectrum=ComputeSpectrum(state.x,context.rank,state.x.cols(),true);
            if(state.valid)
            {
                const auto projected_prepared=PrepareProjectedWidthDerivative(state,context.scale,&context,-1,
                    derivative_tile_rows,true);
                if(projected_prepared.valid)
                {
                    const auto projected=ReduceProjectedWidthDerivative(projected_prepared);
                    if(projected.valid)
                    {
                        const auto width_spectrum=ComputeSpectrum(projected.projected,context.rank,
                            projected.projected.cols(),true);
                        if(width_spectrum.available)
                            destination->projected_width_singular_values.assign(width_spectrum.singular_values.data(),
                                width_spectrum.singular_values.data()+width_spectrum.singular_values.size());
                    }
                }
            }
        }
        for(const auto & [atom,initial_width]:tracked)
        {
            if(destination->tracked_atoms.size()>=capture.maximum_atoms) break;
            FixedNeighborProfileTrial::AtomState value;
            value.atom_id=context.atom_ids[static_cast<std::size_t>(atom)];
            value.first_stage_width=initial_width;
            value.block_entry_width=std::exp(block_entry_eta(atom));
            const double before=trial.lm && trial.lm->accepted_eta.size()==block_entry_eta.size() ?
                trial.lm->accepted_eta(atom) : block_entry_eta(atom);
            const double proposed=trial.endpoint.eta.size()==block_entry_eta.size() ? trial.endpoint.eta(atom) : unavailable;
            const double after=trial.accepted ? proposed : before;
            value.eta_before=before; value.eta_after=after;
            value.width_before=std::exp(before); value.width_after=std::exp(after);
            value.lm_proposed_width=std::isfinite(proposed) ? std::exp(proposed) : unavailable;
            value.accepted_width=std::isfinite(after) ? std::exp(after) : unavailable;
            value.eta_step_inf=trial.lm ? trial.lm->step.lpNorm<Eigen::Infinity>() : 0.0;
            if(trial.lm && trial.lm->step.size()==block_entry_eta.size())
            {
                value.proposed_eta_step=trial.lm->step(atom);
                value.accepted_eta_step=trial.accepted ? value.proposed_eta_step : 0.0;
            }
            if(trial.endpoint.beta.size()==2*block_entry_eta.size())
            {
                value.amplitude_a=trial.endpoint.beta(2*atom);
                value.amplitude_c=trial.endpoint.beta(2*atom+1);
            }
            if(state.valid)
            {
                value.amplitude_active=std::find(state.certificate.active_atoms.begin(),
                    state.certificate.active_atoms.end(),atom)!=state.certificate.active_atoms.end();
                const double gnorm=state.x.col(2*atom).norm(),knorm=state.x.col(2*atom+1).norm();
                value.gaussian_column_norm=gnorm; value.charge_column_norm=knorm;
                value.normalized_a=value.amplitude_a*gnorm; value.normalized_c=value.amplitude_c*knorm;
                value.gaussian_contribution_norm=std::abs(value.normalized_a);
                value.charge_contribution_norm=std::abs(value.normalized_c);
                const Vector atom_contribution=(state.x.col(2*atom)*value.amplitude_a+
                    state.x.col(2*atom+1)*value.amplitude_c).eval();
                const double contribution_norm=atom_contribution.stableNorm();
                const double contribution_sum=value.gaussian_contribution_norm+value.charge_contribution_norm;
                if(contribution_norm>0) value.ac_cancellation_ratio=contribution_sum/contribution_norm;
                else if(contribution_sum==0) value.ac_cancellation_ratio=1.0;
                const double denominator=gnorm*knorm;
                if(denominator>0) value.basis_correlation=SparseColumnDot(state.x,2*atom,2*atom+1)/denominator;
                value.local_kkt=state.certificate.projected_kkt;
                if(state.gradient.size()>atom) value.width_gradient=std::abs(state.gradient(atom));
                value.local_design_rank=spectrum.available ? spectrum.rank : -1;
                value.local_design_condition=spectrum.available ? spectrum.condition : unavailable;
                value.numerical_status=state.certificate.available ? "evaluated" : "certificate-unavailable";
            }
            else value.numerical_status=within_row_limit ? "state-evaluation-unavailable" : "row-limit";
            destination->tracked_atoms.push_back(std::move(value));
        }
        if(!destination->tracked_atoms.empty()) ++capture.events_emitted;
    }
}
#endif
double PredictionDifference(VectorRef a,VectorRef b)
{return (a-b).lpNorm<Eigen::Infinity>();}
void WriteJsonString(std::ostream & out,const std::string & value)
{
    out<<'"';
    for(const char byte:value)
    {
        const auto character=static_cast<unsigned char>(byte);
        if(character=='"' || character=='\\') out<<'\\'<<character;
        else if(character=='\n') out<<"\\n";
        else if(character=='\r') out<<"\\r";
        else if(character=='\t') out<<"\\t";
        else if(character<0x20) out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<static_cast<int>(character)<<std::dec<<std::setfill(' ');
        else out<<character;
    }
    out<<'"';
}
void WriteJsonNumber(std::ostream & out,double value)
{if(std::isfinite(value)) out<<value; else out<<"null";}
void WriteJsonVector(std::ostream & out,VectorRef values)
{
    out<<'['; for(Eigen::Index k=0;k<values.size();++k) {if(k) out<<','; WriteJsonNumber(out,values(k));} out<<']';
}
void WriteJsonStrings(std::ostream & out,const std::vector<std::string> & values)
{
    out<<'['; for(std::size_t k=0;k<values.size();++k) {if(k) out<<','; WriteJsonString(out,values[k]);} out<<']';
}
void WriteJsonIndices(std::ostream & out,const std::vector<Eigen::Index> & values)
{
    out<<'['; for(std::size_t k=0;k<values.size();++k) {if(k) out<<','; out<<values[k];} out<<']';
}
void WriteOptionalJson(std::ostream & out,const std::optional<double> & value)
{if(value) WriteJsonNumber(out,*value); else out<<"null";}
void WriteOptionalJson(std::ostream & out,const std::optional<int> & value)
{if(value) out<<*value; else out<<"null";}
void WriteOptionalJson(std::ostream & out,const std::optional<Eigen::Index> & value)
{if(value) out<<*value; else out<<"null";}
std::string CaptureLocalFailureFixture(const FixedNeighborBlockFailureDiagnostics & d,const Domain & domain,
    VectorRef effective_response,VectorRef eta,const EvaluationContext & parent_context,
    const EvaluationContext & local_context,const Endpoint & initial)
{
    const char * const destination=std::getenv("RHBM_JOINT_LOCAL_FAILURE_CAPTURE");
    if(!destination || !*destination) return "not-requested";
    try
    {
        const std::filesystem::path path(destination);
        if(std::filesystem::exists(path)) return "destination-exists";
        std::ofstream out(path,std::ios::binary);
        if(!out) return "open-failed";
        out<<std::setprecision(std::numeric_limits<double>::max_digits10);
        out<<"{\n  \"format\":\"rhbm-joint-local-ac-kkt-v1\",\n  \"identity\":{\"component_id\":";
        WriteJsonString(out,d.component_id);
        out<<",\"component_index\":"<<d.component_index<<",\"sweep_index\":"<<d.sweep_index
           <<",\"block_index\":"<<d.block_index<<",\"core_atom_count\":"<<d.core_atom_count
           <<",\"core_atom_ids\":"; WriteJsonStrings(out,d.core_atom_ids);
        out<<",\"affected_row_count\":"<<d.affected_row_count<<",\"local_design_rows\":"<<d.local_design_rows
           <<",\"local_design_columns\":"<<d.local_design_columns<<",\"local_design_nonzero_entries\":"
           <<d.local_design_nonzero_entries<<",\"local_design_rank\":";
        WriteOptionalJson(out,d.local_design_rank);
        out<<"},\n  \"atom_ids\":"; std::vector<std::string> atom_ids;
        for(const auto & id:local_context.atom_ids) atom_ids.push_back(id);
        WriteJsonStrings(out,atom_ids);
        out<<",\n  \"row_ids\":"; std::vector<std::string> row_ids;
        for(const auto & id:local_context.row_ids) row_ids.push_back(id);
        WriteJsonStrings(out,row_ids);
        out<<",\n  \"effective_response\":"; WriteJsonVector(out,effective_response);
        out<<",\n  \"old_eta\":"; WriteJsonVector(out,eta);
        out<<",\n  \"old_widths\":"; const Vector widths=eta.array().exp(); WriteJsonVector(out,widths);
        out<<",\n  \"support\":[";
        for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        {
            if(atom) out<<','; out<<'['; std::size_t support_index{};
            for(const auto support:domain.atoms[atom])
            {
                if(support_index++) out<<',';
                out<<"{\"row\":"<<support.row<<",\"squared_distance\":";
                WriteJsonNumber(out,support.square); out<<'}';
            }
            out<<']';
        }
        out<<"],\n  \"numerical_context\":{\"parent_observation_scale\":";
        WriteJsonNumber(out,parent_context.scale);
        out<<",\"parent_response_norm\":"; WriteJsonNumber(out,parent_context.linear.release_response_norm);
        out<<",\"local_effective_response_norm\":"; WriteJsonNumber(out,effective_response.norm());
        out<<",\"local_response_min\":"; WriteJsonNumber(out,effective_response.minCoeff());
        out<<",\"local_response_max\":"; WriteJsonNumber(out,effective_response.maxCoeff());
        out<<",\"parent_rank_relative\":"; WriteJsonNumber(out,parent_context.linear.rank_relative);
        out<<",\"rank_relative\":"; WriteJsonNumber(out,local_context.linear.rank_relative);
        out<<",\"release_factor\":"; WriteJsonNumber(out,local_context.linear.release_factor);
        out<<",\"release_response_norm\":"; WriteJsonNumber(out,local_context.linear.release_response_norm);
        out<<",\"active_set_iteration_factor\":"<<local_context.linear.active_set_iteration_factor<<"},\n  \"production_solver\":{\"valid\":"
           <<(initial.certificate.linear_solver_valid ? "true" : "false")<<",\"reason\":";
        WriteJsonString(out,initial.reason);
        out<<",\"solver_reason\":"; WriteJsonString(out,initial.certificate.linear_solver_reason);
        out<<",\"beta\":"; WriteJsonVector(out,initial.beta);
        out<<",\"linear_rank\":"; WriteOptionalJson(out,initial.certificate.free_rank);
        out<<",\"linear_solve_count\":"; WriteOptionalJson(out,initial.certificate.linear_solves);
        out<<",\"active_set_releases\":"; WriteOptionalJson(out,initial.certificate.linear_releases);
        std::vector<Eigen::Index> active_a_atoms;
        for(const auto amplitude:d.active_amplitude_indices) active_a_atoms.push_back(amplitude/2);
        out<<",\"active_a_atom_indices\":"; WriteJsonIndices(out,active_a_atoms);
        out<<",\"active_amplitude_indices\":"; WriteJsonIndices(out,d.active_amplitude_indices);
        out<<",\"primal_feasibility\":"; WriteOptionalJson(out,d.primal_feasibility);
        out<<",\"scaled_projected_kkt\":"; WriteOptionalJson(out,d.scaled_projected_kkt);
        out<<",\"column_norms\":"; const Sparse design=Design(domain,eta);
        out<<'['; for(Eigen::Index k=0;k<design.cols();++k) {if(k) out<<','; WriteJsonNumber(out,design.col(k).norm());} out<<']';
        out<<",\"raw_gradient\":";
        const Vector residual=design*initial.beta-effective_response;
        const Vector raw_gradient=design.transpose()*residual;
        WriteJsonVector(out,raw_gradient);
        out<<",\"scaled_gradient\":";
        const Vector column_norms=[&] {Vector value(design.cols()); for(Eigen::Index k=0;k<design.cols();++k) value(k)=design.col(k).norm(); return value;}();
        const Vector scaled_gradient=raw_gradient.array()/column_norms.array()/parent_context.scale;
        WriteJsonVector(out,scaled_gradient);
        out<<",\"objective\":"; WriteJsonNumber(out,initial.certificate.objective);
        out<<"}\n}\n";
        out.close();
        return out ? "written" : "write-failed";
    }
    catch(const std::exception &)
    {return "write-failed";}
}
std::string CaptureRejectedTrialFixture(const FixedNeighborBlockFailureDiagnostics & d,const Domain & domain,
    VectorRef effective_response,const EvaluationContext & parent_context,const EvaluationContext & local_context,
    const SearchResult & profile)
{
    const char * const destination=std::getenv("RHBM_JOINT_LOCAL_TRIAL_CAPTURE");
    if(!destination || !*destination) return "not-requested";
    const auto failed=std::find_if(profile.trials.rbegin(),profile.trials.rend(),[](const Trial & trial) {
        return trial.trust && !trial.trust->passed;
    });
    if(failed==profile.trials.rend()) return "no-rejected-trial";
    std::size_t support_count{};
    for(const auto & atom:domain.atoms) support_count+=atom.size();
    if(domain.rows<=0 || domain.rows>4096 || domain.atoms.size()==0 || domain.atoms.size()>256 ||
        support_count>65536 || effective_response.size()!=domain.rows)
        return "fixture-size-limit";
    try
    {
        const auto & endpoint=failed->endpoint;
        const auto primary=EvaluateState(domain,effective_response,endpoint.eta,endpoint.beta,local_context);
        if(!primary.valid) return "primary-reconstruction-failed";
        ReplayDetailsForDiagnostics replay;
        const auto trust=CheckReplay(domain,effective_response,primary,local_context,&replay);
        const Vector primary_prediction=(primary.x*primary.beta).eval();
        const Vector widths=endpoint.eta.array().exp();
        std::vector<std::string> atom_ids,row_ids;
        for(const auto & id:local_context.atom_ids) atom_ids.push_back(id);
        for(const auto & id:local_context.row_ids) row_ids.push_back(id);
        std::ostringstream payload;
        payload<<std::setprecision(std::numeric_limits<double>::max_digits10);
        payload<<"{\n  \"format\":\"rhbm-joint-local-replay-trial-v1\",\n  \"identity\":{\"component_id\":";
        WriteJsonString(payload,d.component_id);
        payload<<",\"component_index\":"<<d.component_index<<",\"sweep_index\":"<<d.sweep_index
            <<",\"block_index\":"<<d.block_index<<",\"trial_index\":"
            <<(profile.trials.size()-1-static_cast<std::size_t>(std::distance(profile.trials.rbegin(),failed)))
            <<",\"profile_evaluation\":"<<failed->evaluation<<",\"local_observation_count\":"<<domain.rows
            <<",\"local_atom_count\":"<<domain.atoms.size()<<"},\n  \"atom_ids\":";
        WriteJsonStrings(payload,atom_ids);
        payload<<",\n  \"row_ids\":"; WriteJsonStrings(payload,row_ids);
        payload<<",\n  \"local_support\":[";
        for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        {
            if(atom) payload<<','; payload<<'[';
            std::size_t item{};
            for(const auto & support:domain.atoms[atom])
            {
                if(item++) payload<<',';
                payload<<"{\"row\":"<<support.row<<",\"squared_distance\":";
                WriteJsonNumber(payload,support.square); payload<<'}';
            }
            payload<<']';
        }
        payload<<"],\n  \"effective_response\":"; WriteJsonVector(payload,effective_response);
        payload<<",\n  \"local_numerical_context\":{\"parent_observation_scale\":";
        WriteJsonNumber(payload,parent_context.scale);
        payload<<",\"local_observation_scale\":"; WriteJsonNumber(payload,local_context.scale);
        payload<<",\"parent_rank_relative\":"; WriteJsonNumber(payload,parent_context.linear.rank_relative);
        payload<<",\"rank_relative\":"; WriteJsonNumber(payload,local_context.linear.rank_relative);
        payload<<",\"release_factor\":"; WriteJsonNumber(payload,local_context.linear.release_factor);
        payload<<",\"release_response_norm\":"; WriteJsonNumber(payload,local_context.linear.release_response_norm);
        payload<<",\"active_set_iteration_factor\":"<<local_context.linear.active_set_iteration_factor<<"},\n";
        payload<<"  \"trial_eta\":"; WriteJsonVector(payload,endpoint.eta);
        payload<<",\n  \"trial_widths\":"; WriteJsonVector(payload,widths);
        payload<<",\n  \"trial_beta\":"; WriteJsonVector(payload,endpoint.beta);
        payload<<",\n  \"primary_prediction\":"; WriteJsonVector(payload,primary_prediction);
        payload<<",\n  \"primary_residual\":"; WriteJsonVector(payload,primary.residual);
        payload<<",\n  \"primary_width_gradient\":"; WriteJsonVector(payload,primary.gradient);
        payload<<",\n  \"primary_objective\":"; WriteJsonNumber(payload,endpoint.certificate.objective);
        payload<<",\n  \"primary_scaled_kkt\":"; WriteJsonNumber(payload,endpoint.certificate.projected_kkt);
        payload<<",\n  \"primary_valid\":"<<(endpoint.valid?"true":"false")
            <<",\n  \"primary_certificate\":{\"available\":"<<(endpoint.certificate.available?"true":"false")
            <<",\"feasible\":"<<(endpoint.certificate.feasible?"true":"false")
            <<",\"kkt_passed\":"<<(endpoint.certificate.kkt_passed?"true":"false")
            <<",\"linear_solver_valid\":"<<(endpoint.certificate.linear_solver_valid?"true":"false")
            <<",\"reason\":";
        WriteJsonString(payload,endpoint.reason);
        payload<<",\"solver_reason\":"; WriteJsonString(payload,endpoint.certificate.linear_solver_reason);
        payload<<",\"free_rank\":"; WriteOptionalJson(payload,endpoint.certificate.free_rank);
        payload<<",\"linear_solves\":"; WriteOptionalJson(payload,endpoint.certificate.linear_solves);
        payload<<",\"linear_releases\":"; WriteOptionalJson(payload,endpoint.certificate.linear_releases);
        payload<<"},\n  \"replay_prediction\":"; WriteJsonVector(payload,replay.prediction);
        payload<<",\n  \"replay_residual\":"; WriteJsonVector(payload,replay.residual);
        payload<<",\n  \"replay_width_gradient\":"; WriteJsonVector(payload,replay.width_gradient);
        payload<<",\n  \"replay_objective\":"; WriteJsonNumber(payload,replay.objective);
        payload<<",\n  \"replay_scaled_kkt\":"; WriteJsonNumber(payload,replay.scaled_kkt);
        payload<<",\n  \"prediction_difference\":"; WriteJsonNumber(payload,trust.prediction_difference);
        payload<<",\n  \"gradient_difference\":"; WriteJsonNumber(payload,trust.gradient_difference);
        payload<<",\n  \"kkt_difference\":"; WriteJsonNumber(payload,trust.kkt_difference);
        payload<<",\n  \"trust\":{\"prediction_passed\":"<<(trust.prediction_passed?"true":"false")
            <<",\"gradient_passed\":"<<(trust.gradient_passed?"true":"false")
            <<",\"kkt_agrees\":"<<(trust.kkt_difference<=1e-13?"true":"false")
            <<",\"trust_passed\":"<<(trust.passed?"true":"false")<<",\"trust_reason\":";
        WriteJsonString(payload,trust.reason);
        payload<<"},\n  \"build_metadata\":{\"compiler\":";
#if defined(__clang__)
        WriteJsonString(payload,std::string("clang ")+__clang_version__);
#elif defined(__GNUC__)
        WriteJsonString(payload,std::string("gcc ")+__VERSION__);
#else
        WriteJsonString(payload,"unknown");
#endif
#if defined(__aarch64__) || defined(__arm64__)
        payload<<",\"architecture\":\"arm64\"";
#elif defined(__x86_64__) || defined(_M_X64)
        payload<<",\"architecture\":\"x86_64\"";
#else
        payload<<",\"architecture\":\"unknown\"";
#endif
        payload<<",\"double_digits\":"<<std::numeric_limits<double>::digits
            <<",\"double_iec559\":"<<(std::numeric_limits<double>::is_iec559?"true":"false");
#ifdef __FAST_MATH__
        payload<<",\"fast_math_macro\":true";
#else
        payload<<",\"fast_math_macro\":false";
#endif
        payload<<"}\n}\n";
        const auto contents=payload.str();
        if(contents.size()>16*1024*1024) return "fixture-size-limit";
        const std::filesystem::path path(destination),temporary(path.string()+".tmp");
        if(std::filesystem::exists(path)) return "destination-exists";
        if(std::filesystem::exists(temporary)) return "temporary-exists";
        {std::ofstream out(temporary,std::ios::binary); out.exceptions(std::ios::failbit|std::ios::badbit); out<<contents;}
        std::filesystem::create_hard_link(temporary,path);
        std::filesystem::remove(temporary);
        return "written";
    }
    catch(const std::exception &)
    {return "write-failed";}
}
FixedNeighborBlockFailureDiagnostics FailureDiagnostics(const FixedNeighborPolicy & policy,
    const FixedNeighborBlockRecord & block,std::size_t sweep_index,const Indices & core_atoms,
    const EvaluationContext & context,const EvaluationContext & local_context,const Domain & domain,VectorRef effective_response,
    const Vector & eta,const SearchResult & profile)
{
    FixedNeighborBlockFailureDiagnostics out;
    out.component_index=policy.component_index; out.component_id=policy.component_id;
    out.sweep_index=sweep_index+1; out.block_index=block.block;
    out.core_atom_count=core_atoms.size(); out.affected_row_count=block.affected_rows;
    out.local_design_rows=domain.rows; out.local_design_columns=2*eta.size();
    out.local_search_stop_reason=profile.stop_reason;
    out.initial_profile_valid=profile.initial.valid;
    out.initial_profile_reason=profile.initial.reason;
    out.initial_profile_trusted=profile.initial_accepted;
    out.local_response_norm=effective_response.norm();
    if(effective_response.size())
    {
        out.local_response_min=effective_response.minCoeff(); out.local_response_max=effective_response.maxCoeff();
    }
    out.parent_observation_scale=context.scale;
    out.parent_response_norm=context.linear.release_response_norm;
    out.parent_rank_relative=context.linear.rank_relative;
    out.release_factor=local_context.linear.release_factor;
    out.rank_relative=local_context.linear.rank_relative;
    out.active_set_iteration_factor=local_context.linear.active_set_iteration_factor;
    const auto design=Design(domain,eta);
    const auto spectrum=ComputeSpectrum(design,{domain.rows,design.cols(),eta.size()},design.cols(),true);
    if(spectrum.available) out.local_design_rank=spectrum.rank;
    out.local_design_nonzero_entries=static_cast<std::size_t>(design.nonZeros());
    if(spectrum.available) out.condition_indicator=spectrum.condition;
    out.local_column_norms.reserve(static_cast<std::size_t>(design.cols()));
    for(Eigen::Index column=0;column<design.cols();++column)
        out.local_column_norms.push_back(design.col(column).norm());
    const auto & certificate=profile.initial.certificate;
    out.solver_valid=certificate.linear_solver_valid;
    out.solver_reason=certificate.linear_solver_reason;
    out.linear_rank=certificate.free_rank; out.linear_solve_count=certificate.linear_solves;
    out.active_set_releases=certificate.linear_releases;
    out.active_amplitude_indices.clear();
    for(const auto atom:certificate.active_atoms) out.active_amplitude_indices.push_back(2*atom);
    out.primal_feasibility=0.0;
    for(Eigen::Index column=0;column<profile.initial.beta.size();column+=2)
        *out.primal_feasibility=std::max(*out.primal_feasibility,std::max(0.0,-profile.initial.beta(column)));
    if(certificate.available) out.scaled_projected_kkt=certificate.projected_kkt;
    if(profile.initial.beta.size()==design.cols() && effective_response.size()==design.rows() &&
        !out.local_column_norms.empty())
    {
        const Vector residual=design*profile.initial.beta-effective_response;
        const Vector raw_gradient=design.transpose()*residual;
        for(Eigen::Index atom=0;atom<eta.size();++atom)
        {
            const Eigen::Index amplitude=2*atom;
            const double scaled_a=raw_gradient(amplitude)/out.local_column_norms[static_cast<std::size_t>(amplitude)]/context.scale;
            if(profile.initial.beta(amplitude)>0)
            {
                const double raw=std::abs(raw_gradient(amplitude));
                const double scaled=std::abs(scaled_a);
                out.maximum_free_a_gradient_raw=out.maximum_free_a_gradient_raw ?
                    std::max(*out.maximum_free_a_gradient_raw,raw) : raw;
                out.maximum_free_a_gradient_scaled=out.maximum_free_a_gradient_scaled ?
                    std::max(*out.maximum_free_a_gradient_scaled,scaled) : scaled;
            }
            else
            {
                const double raw=std::max(0.0,-raw_gradient(amplitude));
                const double scaled=std::max(0.0,-scaled_a);
                out.maximum_active_a_complementarity_violation_raw=out.maximum_active_a_complementarity_violation_raw ?
                    std::max(*out.maximum_active_a_complementarity_violation_raw,raw) : raw;
                out.maximum_active_a_complementarity_violation_scaled=out.maximum_active_a_complementarity_violation_scaled ?
                    std::max(*out.maximum_active_a_complementarity_violation_scaled,scaled) : scaled;
            }
            const Eigen::Index charge=amplitude+1;
            const double raw=std::abs(raw_gradient(charge));
            const double scaled=raw/out.local_column_norms[static_cast<std::size_t>(charge)]/context.scale;
            out.maximum_free_c_gradient_raw=out.maximum_free_c_gradient_raw ?
                std::max(*out.maximum_free_c_gradient_raw,raw) : raw;
            out.maximum_free_c_gradient_scaled=out.maximum_free_c_gradient_scaled ?
                std::max(*out.maximum_free_c_gradient_scaled,scaled) : scaled;
        }
    }
    for(const auto atom:core_atoms)
        if(atom>=0 && static_cast<std::size_t>(atom)<context.atom_ids.size())
            out.core_atom_ids.push_back(context.atom_ids.at(static_cast<std::size_t>(atom)));
    if(policy.halo_census) policy.halo_census(out.halo_count,out.minimum_halo_support_rows);
    for(auto trial=profile.trials.rbegin();trial!=profile.trials.rend();++trial)
        if(trial->trust && !trial->trust->passed)
        {
            out.replay_prediction_discrepancy=trial->trust->prediction_difference;
            out.replay_gradient_discrepancy=trial->trust->gradient_difference;
            out.replay_kkt_discrepancy=trial->trust->kkt_difference;
            break;
        }
    out.failure_capture_status=CaptureLocalFailureFixture(out,domain,effective_response,eta,context,local_context,profile.initial);
    out.replay_trial_capture_status=CaptureRejectedTrialFixture(out,domain,effective_response,context,local_context,profile);
    return out;
}
std::string FormatFailureDiagnostics(const FixedNeighborBlockFailureDiagnostics & d)
{
    std::ostringstream out;
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"[Joint] block-search-failed"
       <<" component_index="<<d.component_index<<" component_id="<<d.component_id
       <<" sweep_index="<<d.sweep_index<<" block_index="<<d.block_index
       <<" core_atom_ids=[";
    for(std::size_t i=0;i<d.core_atom_ids.size();++i) out<<(i ? "," : "")<<d.core_atom_ids[i];
    out<<"] core_atom_count="<<d.core_atom_count<<" affected_row_count="<<d.affected_row_count
       <<" halo_count="<<d.halo_count<<" minimum_halo_support_rows="<<d.minimum_halo_support_rows
       <<" initial_profile_valid="<<(d.initial_profile_valid ? "true" : "false")
       <<" initial_profile_reason="<<d.initial_profile_reason
       <<" initial_profile_trusted="<<(d.initial_profile_trusted ? "true" : "false")
       <<" local_design_rows="<<d.local_design_rows<<" local_design_columns="<<d.local_design_columns
       <<" local_design_nonzero_entries="<<d.local_design_nonzero_entries
       <<" local_design_rank=";
    if(d.local_design_rank) out<<*d.local_design_rank; else out<<"unavailable";
    out<<" condition_indicator="; if(d.condition_indicator) out<<*d.condition_indicator; else out<<"unavailable";
    out<<" local_response_norm="; if(d.local_response_norm) out<<*d.local_response_norm; else out<<"unavailable";
    out<<" local_response_min="; if(d.local_response_min) out<<*d.local_response_min; else out<<"unavailable";
    out<<" local_response_max="; if(d.local_response_max) out<<*d.local_response_max; else out<<"unavailable";
    out<<" parent_observation_scale="; if(d.parent_observation_scale) out<<*d.parent_observation_scale; else out<<"unavailable";
    out<<" parent_response_norm="; if(d.parent_response_norm) out<<*d.parent_response_norm; else out<<"unavailable";
    out<<" rank_relative="; if(d.rank_relative) out<<*d.rank_relative; else out<<"unavailable";
    out<<" parent_rank_relative="; if(d.parent_rank_relative) out<<*d.parent_rank_relative; else out<<"unavailable";
    out<<" release_factor="; if(d.release_factor) out<<*d.release_factor; else out<<"unavailable";
    out<<" solver_valid="<<(d.solver_valid ? "true" : "false")<<" solver_reason="<<d.solver_reason;
    out<<" linear_rank="; if(d.linear_rank) out<<*d.linear_rank; else out<<"unavailable";
    out<<" linear_solve_count="; if(d.linear_solve_count) out<<*d.linear_solve_count; else out<<"unavailable";
    out<<" active_set_releases="; if(d.active_set_releases) out<<*d.active_set_releases; else out<<"unavailable";
    out<<" primal_feasibility="; if(d.primal_feasibility) out<<*d.primal_feasibility; else out<<"unavailable";
    out<<" scaled_projected_kkt="; if(d.scaled_projected_kkt) out<<*d.scaled_projected_kkt; else out<<"unavailable";
    out<<" maximum_free_a_gradient_raw="; if(d.maximum_free_a_gradient_raw) out<<*d.maximum_free_a_gradient_raw; else out<<"unavailable";
    out<<" maximum_free_a_gradient_scaled="; if(d.maximum_free_a_gradient_scaled) out<<*d.maximum_free_a_gradient_scaled; else out<<"unavailable";
    out<<" maximum_free_c_gradient_raw="; if(d.maximum_free_c_gradient_raw) out<<*d.maximum_free_c_gradient_raw; else out<<"unavailable";
    out<<" maximum_free_c_gradient_scaled="; if(d.maximum_free_c_gradient_scaled) out<<*d.maximum_free_c_gradient_scaled; else out<<"unavailable";
    out<<" maximum_active_a_complementarity_violation_raw=";
    if(d.maximum_active_a_complementarity_violation_raw) out<<*d.maximum_active_a_complementarity_violation_raw; else out<<"unavailable";
    out<<" maximum_active_a_complementarity_violation_scaled=";
    if(d.maximum_active_a_complementarity_violation_scaled) out<<*d.maximum_active_a_complementarity_violation_scaled; else out<<"unavailable";
    out<<" failure_capture_status="<<d.failure_capture_status
       <<" replay_trial_capture_status="<<d.replay_trial_capture_status
       <<" local_search_stop_reason="<<d.local_search_stop_reason;
    if(d.replay_prediction_discrepancy)
        out<<" replay_prediction_discrepancy="<<*d.replay_prediction_discrepancy
           <<" replay_gradient_discrepancy="<<*d.replay_gradient_discrepancy
           <<" replay_kkt_discrepancy="<<*d.replay_kkt_discrepancy;
    return out.str();
}
int StabilityProbePriority(const JointBlockConditioning & evidence)
{
    if(evidence.conditioning_classification=="Singular") return 0;
    if(evidence.conditioning_classification=="IllConditioned") return 1;
    if(evidence.conditioning_classification=="ActiveFace") return 2;
    if(evidence.conditioning_classification=="WellConditioned") return 3;
    return -1;
}
double StabilitySensitivityScore(const JointBlockConditioning & evidence)
{
    double score{};
    for(const auto * condition:{&evidence.free_ac_condition,&evidence.projected_width_condition,
        &evidence.profile_condition})
        if(*condition && std::isfinite(**condition) && **condition>=1.0)
            score=std::max(score,std::log10(**condition));
    return score;
}
std::vector<std::size_t> StabilityProbeOrder(const std::vector<JointBlockConditioning> & conditioning)
{
    struct Candidate {std::size_t index{}; int priority{}; double score{}; std::string block_id;};
    std::vector<Candidate> candidates;
    candidates.reserve(conditioning.size());
    for(std::size_t index=0;index<conditioning.size();++index)
    {
        const auto & block=conditioning[index];
        const int priority=StabilityProbePriority(block);
        const bool free_ac_singular=block.free_ac_rank && *block.free_ac_rank<block.free_ac_columns;
        if(priority<0 || (priority==0 && free_ac_singular)) continue;
        candidates.push_back({index,priority,StabilitySensitivityScore(block),block.block_id});
    }
    std::sort(candidates.begin(),candidates.end(),[](const auto & left,const auto & right) {
        if(left.priority!=right.priority) return left.priority<right.priority;
        if(left.score!=right.score) return left.score>right.score;
        return left.block_id<right.block_id;
    });
    std::vector<std::size_t> order; order.reserve(candidates.size());
    for(const auto & candidate:candidates) order.push_back(candidate.index);
    return order;
}
double ScaledCoefficientDifference(VectorRef lhs,VectorRef rhs)
{
    if(lhs.size()!=rhs.size() || lhs.size()==0) return std::numeric_limits<double>::infinity();
    return ((lhs-rhs).array().abs()/(1.0+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
}
}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
std::vector<std::size_t> RankFixedNeighborStabilityProbesForTesting(
    const std::vector<JointBlockConditioning> & conditioning)
{return StabilityProbeOrder(conditioning);}
#endif
bool IsFixedNeighborEtaChangeConfirmed(double eta_change_inf,bool has_previous_complete_sweep,
    const JointConvergenceCriteria & criteria)
{
    return IsValidJointConvergenceCriteria(criteria) && has_previous_complete_sweep &&
        std::isfinite(eta_change_inf) && eta_change_inf>=0.0 &&
        eta_change_inf<=criteria.confirmed_log_width_change;
}
FixedNeighborStationarityState ClassifyFixedNeighborStationarity(
    double ac_kkt,double width_gradient_inf,double eta_change_inf,bool has_previous_complete_sweep,
    const JointConvergenceCriteria & criteria)
{
    if(!IsValidJointConvergenceCriteria(criteria) || !std::isfinite(ac_kkt) || ac_kkt<0.0 ||
        ac_kkt>criteria.global_ac_kkt || !std::isfinite(width_gradient_inf) || width_gradient_inf<0.0 ||
        width_gradient_inf>criteria.global_width_gradient)
        return FixedNeighborStationarityState::NotStationary;
    return IsFixedNeighborEtaChangeConfirmed(eta_change_inf,has_previous_complete_sweep,criteria) ?
        FixedNeighborStationarityState::ConfirmedStationary : FixedNeighborStationarityState::CandidateStationary;
}
bool ShouldAttemptFixedNeighborCertification(FixedNeighborStationarityState stationarity,std::size_t sweep,
    std::size_t & next_retry_sweep)
{
    if(stationarity!=FixedNeighborStationarityState::ConfirmedStationary)
    {next_retry_sweep=0; return false;}
    if(next_retry_sweep!=0 && sweep<next_retry_sweep) return false;
    next_retry_sweep=sweep+EndpointCertificationRetryInterval;
    return true;
}
FixedNeighborSearchResult SearchFixedNeighborComponent(
    const Domain & domain,const JointParameterLayout & layout,VectorRef y,VectorRef initial_eta,const EvaluationContext & context,
    const FixedNeighborPolicy & policy,Vector initial_beta)
{
    FixedNeighborSearchResult out;
    const bool has_initial_beta=initial_beta.size()!=0;
    if(initial_eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) || !initial_eta.allFinite() ||
        y.size()!=domain.rows || !y.allFinite() || !(context.scale>0) || !std::isfinite(context.scale) ||
        (has_initial_beta && (initial_beta.size()!=2*static_cast<Eigen::Index>(domain.atoms.size()) ||
            !initial_beta.allFinite())) ||
        policy.core_atoms==0 || layout.full_atoms.size()!=domain.atoms.size() ||
        layout.informative_rows.size()!=static_cast<std::size_t>(domain.rows) ||
        context.atom_ids.size()!=domain.atoms.size() || context.row_ids.size()!=static_cast<std::size_t>(domain.rows))
    {out.reason="block-invalid-partition"; return out;}
    for(Eigen::Index a=0;a<initial_eta.size();++a)
        if(!(std::exp(initial_eta(a))>0) || !std::isfinite(std::exp(initial_eta(a))))
        {out.reason="block-invalid-partition"; return out;}

    auto structural_timer=Measure(policy.collect_telemetry,out.work.structural_preparation_seconds);
    MemoryShape memory_shape; memory_shape.component_atoms=domain.atoms.size();
    memory_shape.component_observations=static_cast<std::size_t>(domain.rows);
    if(!CheckedMultiplyBytes(memory_shape.component_atoms,2,memory_shape.component_free_columns))
    {out.memory_estimate.stage=MemoryStage::StructuralConstruction; out.memory_estimate.status=MemoryStatus::Overflow;
        out.memory_estimate.reason="size-overflow"; out.reason="resource-limit:structural-construction:size-overflow"; return out;}
    for(std::size_t a=0;a<domain.atoms.size();++a)
        if(!CheckedAddBytes(memory_shape.component_support_entries,domain.atoms[a].size(),memory_shape.component_support_entries))
        {out.memory_estimate.stage=MemoryStage::StructuralConstruction; out.memory_estimate.status=MemoryStatus::Overflow;
            out.memory_estimate.reason="size-overflow"; out.reason="resource-limit:structural-construction:size-overflow"; return out;}
    const auto * input=domain.atoms.Storage();
    const auto snapshot_bytes=input ? EstimateSnapshotBytes(*input) : std::optional<std::size_t>{0};
    std::size_t base_resident=snapshot_bytes.value_or(0),root_mappings{},root_entries{},selector_entries{},twice_rows{};
    const auto local_support_bytes=domain.atoms.LocalSupportIndexBytes();
    if(!snapshot_bytes || local_support_bytes==std::numeric_limits<std::size_t>::max() ||
        (input && (!CheckedAddBytes(input->atom_ids.size(),input->observations.size(),root_entries) ||
            !AddMemoryBytes(root_mappings,root_entries,4*sizeof(Eigen::Index)) ||
            !CheckedAddBytes(base_resident,root_mappings,base_resident))) ||
        !AddMemoryBytes(base_resident,local_support_bytes,1) ||
        !CheckedMultiplyBytes(memory_shape.component_observations,2,twice_rows) ||
        !CheckedAddBytes(memory_shape.component_atoms,twice_rows,selector_entries) ||
        !AddMemoryBytes(base_resident,selector_entries,sizeof(Eigen::Index)))
    {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::StructuralConstruction;
        out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
        out.reason="resource-limit:structural-construction:size-overflow"; return out;}
    memory_shape.resident_bytes=base_resident;
    memory_shape.component_blocks=memory_shape.component_atoms;
    out.memory_estimate=PlanMemory(memory_shape,MemoryStage::StructuralConstruction);
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {LogMemoryEstimate(out.memory_estimate); out.reason="resource-limit:structural-construction:"+out.memory_estimate.reason; return out;}

    StructuralBlockPartition partition;
    StructuralBlockCensus census;
    std::size_t partition_bytes{};
    try
    {
        StructuralTopology topology(domain,context.atom_ids,layout);
        partition=topology.BuildCores(policy.core_atoms,false);
        if(partition.cores.empty()) {out.reason="block-invalid-partition"; return out;}
        census=topology.Census(partition);
        memory_shape.component_blocks=partition.cores.size();
        memory_shape.affected_row_entries=census.affected_row_entries;
        memory_shape.component_neighbor_entries=census.neighbor_atom_entries;
        memory_shape.maximum_block_rows=census.maximum_block_rows;
        memory_shape.maximum_block_support_entries=census.maximum_block_support_entries;
        memory_shape.total_block_support_entries=census.total_block_support_entries;
        memory_shape.maximum_neighbor_atoms=census.maximum_neighbor_atoms;
        std::size_t topology_bytes{},core_bytes{},planned_metadata{},planned_entries{};
        if(!CheckedAddBytes(census.affected_row_entries,census.neighbor_atom_entries,planned_entries) ||
            !CheckedMultiplyBytes(planned_entries,2*sizeof(Eigen::Index),planned_metadata) ||
            !AddMemoryBytes(planned_metadata,partition.cores.size(),4*32) ||
            !PartitionStorageBytes(partition,false,core_bytes) ||
            (topology_bytes=topology.StorageBytes())==std::numeric_limits<std::size_t>::max() ||
            !CheckedAddBytes(base_resident,topology_bytes,memory_shape.resident_bytes) ||
            !CheckedAddBytes(memory_shape.resident_bytes,core_bytes,memory_shape.resident_bytes) ||
            !CheckedAddBytes(memory_shape.resident_bytes,planned_metadata,memory_shape.resident_bytes))
        {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
            out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
            out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
        const auto structural_estimate=PlanMemory(memory_shape,MemoryStage::FixedNeighborPreparation);
        if(structural_estimate.status!=MemoryStatus::Pass)
        {out.memory_estimate=structural_estimate; LogMemoryEstimate(out.memory_estimate);
            out.reason="resource-limit:fixed-neighbor-preparation:"+out.memory_estimate.reason; return out;}
        topology.PopulateMetadata(partition);
        if(!PartitionStorageBytes(partition,true,partition_bytes))
        {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
            out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
            out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
    }
    catch(const std::overflow_error &)
    {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
        out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
        out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
    catch(const std::exception &) {out.reason="block-invalid-partition"; return out;}
    std::size_t resident{};
    if(!CheckedAddBytes(base_resident,partition_bytes,resident))
    {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
        out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
        out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
    memory_shape.resident_bytes=resident;
    out.memory_estimate=PlanMemory(memory_shape,MemoryStage::FixedNeighborPreparation);
    if(out.memory_estimate.status!=MemoryStatus::Pass ||
        out.memory_estimate.estimated_peak_bytes>out.memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.memory_estimate);
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {out.reason="resource-limit:fixed-neighbor-preparation:"+out.memory_estimate.reason; return out;}
    if(structural_timer.target)
    {
        out.work.structural_preparation_seconds+=Seconds(structural_timer.started);
        structural_timer.target=nullptr;
    }
    auto global_shape=memory_shape;
    auto global_resident=resident;
    if(!AddMemoryBytes(global_resident,memory_shape.component_support_entries,sizeof(LocalSupportIndex)) ||
        !AddMemoryBytes(global_resident,memory_shape.affected_row_entries,5*sizeof(Eigen::Index)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_atoms,4*sizeof(Eigen::Index)+3*sizeof(double)+
            sizeof(std::vector<LocalSupportIndex>)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_blocks,
            sizeof(PreparedFixedNeighborBlock)+sizeof(CompactRowMapping)+sizeof(std::size_t)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_observations,3*sizeof(double)+sizeof(Eigen::Index)))
    {
        out.global_state_memory_estimate.shape=memory_shape;
        out.global_state_memory_estimate.stage=MemoryStage::GlobalStateEvaluation;
        out.global_state_memory_estimate.status=MemoryStatus::Overflow;
        out.global_state_memory_estimate.reason="size-overflow";
        out.memory_estimate=out.global_state_memory_estimate;
        out.reason="resource-limit:global-state-evaluation:size-overflow"; return out;
    }
    global_shape.resident_bytes=global_resident;
    out.global_state_memory_estimate=PlanMemory(global_shape,MemoryStage::GlobalStateEvaluation);
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass ||
        out.global_state_memory_estimate.estimated_peak_bytes>out.global_state_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.global_state_memory_estimate);
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass)
    {
        out.memory_estimate=out.global_state_memory_estimate;
        out.reason="resource-limit:global-state-evaluation:"+out.memory_estimate.reason; return out;
    }
    auto * const profile_work=policy.collect_telemetry ? &out.work.local_profile_work : nullptr;
    const auto certify_endpoint=[&] {
        const auto assessment_timer=Measure(policy.collect_telemetry,out.work.endpoint_verification_seconds);
        const auto record_assessment_work=[&](FixedNeighborEndpointCertification certification) {
            if(policy.collect_telemetry)
            {
                out.work.primary_endpoint_seconds+=certification.assessment.primary_endpoint_seconds;
                out.work.independent_replay_seconds+=certification.assessment.independent_replay_seconds;
            }
            return certification;
        };
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.certification_for_testing)
            return record_assessment_work(policy.certification_for_testing(domain,y,out.state,context,profile_work));
#endif
        auto certification_resident=out.global_state_memory_estimate.shape.resident_bytes;
        if(!AddMemoryBytes(certification_resident,memory_shape.component_observations,3*sizeof(double)) ||
            !AddMemoryBytes(certification_resident,memory_shape.component_atoms,8*sizeof(double)))
            certification_resident=std::numeric_limits<std::size_t>::max();
        return record_assessment_work(CertifyFixedNeighborEndpoint(domain,y,out.state,context,profile_work,
            certification_resident,policy.convergence_criteria));
    };

    out.state.eta=initial_eta;
    out.state.beta=has_initial_beta ? std::move(initial_beta) :
        Vector::Zero(2*static_cast<Eigen::Index>(domain.atoms.size()));
    auto initial=Replay(domain,layout,y,out.state.eta,out.state.beta,context.scale);
    out.state.prediction=std::move(initial.prediction); out.state.residual=std::move(initial.residual);
    out.state.objective=initial.objective;
    auto mapping_timer=Measure(policy.collect_telemetry,out.work.structural_preparation_seconds);
    Indices atom_position(domain.atoms.size(),-1),row_position(static_cast<std::size_t>(domain.rows),-1);
    for(std::size_t k=0;k<layout.full_atoms.size();++k) atom_position.at(layout.full_atoms[k])=static_cast<Eigen::Index>(k);
    for(std::size_t k=0;k<layout.informative_rows.size();++k) row_position.at(layout.informative_rows[k])=static_cast<Eigen::Index>(k);
    std::vector<PreparedFixedNeighborBlock> prepared_blocks;
    prepared_blocks.reserve(partition.cores.size());
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
        prepared.row_mapping=CompactRowMapping::FromParentRows(prepared.profile_rows);
        prepared.domain=domain.Select(prepared.profile_atoms,
            static_cast<Eigen::Index>(prepared.profile_rows.size()),prepared.row_mapping);
        std::size_t mapping_entries{},mapping_bytes{};
        const auto support_entries=prepared.domain.atoms.LocalSupportIndexEntries();
        const auto support_bytes=prepared.domain.atoms.LocalSupportIndexBytes();
        if(!CheckedAddBytes(prepared.row_mapping->parent_rows.size(),prepared.row_mapping->local_rows.size(),mapping_entries) ||
            !CheckedAddBytes(mapping_entries,support_entries,mapping_entries) ||
            !CheckedAddBytes(prepared.row_mapping->StorageBytes(),support_bytes,mapping_bytes) ||
            !CheckedAddBytes(out.prepared_mapping_entries,mapping_entries,out.prepared_mapping_entries) ||
            !CheckedAddBytes(out.prepared_mapping_bytes,mapping_bytes,out.prepared_mapping_bytes))
        {out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
            out.reason="resource-limit:fixed-neighbor-mapping:size-overflow"; return out;}
        prepared.context=context;
        prepared.context.atom_ids=context.atom_ids.Select(prepared.profile_atoms);
        prepared.context.row_ids=context.row_ids.Select(prepared.profile_rows);
        prepared.context.independent_search=true;
        prepared.local_atoms=static_cast<Eigen::Index>(prepared.profile_atoms.size());
        prepared.context.rank={prepared.domain.rows,2*prepared.local_atoms,prepared.local_atoms};
        prepared.context.linear.rank_relative=prepared.context.rank.Relative(2*prepared.local_atoms);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(!policy.fixed_width_atom_ids.empty())
        {
            prepared.fixed_width_mask.resize(prepared.profile_atoms.size());
            for(std::size_t k=0;k<prepared.profile_atoms.size();++k)
                prepared.fixed_width_mask[k]=std::find(policy.fixed_width_atom_ids.begin(),
                    policy.fixed_width_atom_ids.end(),prepared.context.atom_ids[k])!=policy.fixed_width_atom_ids.end();
        }
#endif
        // FixedNeighbor performs at most one trusted accepted local update per block visit.
        prepared.context.update_budget=1;
        prepared_blocks.push_back(std::move(prepared));
    }
    out.block_preparations=prepared_blocks.size();
    out.domain_preparations=prepared_blocks.size();
    out.mapping_preparations=prepared_blocks.size();
    std::vector<std::size_t> order(partition.cores.size()); std::iota(order.begin(),order.end(),0);
    if(policy.order==FixedNeighborBlockOrder::Reverse) std::reverse(order.begin(),order.end());
    if(mapping_timer.target)
    {
        out.work.structural_preparation_seconds+=Seconds(mapping_timer.started);
        mapping_timer.target=nullptr;
    }
    if(policy.preparation_observer) policy.preparation_observer();
    const bool retain_diagnostics=policy.collect_diagnostics || policy.capture_local_trajectory ||
        static_cast<bool>(policy.state_observer);
    bool stop=false;
    std::size_t next_certification_retry_sweep{};
    std::vector<std::size_t> block_failure_visits(partition.cores.size());
    for(std::size_t sweep_index=0;sweep_index<policy.maximum_sweeps && !stop;++sweep_index)
    {
        const auto sweep_started=Clock::now(); FixedNeighborBlockSweep sweep; sweep.sweep=sweep_index+1;
        const Vector eta_before=out.state.eta,beta_before=out.state.beta;
        sweep.objective_before=out.state.objective;
        bool sweep_had_block_failure=false;
        auto finish_record=[&](FixedNeighborBlockRecord record,Clock::time_point block_started) {
            record.block_seconds=Seconds(block_started);
            if(policy.collect_telemetry)
            {
                ++out.work.block_visits;
                out.work.maximum_block_seconds=std::max(out.work.maximum_block_seconds,record.block_seconds);
            }
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
            const auto local_atoms=prepared.local_atoms; const auto block_started=Clock::now(); FixedNeighborBlockRecord record;
            record.sweep=sweep_index+1; record.block=block_index+1; record.core_atoms=core_atoms;
            record.affected_rows=affected_rows.size(); record.objective_before=out.state.objective;
            record.recovery_retry=block_failure_visits[block_index]>0;
            if(record.recovery_retry)
            {
                record.recovery_attempt=block_failure_visits[block_index];
                ++out.recovery_attempts;
                ++sweep.recovery_attempts;
            }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            record.component_index=policy.component_index; record.component_id=policy.component_id;
#endif

            Vector old_eta(local_atoms),old_beta(2*local_atoms),old_core;
            {
                const auto old_core_timer=Measure(policy.collect_telemetry,out.work.old_core_seconds);
                if(policy.collect_telemetry) ++out.work.old_core_basis_builds;
                for(std::size_t k=0;k<core_atoms.size();++k)
                {
                    const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                    old_eta(static_cast<Eigen::Index>(k))=out.state.eta(atom);
                    old_beta.segment<2>(2*static_cast<Eigen::Index>(k))=out.state.beta.segment<2>(2*atom);
                }
                const Sparse old_design=Design(local_domain,old_eta);
                old_core=(old_design*old_beta).eval();
                (void)old_core_timer;
            }

            Vector local_y(static_cast<Eigen::Index>(affected_rows.size()));
            {
                const auto response_timer=Measure(policy.collect_telemetry,out.work.effective_response_seconds);
                double local_before_squared{};
                for(std::size_t k=0;k<affected_rows.size();++k)
                {
                    const auto row=static_cast<Eigen::Index>(affected_rows[k]);
                    local_y(static_cast<Eigen::Index>(k))=old_core(static_cast<Eigen::Index>(k))-out.state.residual(row);
                    local_before_squared+=out.state.residual(row)*out.state.residual(row);
                }
                record.local_objective_before=.5*local_before_squared/(context.scale*context.scale);
                (void)response_timer;
            }

            const auto search_started=Clock::now();
            if(policy.collect_telemetry)
            {
                const auto atoms=static_cast<std::size_t>(local_atoms);
                ++out.work.local_problem_count;
                out.work.total_local_atoms+=static_cast<double>(atoms);
                if(out.work.local_problem_count==1)
                    out.work.minimum_local_atoms=out.work.maximum_local_atoms=atoms;
                else
                {
                    out.work.minimum_local_atoms=std::min(out.work.minimum_local_atoms,atoms);
                    out.work.maximum_local_atoms=std::max(out.work.maximum_local_atoms,atoms);
                }
            }
            const auto old_widths=old_eta.array().exp().eval();
            const auto profile_search=SearchProfile(local_domain,local_y,old_widths,local_context,{},nullptr,profile_work,
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                prepared.fixed_width_mask.empty() ? nullptr : &prepared.fixed_width_mask
#else
                nullptr
#endif
            );
            if(policy.collect_telemetry) out.work.local_profile_work.Merge(profile_search.profile_work);
            record.search_seconds=Seconds(search_started); record.profile_evaluations=profile_search.evaluations;
            if(policy.collect_telemetry) out.work.local_search_seconds+=record.search_seconds;
            record.accepted_updates=profile_search.accepted; record.local_search_stop_reason=profile_search.stop_reason;
            if(policy.capture_local_trajectory)
            {
                record.profile_trials.reserve(profile_search.trials.size());
                double accepted_objective=record.local_objective_before;
                Vector accepted_eta=old_eta; double cumulative_factor_seconds{};
                for(std::size_t trial_index=0;trial_index<profile_search.trials.size();++trial_index)
                {
                    const auto & trial=profile_search.trials[trial_index];
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(policy.width_trajectory)
                CaptureWidthTrajectory(record,profile_search,local_domain,local_y,old_eta,local_context,
                    *policy.width_trajectory);
#endif
            ++sweep.block_solves; ++out.total_block_solves;
            if(policy.collect_telemetry) ++out.work.block_solves;
            sweep.profile_evaluations+=static_cast<std::size_t>(profile_search.evaluations);
            out.total_profile_evaluations+=static_cast<std::size_t>(profile_search.evaluations);
            sweep.maximum_block_rows=std::max(sweep.maximum_block_rows,static_cast<std::size_t>(local_domain.rows));
            sweep.maximum_block_columns=std::max(sweep.maximum_block_columns,static_cast<std::size_t>(2*local_atoms));
            const Trial * accepted=nullptr;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const bool force_no_trusted_profile=policy.force_no_trusted_profile_for_testing &&
                policy.force_no_trusted_profile_for_testing(block_index,sweep_index);
#else
            const bool force_no_trusted_profile=false;
#endif
            if(!force_no_trusted_profile)
                for(auto it=profile_search.trials.rbegin();it!=profile_search.trials.rend();++it)
                    if(it->accepted && it->endpoint.valid && it->trust && it->trust->passed) {accepted=&*it; break;}
            if(!accepted)
            {
                record.failure_diagnostics=FailureDiagnostics(policy,record,sweep_index,core_atoms,
                    context,local_context,local_domain,local_y,old_eta,profile_search);
                Logger::Log(LogLevel::Warning,FormatFailureDiagnostics(*record.failure_diagnostics));
                ++block_failure_visits[block_index]; sweep_had_block_failure=true;
                if(block_failure_visits[block_index]>=MaximumBlockFailureVisits)
                {
                    record.status="failed"; record.reason="block-recovery-exhausted";
                    out.reason=record.reason; stop=true;
                }
                else
                {
                    record.status="retry-pending"; record.reason="block-retry-pending";
                }
                finish_record(std::move(record),block_started);
                if(stop) break;
                continue;
            }
            const auto local_state_started=Clock::now();
            const auto local_state=EvaluateState(local_domain,local_y,accepted->endpoint.eta,accepted->endpoint.beta,local_context,
                ProfileEvaluationRole::AcceptedEndpoint,profile_work);
            if(policy.collect_telemetry) out.work.local_state_seconds+=Seconds(local_state_started);
            if(!local_state.valid || !local_state.certificate.available || !local_state.certificate.feasible ||
                !local_state.certificate.kkt_passed)
            {
                record.status="failed"; record.reason="block-inner-invalid"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            record.local_final_gradient_inf_norm=local_state.gradient.lpNorm<Eigen::Infinity>();
            record.local_final_ac_kkt=local_state.certificate.projected_kkt;
            record.local_objective_after=local_state.certificate.objective/(context.scale*context.scale);
            if(record.local_objective_after>record.local_objective_before)
            {
                record.status=record.recovery_retry ? "failed" : "unchanged";
                record.reason=record.recovery_retry ? "block-recovery-exhausted" : "block-objective-increase";
                record.objective_after=record.objective_before;
                if(record.recovery_retry) {out.reason=record.reason; stop=true;}
                finish_record(std::move(record),block_started);
                if(stop) break;
                continue;
            }
            Vector candidate_eta,candidate_beta;
            {
                const auto candidate_copy_timer=Measure(policy.collect_telemetry,out.work.candidate_copy_seconds);
                if(policy.collect_telemetry) ++out.work.candidate_state_full_copies;
                candidate_eta=out.state.eta; candidate_beta=out.state.beta;
                for(std::size_t k=0;k<core_atoms.size();++k)
                {
                    const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                    candidate_eta(atom)=accepted->endpoint.eta(static_cast<Eigen::Index>(k));
                    candidate_beta.segment<2>(2*atom)=accepted->endpoint.beta.segment<2>(2*static_cast<Eigen::Index>(k));
                }
                (void)candidate_copy_timer;
            }
            std::optional<IncrementalCandidateReplay> incremental;
            if(policy.compare_incremental_replay)
            {
                const auto incremental_started=Clock::now();
                const Vector new_core=(local_state.x*local_state.beta).eval();
                incremental=EvaluateCandidateDelta(out.state.residual,affected_rows,old_core,new_core,
                    out.state.objective,context.scale);
                out.work.incremental_candidate_replay_seconds+=Seconds(incremental_started);
                ++out.work.incremental_candidate_replays;
            }
            ReplayState replay;
            {
                const auto candidate_replay_timer=Measure(policy.collect_telemetry,out.work.candidate_replay_seconds);
                if(policy.collect_telemetry) ++out.work.full_candidate_replays;
                replay=Replay(domain,layout,y,candidate_eta,candidate_beta,context.scale);
                (void)candidate_replay_timer;
            }
            record.objective_after=replay.objective; record.global_replay_delta=replay.objective-out.state.objective;
            const auto replay_reference=std::max(std::abs(replay.objective),std::abs(out.state.objective));
            record.objective_replay_enclosure=BlockObjectiveReplayEnclosure(replay_reference);
            record.local_global_delta_error=std::abs(record.global_replay_delta+
                record.local_objective_before-record.local_objective_after);
            const bool delta_enclosed=WithinBlockObjectiveReplay(record.local_global_delta_error,replay_reference);
            const bool candidate_unchanged=(out.state.eta.array()==candidate_eta.array()).all() &&
                (out.state.beta.array()==candidate_beta.array()).all();
            if(incremental)
            {
                const double objective_difference=std::abs(incremental->objective-replay.objective);
                const double relative_scale=std::max({std::abs(incremental->objective),std::abs(replay.objective),
                    std::numeric_limits<double>::min()});
                std::size_t affected_index{};
                double residual_difference{};
                for(auto row:layout.informative_rows)
                {
                    double change{};
                    if(affected_index<affected_rows.size() &&
                        affected_rows[affected_index]==static_cast<Eigen::Index>(row))
                    {
                        change=incremental->residual_delta(static_cast<Eigen::Index>(affected_index));
                        ++affected_index;
                    }
                    residual_difference=std::max(residual_difference,
                        std::abs(out.state.residual(static_cast<Eigen::Index>(row))+change-
                            replay.residual(static_cast<Eigen::Index>(row))));
                }
                const double incremental_reference=std::max(std::abs(incremental->objective),std::abs(out.state.objective));
                const double incremental_error=std::abs(incremental->objective_delta+
                    record.local_objective_before-record.local_objective_after);
                const double incremental_enclosure=BlockObjectiveReplayEnclosure(incremental_reference);
                const bool incremental_enclosed=WithinBlockObjectiveReplay(incremental_error,incremental_reference);
                const bool full_accept=delta_enclosed && !candidate_unchanged &&
                    !(record.global_replay_delta>0 && !(record.local_objective_after<=record.local_objective_before &&
                        record.global_replay_delta<=record.objective_replay_enclosure));
                const bool incremental_accept=incremental_enclosed && !candidate_unchanged &&
                    !(incremental->objective_delta>0 && !(record.local_objective_after<=record.local_objective_before &&
                        incremental->objective_delta<=incremental_enclosure));
                if(full_accept!=incremental_accept) ++out.work.incremental_acceptance_mismatches;
                out.work.maximum_incremental_objective_difference=std::max(
                    out.work.maximum_incremental_objective_difference,objective_difference);
                out.work.maximum_incremental_relative_objective_difference=std::max(
                    out.work.maximum_incremental_relative_objective_difference,objective_difference/relative_scale);
                out.work.maximum_incremental_residual_difference=std::max(
                    out.work.maximum_incremental_residual_difference,residual_difference);
                out.work.maximum_incremental_delta_condition=std::max(
                    out.work.maximum_incremental_delta_condition,incremental->delta_condition);
                out.work.last_incremental_candidate_objective=incremental->objective;
                out.work.last_full_candidate_objective=replay.objective;
            }
            if(!delta_enclosed)
            {
                record.status="failed"; record.reason="block-objective-replay-failed"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            if(record.global_replay_delta>0 && !(record.local_objective_after<=record.local_objective_before &&
                record.global_replay_delta<=record.objective_replay_enclosure))
            {
                record.status=record.recovery_retry ? "failed" : "unchanged";
                record.reason=record.recovery_retry ? "block-recovery-exhausted" : "block-objective-increase";
                if(record.recovery_retry) {out.reason=record.reason; stop=true;}
                finish_record(std::move(record),block_started);
                if(stop) break;
                continue;
            }
            const bool unchanged=candidate_unchanged;
            if(!unchanged)
                for(const auto atom:core_atoms)
                    sweep.maximum_accepted_log_width_step=std::max(sweep.maximum_accepted_log_width_step,
                        std::abs(candidate_eta(static_cast<Eigen::Index>(atom))-
                            out.state.eta(static_cast<Eigen::Index>(atom))));
            const auto cache_started=Clock::now();
            out.state.prediction=std::move(replay.prediction);
            out.state.residual=std::move(replay.residual);
            out.state.objective=replay.objective;
            if(policy.collect_telemetry) out.work.affected_row_updates+=affected_rows.size();
            out.state.eta=std::move(candidate_eta); out.state.beta=std::move(candidate_beta);
            if(policy.collect_telemetry) out.work.cache_update_seconds+=Seconds(cache_started);
            record.accepted=!unchanged; record.status=unchanged ? "unchanged" : "accepted";
            record.reason=unchanged ? "block-unchanged" : "";
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(policy.width_trajectory && !unchanged)
                for(auto & trial:record.profile_trials) if(trial.profile_evaluation==accepted->evaluation)
                    trial.accepted_global_objective_improvement=record.objective_before-record.objective_after;
#endif
            if(record.recovery_retry)
            {
                ++out.recovery_successes;
                ++sweep.recovery_successes;
                block_failure_visits[block_index]=0;
                record.recovery_success=true;
            }
            finish_record(std::move(record),block_started);
        }
        if(stop) break;
        sweep.objective_after=out.state.objective;
        ReplayState replay;
        {
            const auto timer=Measure(policy.collect_telemetry,out.work.sweep_replay_seconds);
            replay=Replay(domain,layout,y,out.state.eta,out.state.beta,context.scale);
            (void)timer;
        }
        sweep.cache_replay_error=PredictionDifference(out.state.prediction,replay.prediction);
        sweep.objective_replay_error=std::abs(out.state.objective-replay.objective);
        if(policy.collect_telemetry)
            out.work.maximum_cache_replay_error=std::max(out.work.maximum_cache_replay_error,sweep.cache_replay_error);
        const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
        const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);
        const auto global_started=Clock::now();
        const auto global=EvaluateStateStreaming(domain,y,eta,beta,context,ProfileEvaluationRole::AcceptedEndpoint,profile_work);
        if(policy.collect_telemetry) out.work.sweep_global_state_seconds+=Seconds(global_started);
        if(!global.valid || !global.certificate.available)
        {out.reason="block-inner-invalid"; stop=true; break;}
        sweep.global_a_feasibility=global.certificate.feasible ? 0.0 : 1.0;
        sweep.global_ac_kkt=global.certificate.projected_kkt;
        sweep.global_width_gradient_inf_norm=global.gradient.lpNorm<Eigen::Infinity>();
        sweep.eta_change_inf=(out.state.eta-eta_before).lpNorm<Eigen::Infinity>();
        sweep.beta_scaled_change=ScaledCoefficientDifference(out.state.beta,beta_before);
        sweep.maximum_target_width=out.state.eta.size() ? std::exp(out.state.eta.maxCoeff()) : 0.0;
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.width_trajectory && policy.width_trajectory->state_observer)
            policy.width_trajectory->state_observer(sweep_index+1,out.state,sweep);
#endif
        if(sweep.cache_replay_error>2e-12+2e-13*std::max(1.0,replay.prediction.cwiseAbs().maxCoeff()) ||
            !WithinBlockObjectiveReplay(sweep.objective_replay_error,replay.objective))
        {out.reason="block-cache-replay-failed"; stop=true; break;}
        out.recent_sweep_changes.emplace_back(sweep.eta_change_inf,sweep.beta_scaled_change);
        if(out.recent_sweep_changes.size()>3) out.recent_sweep_changes.erase(out.recent_sweep_changes.begin());
        double stationarity_width_gradient=sweep.global_width_gradient_inf_norm;
        double stationarity_eta_change=sweep.eta_change_inf;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(!policy.fixed_width_atom_ids.empty())
        {
            stationarity_width_gradient=0.0; stationarity_eta_change=0.0;
            const auto count=std::min({static_cast<std::size_t>(global.gradient.size()),
                static_cast<std::size_t>(out.state.eta.size()),context.atom_ids.size()});
            for(std::size_t atom=0;atom<count;++atom)
            {
                const auto & id=context.atom_ids[atom];
                if(std::find(policy.fixed_width_atom_ids.begin(),policy.fixed_width_atom_ids.end(),id)!=
                    policy.fixed_width_atom_ids.end()) continue;
                stationarity_width_gradient=std::max(stationarity_width_gradient,
                    std::abs(global.gradient(static_cast<Eigen::Index>(atom))));
                stationarity_eta_change=std::max(stationarity_eta_change,
                    std::abs(out.state.eta(static_cast<Eigen::Index>(atom))-
                        eta_before(static_cast<Eigen::Index>(atom))));
            }
        }
#endif
        const auto stationarity=sweep_had_block_failure ? FixedNeighborStationarityState::NotStationary :
            ClassifyFixedNeighborStationarity(sweep.global_ac_kkt,stationarity_width_gradient,
                stationarity_eta_change,sweep.coordinate_confirmation_available,policy.convergence_criteria);
        if(stationarity!=FixedNeighborStationarityState::NotStationary && out.first_order_stationarity_sweep==0)
            out.first_order_stationarity_sweep=sweep_index+1;
        const bool certification_retry_due=ShouldAttemptFixedNeighborCertification(
            stationarity,sweep_index+1,next_certification_retry_sweep);
        if(stationarity==FixedNeighborStationarityState::ConfirmedStationary)
        {
            if(out.confirmed_stationarity_sweep==0) out.confirmed_stationarity_sweep=sweep_index+1;
            if(!policy.stop_after_stationarity || !policy.certify_endpoint)
            {
                out.search_converged=true; out.reason="block-stationary";
                if(policy.stop_after_stationarity) break;
                continue;
            }
            const bool retry_due=!out.certification || certification_retry_due;
            if(!retry_due) continue;
            out.certification=certify_endpoint();
            out.certification_sweep=sweep_index+1;
            if(out.certification->qualified())
            {
                out.search_converged=true; out.reason="block-stationary";
                break;
            }
            if(out.certification->resource_limited())
            {
                out.reason=out.certification->resource_limit_reason();
                stop=true;
            }
            else if(!out.certification->trust.primary_valid || !out.certification->assessment.primary.valid)
            {
                out.reason=!out.certification->trust.reason.empty() ? out.certification->trust.reason :
                    out.certification->assessment.failure.empty() ? "endpoint-certification-invalid-primary" :
                        out.certification->assessment.failure;
                stop=true;
            }
        }
    }
    if(!out.search_converged && !stop)
    {
        const bool recovery_pending=std::any_of(block_failure_visits.begin(),block_failure_visits.end(),
            [](std::size_t visits){return visits>0;});
        out.reason=recovery_pending ? "block-recovery-sweep-budget" : "block-sweep-budget";
        if(!recovery_pending && policy.stop_after_stationarity && policy.certify_endpoint && out.sweep_count &&
            out.certification_sweep!=out.sweep_count)
        {
            out.certification=certify_endpoint();
            out.certification_sweep=out.sweep_count;
            if(out.certification->qualified() && out.final_sweep &&
                ClassifyFixedNeighborStationarity(out.final_sweep->global_ac_kkt,
                    out.final_sweep->global_width_gradient_inf_norm,out.final_sweep->eta_change_inf,
                    out.final_sweep->coordinate_confirmation_available,policy.convergence_criteria)==
                        FixedNeighborStationarityState::ConfirmedStationary)
            {
                out.search_converged=true; out.reason="block-stationary";
            }
            else if(out.certification->resource_limited()) out.reason=out.certification->resource_limit_reason();
            else if(!out.certification->trust.primary_valid || !out.certification->assessment.primary.valid)
                out.reason=!out.certification->trust.reason.empty() ? out.certification->trust.reason :
                    out.certification->assessment.failure.empty() ? "endpoint-certification-invalid-primary" :
                        out.certification->assessment.failure;
        }
    }
    if(policy.assess_local_conditioning && out.search_converged && out.certification &&
        out.certification->qualified())
    {
        out.local_conditioning.reserve(prepared_blocks.size());
        for(const auto & prepared:prepared_blocks)
        {
            const auto local_atoms=prepared.local_atoms;
            Vector eta(local_atoms),beta(2*local_atoms);
            for(std::size_t k=0;k<prepared.core_atoms.size();++k)
            {
                const auto atom=static_cast<Eigen::Index>(prepared.core_atoms[k]);
                eta(static_cast<Eigen::Index>(k))=out.state.eta(atom);
                beta.segment<2>(2*static_cast<Eigen::Index>(k))=out.state.beta.segment<2>(2*atom);
            }
            Vector local_y(static_cast<Eigen::Index>(prepared.affected_rows.size()));
            {
                const Sparse design=Design(prepared.domain,eta);
                const Vector local_core=(design*beta).eval();
                for(std::size_t k=0;k<prepared.affected_rows.size();++k)
                {
                    const auto row=static_cast<Eigen::Index>(prepared.affected_rows[k]);
                    local_y(static_cast<Eigen::Index>(k))=local_core(static_cast<Eigen::Index>(k))-out.state.residual(row);
                }
            }
            const auto block_id=(policy.component_id.empty() ? std::string("component") : policy.component_id)+
                "/block-"+std::to_string(prepared.block+1);
            auto conditioning=AssessFixedNeighborLocalConditioning(prepared.domain,local_y,eta,beta,
                prepared.context,block_id,out.global_state_memory_estimate.shape.resident_bytes);
            conditioning.evidence.core_atoms.assign(prepared.core_atoms.begin(),prepared.core_atoms.end());
            out.work.local_conditioning_seconds+=conditioning.seconds;
            if(conditioning.evidence.workspace_peak_memory_estimate_bytes)
                out.work.local_conditioning_peak_memory_estimate_bytes=std::max(
                    out.work.local_conditioning_peak_memory_estimate_bytes,
                    *conditioning.evidence.workspace_peak_memory_estimate_bytes);
            out.local_conditioning.push_back(std::move(conditioning.evidence));
        }
    }
    if(policy.assess_parameter_stability && out.search_converged && out.certification &&
        out.certification->qualified())
    {
        const auto started=Clock::now(); JointParameterStability stability;
        stability.recent_sweep_count=out.recent_sweep_changes.size();
        if(stability.recent_sweep_count)
        {
            double maximum_eta{},maximum_ac{};
            for(const auto & [eta_change,ac_change]:out.recent_sweep_changes)
            {maximum_eta=std::max(maximum_eta,eta_change); maximum_ac=std::max(maximum_ac,ac_change);}
            stability.maximum_eta_change=maximum_eta;
            stability.maximum_scaled_ac_change=maximum_ac;
        }
        stability.trajectory_stability_status=stability.recent_sweep_count>=3 ? "observed" : "insufficient-history";
        std::vector<JointBlockParameterStability> block_evidence(prepared_blocks.size());
        for(std::size_t index=0;index<prepared_blocks.size();++index)
        {
            const auto & prepared=prepared_blocks[index];
            auto & evidence=block_evidence[index];
            evidence.block_id=(policy.component_id.empty() ? std::string("component") : policy.component_id)+
                "/block-"+std::to_string(prepared.block+1);
            evidence.core_atoms.assign(prepared.core_atoms.begin(),prepared.core_atoms.end());
            if(index>=out.local_conditioning.size())
                evidence.reason="condition-data-unavailable";
            else
            {
                const auto & condition=out.local_conditioning[index];
                const bool free_ac_singular=condition.free_ac_rank && *condition.free_ac_rank<condition.free_ac_columns;
                if(condition.conditioning_classification=="Singular" && free_ac_singular)
                    evidence.reason="free-ac-rank-deficient-no-width-probe";
                else if(StabilityProbePriority(condition)<0)
                    evidence.reason="condition-data-unavailable";
                else evidence.reason="not-selected-lower-sensitivity";
            }
        }
        std::vector<std::size_t> selected;
        std::size_t representative_count{};
        for(const auto index:StabilityProbeOrder(out.local_conditioning))
        {
            const auto & condition=out.local_conditioning[index];
            const auto priority=StabilityProbePriority(condition);
            if(priority==3 && representative_count>=MaximumRepresentativeLocalStabilityProbesPerComponent)
            {
                block_evidence[index].reason="representative-probe-limit";
                continue;
            }
            if(selected.size()>=MaximumPerturbedBlocksPerFit)
            {
                block_evidence[index].reason="probe-block-budget-exhausted";
                continue;
            }
            selected.push_back(index);
            if(priority==3) ++representative_count;
        }
        bool budget_exhausted=false;
        for(const auto index:selected)
        {
            auto & evidence=block_evidence[index];
            if(budget_exhausted || !policy.claim_parameter_stability_block || !policy.claim_parameter_stability_block())
            {
                evidence.reason="perturbation-budget-exhausted";
                budget_exhausted=true;
                continue;
            }
            const auto & prepared=prepared_blocks[index];
            const auto block_id=evidence.block_id;
            const auto local_atoms=prepared.local_atoms;
            Vector eta(local_atoms),beta(2*local_atoms);
            for(std::size_t k=0;k<prepared.core_atoms.size();++k)
            {
                const auto atom=static_cast<Eigen::Index>(prepared.core_atoms[k]);
                eta(static_cast<Eigen::Index>(k))=out.state.eta(atom);
                beta.segment<2>(2*static_cast<Eigen::Index>(k))=out.state.beta.segment<2>(2*atom);
            }
            Vector local_y(static_cast<Eigen::Index>(prepared.affected_rows.size()));
            {
                const Sparse design=Design(prepared.domain,eta);
                const Vector local_core=(design*beta).eval();
                for(std::size_t k=0;k<prepared.affected_rows.size();++k)
                {
                    const auto row=static_cast<Eigen::Index>(prepared.affected_rows[k]);
                    local_y(static_cast<Eigen::Index>(k))=local_core(static_cast<Eigen::Index>(k))-out.state.residual(row);
                }
            }
            auto assessment=AssessLocalParameterStability(prepared.domain,local_y,eta,beta,prepared.context,
                block_id,out.global_state_memory_estimate.shape.resident_bytes);
            evidence=std::move(assessment.evidence);
            evidence.core_atoms.assign(prepared.core_atoms.begin(),prepared.core_atoms.end());
            stability.perturbation_directions+=evidence.perturbation_directions;
            stability.perturbation_seeds+=evidence.perturbation_seeds;
            stability.perturbation_solves+=evidence.perturbation_solves;
            stability.perturbation_evaluations+=evidence.perturbation_evaluations;
            stability.cpu_seconds+=evidence.cpu_seconds;
            if(assessment.memory_estimate.status!=MemoryStatus::Overflow)
                stability.workspace_peak_memory_estimate_bytes=std::max(
                    stability.workspace_peak_memory_estimate_bytes.value_or(0),
                    assessment.memory_estimate.estimated_peak_bytes);
            out.work.perturbation_solves+=evidence.perturbation_solves;
            out.work.perturbation_evaluations+=evidence.perturbation_evaluations;
            if(assessment.memory_estimate.status!=MemoryStatus::Overflow)
                out.work.parameter_stability_peak_memory_estimate_bytes=std::max(
                    out.work.parameter_stability_peak_memory_estimate_bytes,
                    assessment.memory_estimate.estimated_peak_bytes);
            evidence.block_id=block_id;
        }
        stability.blocks=std::move(block_evidence);
        int worst{}; bool all_assessed=true;
        for(const auto & block:stability.blocks)
        {
            all_assessed&=block.status!="NotAssessed";
            const int severity=block.status=="Unstable" ? 4 : block.status=="Sensitive" ? 3 :
                block.status=="Unavailable" ? 2 : block.status=="Stable" ? 1 : 0;
            worst=std::max(worst,severity);
        }
        stability.perturbation_status=worst==4 ? "Unstable" : worst==3 ? "Sensitive" :
            worst==2 ? "Unavailable" : worst==1 && all_assessed ? "Stable" : "NotAssessed";
        stability.seconds=Seconds(started);
        out.work.parameter_stability_seconds+=stability.seconds;
        out.parameter_stability=std::move(stability);
    }
    return out;
}

FixedNeighborEndpointCertification CertifyFixedNeighborEndpoint(
    const Domain & domain,VectorRef y,const BlockCoordinateState & state,const EvaluationContext & context,
    ProfileSearchWork * profile_work,std::size_t resident_override,
    const JointConvergenceCriteria & convergence_criteria)
{
    FixedNeighborEndpointCertification out;
    const auto atom_count=static_cast<Eigen::Index>(domain.atoms.size());
    if(state.eta.size()!=atom_count || state.beta.size()!=2*atom_count || y.size()!=domain.rows ||
        !state.eta.allFinite() || !state.beta.allFinite())
    {
        out.trust.reason="invalid-primary";
        return out;
    }
    for(Eigen::Index atom=0;atom<atom_count;++atom)
    {
        const double width=std::exp(state.eta(atom));
        if(state.beta(2*atom)<0 || !std::isfinite(width) || !(width>0))
        {
            out.trust.reason="invalid-primary";
            out.assessment.failure="invalid-accepted-parameter-state";
            out.assessment.primary.reason=out.assessment.failure;
            return out;
        }
    }
    out.verification_memory_estimate=EndpointVerificationMemoryEstimate(domain,y,context,resident_override);
    if(out.verification_memory_estimate.status!=MemoryStatus::Pass)
    {
        LogMemoryEstimate(out.verification_memory_estimate);
        const auto reason=std::string("resource-limit:endpoint-verification:")+out.verification_memory_estimate.reason;
        out.assessment.primary.eta=state.eta; out.assessment.primary.reason=reason;
        out.assessment.failure=reason;
        out.assessment.verification_memory_estimate=out.verification_memory_estimate;
        out.trust.reason=reason;
        return out;
    }
    if(out.verification_memory_estimate.estimated_peak_bytes>out.verification_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.verification_memory_estimate);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting primary_stage("endpoint-primary-evaluation",y.size(),state.eta.size(),
        &out.verification_memory_estimate);
#endif
    const auto primary_started=Clock::now();
    const auto endpoint=EvaluateStateStreaming(domain,y,state.eta,state.beta,context,
        ProfileEvaluationRole::AcceptedEndpoint,profile_work);
    const double primary_seconds=Seconds(primary_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    primary_stage.Finish();
    AssessmentStageTimerForTesting replay_stage("independent-streaming-replay",y.size(),state.eta.size(),
        &out.verification_memory_estimate);
#endif
    const auto replay_started=Clock::now();
    const auto replay=CheckIndependentReplay(domain,y,endpoint,context,nullptr,nullptr,convergence_criteria);
    const double replay_seconds=Seconds(replay_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    replay_stage.Finish();
#endif
    out.assessment=VerifyJointEndpoint(domain,y,endpoint,replay,context,convergence_criteria);
    out.assessment.verification_memory_estimate=out.verification_memory_estimate;
    out.assessment.primary_endpoint_seconds=primary_seconds;
    out.assessment.independent_replay_seconds=replay_seconds;
    out.trust=CheckTrust(domain,y,endpoint,context,replay);
    return out;
}

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

FixedNeighborComponentResult SolveFixedNeighborComponent(
    const ComponentProblemView & prepared,VectorRef initial_b,const EvaluationContext & parent_context,
    const FixedNeighborSearchPolicy & production_policy,const JointProgressObserver & observer,
    const JointProgressComponent * progress_component)
{
    FixedNeighborComponentResult out;
    const auto & local_layout=prepared.local_layout;
    if(local_layout.full_atoms.empty() || local_layout.informative_rows.empty())
    {
        out.stop_reason="unobserved-full-parameters";
        return out;
    }
    for(const auto atom:prepared.local_to_parent_atoms)
        if(atom<0 || atom>=initial_b.size() || !std::isfinite(initial_b(atom)) || !(initial_b(atom)>0))
        {
            out.stop_reason="invalid-initial-widths"; return out;
        }

    Vector local_y(static_cast<Eigen::Index>(prepared.local_to_parent_rows.size()));
    for(std::size_t k=0;k<prepared.local_to_parent_rows.size();++k)
        local_y(static_cast<Eigen::Index>(k))=(*parent_context.observations)(prepared.local_to_parent_rows[k]);
    auto context=parent_context;
    context.atom_ids=prepared.atom_ids; context.row_ids=prepared.row_ids;
    context.independent_search=true;
    const auto atom_count=static_cast<Eigen::Index>(local_layout.full_atoms.size());
    context.rank={prepared.parent_row_count,2*atom_count,atom_count};
    context.linear.rank_relative=context.rank.Relative(2*atom_count);

    Vector initial_eta(atom_count);
    for(std::size_t k=0;k<prepared.local_to_parent_atoms.size();++k)
        initial_eta(static_cast<Eigen::Index>(k))=std::log(initial_b(prepared.local_to_parent_atoms[k]));

    FixedNeighborPolicy policy;
    policy.core_atoms=production_policy.core_atoms;
    policy.maximum_sweeps=production_policy.maximum_sweeps;
    policy.convergence_criteria=production_policy.convergence_criteria;
    policy.order=production_policy.order;
    policy.collect_telemetry=production_policy.collect_telemetry;
    policy.component_index=production_policy.component_index;
    policy.component_id=production_policy.component_id;
    policy.halo_census=production_policy.halo_census;
    policy.assess_local_conditioning=true;
    policy.assess_parameter_stability=true;
    policy.claim_parameter_stability_block=production_policy.claim_parameter_stability_block;
    if(production_policy.state_observer)
        policy.state_observer=[&production_policy](std::size_t sweep,const BlockCoordinateState & state,
            const FixedNeighborBlockSweep & record,const std::vector<FixedNeighborBlockRecord> &) {
            production_policy.state_observer(production_policy.component_index,sweep,state,record);
        };
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(production_policy.configure_search_for_testing) production_policy.configure_search_for_testing(policy);
#endif
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
            JointFixedNeighborProgress progress;
            progress.sweep=sweep.sweep;
            progress.maximum_sweeps=production_policy.maximum_sweeps;
            progress.block_solves=sweep.block_solves;
            progress.accepted_blocks=sweep.accepted_blocks;
            progress.local_accepted_updates=sweep.accepted_local_updates;
            progress.profile_evaluations=sweep.profile_evaluations;
            progress.unchanged_blocks=sweep.unchanged_blocks;
            progress.failed_blocks=sweep.block_solves-sweep.accepted_blocks-sweep.unchanged_blocks;
            progress.recovery_attempts=sweep.recovery_attempts;
            progress.recovery_successes=sweep.recovery_successes;
            progress.objective=sweep.objective_after;
            progress.global_ac_kkt=sweep.global_ac_kkt;
            progress.width_gradient_inf_norm=sweep.global_width_gradient_inf_norm;
            progress.eta_change_inf=sweep.eta_change_inf;
            progress.maximum_accepted_log_width_step=sweep.maximum_accepted_log_width_step;
            progress.maximum_target_width=sweep.maximum_target_width;
            event.fixed_neighbor=std::move(progress);
            observer(event);
        };
    }
    out.search=SearchFixedNeighborComponent(prepared.domain,local_layout,local_y,initial_eta,context,policy);
    const bool qualified=out.search.certification && out.search.certification->qualified();
    out.succeeded=out.search.search_converged && qualified;
    out.stop_reason=out.search.reason;
    if(progress_component)
        NotifyJointProgress(observer,JointProgressPhase::CertificationStarted,*progress_component,
            static_cast<int>(out.search.total_profile_evaluations),parent_context.profile_budget,
            static_cast<int>(out.search.total_accepted_local_updates),parent_context.update_budget,
            out.search.search_seconds,out.stop_reason,qualified,
            LastFixedNeighborObjective(out.search),LastFixedNeighborGradient(out.search));
    return out;
}

}
