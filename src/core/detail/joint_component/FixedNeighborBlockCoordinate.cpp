#include "FixedNeighborBlockCoordinate.hpp"
#include "LocalConditioning.hpp"
#include "ParameterStability.hpp"
#include "Problem.hpp"
#include "ResourceWork.hpp"
#include "StructuralPartition.hpp"
#include "TiledDerivative.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>

namespace rhbm_gem::core::joint_component {
double BlockObjectiveReplayEnclosure(double reference)
{return 1e-12+2e-12*std::abs(reference);}
bool WithinBlockObjectiveReplay(double error,double reference)
{return std::isfinite(error) && std::isfinite(reference) && error<=BlockObjectiveReplayEnclosure(reference);}
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
            !AddMemoryBytes(bytes,core.neighbor_atoms.capacity(),sizeof(Eigen::Index)))) return false;
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
struct ReplayStructure
{
    std::vector<std::uint8_t> informative_rows;
    ReplayStructure(Eigen::Index row_count,const JointParameterLayout & layout)
        :informative_rows(static_cast<std::size_t>(row_count))
    {
        for(const auto row:layout.informative_rows) informative_rows.at(row)=1;
    }
};
void ReplayInto(const Domain & domain,const JointParameterLayout & layout,VectorRef observations,
    VectorRef eta,VectorRef beta,double scale,const ReplayStructure & structure,ReplayState & out)
{
    if(out.prediction.size()!=observations.size()) out.prediction.resize(observations.size());
    if(out.residual.size()!=observations.size()) out.residual.resize(observations.size());
    out.prediction.setZero(); out.residual.setZero();
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); const double width=std::exp(eta(a));
        for(const auto & support:domain.atoms[atom]) if(structure.informative_rows.at(static_cast<std::size_t>(support.row)))
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
    out.objective=.5*squared/(scale*scale);
}
struct AffectedReplayContributor
{
    Eigen::Index atom_position{};
    double square{};
};
static_assert(sizeof(AffectedReplayContributor)==sizeof(Eigen::Index)+sizeof(double));
struct AffectedReplayIncidence
{
    Indices atoms;
    std::vector<std::size_t> offsets;
    std::vector<AffectedReplayContributor> contributors;
    std::size_t build_scratch_bytes{};
};
struct AffectedReplayResult
{
    Vector prediction,residual;
    double objective{};
    std::size_t scratch_bytes{};
};
Eigen::Index AffectedReplayPosition(const Indices & row_positions,Eigen::Index row)
{
    if(row<0 || static_cast<std::size_t>(row)>=row_positions.size()) return -1;
    return row_positions[static_cast<std::size_t>(row)];
}
AffectedReplayIncidence PrepareAffectedReplayIncidence(const Domain & domain,const Indices & core_atoms,
    const Indices & neighbor_atoms,const Indices & affected_rows,const Indices & atom_positions,
    const Indices & row_positions)
{
    if(neighbor_atoms.size()>std::numeric_limits<std::size_t>::max()-core_atoms.size())
        throw std::overflow_error("Affected-row replay atom count overflow.");
    AffectedReplayIncidence out;
    out.atoms.reserve(core_atoms.size()+neighbor_atoms.size());
    out.atoms.insert(out.atoms.end(),core_atoms.begin(),core_atoms.end());
    out.atoms.insert(out.atoms.end(),neighbor_atoms.begin(),neighbor_atoms.end());
    std::sort(out.atoms.begin(),out.atoms.end(),[&](Eigen::Index left,Eigen::Index right) {
        return atom_positions.at(static_cast<std::size_t>(left))<atom_positions.at(static_cast<std::size_t>(right));
    });
    out.atoms.erase(std::unique(out.atoms.begin(),out.atoms.end()),out.atoms.end());
    if(affected_rows.size()==std::numeric_limits<std::size_t>::max())
        throw std::overflow_error("Affected-row replay row count overflow.");
    out.offsets.assign(affected_rows.size()+1,0);
    for(const auto atom_value:out.atoms)
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom_value)])
        {
            const auto local=AffectedReplayPosition(row_positions,support.row);
            if(local>=0 && !CheckedAddBytes(out.offsets[static_cast<std::size_t>(local)+1],1,
                out.offsets[static_cast<std::size_t>(local)+1]))
                throw std::overflow_error("Affected-row replay contributor count overflow.");
        }
    for(std::size_t row=1;row<out.offsets.size();++row)
        if(!CheckedAddBytes(out.offsets[row-1],out.offsets[row],out.offsets[row]))
            throw std::overflow_error("Affected-row replay contributor prefix overflow.");
    std::vector<std::size_t> cursor(out.offsets.begin(),out.offsets.end()-1);
    out.contributors.resize(out.offsets.back());
    for(std::size_t atom_position=0;atom_position<out.atoms.size();++atom_position)
    {
        const auto atom=static_cast<std::size_t>(out.atoms[atom_position]);
        for(const auto & support:domain.atoms[atom])
        {
            const auto local=AffectedReplayPosition(row_positions,support.row);
            if(local<0) continue;
            const auto row=static_cast<std::size_t>(local);
            out.contributors[cursor[row]++]={static_cast<Eigen::Index>(atom_position),support.square};
        }
    }
    std::size_t bytes{};
    if(!AddMemoryBytes(bytes,out.atoms.capacity(),sizeof(Eigen::Index)) ||
        !AddMemoryBytes(bytes,out.offsets.capacity(),sizeof(std::size_t)) ||
        !AddMemoryBytes(bytes,cursor.capacity(),sizeof(std::size_t)) ||
        !AddMemoryBytes(bytes,out.contributors.capacity(),sizeof(AffectedReplayContributor)))
        throw std::overflow_error("Affected-row replay memory estimate overflow.");
    out.build_scratch_bytes=bytes;
    return out;
}
AffectedReplayResult ReplayAffectedRows(const Domain & domain,VectorRef observations,
    const BlockCoordinateState & accepted,VectorRef candidate_eta,
    VectorRef candidate_beta,double scale,const Indices & core_atoms,const Indices & neighbor_atoms,
    const Indices & affected_rows,const Indices & atom_position,const Indices & row_positions)
{
    const auto incidence=PrepareAffectedReplayIncidence(domain,core_atoms,neighbor_atoms,affected_rows,
        atom_position,row_positions);
    Vector widths(static_cast<Eigen::Index>(incidence.atoms.size()));
    for(std::size_t k=0;k<incidence.atoms.size();++k)
        widths(static_cast<Eigen::Index>(k))=std::exp(candidate_eta(incidence.atoms[k]));
    AffectedReplayResult out;
    out.prediction.resize(static_cast<Eigen::Index>(affected_rows.size()));
    out.residual.resize(static_cast<Eigen::Index>(affected_rows.size()));
    for(std::size_t row=0;row<affected_rows.size();++row)
    {
        double prediction{};
        for(auto entry=incidence.offsets[row];entry<incidence.offsets[row+1];++entry)
        {
            const auto & contributor=incidence.contributors[entry];
            const auto atom=incidence.atoms[static_cast<std::size_t>(contributor.atom_position)];
            const auto basis=EvaluateKernel(contributor.square,widths(contributor.atom_position),2.5);
            prediction+=candidate_beta(2*atom)*basis.gaussian+candidate_beta(2*atom+1)*basis.charge;
        }
        out.prediction(static_cast<Eigen::Index>(row))=prediction;
        out.residual(static_cast<Eigen::Index>(row))=prediction-observations(affected_rows[row]);
    }
    {
        double delta_squared{};
        for(std::size_t k=0;k<affected_rows.size();++k)
        {
            const double old_residual=accepted.residual(affected_rows[k]);
            const double new_residual=out.residual(static_cast<Eigen::Index>(k));
            delta_squared+=new_residual*new_residual-old_residual*old_residual;
        }
        out.objective=accepted.objective+.5*delta_squared/(scale*scale);
    }
    std::size_t evaluation_bytes{};
    if(!AddMemoryBytes(evaluation_bytes,incidence.atoms.capacity(),sizeof(Eigen::Index)) ||
        !AddMemoryBytes(evaluation_bytes,incidence.offsets.capacity(),sizeof(std::size_t)) ||
        !AddMemoryBytes(evaluation_bytes,incidence.contributors.capacity(),sizeof(AffectedReplayContributor)) ||
        !AddMemoryBytes(evaluation_bytes,static_cast<std::size_t>(widths.size()),sizeof(double)) ||
        !AddMemoryBytes(evaluation_bytes,static_cast<std::size_t>(out.prediction.size()),sizeof(double)) ||
        !AddMemoryBytes(evaluation_bytes,static_cast<std::size_t>(out.residual.size()),sizeof(double)))
        throw std::overflow_error("Affected-row replay memory estimate overflow.");
    out.scratch_bytes=std::max(incidence.build_scratch_bytes,evaluation_bytes);
    return out;
}
bool BitwiseEqual(double left,double right)
{return std::bit_cast<std::uint64_t>(left)==std::bit_cast<std::uint64_t>(right);}
bool IsBitwiseIdenticalCoreCandidate(const BlockCoordinateState & state,
    const Indices & core_atoms,const Endpoint & endpoint)
{
    const auto core_count=static_cast<Eigen::Index>(core_atoms.size());
    if(state.beta.size()!=2*state.eta.size() || endpoint.eta.size()!=core_count || endpoint.beta.size()!=2*core_count)
        return false;
    for(std::size_t k=0;k<core_atoms.size();++k)
    {
        const auto atom=core_atoms[k];
        const auto local=static_cast<Eigen::Index>(k);
        const double current_eta=state.eta(atom),candidate_eta=endpoint.eta(local);
        const double current_a=state.beta(2*atom),candidate_a=endpoint.beta(2*local);
        const double current_c=state.beta(2*atom+1),candidate_c=endpoint.beta(2*local+1);
        if(!std::isfinite(current_eta) || !std::isfinite(candidate_eta) ||
            !std::isfinite(current_a) || !std::isfinite(candidate_a) ||
            !std::isfinite(current_c) || !std::isfinite(candidate_c) ||
            !BitwiseEqual(current_eta,candidate_eta) || !BitwiseEqual(current_a,candidate_a) ||
            !BitwiseEqual(current_c,candidate_c)) return false;
    }
    return true;
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
struct LocalEndpointWorkspace
{
    Vector eta,beta,effective_response;
    EvaluationContext context;
};
LocalEndpointWorkspace PrepareLocalEndpointWorkspace(const PreparedFixedNeighborBlock & prepared,
    const BlockCoordinateState & state,const std::optional<double> & eta_upper_bound)
{
    LocalEndpointWorkspace out;
    out.eta.resize(prepared.local_atoms); out.beta.resize(2*prepared.local_atoms);
    for(std::size_t k=0;k<prepared.core_atoms.size();++k)
    {
        const auto atom=static_cast<Eigen::Index>(prepared.core_atoms[k]);
        out.eta(static_cast<Eigen::Index>(k))=state.eta(atom);
        out.beta.segment<2>(2*static_cast<Eigen::Index>(k))=state.beta.segment<2>(2*atom);
    }
    const Sparse design=Design(prepared.domain,out.eta);
    const Vector local_core=(design*out.beta).eval();
    out.effective_response.resize(static_cast<Eigen::Index>(prepared.affected_rows.size()));
    for(std::size_t k=0;k<prepared.affected_rows.size();++k)
    {
        const auto row=static_cast<Eigen::Index>(prepared.affected_rows[k]);
        out.effective_response(static_cast<Eigen::Index>(k))=local_core(static_cast<Eigen::Index>(k))-
            state.residual(row);
    }
    out.context=prepared.context;
    if(eta_upper_bound) out.context.eta_upper_bound=eta_upper_bound;
    return out;
}
double PredictionDifference(VectorRef a,VectorRef b)
{return (a-b).lpNorm<Eigen::Infinity>();}
FixedNeighborBlockFailureDiagnostics FailureDiagnostics(const FixedNeighborPolicy & policy,
    const FixedNeighborBlockRecord & block,std::size_t sweep_index,const Indices & core_atoms,
    const EvaluationContext & context,const SearchResult & profile)
{
    FixedNeighborBlockFailureDiagnostics out;
    out.component_index=policy.component_index; out.component_id=policy.component_id;
    out.sweep_index=sweep_index+1; out.block_index=block.block;
    out.core_atom_count=core_atoms.size(); out.affected_row_count=block.affected_rows;
    out.local_search_stop_reason=profile.stop_reason;
    out.initial_profile_valid=profile.initial.valid;
    out.initial_profile_reason=profile.initial.reason;
    out.initial_profile_trusted=profile.initial_accepted;
    out.width_upper_bound=context.eta_upper_bound ? std::optional<double>(std::exp(*context.eta_upper_bound)) : std::nullopt;
    const auto & certificate=profile.initial.certificate;
    out.solver_valid=certificate.linear_solver_valid;
    out.solver_reason=certificate.linear_solver_reason;
    out.linear_rank=certificate.free_rank;
    out.linear_solve_count=certificate.linear_solves;
    out.active_set_releases=certificate.linear_releases;
    if(profile.initial.beta.size())
    {
        out.primal_feasibility=0.0;
        for(Eigen::Index column=0;column<profile.initial.beta.size();column+=2)
            *out.primal_feasibility=std::max(*out.primal_feasibility,
                std::max(0.0,-profile.initial.beta(column)));
    }
    if(certificate.available) out.scaled_projected_kkt=certificate.projected_kkt;
    for(Eigen::Index column=0;column<profile.initial.beta.size();++column)
        if(!std::isfinite(profile.initial.beta(column)) ||
            (column%2==0 && profile.initial.beta(column)<0))
        {
            out.first_invalid_coefficient=column;
            out.first_invalid_coefficient_value=profile.initial.beta(column);
            const auto atom=column/2;
            if(atom>=0 && static_cast<std::size_t>(atom)<context.atom_ids.size())
            {
                out.first_invalid_atom_id=context.atom_ids[static_cast<std::size_t>(atom)];
                out.first_invalid_parameter=column%2==0 ? "A" : "C";
            }
            break;
        }
    for(const auto atom:core_atoms)
        if(atom>=0 && static_cast<std::size_t>(atom)<context.atom_ids.size())
            out.core_atom_ids.push_back(context.atom_ids.at(static_cast<std::size_t>(atom)));
    for(auto trial=profile.trials.rbegin();trial!=profile.trials.rend();++trial)
        if(trial->lm)
        {
            out.lm_subproblem_status=trial->lm->subproblem_status;
            out.lm_actual_decrease=trial->lm->actual_decrease;
            out.lm_predicted_decrease=trial->lm->predicted_decrease;
            out.lm_ratio=trial->lm->ratio;
            out.width_step_inf=trial->lm->step.size() ?
                trial->lm->step.lpNorm<Eigen::Infinity>() : 0.0;
            out.weak_width_step_fraction=trial->weak_width_step_fraction;
            out.width_bound_active=context.eta_upper_bound && trial->endpoint.eta.size()>0 &&
                (trial->endpoint.eta.array()>=*context.eta_upper_bound-
                    64.0*std::numeric_limits<double>::epsilon()*
                        std::max(1.0,std::abs(*context.eta_upper_bound))).any();
            break;
        }
    for(auto trial=profile.trials.rbegin();trial!=profile.trials.rend();++trial)
        if(trial->trust && !trial->trust->passed)
        {
            out.replay_prediction_discrepancy=trial->trust->prediction_difference;
            out.replay_gradient_discrepancy=trial->trust->gradient_difference;
            out.replay_kkt_discrepancy=trial->trust->kkt_difference;
            break;
        }
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
       <<" initial_profile_valid="<<(d.initial_profile_valid ? "true" : "false")
       <<" initial_profile_reason="<<d.initial_profile_reason
       <<" initial_profile_trusted="<<(d.initial_profile_trusted ? "true" : "false")
       <<" solver_valid="<<(d.solver_valid ? "true" : "false")<<" solver_reason="<<d.solver_reason;
    out<<" linear_rank="; if(d.linear_rank) out<<*d.linear_rank; else out<<"unavailable";
    out<<" linear_solve_count="; if(d.linear_solve_count) out<<*d.linear_solve_count; else out<<"unavailable";
    out<<" active_set_releases="; if(d.active_set_releases) out<<*d.active_set_releases; else out<<"unavailable";
    out<<" primal_feasibility="; if(d.primal_feasibility) out<<*d.primal_feasibility; else out<<"unavailable";
    out<<" scaled_projected_kkt="; if(d.scaled_projected_kkt) out<<*d.scaled_projected_kkt; else out<<"unavailable";
    out<<" width_upper_bound="; if(d.width_upper_bound) out<<*d.width_upper_bound; else out<<"unavailable";
    out<<" width_bound_active="<<(d.width_bound_active ? "true" : "false")
       <<" width_step_inf="; if(d.width_step_inf) out<<*d.width_step_inf; else out<<"unavailable";
    out<<" lm_subproblem_status="<<(d.lm_subproblem_status.empty() ? "unavailable" : d.lm_subproblem_status)
       <<" lm_actual_decrease="; if(d.lm_actual_decrease) out<<*d.lm_actual_decrease; else out<<"unavailable";
    out<<" lm_predicted_decrease="; if(d.lm_predicted_decrease) out<<*d.lm_predicted_decrease; else out<<"unavailable";
    out<<" lm_ratio="; if(d.lm_ratio) out<<*d.lm_ratio; else out<<"unavailable";
    out<<" weak_width_step_fraction="; if(d.weak_width_step_fraction) out<<*d.weak_width_step_fraction; else out<<"unavailable";
    out<<" local_search_stop_reason="<<d.local_search_stop_reason;
    if(d.first_invalid_coefficient)
        out<<" first_invalid_coefficient="<<*d.first_invalid_coefficient
           <<" first_invalid_coefficient_value="<<*d.first_invalid_coefficient_value
           <<" first_invalid_atom_id="<<d.first_invalid_atom_id
           <<" first_invalid_parameter="<<d.first_invalid_parameter;
    if(d.replay_prediction_discrepancy)
        out<<" replay_prediction_discrepancy="<<*d.replay_prediction_discrepancy
           <<" replay_gradient_discrepancy="<<*d.replay_gradient_discrepancy
           <<" replay_kkt_discrepancy="<<*d.replay_kkt_discrepancy;
    return out.str();
}
int StabilityProbePriority(const JointBlockConditioning & evidence)
{
    if(evidence.active_width_bound_targets.value_or(0)>0) return 0;
    const bool free_ac_singular=evidence.free_ac_rank && *evidence.free_ac_rank<evidence.free_ac_columns;
    if(free_ac_singular) return -1;
    double practical_kkt_ratio{};
    if(evidence.local_ac_projected_kkt && std::isfinite(*evidence.local_ac_projected_kkt))
        practical_kkt_ratio=std::max(practical_kkt_ratio,
            *evidence.local_ac_projected_kkt/PracticalJointConvergenceCriteria.global_ac_kkt);
    if(evidence.local_projected_width_kkt && std::isfinite(*evidence.local_projected_width_kkt))
        practical_kkt_ratio=std::max(practical_kkt_ratio,
            *evidence.local_projected_width_kkt/PracticalJointConvergenceCriteria.global_width_gradient);
    if(practical_kkt_ratio>=.5) return 1;
    if(evidence.conditioning_classification=="Singular") return 2;
    if(evidence.conditioning_classification=="IllConditioned") return 3;
    if(evidence.conditioning_classification=="ActiveFace" || evidence.active_amplitude_count>0) return 4;
    if(evidence.conditioning_classification=="WellConditioned") return 5;
    return -1;
}
double StabilitySensitivityScore(const JointBlockConditioning & evidence)
{
    double score{};
    double conditioning_score{};
    if(evidence.local_ac_projected_kkt && std::isfinite(*evidence.local_ac_projected_kkt))
        score=std::max(score,*evidence.local_ac_projected_kkt/
            PracticalJointConvergenceCriteria.global_ac_kkt);
    if(evidence.local_projected_width_kkt && std::isfinite(*evidence.local_projected_width_kkt))
        score=std::max(score,*evidence.local_projected_width_kkt/
            PracticalJointConvergenceCriteria.global_width_gradient);
    for(const auto * condition:{&evidence.free_ac_condition,&evidence.projected_width_condition,
        &evidence.profile_condition})
        if(*condition && std::isfinite(**condition) && **condition>=1.0)
            conditioning_score=std::max(conditioning_score,std::log10(**condition));
    return score+conditioning_score;
}
double StabilityBlockCenter(const JointBlockConditioning & evidence)
{
    if(evidence.core_atoms.empty()) return std::numeric_limits<double>::infinity();
    const auto [minimum,maximum]=std::minmax_element(evidence.core_atoms.begin(),evidence.core_atoms.end());
    return .5*(static_cast<double>(*minimum)+static_cast<double>(*maximum));
}
std::vector<std::size_t> RepresentativeStabilityProbeSelection(
    const std::vector<JointBlockConditioning> & conditioning,std::vector<std::size_t> candidates,std::size_t limit)
{
    std::sort(candidates.begin(),candidates.end(),[&](std::size_t left,std::size_t right) {
        const auto & a=conditioning.at(left); const auto & b=conditioning.at(right);
        const double center_a=StabilityBlockCenter(a),center_b=StabilityBlockCenter(b);
        if(center_a!=center_b) return center_a<center_b;
        const double support_a=static_cast<double>(a.affected_row_count)/
            static_cast<double>(std::max<std::size_t>(a.core_atom_count,1));
        const double support_b=static_cast<double>(b.affected_row_count)/
            static_cast<double>(std::max<std::size_t>(b.core_atom_count,1));
        if(support_a!=support_b) return support_a>support_b;
        return a.block_id<b.block_id;
    });
    if(limit==0 || candidates.empty()) return {};
    if(candidates.size()<=limit) return candidates;
    if(limit==1) return {candidates[candidates.size()/2]};
    std::vector<std::size_t> selected; selected.reserve(limit);
    for(std::size_t slot=0;slot<limit;++slot)
        selected.push_back(candidates[slot*(candidates.size()-1)/(limit-1)]);
    return selected;
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
        if(priority<0) continue;
        candidates.push_back({index,priority,StabilitySensitivityScore(block),block.block_id});
    }
    std::sort(candidates.begin(),candidates.end(),[&conditioning](const auto & left,const auto & right) {
        if(left.priority!=right.priority) return left.priority<right.priority;
        if(left.score!=right.score) return left.score>right.score;
        const auto center_left=StabilityBlockCenter(conditioning[left.index]);
        const auto center_right=StabilityBlockCenter(conditioning[right.index]);
        if(center_left!=center_right) return center_left<center_right;
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
    const Domain & domain,const JointParameterLayout & layout,VectorRef y,VectorRef initial_eta,const EvaluationContext & parent_context,
    const FixedNeighborPolicy & policy)
{
    FixedNeighborSearchResult out;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool stop_after_stationarity=policy.stop_after_stationarity;
    const bool assess_local_conditioning=policy.assess_local_conditioning;
    const bool assess_parameter_stability=policy.assess_parameter_stability;
    const bool collect_telemetry=policy.collect_telemetry;
    const bool collect_diagnostics=policy.collect_diagnostics;
#else
    constexpr bool stop_after_stationarity=true;
    constexpr bool assess_local_conditioning=true;
    constexpr bool assess_parameter_stability=true;
    constexpr bool collect_telemetry=true;
    constexpr bool collect_diagnostics=false;
#endif
    std::array<std::pair<double,double>,3> recent_sweep_changes{};
    std::size_t next_recent_sweep_change{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool endpoint_certification_enabled=policy.certify_endpoint;
#else
    constexpr bool endpoint_certification_enabled=true;
#endif
    EvaluationContext context=parent_context;
    if(policy.eta_upper_bound)
    {
        if(!std::isfinite(*policy.eta_upper_bound) ||
            (context.eta_upper_bound &&
                *context.eta_upper_bound!=*policy.eta_upper_bound))
        {out.reason="invalid-global-width-bound"; return out;}
        context.eta_upper_bound=policy.eta_upper_bound;
    }
    if(initial_eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) || !initial_eta.allFinite() ||
        y.size()!=domain.rows || !y.allFinite() || !(context.scale>0) || !std::isfinite(context.scale) ||
        policy.core_atoms==0 || layout.full_atoms.size()!=domain.atoms.size() ||
        layout.informative_rows.size()!=static_cast<std::size_t>(domain.rows) ||
        context.atom_ids.size()!=domain.atoms.size() || context.row_ids.size()!=static_cast<std::size_t>(domain.rows))
    {out.reason="block-invalid-partition"; return out;}
    for(Eigen::Index a=0;a<initial_eta.size();++a)
        if(!(std::exp(initial_eta(a))>0) || !std::isfinite(std::exp(initial_eta(a))))
        {out.reason="block-invalid-partition"; return out;}
    if(!context.eta_upper_bound || !std::isfinite(*context.eta_upper_bound))
    {out.reason="invalid-width-upper-bound"; return out;}
    if((initial_eta.array()>*context.eta_upper_bound).any())
    {out.reason="bounded-initial-width-infeasible"; return out;}

    auto structural_timer=Measure(collect_telemetry,out.work.structural_preparation_seconds);
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
        !AddMemoryBytes(base_resident,memory_shape.component_observations,sizeof(std::uint8_t)) ||
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
        partition=topology.BuildCores(policy.core_atoms);
        if(partition.cores.empty()) {out.reason="block-invalid-partition"; return out;}
        census=topology.Census(partition);
        memory_shape.component_blocks=partition.cores.size();
        memory_shape.affected_row_entries=census.affected_row_entries;
        memory_shape.component_neighbor_entries=census.neighbor_atom_entries;
        memory_shape.maximum_block_rows=census.maximum_block_rows;
        memory_shape.maximum_block_support_entries=census.maximum_block_support_entries;
        memory_shape.total_block_support_entries=census.total_block_support_entries;
        memory_shape.maximum_neighbor_atoms=census.maximum_neighbor_atoms;
        memory_shape.maximum_block_replay_entries=census.maximum_block_replay_entries;
        memory_shape.maximum_block_atoms=census.maximum_block_atoms;
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
    global_shape.affected_replay_scratch=true;
    auto global_resident=resident;
    if(!AddMemoryBytes(global_resident,memory_shape.component_support_entries,sizeof(LocalSupportIndex)) ||
        !AddMemoryBytes(global_resident,memory_shape.affected_row_entries,5*sizeof(Eigen::Index)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_atoms,4*sizeof(Eigen::Index)+3*sizeof(double)+
            sizeof(std::vector<LocalSupportIndex>)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_blocks,
            sizeof(PreparedFixedNeighborBlock)+sizeof(CompactRowMapping)+sizeof(std::size_t)) ||
        // The reusable replay vectors remain allocated between block visits.
        !AddMemoryBytes(global_resident,memory_shape.component_observations,5*sizeof(double)+sizeof(Eigen::Index)))
    {
        out.global_state_memory_estimate.shape=memory_shape;
        out.global_state_memory_estimate.stage=MemoryStage::GlobalStateEvaluation;
        out.global_state_memory_estimate.status=MemoryStatus::Overflow;
        out.global_state_memory_estimate.reason="size-overflow";
        out.memory_estimate=out.global_state_memory_estimate;
        out.reason="resource-limit:global-state-evaluation:size-overflow"; return out;
    }
    global_shape.resident_bytes=global_resident;
    std::size_t global_memory_budget=JointMemoryPlanningBudgetBytes;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(policy.global_state_memory_budget_for_testing)
        global_memory_budget=*policy.global_state_memory_budget_for_testing;
#endif
    out.global_state_memory_estimate=PlanMemory(global_shape,MemoryStage::GlobalStateEvaluation,global_memory_budget);
    bool use_affected_row_replay=true;
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass)
    {
        auto full_replay_shape=global_shape;
        full_replay_shape.affected_replay_scratch=false;
        const auto full_replay_estimate=PlanMemory(full_replay_shape,MemoryStage::GlobalStateEvaluation,
            global_memory_budget);
        if(full_replay_estimate.status==MemoryStatus::Pass)
        {
            use_affected_row_replay=false;
            global_shape=full_replay_shape;
            out.global_state_memory_estimate=full_replay_estimate;
            out.work.full_replay_memory_fallback=true;
        }
        else
        {
            out.global_state_memory_estimate=full_replay_estimate;
        }
    }
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass ||
        out.global_state_memory_estimate.estimated_peak_bytes>out.global_state_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.global_state_memory_estimate);
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass)
    {
        out.memory_estimate=out.global_state_memory_estimate;
        out.reason="resource-limit:global-state-evaluation:"+out.memory_estimate.reason; return out;
    }
    auto * const profile_work=&out.work.local_profile_work;
    const auto replay_structure_started=Clock::now();
    const ReplayStructure replay_structure(domain.rows,layout);
    ReplayState replay_workspace;
    if(collect_telemetry)
        out.work.structural_preparation_seconds+=Seconds(replay_structure_started);
    const auto certify_endpoint=[&] {
        const auto assessment_timer=Measure(collect_telemetry,out.work.endpoint_verification_seconds);
        const auto record_assessment_work=[&](FixedNeighborEndpointCertification certification) {
            if(collect_telemetry)
            {
                out.work.primary_endpoint_seconds+=certification.assessment.primary_endpoint_seconds;
                out.work.independent_replay_seconds+=certification.assessment.independent_replay_seconds;
            }
            return certification;
        };
        auto certification_resident=out.global_state_memory_estimate.shape.resident_bytes;
        if(!AddMemoryBytes(certification_resident,memory_shape.component_observations,3*sizeof(double)) ||
            !AddMemoryBytes(certification_resident,memory_shape.component_atoms,8*sizeof(double)))
            certification_resident=std::numeric_limits<std::size_t>::max();
        return record_assessment_work(CertifyFixedNeighborEndpoint(domain,y,out.state,context,profile_work,
            certification_resident,policy.convergence_criteria));
    };

    out.state.eta=initial_eta;
    out.state.beta=Vector::Zero(2*static_cast<Eigen::Index>(domain.atoms.size()));
    ReplayState initial;
    ReplayInto(domain,layout,y,out.state.eta,out.state.beta,context.scale,replay_structure,initial);
    out.state.prediction=std::move(initial.prediction); out.state.residual=std::move(initial.residual);
    out.state.objective=initial.objective;
    bool state_cache_valid=out.state.prediction.size()==y.size() && out.state.residual.size()==y.size() &&
        out.state.eta.allFinite() && out.state.beta.allFinite() && out.state.prediction.allFinite() &&
        out.state.residual.allFinite() && std::isfinite(out.state.objective);
    auto mapping_timer=Measure(collect_telemetry,out.work.structural_preparation_seconds);
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
        prepared.local_atoms=static_cast<Eigen::Index>(prepared.profile_atoms.size());
        prepared.context.rank={prepared.domain.rows,2*prepared.local_atoms,prepared.local_atoms};
        prepared.context.linear.rank_relative=prepared.context.rank.Relative(2*prepared.local_atoms);
        // FixedNeighbor performs at most one trusted accepted local update per block visit.
        prepared.context.update_budget=1;
        prepared_blocks.push_back(std::move(prepared));
    }
    out.block_preparations=prepared_blocks.size();
    out.domain_preparations=prepared_blocks.size();
    out.mapping_preparations=prepared_blocks.size();
    if(mapping_timer.target)
    {
        out.work.structural_preparation_seconds+=Seconds(mapping_timer.started);
        mapping_timer.target=nullptr;
    }
    for(const auto row:layout.informative_rows)
        row_position.at(static_cast<std::size_t>(row))=-1;
    std::size_t mapped_replay_block=std::numeric_limits<std::size_t>::max();
    const bool retain_diagnostics=collect_diagnostics;
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
            if(collect_telemetry) out.work.maximum_block_seconds=std::max(
                out.work.maximum_block_seconds,Seconds(block_started));
            if(record.accepted)
            {
                ++sweep.accepted_blocks;
                sweep.accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
                out.total_accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
            }
            else if(record.status=="unchanged") ++sweep.unchanged_blocks;
            if(retain_diagnostics) out.blocks.push_back(std::move(record));
        };
        for(std::size_t block_index=0;block_index<partition.cores.size();++block_index)
        {
            const auto & prepared=prepared_blocks.at(block_index);
            const auto & core_atoms=prepared.core_atoms; const auto & affected_rows=prepared.affected_rows;
            const auto & local_domain=prepared.domain; const auto & local_context=prepared.context;
            const auto local_atoms=prepared.local_atoms;
            const auto block_started=collect_telemetry ? Clock::now() : Clock::time_point{};
            FixedNeighborBlockRecord record;
            record.sweep=sweep_index+1; record.block=block_index+1;
            record.affected_rows=affected_rows.size(); record.objective_before=out.state.objective;
            record.recovery_retry=block_failure_visits[block_index]>0;
            if(record.recovery_retry)
            {
                record.recovery_attempt=block_failure_visits[block_index];
                ++out.recovery_attempts;
                ++sweep.recovery_attempts;
            }
            Vector old_eta(local_atoms),old_beta(2*local_atoms),old_core;
            {
                const auto old_core_timer=Measure(collect_telemetry,out.work.old_core_seconds);
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
                const auto response_timer=Measure(collect_telemetry,out.work.effective_response_seconds);
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

            const auto search_started=collect_telemetry ? Clock::now() :
                Clock::time_point{};
            if(collect_telemetry)
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
            const auto profile_search=SearchProfile(local_domain,local_y,old_widths,local_context,profile_work);
            out.work.local_profile_work.Merge(profile_search.profile_work);
            if(collect_telemetry) out.work.local_search_seconds+=Seconds(search_started);
            record.accepted_updates=profile_search.accepted;
            ++sweep.block_solves; ++out.total_block_solves;
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
                {
                    if(!it->accepted || !it->endpoint.valid || !it->trust || !it->trust->passed) continue;
                    accepted=&*it; break;
                }
            if(!accepted)
            {
                record.failure_diagnostics=FailureDiagnostics(policy,record,sweep_index,core_atoms,
                    context,profile_search);
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
            const Endpoint & local_state=accepted->endpoint;
            if(!local_state.valid || !local_state.certificate.available || !local_state.certificate.feasible ||
                !local_state.certificate.kkt_passed)
            {
                record.status="failed"; record.reason="block-inner-invalid"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            record.local_objective_after=local_state.certificate.objective/(context.scale*context.scale);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto candidate_fault=policy.candidate_fault_for_testing ?
                policy.candidate_fault_for_testing(sweep_index+1,block_index+1,out.state) :
                    FixedNeighborCandidateFaultForTesting::None;
            if(candidate_fault==FixedNeighborCandidateFaultForTesting::LocalObjectiveIncrease)
                record.local_objective_after+=std::max(1.0,std::abs(record.local_objective_before));
#endif
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
            // The candidate path only overwrites core parameters in copies of the current state.
            // Once the accepted endpoint is qualified, bitwise-equal core values leave the full
            // parameter state unchanged and the current full-replay cache remains valid.
            bool candidate_result_affected_rows_only{};
            Vector candidate_affected_prediction,candidate_affected_residual;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            double incremental_candidate_objective=std::numeric_limits<double>::quiet_NaN();
#endif
            const bool candidate_fast_unchanged=state_cache_valid && std::isfinite(out.state.objective) &&
                IsBitwiseIdenticalCoreCandidate(out.state,core_atoms,accepted->endpoint)
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                && candidate_fault==FixedNeighborCandidateFaultForTesting::None
#endif
                ;
            Vector candidate_eta,candidate_beta;
            if(candidate_fast_unchanged)
            {
                if(collect_telemetry) ++out.work.candidate_replay_skipped_unchanged;
            }
            else
            {
                {
                    const auto candidate_copy_timer=Measure(collect_telemetry,out.work.candidate_copy_seconds);
                    candidate_eta=out.state.eta; candidate_beta=out.state.beta;
                    for(std::size_t k=0;k<core_atoms.size();++k)
                    {
                        const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                        candidate_eta(atom)=accepted->endpoint.eta(static_cast<Eigen::Index>(k));
                        candidate_beta.segment<2>(2*atom)=accepted->endpoint.beta.segment<2>(2*static_cast<Eigen::Index>(k));
                    }
                    (void)candidate_copy_timer;
                }
                const auto candidate_replay_timer=Measure(collect_telemetry,out.work.candidate_replay_seconds);
                if(use_affected_row_replay)
                {
                    if(mapped_replay_block!=std::numeric_limits<std::size_t>::max())
                        for(const auto row:prepared_blocks[mapped_replay_block].affected_rows)
                            row_position.at(static_cast<std::size_t>(row))=-1;
                    for(std::size_t k=0;k<affected_rows.size();++k)
                        row_position.at(static_cast<std::size_t>(affected_rows[k]))=static_cast<Eigen::Index>(k);
                    mapped_replay_block=block_index;
                    auto candidate=ReplayAffectedRows(domain
                        ,y,out.state,candidate_eta,candidate_beta,
                        context.scale,core_atoms,partition.cores[block_index].neighbor_atoms,affected_rows,
                        atom_position,row_position);
                    candidate_result_affected_rows_only=true;
                    candidate_affected_prediction=std::move(candidate.prediction);
                    candidate_affected_residual=std::move(candidate.residual);
                    replay_workspace.objective=candidate.objective;
                    if(collect_telemetry)
                    {
                        ++out.work.affected_candidate_replays;
                        ++out.work.incremental_objective_evaluations;
                        out.work.affected_replay_peak_scratch_bytes=std::max(
                            out.work.affected_replay_peak_scratch_bytes,candidate.scratch_bytes);
                    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    incremental_candidate_objective=candidate.objective;
#endif
                }
                else
                {
                    if(collect_telemetry) ++out.work.full_candidate_replays;
                    ReplayInto(domain,layout,y,candidate_eta,candidate_beta,context.scale,replay_structure,
                        replay_workspace);
                }
                (void)candidate_replay_timer;
            }
            double candidate_objective=candidate_fast_unchanged ? out.state.objective : replay_workspace.objective;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            switch(candidate_fault)
            {
            case FixedNeighborCandidateFaultForTesting::ObjectiveConsistencyFailure:
                candidate_objective+=std::max(1.0,std::abs(candidate_objective));
                replay_workspace.objective=candidate_objective;
                break;
            case FixedNeighborCandidateFaultForTesting::NaNObjective:
                candidate_objective=std::numeric_limits<double>::quiet_NaN();
                replay_workspace.objective=candidate_objective;
                break;
            case FixedNeighborCandidateFaultForTesting::InfiniteObjective:
                candidate_objective=std::numeric_limits<double>::infinity();
                replay_workspace.objective=candidate_objective;
                break;
            case FixedNeighborCandidateFaultForTesting::NaNParameter:
                candidate_eta(core_atoms.front())=std::numeric_limits<double>::quiet_NaN();
                break;
            case FixedNeighborCandidateFaultForTesting::InfiniteParameter:
                candidate_eta(core_atoms.front())=std::numeric_limits<double>::infinity();
                break;
            case FixedNeighborCandidateFaultForTesting::InfinitePrediction:
                if(candidate_result_affected_rows_only)
                    candidate_affected_prediction(0)=std::numeric_limits<double>::infinity();
                else replay_workspace.prediction(affected_rows.front())=std::numeric_limits<double>::infinity();
                break;
            case FixedNeighborCandidateFaultForTesting::InfiniteResidual:
                if(candidate_result_affected_rows_only)
                    candidate_affected_residual(0)=std::numeric_limits<double>::infinity();
                else replay_workspace.residual(affected_rows.front())=std::numeric_limits<double>::infinity();
                break;
            default: break;
            }
#endif
            record.objective_after=candidate_objective;
            record.global_replay_delta=candidate_fast_unchanged ? 0.0 : candidate_objective-out.state.objective;
            const auto replay_reference=std::max(std::abs(candidate_objective),std::abs(out.state.objective));
            record.objective_replay_enclosure=BlockObjectiveReplayEnclosure(replay_reference);
            record.local_global_delta_error=std::abs(record.global_replay_delta+
                record.local_objective_before-record.local_objective_after);
            const bool delta_enclosed=WithinBlockObjectiveReplay(record.local_global_delta_error,replay_reference);
            const bool candidate_unchanged=candidate_fast_unchanged ||
                ((out.state.eta.array()==candidate_eta.array()).all() &&
                    (out.state.beta.array()==candidate_beta.array()).all());
            const bool candidate_values_finite=candidate_fast_unchanged ?
                out.state.eta.allFinite() && out.state.beta.allFinite() && out.state.prediction.allFinite() &&
                    out.state.residual.allFinite() && std::isfinite(candidate_objective) :
                candidate_eta.allFinite() && candidate_beta.allFinite() && std::isfinite(candidate_objective) &&
                    (candidate_result_affected_rows_only ?
                        candidate_affected_prediction.allFinite() && candidate_affected_residual.allFinite() :
                        replay_workspace.prediction.allFinite() && replay_workspace.residual.allFinite());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto observe_candidate=[&](bool accepted_candidate,bool unchanged_candidate) {
                if(!policy.candidate_observer) return;
                const auto & observed_eta=candidate_fast_unchanged ? out.state.eta : candidate_eta;
                const auto & observed_beta=candidate_fast_unchanged ? out.state.beta : candidate_beta;
                static const Vector empty_vector;
                const auto & observed_affected_prediction=candidate_result_affected_rows_only ?
                    candidate_affected_prediction : empty_vector;
                const auto & observed_affected_residual=candidate_result_affected_rows_only ?
                    candidate_affected_residual : empty_vector;
                policy.candidate_observer({sweep_index+1,block_index+1,out.state,core_atoms,affected_rows,
                    observed_eta,observed_beta,observed_affected_prediction,observed_affected_residual,
                    candidate_objective,
                    record.local_objective_before,record.local_objective_after,record.global_replay_delta,
                    record.local_global_delta_error,record.objective_replay_enclosure,incremental_candidate_objective,
                    accepted_candidate,
                    unchanged_candidate,!candidate_fast_unchanged});
            };
#endif
            if(!candidate_values_finite)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                observe_candidate(false,candidate_unchanged);
#endif
                record.status="failed"; record.reason="block-objective-replay-failed";
                out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            if(!delta_enclosed)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                observe_candidate(false,candidate_unchanged);
#endif
                record.status="failed"; record.reason="block-objective-replay-failed"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            if(record.global_replay_delta>0 && !(record.local_objective_after<=record.local_objective_before &&
                record.global_replay_delta<=record.objective_replay_enclosure))
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                observe_candidate(false,candidate_unchanged);
#endif
                record.status=record.recovery_retry ? "failed" : "unchanged";
                record.reason=record.recovery_retry ? "block-recovery-exhausted" : "block-objective-increase";
                if(record.recovery_retry) {out.reason=record.reason; stop=true;}
                finish_record(std::move(record),block_started);
                if(stop) break;
                continue;
            }
            const bool unchanged=candidate_unchanged;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            observe_candidate(!unchanged,unchanged);
#endif
            if(!unchanged)
                for(const auto atom:core_atoms)
                    sweep.maximum_accepted_log_width_step=std::max(sweep.maximum_accepted_log_width_step,
                        std::abs(candidate_eta(static_cast<Eigen::Index>(atom))-
                            out.state.eta(static_cast<Eigen::Index>(atom))));
            if(!candidate_fast_unchanged)
            {
                const auto cache_started=Clock::now();
                if(candidate_result_affected_rows_only)
                {
                    for(std::size_t k=0;k<affected_rows.size();++k)
                    {
                        const auto row=affected_rows[k];
                        out.state.prediction(row)=candidate_affected_prediction(static_cast<Eigen::Index>(k));
                        out.state.residual(row)=candidate_affected_residual(static_cast<Eigen::Index>(k));
                    }
                    out.state.objective=replay_workspace.objective;
                }
                else
                {
                    out.state.prediction.swap(replay_workspace.prediction);
                    out.state.residual.swap(replay_workspace.residual);
                    out.state.objective=replay_workspace.objective;
                }
                if(collect_telemetry) out.work.affected_row_updates+=affected_rows.size();
                out.state.eta.swap(candidate_eta); out.state.beta.swap(candidate_beta);
                state_cache_valid=out.state.prediction.size()==y.size() && out.state.residual.size()==y.size() &&
                    std::isfinite(out.state.objective);
                if(collect_telemetry) out.work.cache_update_seconds+=Seconds(cache_started);
            }
            record.accepted=!unchanged; record.status=unchanged ? "unchanged" : "accepted";
            record.reason=unchanged ? "block-unchanged" : "";
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
        {
            const auto timer=Measure(collect_telemetry,out.work.sweep_replay_seconds);
            if(collect_telemetry) ++out.work.sweep_replays;
            ReplayInto(domain,layout,y,out.state.eta,out.state.beta,context.scale,replay_structure,replay_workspace);
            (void)timer;
        }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.perturb_sweep_replay_objective_for_testing)
            replay_workspace.objective+=1.0;
#endif
        const auto & replay=replay_workspace;
        const auto replay_audit_started=collect_telemetry ? Clock::now() : Clock::time_point{};
        sweep.cache_replay_error=PredictionDifference(out.state.prediction,replay.prediction);
        sweep.objective_replay_error=std::abs(out.state.objective-replay.objective);
        const bool sweep_replay_valid=
            sweep.cache_replay_error<=2e-12+2e-13*std::max(1.0,replay.prediction.cwiseAbs().maxCoeff()) &&
            WithinBlockObjectiveReplay(sweep.objective_replay_error,replay.objective);
        if(collect_telemetry)
        {
            out.work.maximum_cache_replay_error=std::max(out.work.maximum_cache_replay_error,sweep.cache_replay_error);
            out.work.sweep_replay_audit_seconds+=Seconds(replay_audit_started);
        }
        const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
        const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);
        const auto global_started=Clock::now();
        const auto global=EvaluateStateStreaming(domain,y,eta,beta,context,ProfileEvaluationRole::AcceptedEndpoint,profile_work);
        if(collect_telemetry) out.work.sweep_global_state_seconds+=Seconds(global_started);
        if(!global.valid || !global.certificate.available)
        {out.reason="block-inner-invalid"; stop=true; break;}
        sweep.global_a_feasibility=global.certificate.feasible ? 0.0 : 1.0;
        sweep.global_ac_kkt=global.certificate.projected_kkt;
        sweep.global_width_gradient_inf_norm=global.gradient.lpNorm<Eigen::Infinity>();
        sweep.width_bound_applied=global.certificate.width_bound_applied;
        sweep.projected_width_gradient_inf_norm=global.certificate.width_bound_applied ?
            global.certificate.width_bound_projected_kkt : sweep.global_width_gradient_inf_norm;
        sweep.maximum_width_bound_violation=global.certificate.maximum_width_bound_violation;
        sweep.boundary_gradient_max=global.certificate.boundary_gradient_max;
        sweep.boundary_multiplier_max=global.certificate.boundary_multiplier_max;
        sweep.width_bound_complementarity=global.certificate.width_bound_complementarity;
        sweep.active_width_bound_count=global.certificate.active_width_bound_count;
        sweep.eta_change_inf=(out.state.eta-eta_before).lpNorm<Eigen::Infinity>();
        sweep.beta_scaled_change=ScaledCoefficientDifference(out.state.beta,beta_before);
        sweep.maximum_target_width=out.state.eta.size() ? std::exp(out.state.eta.maxCoeff()) : 0.0;
        sweep.coordinate_confirmation_available=sweep_index>0;
        double stationarity_width_gradient=sweep.projected_width_gradient_inf_norm;
        double stationarity_eta_change=sweep.eta_change_inf;
        if(sweep.width_bound_applied && !global.certificate.width_bound_feasible)
            stationarity_width_gradient=std::numeric_limits<double>::infinity();
        sweep.stationarity_width_gradient_inf_norm=stationarity_width_gradient;
        sweep.wall_seconds=Seconds(sweep_started); out.search_seconds+=sweep.wall_seconds;
        ++out.sweep_count; out.final_sweep=sweep;
        if(retain_diagnostics) out.sweeps.push_back(sweep);
        if(!sweep_replay_valid)
        {out.reason="block-cache-replay-failed"; stop=true; break;}
        // Reconciliation follows the successful audit and global-state evaluation.
        const auto objective_reconciliation_started=collect_telemetry ? Clock::now() : Clock::time_point{};
        ++out.work.objective_reconciliations;
        out.work.maximum_objective_reconciliation_error=std::max(
            out.work.maximum_objective_reconciliation_error,sweep.objective_replay_error);
        out.state.objective=replay.objective;
        sweep.objective_after=replay.objective;
        sweep.objective_cache_refreshed=true;
        if(collect_telemetry)
            out.work.objective_reconciliation_seconds+=Seconds(objective_reconciliation_started);
        out.final_sweep=sweep;
        if(retain_diagnostics) out.sweeps.back()=sweep;
        if(policy.production_progress) policy.production_progress(sweep);
        if(assess_parameter_stability)
        {
            recent_sweep_changes[next_recent_sweep_change]={sweep.eta_change_inf,sweep.beta_scaled_change};
            next_recent_sweep_change=(next_recent_sweep_change+1)%recent_sweep_changes.size();
            out.recent_sweep_count=std::min(out.recent_sweep_count+1,recent_sweep_changes.size());
            out.maximum_eta_change=0; out.maximum_scaled_ac_change=0;
            for(std::size_t k=0;k<out.recent_sweep_count;++k)
            {
                out.maximum_eta_change=std::max(out.maximum_eta_change,recent_sweep_changes[k].first);
                out.maximum_scaled_ac_change=std::max(out.maximum_scaled_ac_change,recent_sweep_changes[k].second);
            }
        }
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
            if(!stop_after_stationarity || !endpoint_certification_enabled)
            {
                out.search_converged=true; out.reason="block-stationary";
                if(stop_after_stationarity) break;
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
        if(!recovery_pending && stop_after_stationarity && endpoint_certification_enabled && out.sweep_count &&
            out.certification_sweep!=out.sweep_count)
        {
            out.certification=certify_endpoint();
            out.certification_sweep=out.sweep_count;
            if(out.certification->qualified() && out.final_sweep &&
                ClassifyFixedNeighborStationarity(out.final_sweep->global_ac_kkt,
                    out.final_sweep->stationarity_width_gradient_inf_norm,out.final_sweep->eta_change_inf,
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
    const bool assess_failed_diagnostic=collect_diagnostics && out.certification.has_value();
    if(assess_local_conditioning && out.certification &&
        ((out.search_converged && out.certification->qualified()) || assess_failed_diagnostic))
    {
        out.local_conditioning.reserve(prepared_blocks.size());
        for(const auto & prepared:prepared_blocks)
        {
            auto workspace=PrepareLocalEndpointWorkspace(prepared,out.state,policy.eta_upper_bound);
            const auto block_id=(policy.component_id.empty() ? std::string("component") : policy.component_id)+
                "/block-"+std::to_string(prepared.block+1);
            auto conditioning=AssessFixedNeighborLocalConditioning(prepared.domain,workspace.effective_response,
                workspace.eta,workspace.beta,workspace.context,block_id,
                out.global_state_memory_estimate.shape.resident_bytes);
            conditioning.evidence.core_atoms.assign(prepared.core_atoms.begin(),prepared.core_atoms.end());
            out.work.local_conditioning_seconds+=conditioning.seconds;
            if(conditioning.evidence.workspace_peak_memory_estimate_bytes)
                out.work.local_conditioning_peak_memory_estimate_bytes=std::max(
                    out.work.local_conditioning_peak_memory_estimate_bytes,
                    *conditioning.evidence.workspace_peak_memory_estimate_bytes);
            out.local_conditioning.push_back(std::move(conditioning.evidence));
        }
    }
    if(assess_parameter_stability && out.search_converged && out.certification &&
        out.certification->qualified())
    {
        const auto started=Clock::now(); JointParameterStability stability;
        stability.recent_sweep_count=out.recent_sweep_count;
        stability.maximum_eta_change=out.maximum_eta_change;
        stability.maximum_scaled_ac_change=out.maximum_scaled_ac_change;
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
        const auto probe_order=StabilityProbeOrder(out.local_conditioning);
        std::vector<std::size_t> representative_candidates;
        for(const auto index:probe_order)
        {
            const auto & condition=out.local_conditioning[index];
            const auto priority=StabilityProbePriority(condition);
            if(priority==5)
            {
                representative_candidates.push_back(index);
                continue;
            }
            if(selected.size()>=MaximumPerturbedBlocksPerFit)
            {
                block_evidence[index].reason="probe-block-budget-exhausted";
                continue;
            }
            selected.push_back(index);
        }
        const auto representatives=RepresentativeStabilityProbeSelection(out.local_conditioning,
            representative_candidates,MaximumRepresentativeLocalStabilityProbesPerComponent);
        for(const auto index:representative_candidates)
        {
            if(std::find(representatives.begin(),representatives.end(),index)==representatives.end())
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
        }
        bool budget_exhausted=false;
        for(const auto index:selected)
        {
            auto & evidence=block_evidence[index];
            const int priority=StabilityProbePriority(out.local_conditioning[index]);
            if(budget_exhausted || !policy.claim_parameter_stability_block ||
                !policy.claim_parameter_stability_block(priority))
            {
                evidence.reason="perturbation-budget-exhausted";
                budget_exhausted=true;
                continue;
            }
            const auto & prepared=prepared_blocks[index];
            const auto block_id=evidence.block_id;
            auto workspace=PrepareLocalEndpointWorkspace(prepared,out.state,policy.eta_upper_bound);
            auto assessment=AssessLocalParameterStability(prepared.domain,workspace.effective_response,
                workspace.eta,workspace.beta,workspace.context,block_id,
                out.global_state_memory_estimate.shape.resident_bytes);
            evidence=std::move(assessment.evidence);
            evidence.probe_selected=true;
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
    const auto primary_started=Clock::now();
    const auto endpoint=EvaluateStateStreaming(domain,y,state.eta,state.beta,context,
        ProfileEvaluationRole::AcceptedEndpoint,profile_work);
    const double primary_seconds=Seconds(primary_started);
    const auto replay_started=Clock::now();
    const auto replay=CheckIndependentReplay(domain,y,endpoint,context,nullptr,nullptr,convergence_criteria);
    const double replay_seconds=Seconds(replay_started);
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
    context.eta_upper_bound=production_policy.eta_upper_bound;
    const auto atom_count=static_cast<Eigen::Index>(local_layout.full_atoms.size());
    context.rank={prepared.parent_row_count,2*atom_count,atom_count};
    context.linear.rank_relative=context.rank.Relative(2*atom_count);

    Vector initial_eta(atom_count);
    for(std::size_t k=0;k<prepared.local_to_parent_atoms.size();++k)
        initial_eta(static_cast<Eigen::Index>(k))=std::log(initial_b(prepared.local_to_parent_atoms[k]));

    FixedNeighborPolicy policy;
    policy.core_atoms=ProductionFixedNeighborCoreAtoms;
    policy.maximum_sweeps=ProductionFixedNeighborMaximumSweeps;
    policy.convergence_criteria=PracticalJointConvergenceCriteria;
    policy.component_index=production_policy.component_index;
    policy.component_id=production_policy.component_id;
    policy.eta_upper_bound=production_policy.eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    policy.assess_local_conditioning=true;
    policy.assess_parameter_stability=true;
#endif
    policy.claim_parameter_stability_block=production_policy.claim_parameter_stability_block;
    std::optional<JointFixedNeighborProgress> last_fixed_neighbor_progress;
    if(progress_component && observer)
    {
        const auto maximum_sweeps=policy.maximum_sweeps;
        policy.production_progress=[&observer,&last_fixed_neighbor_progress,progress_component,maximum_sweeps,context](
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
            progress.maximum_sweeps=maximum_sweeps;
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
            last_fixed_neighbor_progress=progress;
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
            LastFixedNeighborObjective(out.search),LastFixedNeighborGradient(out.search),last_fixed_neighbor_progress);
    return out;
}

}
