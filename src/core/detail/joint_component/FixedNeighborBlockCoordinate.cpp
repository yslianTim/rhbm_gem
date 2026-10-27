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
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace rhbm_gem::core::joint_component {
double BlockObjectiveReplayEnclosure(double reference)
{return 1e-12+2e-12*std::abs(reference);}
bool WithinBlockObjectiveReplay(double error,double reference)
{return error<=BlockObjectiveReplayEnclosure(reference);}
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
Eigen::Index AffectedReplayPosition(const Indices & affected_rows,Eigen::Index row)
{
    const auto found=std::lower_bound(affected_rows.begin(),affected_rows.end(),row);
    return found==affected_rows.end() || *found!=row ? -1 : static_cast<Eigen::Index>(found-affected_rows.begin());
}
AffectedReplayIncidence PrepareAffectedReplayIncidence(const Domain & domain,const Indices & core_atoms,
    const Indices & neighbor_atoms,const Indices & affected_rows,const Indices & atom_positions)
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
            const auto local=AffectedReplayPosition(affected_rows,support.row);
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
            const auto local=AffectedReplayPosition(affected_rows,support.row);
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
AffectedReplayResult ReplayAffectedRows(const Domain & domain,const JointParameterLayout & layout,
    VectorRef observations,const BlockCoordinateState & accepted,VectorRef candidate_eta,
    VectorRef candidate_beta,double scale,const Indices & core_atoms,const Indices & neighbor_atoms,
    const Indices & affected_rows,const Indices & atom_position
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ,bool incremental_objective=false
#endif
    )
{
    const auto incidence=PrepareAffectedReplayIncidence(domain,core_atoms,neighbor_atoms,affected_rows,atom_position);
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(incremental_objective)
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
    else
#endif
    {
        double squared{};
        for(const auto row_value:layout.informative_rows)
        {
            const auto row=static_cast<Eigen::Index>(row_value);
            const auto local=AffectedReplayPosition(affected_rows,row);
            const double residual=local>=0 ? out.residual(local) : accepted.residual(row);
            squared+=residual*residual;
        }
        out.objective=.5*squared/(scale*scale);
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
double FullObjectiveFromAffectedRows(const JointParameterLayout & layout,VectorRef accepted_residual,
    const Indices & affected_rows,VectorRef candidate_residual,double scale)
{
    double squared{};
    for(const auto row_value:layout.informative_rows)
    {
        const auto row=static_cast<Eigen::Index>(row_value);
        const auto local=AffectedReplayPosition(affected_rows,row);
        const double residual=local>=0 ? candidate_residual(local) : accepted_residual(row);
        squared+=residual*residual;
    }
    return .5*squared/(scale*scale);
}
#endif
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
    if(spectrum.available)
    {
        out.local_design_rank=spectrum.rank;
        out.local_design_rank_threshold=spectrum.threshold;
        out.local_design_singular_values.assign(spectrum.singular_values.data(),
            spectrum.singular_values.data()+spectrum.singular_values.size());
    }
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
    out.initial_widths.reserve(static_cast<std::size_t>(eta.size()));
    for(Eigen::Index atom=0;atom<eta.size();++atom) out.initial_widths.push_back(std::exp(eta(atom)));
    if(context.eta_upper_bound)
        out.width_upper_bound=std::exp(*context.eta_upper_bound);
    if(profile.initial.beta.size()==design.cols())
        for(Eigen::Index column=0;column<profile.initial.beta.size();++column)
            if(!std::isfinite(profile.initial.beta(column)) ||
                (column%2==0 && profile.initial.beta(column)<0))
            {
                out.first_invalid_coefficient=column;
                out.first_invalid_coefficient_value=profile.initial.beta(column);
                const auto atom=column/2;
                if(atom>=0 && static_cast<std::size_t>(atom)<local_context.atom_ids.size())
                {
                    out.first_invalid_atom_id=local_context.atom_ids[static_cast<std::size_t>(atom)];
                    out.first_invalid_parameter=column%2==0 ? "A" : "C";
                }
                break;
            }
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(policy.collect_diagnostics && profile.initial.beta.size()==design.cols())
    {
        const auto initial_state=EvaluateState(domain,effective_response,eta,profile.initial.beta,context);
        if(initial_state.valid && initial_state.derivative.cols()==eta.size())
        {
            out.width_derivative_column_norms.reserve(static_cast<std::size_t>(eta.size()));
            for(Eigen::Index column=0;column<initial_state.derivative.cols();++column)
                out.width_derivative_column_norms.push_back(initial_state.derivative.col(column).norm());
            if(initial_state.gradient.size())
                out.width_gradient_inf=initial_state.gradient.lpNorm<Eigen::Infinity>();
            const auto projected_prepared=PrepareProjectedWidthDerivative(
                initial_state,context.scale,&context);
            if(projected_prepared.valid)
            {
                const auto projected=ReduceProjectedWidthDerivative(projected_prepared);
                if(projected.valid)
                {
                    out.projected_width_column_norms.assign(projected.projected_norms.data(),
                        projected.projected_norms.data()+projected.projected_norms.size());
                    const Sparse projected_sparse=projected.projected.sparseView();
                    const auto width_spectrum=ComputeSpectrum(projected_sparse,context.rank,
                        projected.projected.cols(),true);
                    if(width_spectrum.available)
                    {
                        out.projected_width_rank_threshold=width_spectrum.threshold;
                        out.projected_width_singular_values.assign(width_spectrum.singular_values.data(),
                            width_spectrum.singular_values.data()+width_spectrum.singular_values.size());
                    }
                }
            }
        }
    }
#endif
    for(auto trial=profile.trials.rbegin();trial!=profile.trials.rend();++trial)
        if(trial->lm)
        {
            out.lm_subproblem_status=trial->lm->subproblem_status;
            out.lm_actual_decrease=trial->lm->actual_decrease;
            out.lm_predicted_decrease=trial->lm->predicted_decrease;
            out.lm_ratio=trial->lm->ratio;
            out.width_step_inf=trial->lm->step.size() ?
                trial->lm->step.lpNorm<Eigen::Infinity>() : 0.0;
            out.proposed_eta_step.assign(trial->lm->step.data(),
                trial->lm->step.data()+trial->lm->step.size());
            out.weak_width_step_fraction=trial->weak_width_step_fraction;
            out.width_bound_active=context.eta_upper_bound && trial->endpoint.eta.size()==eta.size() &&
                (trial->endpoint.eta.array()>=*context.eta_upper_bound-
                    64.0*std::numeric_limits<double>::epsilon()*
                        std::max(1.0,std::abs(*context.eta_upper_bound))).any();
            break;
        }
    for(const auto atom:core_atoms)
        if(atom>=0 && static_cast<std::size_t>(atom)<context.atom_ids.size())
            out.core_atom_ids.push_back(context.atom_ids.at(static_cast<std::size_t>(atom)));
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
FixedNeighborReplayResult FullReplayForTesting(const Domain & domain,const JointParameterLayout & layout,
    VectorRef observations,VectorRef eta,VectorRef beta,double scale)
{
    const ReplayStructure structure(domain.rows,layout);
    ReplayState replay;
    ReplayInto(domain,layout,observations,eta,beta,scale,structure,replay);
    return {std::move(replay.prediction),std::move(replay.residual),replay.objective};
}
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
    constexpr bool collect_telemetry=false;
    constexpr bool collect_diagnostics=false;
#endif
    std::array<std::pair<double,double>,3> recent_sweep_changes{};
    std::size_t next_recent_sweep_change{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const bool endpoint_certification_enabled=policy.certify_endpoint;
    const bool capture_local_trajectory=policy.capture_local_trajectory;
#else
    constexpr bool endpoint_certification_enabled=true;
    constexpr bool capture_local_trajectory=false;
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(!std::isfinite(policy.objective_guard_band_fraction_for_testing) ||
        policy.objective_guard_band_fraction_for_testing<0.0)
    {out.reason="invalid-objective-guard-policy"; return out;}
    if(policy.objective_strategy_for_testing!=FixedNeighborObjectiveStrategyForTesting::FullScan &&
        policy.candidate_replay_for_testing)
    {out.reason="invalid-objective-test-strategy"; return out;}
#endif
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
        std::isfinite(out.state.objective);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    bool objective_cache_trusted=true;
    double last_objective_cache_error{};
    std::size_t accepted_updates_since_objective_refresh{};
#endif
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
    const bool retain_diagnostics=collect_diagnostics || capture_local_trajectory
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        || static_cast<bool>(policy.state_observer)
#endif
        ;
    const bool capture_block_records=retain_diagnostics
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        || (policy.diagnostic_sink && static_cast<bool>(policy.diagnostic_sink->block))
#endif
        ;
    bool stop=false;
    std::size_t next_certification_retry_sweep{};
    std::vector<std::size_t> block_failure_visits(partition.cores.size());
    for(std::size_t sweep_index=0;sweep_index<policy.maximum_sweeps && !stop;++sweep_index)
    {
        const auto sweep_started=Clock::now(); FixedNeighborBlockSweep sweep; sweep.sweep=sweep_index+1;
        const Vector eta_before=out.state.eta,beta_before=out.state.beta;
        sweep.objective_before=out.state.objective;
        bool sweep_had_block_failure=false;
        auto finish_record=[&](FixedNeighborBlockRecord record,const Indices & core_atoms,
            Clock::time_point block_started) {
            double block_seconds{};
            if(capture_block_records || collect_telemetry) block_seconds=Seconds(block_started);
            if(capture_block_records)
            {
                record.core_atoms=core_atoms;
                record.block_seconds=block_seconds;
            }
            if(collect_telemetry)
            {
                out.work.maximum_block_seconds=std::max(out.work.maximum_block_seconds,block_seconds);
            }
            if(record.accepted)
            {
                ++sweep.accepted_blocks;
                sweep.accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
                out.total_accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
            }
            else if(record.status=="unchanged") ++sweep.unchanged_blocks;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(policy.diagnostic_sink && policy.diagnostic_sink->block) policy.diagnostic_sink->block(record);
#endif
            if(retain_diagnostics) out.blocks.push_back(std::move(record));
        };
        for(std::size_t block_index=0;block_index<partition.cores.size();++block_index)
        {
            const auto & prepared=prepared_blocks.at(block_index);
            const auto & core_atoms=prepared.core_atoms; const auto & affected_rows=prepared.affected_rows;
            const auto & local_domain=prepared.domain; const auto & local_context=prepared.context;
            const auto local_atoms=prepared.local_atoms;
            const auto block_started=capture_block_records || collect_telemetry ? Clock::now() : Clock::time_point{};
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(capture_block_records)
            {record.component_index=policy.component_index; record.component_id=policy.component_id;}
#endif

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

            const auto search_started=capture_block_records || collect_telemetry ? Clock::now() :
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
            const double search_seconds=capture_block_records || collect_telemetry ?
                Seconds(search_started) : 0.0;
            if(capture_block_records)
            {
                record.search_seconds=search_seconds; record.profile_evaluations=profile_search.evaluations;
                record.local_search_stop_reason=profile_search.stop_reason;
            }
            if(collect_telemetry) out.work.local_search_seconds+=search_seconds;
            record.accepted_updates=profile_search.accepted;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(capture_local_trajectory)
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
                        (trial.endpoint.eta-accepted_eta).cwiseAbs().maxCoeff() : unavailable;
                    telemetry.gradient_inf_norm=trial.endpoint.gradient.size() ?
                        trial.endpoint.gradient.cwiseAbs().maxCoeff() : unavailable;
                    telemetry.profile_seconds=trial.seconds; telemetry.factor_seconds=trial.factor_seconds;
                    cumulative_factor_seconds+=trial.factor_seconds;
                    telemetry.cumulative_factor_seconds=cumulative_factor_seconds;
                    if(trial.lm)
                    {
                        telemetry.lm_actual_decrease=trial.lm->actual_decrease;
                        telemetry.lm_predicted_decrease=trial.lm->predicted_decrease;
                        telemetry.lm_ratio=trial.lm->ratio;
                        telemetry.lm_proposed_acceptance=trial.lm->proposed_acceptance;
                        telemetry.bounded_subproblem=trial.lm->bounded_subproblem;
                        telemetry.bounded_subproblem_status=trial.lm->subproblem_status;
                        telemetry.lm_subproblem_kkt_inf_norm=trial.lm->subproblem_kkt_inf_norm;
                        telemetry.active_set_iterations=trial.lm->active_set_iterations;
                        telemetry.active_set_releases=trial.lm->active_set_releases;
                        telemetry.active_bound_indices=trial.lm->active_bound_indices;
                        telemetry.active_bound_multipliers=trial.lm->active_bound_multipliers;
                    }
                    record.profile_trials.push_back(std::move(telemetry));
                    if(trial.accepted && std::isfinite(objective))
                    {
                        accepted_objective=objective;
                        if(trial.endpoint.eta.size()==accepted_eta.size()) accepted_eta=trial.endpoint.eta;
                    }
                }
                record.profile_factor_seconds=cumulative_factor_seconds;
            }
#endif
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
                finish_record(std::move(record),core_atoms,block_started);
                if(stop) break;
                continue;
            }
            const Endpoint & local_state=accepted->endpoint;
            if(!local_state.valid || !local_state.certificate.available || !local_state.certificate.feasible ||
                !local_state.certificate.kkt_passed)
            {
                record.status="failed"; record.reason="block-inner-invalid"; out.reason=record.reason;
                finish_record(std::move(record),core_atoms,block_started); stop=true; break;
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
                finish_record(std::move(record),core_atoms,block_started);
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
            bool candidate_objective_full_scan_fallback{};
#endif
            const bool candidate_fast_unchanged=state_cache_valid && std::isfinite(out.state.objective) &&
                IsBitwiseIdenticalCoreCandidate(out.state,core_atoms,accepted->endpoint);
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if(policy.candidate_replay_for_testing)
                {
                    auto candidate=policy.candidate_replay_for_testing(out.state,core_atoms,affected_rows,
                        candidate_eta,candidate_beta);
                    replay_workspace.prediction=std::move(candidate.prediction);
                    replay_workspace.residual=std::move(candidate.residual);
                    replay_workspace.objective=candidate.objective;
                    if(collect_telemetry) ++out.work.full_candidate_replays;
                }
                else
#endif
                if(use_affected_row_replay)
                {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    const bool incremental_objective=policy.objective_strategy_for_testing!=
                        FixedNeighborObjectiveStrategyForTesting::FullScan;
#endif
                    auto candidate=ReplayAffectedRows(domain,layout,y,out.state,candidate_eta,candidate_beta,
                        context.scale,core_atoms,partition.cores[block_index].neighbor_atoms,affected_rows,atom_position
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                        ,incremental_objective
#endif
                        );
                    candidate_result_affected_rows_only=true;
                    candidate_affected_prediction=std::move(candidate.prediction);
                    candidate_affected_residual=std::move(candidate.residual);
                    replay_workspace.objective=candidate.objective;
                    if(collect_telemetry)
                    {
                        ++out.work.affected_candidate_replays;
                        out.work.affected_replay_peak_scratch_bytes=std::max(
                            out.work.affected_replay_peak_scratch_bytes,candidate.scratch_bytes);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                        if(incremental_objective) ++out.work.incremental_objective_evaluations;
#endif
                    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                    if(incremental_objective) incremental_candidate_objective=candidate.objective;
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
            if(!candidate_fast_unchanged && candidate_result_affected_rows_only &&
                policy.objective_strategy_for_testing==FixedNeighborObjectiveStrategyForTesting::IncrementalGuarded)
            {
                const double estimated_delta=candidate_objective-out.state.objective;
                const double reference=std::max(std::abs(candidate_objective),std::abs(out.state.objective));
                const double enclosure=BlockObjectiveReplayEnclosure(reference);
                const double consistency=std::abs(estimated_delta+record.local_objective_before-
                    record.local_objective_after);
                const double boundary_distance=std::min({std::abs(estimated_delta),
                    std::abs(estimated_delta-enclosure),
                    std::abs(record.local_objective_after-record.local_objective_before),
                    std::abs(consistency-enclosure)});
                const double guard=policy.objective_guard_band_fraction_for_testing*enclosure;
                const bool cache_drifted=last_objective_cache_error>
                    policy.objective_guard_band_fraction_for_testing*BlockObjectiveReplayEnclosure(out.state.objective);
                candidate_objective_full_scan_fallback=!objective_cache_trusted || cache_drifted ||
                    boundary_distance<=guard;
                if(candidate_objective_full_scan_fallback)
                {
                    const auto fallback_started=collect_telemetry ? Clock::now() : Clock::time_point{};
                    candidate_objective=FullObjectiveFromAffectedRows(layout,out.state.residual,affected_rows,
                        candidate_affected_residual,context.scale);
                    replay_workspace.objective=candidate_objective;
                    if(collect_telemetry)
                    {
                        ++out.work.objective_full_scan_fallbacks;
                        out.work.objective_full_scan_fallback_seconds+=Seconds(fallback_started);
                    }
                }
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto observe_candidate=[&](bool accepted_candidate,bool unchanged_candidate) {
                if(!policy.candidate_observer) return;
                const auto & observed_eta=candidate_fast_unchanged ? out.state.eta : candidate_eta;
                const auto & observed_beta=candidate_fast_unchanged ? out.state.beta : candidate_beta;
                Vector observed_candidate_prediction,observed_candidate_residual;
                const Vector * observed_prediction=&replay_workspace.prediction;
                const Vector * observed_residual=&replay_workspace.residual;
                if(candidate_fast_unchanged)
                {
                    observed_prediction=&out.state.prediction; observed_residual=&out.state.residual;
                }
                else if(candidate_result_affected_rows_only)
                {
                    observed_candidate_prediction=out.state.prediction;
                    observed_candidate_residual=out.state.residual;
                    for(std::size_t k=0;k<affected_rows.size();++k)
                    {
                        const auto row=affected_rows[k];
                        observed_candidate_prediction(row)=candidate_affected_prediction(static_cast<Eigen::Index>(k));
                        observed_candidate_residual(row)=candidate_affected_residual(static_cast<Eigen::Index>(k));
                    }
                    observed_prediction=&observed_candidate_prediction;
                    observed_residual=&observed_candidate_residual;
                }
                policy.candidate_observer({sweep_index+1,block_index+1,out.state,core_atoms,affected_rows,
                    observed_eta,observed_beta,*observed_prediction,*observed_residual,candidate_objective,
                    record.local_objective_before,record.local_objective_after,record.global_replay_delta,
                    record.local_global_delta_error,record.objective_replay_enclosure,incremental_candidate_objective,
                    candidate_objective_full_scan_fallback,accepted_candidate,
                    unchanged_candidate,!candidate_fast_unchanged});
            };
#endif
            if(!delta_enclosed)
            {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                observe_candidate(false,candidate_unchanged);
#endif
                record.status="failed"; record.reason="block-objective-replay-failed"; out.reason=record.reason;
                finish_record(std::move(record),core_atoms,block_started); stop=true; break;
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
                finish_record(std::move(record),core_atoms,block_started);
                if(stop) break;
                continue;
            }
            const bool unchanged=candidate_unchanged;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const bool periodic_objective_refresh_due=!candidate_fast_unchanged && !unchanged &&
                policy.objective_refresh_interval_for_testing>0 &&
                accepted_updates_since_objective_refresh+1>=policy.objective_refresh_interval_for_testing;
            const auto objective_refresh_started=periodic_objective_refresh_due && collect_telemetry ?
                Clock::now() : Clock::time_point{};
            double old_objective_reference{};
            if(periodic_objective_refresh_due)
            {
                const Indices no_affected_rows;
                const Vector no_candidate_residual;
                old_objective_reference=FullObjectiveFromAffectedRows(layout,out.state.residual,
                    no_affected_rows,no_candidate_residual,context.scale);
            }
#endif
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
                out.state.eta=std::move(candidate_eta); out.state.beta=std::move(candidate_beta);
                state_cache_valid=out.state.prediction.size()==y.size() && out.state.residual.size()==y.size() &&
                    std::isfinite(out.state.objective);
                if(collect_telemetry) out.work.cache_update_seconds+=Seconds(cache_started);
            }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(!candidate_fast_unchanged && !unchanged)
            {
                ++accepted_updates_since_objective_refresh;
                if(periodic_objective_refresh_due)
                {
                    ReplayInto(domain,layout,y,out.state.eta,out.state.beta,context.scale,replay_structure,replay_workspace);
                    const double prediction_error=PredictionDifference(out.state.prediction,replay_workspace.prediction);
                    const double residual_error=PredictionDifference(out.state.residual,replay_workspace.residual);
                    const double objective_error=std::abs(out.state.objective-replay_workspace.objective);
                    const double reference=std::max(std::abs(old_objective_reference),std::abs(replay_workspace.objective));
                    const double true_delta=replay_workspace.objective-old_objective_reference;
                    const double prediction_limit=2e-12+2e-13*replay_workspace.prediction.cwiseAbs().maxCoeff();
                    const double residual_limit=2e-12+2e-13*replay_workspace.residual.cwiseAbs().maxCoeff();
                    const bool state_valid=prediction_error<=prediction_limit && residual_error<=residual_limit &&
                        WithinBlockObjectiveReplay(objective_error,replay_workspace.objective);
                    const bool true_increase_allowed=true_delta<=0 ||
                        (record.local_objective_after<=record.local_objective_before &&
                            true_delta<=BlockObjectiveReplayEnclosure(reference));
                    if(collect_telemetry)
                    {
                        ++out.work.objective_cache_refreshes;
                        out.work.maximum_accepted_updates_between_objective_refreshes=std::max(
                            out.work.maximum_accepted_updates_between_objective_refreshes,
                            accepted_updates_since_objective_refresh);
                        out.work.objective_cache_refresh_seconds+=Seconds(objective_refresh_started);
                        out.work.maximum_objective_refresh_error=std::max(
                            out.work.maximum_objective_refresh_error,objective_error);
                    }
                    ++sweep.objective_cache_refreshes;
                    sweep.maximum_objective_refresh_error=std::max(sweep.maximum_objective_refresh_error,objective_error);
                    if(!state_valid || !true_increase_allowed)
                    {
                        if(collect_telemetry) ++out.work.objective_refresh_failures;
                        out.reason=!state_valid ? "block-cache-replay-failed" : "block-objective-increase";
                        stop=true;
                    }
                    else
                    {
                        if(true_delta>0 && collect_telemetry) ++out.work.objective_refresh_positive_deltas;
                        out.state.objective=replay_workspace.objective;
                        sweep.objective_cache_refreshed=true;
                        objective_cache_trusted=true;
                        last_objective_cache_error=0;
                        accepted_updates_since_objective_refresh=0;
                    }
                }
            }
#endif
            record.accepted=!unchanged; record.status=unchanged ? "unchanged" : "accepted";
            record.reason=unchanged ? "block-unchanged" : "";
            if(record.recovery_retry)
            {
                ++out.recovery_successes;
                ++sweep.recovery_successes;
                block_failure_visits[block_index]=0;
                record.recovery_success=true;
            }
            finish_record(std::move(record),core_atoms,block_started);
        }
        if(stop) break;
        sweep.objective_after=out.state.objective;
        {
            const auto timer=Measure(collect_telemetry,out.work.sweep_replay_seconds);
            if(collect_telemetry) ++out.work.sweep_replays;
            ReplayInto(domain,layout,y,out.state.eta,out.state.beta,context.scale,replay_structure,replay_workspace);
            (void)timer;
        }
        const auto & replay=replay_workspace;
        sweep.cache_replay_error=PredictionDifference(out.state.prediction,replay.prediction);
        sweep.objective_replay_error=std::abs(out.state.objective-replay.objective);
        if(collect_telemetry)
            out.work.maximum_cache_replay_error=std::max(out.work.maximum_cache_replay_error,sweep.cache_replay_error);
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.diagnostic_sink && policy.diagnostic_sink->sweep) policy.diagnostic_sink->sweep(sweep);
#endif
        if(policy.production_progress) policy.production_progress(sweep);
        if(retain_diagnostics) out.sweeps.push_back(sweep);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.sweep_observer)
        {
            const auto & observed=retain_diagnostics ? out.sweeps.back() : sweep;
            policy.sweep_observer(observed);
        }
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.state_observer) policy.state_observer(sweep_index+1,out.state,out.sweeps.back(),out.blocks);
#endif
        const bool sweep_replay_valid=
            sweep.cache_replay_error<=2e-12+2e-13*std::max(1.0,replay.prediction.cwiseAbs().maxCoeff()) &&
            WithinBlockObjectiveReplay(sweep.objective_replay_error,replay.objective);
        if(!sweep_replay_valid)
        {out.reason="block-cache-replay-failed"; stop=true; break;}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        last_objective_cache_error=sweep.objective_replay_error;
        objective_cache_trusted=sweep.objective_replay_error<=
            policy.objective_guard_band_fraction_for_testing*BlockObjectiveReplayEnclosure(replay.objective);
        if(policy.refresh_objective_after_sweep_for_testing)
        {
            // The replay audit above must pass before objective-only refresh is allowed.
            ++sweep.objective_cache_refreshes;
            sweep.maximum_objective_refresh_error=sweep.objective_replay_error;
            if(collect_telemetry)
            {
                ++out.work.objective_cache_refreshes;
                out.work.maximum_objective_refresh_error=std::max(
                    out.work.maximum_objective_refresh_error,sweep.objective_replay_error);
                out.work.maximum_accepted_updates_between_objective_refreshes=std::max(
                    out.work.maximum_accepted_updates_between_objective_refreshes,
                    accepted_updates_since_objective_refresh);
            }
            accepted_updates_since_objective_refresh=0;
            out.state.objective=replay.objective;
            sweep.objective_after=replay.objective;
            sweep.objective_cache_refreshed=true;
            objective_cache_trusted=true;
        }
        out.final_sweep=sweep;
        if(retain_diagnostics) out.sweeps.back()=sweep;
#endif
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(production_policy.configure_search_for_testing) production_policy.configure_search_for_testing(policy);
#endif
    if(progress_component && observer)
    {
        const auto maximum_sweeps=policy.maximum_sweeps;
        policy.production_progress=[&observer,progress_component,maximum_sweeps,context](
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
