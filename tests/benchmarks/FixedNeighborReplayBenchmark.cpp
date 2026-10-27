#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/JointReplayExperiment.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;
using Vector=Eigen::VectorXd;

struct ReplayCandidate
{
    joint::Indices core_atoms,neighbor_atoms,affected_rows;
    Vector eta,beta;
};
struct GlobalReplayCsr
{
    std::vector<std::size_t> offsets;
    std::vector<Eigen::Index> atoms;
    std::vector<double> squares;
    std::size_t StorageBytes() const
    {
        return offsets.capacity()*sizeof(std::size_t)+atoms.capacity()*sizeof(Eigen::Index)+
            squares.capacity()*sizeof(double);
    }
};
struct GlobalReplayCsrEntry
{
    std::uint32_t atom_position{},support_position{};
};
static_assert(sizeof(GlobalReplayCsrEntry)==2*sizeof(std::uint32_t));
struct GlobalReplayCsrReference
{
    std::vector<std::size_t> offsets;
    std::vector<GlobalReplayCsrEntry> entries;
    std::size_t build_temporary_bytes{};
    std::size_t StorageBytes() const
    {
        std::size_t offsets_bytes{},entries_bytes{},total{};
        if(!joint::CheckedMultiplyBytes(offsets.capacity(),sizeof(std::size_t),offsets_bytes) ||
            !joint::CheckedMultiplyBytes(entries.capacity(),sizeof(GlobalReplayCsrEntry),entries_bytes) ||
            !joint::CheckedAddBytes(offsets_bytes,entries_bytes,total))
            throw std::overflow_error("E2 CSR storage estimate overflow.");
        return total;
    }
};
struct ErrorSummary
{
    std::size_t bitwise_mismatches{},normal_pairs{},nonfinite_pairs{};
    double maximum_absolute{},maximum_relative{};
    std::uint64_t maximum_ulp{};
};
struct MethodAccuracy
{
    ErrorSummary prediction,residual,objective;
    std::size_t objective_bitwise_mismatches{};
    Eigen::Index first_prediction_mismatch{-1};
    Eigen::Index first_residual_mismatch{-1};
};
struct MethodTiming
{
    std::vector<double> preparation,candidate,total;
    std::size_t maximum_temporary_bytes{};
    double checksum{};
};
struct Method
{
    const char * name;
    second_stage_test::ReplayPredictionStrategy prediction;
    second_stage_test::ReplayObjectiveStrategy objective;
    bool exact,global_csr,global_csr_reference{};
};
constexpr std::array<Method,7> Methods{{
    {"B1",second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false,false},
    {"B2",second_stage_test::ReplayPredictionStrategy::DeltaSubtractThenAdd,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false,false},
    {"C",second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayObjectiveStrategy::Incremental,false,false},
    {"D2",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false,true},
    {"D1",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::FullScan,true,false},
    {"E",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::Incremental,true,false},
    {"E2",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::Incremental,false,false,true}
}};
double Seconds(Clock::time_point started)
{return std::chrono::duration<double>(Clock::now()-started).count();}
bool BitwiseEqual(double left,double right)
{return std::bit_cast<std::uint64_t>(left)==std::bit_cast<std::uint64_t>(right);}
std::uint64_t OrderedBits(double value)
{
    const auto bits=std::bit_cast<std::uint64_t>(value);
    return bits&(std::uint64_t{1}<<63) ? ~bits : bits|(std::uint64_t{1}<<63);
}
std::uint64_t UlpDistance(double left,double right)
{
    const auto a=OrderedBits(left),b=OrderedBits(right);
    return a>b ? a-b : b-a;
}
void AddError(ErrorSummary & summary,double expected,double actual)
{
    if(!std::isfinite(expected) || !std::isfinite(actual))
    {++summary.nonfinite_pairs; return;}
    if(!BitwiseEqual(expected,actual)) ++summary.bitwise_mismatches;
    const double absolute=std::abs(actual-expected);
    const double relative=absolute/std::max(std::abs(expected),std::numeric_limits<double>::min());
    summary.maximum_absolute=std::max(summary.maximum_absolute,absolute);
    summary.maximum_relative=std::max(summary.maximum_relative,relative);
    if(std::isnormal(expected) && std::isnormal(actual))
    {
        ++summary.normal_pairs;
        summary.maximum_ulp=std::max(summary.maximum_ulp,UlpDistance(expected,actual));
    }
}
Eigen::Index AffectedPosition(const joint::Indices & rows,Eigen::Index row)
{
    const auto found=std::lower_bound(rows.begin(),rows.end(),row);
    return found==rows.end() || *found!=row ? -1 : static_cast<Eigen::Index>(found-rows.begin());
}
GlobalReplayCsr BuildGlobalReplayCsr(const joint::Domain & domain,
    const rhbm_gem::JointParameterLayout & layout)
{
    const auto row_count=static_cast<std::size_t>(domain.rows);
    std::vector<std::uint8_t> informative(row_count);
    std::vector<std::size_t> counts(row_count);
    for(const auto row:layout.informative_rows) informative.at(row)=1;
    for(const auto atom_value:layout.full_atoms)
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom_value)])
            if(informative.at(static_cast<std::size_t>(support.row)))
                ++counts[static_cast<std::size_t>(support.row)];
    GlobalReplayCsr out; out.offsets.resize(row_count+1);
    for(std::size_t row=0;row<row_count;++row) out.offsets[row+1]=out.offsets[row]+counts[row];
    out.atoms.resize(out.offsets.back()); out.squares.resize(out.offsets.back());
    std::vector<std::size_t> cursor(out.offsets.begin(),out.offsets.end()-1);
    for(const auto atom_value:layout.full_atoms)
    {
        const auto atom=static_cast<Eigen::Index>(atom_value);
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto row=static_cast<std::size_t>(support.row);
            if(!informative.at(row)) continue;
            const auto position=cursor[row]++;
            out.atoms[position]=atom; out.squares[position]=support.square;
        }
    }
    return out;
}
GlobalReplayCsrReference BuildGlobalReplayCsrReference(const joint::Domain & domain,
    const rhbm_gem::JointParameterLayout & layout)
{
    const auto row_count=static_cast<std::size_t>(domain.rows);
    if(row_count==std::numeric_limits<std::size_t>::max())
        throw std::overflow_error("Global CSR row count overflow.");
    if(layout.full_atoms.size()>std::numeric_limits<std::uint32_t>::max())
        throw std::overflow_error("Global CSR atom index does not fit its entry.");
    std::vector<std::uint8_t> informative(row_count);
    for(const auto row:layout.informative_rows) informative.at(row)=1;
    std::vector<std::size_t> counts(row_count);
    for(const auto atom_value:layout.full_atoms)
    {
        const auto & supports=domain.atoms[atom_value];
        if(supports.size()>std::numeric_limits<std::uint32_t>::max())
            throw std::overflow_error("Global CSR support index does not fit its entry.");
        for(const auto & support:supports)
        {
            const auto row=static_cast<std::size_t>(support.row);
            if(informative.at(row))
            {
                if(counts[row]==std::numeric_limits<std::size_t>::max())
                    throw std::overflow_error("Global CSR contributor count overflow.");
                ++counts[row];
            }
        }
    }
    GlobalReplayCsrReference out;
    out.offsets.resize(row_count+1);
    for(std::size_t row=0;row<row_count;++row)
    {
        if(counts[row]>std::numeric_limits<std::size_t>::max()-out.offsets[row])
            throw std::overflow_error("Global CSR row offsets overflow.");
        out.offsets[row+1]=out.offsets[row]+counts[row];
    }
    out.entries.resize(out.offsets.back());
    std::vector<std::size_t> cursor(out.offsets.begin(),out.offsets.end()-1);
    for(std::size_t atom_position=0;atom_position<layout.full_atoms.size();++atom_position)
    {
        const auto atom=layout.full_atoms[atom_position];
        const auto & supports=domain.atoms[atom];
        for(std::size_t support_position=0;support_position<supports.size();++support_position)
        {
            const auto row=static_cast<std::size_t>(supports[support_position].row);
            if(!informative.at(row)) continue;
            const auto entry_position=cursor[row]++;
            out.entries[entry_position]={static_cast<std::uint32_t>(atom_position),
                static_cast<std::uint32_t>(support_position)};
        }
    }
    std::size_t counts_bytes{},cursor_bytes{},informative_bytes{};
    if(!joint::CheckedMultiplyBytes(counts.capacity(),sizeof(std::size_t),counts_bytes) ||
        !joint::CheckedMultiplyBytes(cursor.capacity(),sizeof(std::size_t),cursor_bytes) ||
        !joint::CheckedMultiplyBytes(informative.capacity(),sizeof(std::uint8_t),informative_bytes) ||
        !joint::CheckedAddBytes(counts_bytes,cursor_bytes,counts_bytes) ||
        !joint::CheckedAddBytes(counts_bytes,informative_bytes,out.build_temporary_bytes))
        throw std::overflow_error("Global CSR build memory estimate overflow.");
    return out;
}
second_stage_test::ReplayCandidateResult EvaluateGlobalCsrCandidate(
    const rhbm_gem::JointParameterLayout & layout,
    joint::VectorRef observations,const joint::BlockCoordinateState & accepted,
    joint::VectorRef candidate_eta,joint::VectorRef candidate_beta,double scale,
    const joint::Indices & affected_rows,const GlobalReplayCsr & incidence)
{
    second_stage_test::ReplayCandidateResult out;
    out.prediction.resize(static_cast<Eigen::Index>(affected_rows.size()));
    out.residual.resize(static_cast<Eigen::Index>(affected_rows.size()));
    Vector widths(candidate_eta.size());
    for(Eigen::Index atom=0;atom<candidate_eta.size();++atom) widths(atom)=std::exp(candidate_eta(atom));
    for(std::size_t local=0;local<affected_rows.size();++local)
    {
        const auto row=static_cast<std::size_t>(affected_rows[local]);
        double prediction{};
        for(std::size_t entry=incidence.offsets[row];entry<incidence.offsets[row+1];++entry)
        {
            const auto atom=incidence.atoms[entry];
            const auto basis=joint::EvaluateKernel(incidence.squares[entry],widths(atom),2.5);
            prediction+=candidate_beta(2*atom)*basis.gaussian+candidate_beta(2*atom+1)*basis.charge;
        }
        out.prediction(static_cast<Eigen::Index>(local))=prediction;
        out.residual(static_cast<Eigen::Index>(local))=prediction-observations(affected_rows[local]);
    }
    double squared{};
    for(const auto row_value:layout.informative_rows)
    {
        const auto row=static_cast<Eigen::Index>(row_value);
        const auto local=AffectedPosition(affected_rows,row);
        const double residual=local>=0 ? out.residual(local) : accepted.residual(row);
        squared+=residual*residual;
    }
    out.objective=.5*squared/(scale*scale);
    out.scratch_bytes=static_cast<std::size_t>(2*out.prediction.size()+widths.size())*sizeof(double);
    return out;
}
second_stage_test::ReplayCandidateResult EvaluateGlobalCsrReferenceCandidate(
    const joint::Domain & domain,const rhbm_gem::JointParameterLayout & layout,
    joint::VectorRef observations,const joint::BlockCoordinateState & accepted,
    joint::VectorRef candidate_eta,joint::VectorRef candidate_beta,double scale,
    const joint::Indices & core_atoms,const joint::Indices & neighbor_atoms,
    const joint::Indices & affected_rows,const GlobalReplayCsrReference & incidence)
{
    std::vector<std::uint8_t> active(layout.full_atoms.size());
    for(const auto atom:core_atoms)
    {
        const auto found=std::find(layout.full_atoms.begin(),layout.full_atoms.end(),static_cast<std::size_t>(atom));
        if(found==layout.full_atoms.end()) throw std::invalid_argument("E2 core atom is absent from component layout.");
        active[static_cast<std::size_t>(found-layout.full_atoms.begin())]=1;
    }
    for(const auto atom:neighbor_atoms)
    {
        const auto found=std::find(layout.full_atoms.begin(),layout.full_atoms.end(),static_cast<std::size_t>(atom));
        if(found==layout.full_atoms.end()) throw std::invalid_argument("E2 neighbor atom is absent from component layout.");
        active[static_cast<std::size_t>(found-layout.full_atoms.begin())]=1;
    }
    Vector widths(static_cast<Eigen::Index>(layout.full_atoms.size()));
    for(std::size_t position=0;position<layout.full_atoms.size();++position)
        if(active[position]) widths(static_cast<Eigen::Index>(position))=
            std::exp(candidate_eta(static_cast<Eigen::Index>(layout.full_atoms[position])));
    second_stage_test::ReplayCandidateResult out;
    out.prediction.resize(static_cast<Eigen::Index>(affected_rows.size()));
    out.residual.resize(static_cast<Eigen::Index>(affected_rows.size()));
    for(std::size_t local=0;local<affected_rows.size();++local)
    {
        const auto row=static_cast<std::size_t>(affected_rows[local]);
        double prediction{};
        for(std::size_t entry=incidence.offsets[row];entry<incidence.offsets[row+1];++entry)
        {
            const auto & contributor=incidence.entries[entry];
            const auto atom_position=static_cast<std::size_t>(contributor.atom_position);
            if(!active[atom_position]) continue;
            const auto atom=static_cast<std::size_t>(layout.full_atoms[atom_position]);
            const auto & support=domain.atoms[atom][static_cast<std::size_t>(contributor.support_position)];
            const auto basis=joint::EvaluateKernel(support.square,widths(static_cast<Eigen::Index>(atom_position)),2.5);
            prediction+=candidate_beta(2*static_cast<Eigen::Index>(atom))*basis.gaussian+
                candidate_beta(2*static_cast<Eigen::Index>(atom)+1)*basis.charge;
        }
        out.prediction(static_cast<Eigen::Index>(local))=prediction;
        out.residual(static_cast<Eigen::Index>(local))=prediction-observations(affected_rows[local]);
    }
    double delta_squared{};
    for(std::size_t local=0;local<affected_rows.size();++local)
    {
        const double old_residual=accepted.residual(affected_rows[local]);
        const double new_residual=out.residual(static_cast<Eigen::Index>(local));
        delta_squared+=new_residual*new_residual-old_residual*old_residual;
    }
    out.objective=accepted.objective+.5*delta_squared/(scale*scale);
    std::size_t scratch{};
    if(!joint::CheckedMultiplyBytes(active.capacity(),sizeof(std::uint8_t),scratch))
        throw std::overflow_error("E2 candidate scratch estimate overflow.");
    std::size_t values{};
    if(!joint::CheckedMultiplyBytes(static_cast<std::size_t>(widths.size()+out.prediction.size()+out.residual.size()),
        sizeof(double),values) || !joint::CheckedAddBytes(scratch,values,out.scratch_bytes))
        throw std::overflow_error("E2 candidate scratch estimate overflow.");
    return out;
}
void CompareMethod(const rhbm_gem::JointParameterLayout & layout,
    const joint::BlockCoordinateState & accepted,const ReplayCandidate & candidate,
    const joint::FixedNeighborReplayResult & reference,
    const second_stage_test::ReplayCandidateResult & result,MethodAccuracy & accuracy)
{
    for(const auto row_value:layout.informative_rows)
    {
        const auto row=static_cast<Eigen::Index>(row_value);
        const auto local=AffectedPosition(candidate.affected_rows,row);
        const double prediction=local>=0 ? result.prediction(local) : accepted.prediction(row);
        const double residual=local>=0 ? result.residual(local) : accepted.residual(row);
        AddError(accuracy.prediction,reference.prediction(row),prediction);
        AddError(accuracy.residual,reference.residual(row),residual);
        if(accuracy.first_prediction_mismatch<0 && !BitwiseEqual(reference.prediction(row),prediction))
            accuracy.first_prediction_mismatch=row;
        if(accuracy.first_residual_mismatch<0 && !BitwiseEqual(reference.residual(row),residual))
            accuracy.first_residual_mismatch=row;
    }
    AddError(accuracy.objective,reference.objective,result.objective);
    accuracy.objective_bitwise_mismatches+=!BitwiseEqual(reference.objective,result.objective);
}
std::size_t CandidateSequenceBytes(const std::vector<ReplayCandidate> & candidates)
{
    std::size_t bytes=candidates.capacity()*sizeof(ReplayCandidate);
    for(const auto & candidate:candidates)
    {
        bytes+=candidate.core_atoms.capacity()*sizeof(Eigen::Index);
        bytes+=candidate.neighbor_atoms.capacity()*sizeof(Eigen::Index);
        bytes+=candidate.affected_rows.capacity()*sizeof(Eigen::Index);
        bytes+=static_cast<std::size_t>(candidate.eta.size()+candidate.beta.size())*sizeof(double);
    }
    return bytes;
}
void ReportAccuracy(const char * method,const MethodAccuracy & accuracy)
{
    const auto report=[](const char * name,const ErrorSummary & value) {
        std::cout<<' '<<name<<"_mismatch="<<value.bitwise_mismatches
            <<"_max_abs="<<value.maximum_absolute<<"_max_rel="<<value.maximum_relative
            <<"_max_ulp_normal="<<value.maximum_ulp<<"_normal_pairs="<<value.normal_pairs
            <<"_nonfinite_pairs="<<value.nonfinite_pairs;
    };
    std::cout<<"accuracy_"<<method;
    report("prediction",accuracy.prediction); report("residual",accuracy.residual);
    report("objective",accuracy.objective);
    std::cout<<"_objective_bitwise_mismatches="<<accuracy.objective_bitwise_mismatches
        <<"_first_prediction_mismatch_row="<<accuracy.first_prediction_mismatch
        <<"_first_residual_mismatch_row="<<accuracy.first_residual_mismatch<<'\n';
}
void ReportTimes(const char * method,const MethodTiming & timing,std::size_t candidates)
{
    const auto median=[](std::vector<double> values) {
        std::sort(values.begin(),values.end()); return values[values.size()/2];
    };
    const auto minmax=[](std::vector<double> values) {
        std::sort(values.begin(),values.end()); return std::pair{values.front(),values.back()};
    };
    const auto prep=median(timing.preparation)/static_cast<double>(candidates)*1000.0;
    const auto candidate=median(timing.candidate)/static_cast<double>(candidates)*1000.0;
    const auto total=median(timing.total)/static_cast<double>(candidates)*1000.0;
    const auto [minimum,maximum]=minmax(timing.total);
    std::cout<<"timing_"<<method<<"_median_ms_per_candidate="<<total
        <<" range_ms_per_candidate="<<minimum*1000.0/static_cast<double>(candidates)
        <<".."<<maximum*1000.0/static_cast<double>(candidates)
        <<" prep_ms_per_candidate="<<prep<<" replay_ms_per_candidate="<<candidate
        <<" temporary_bytes_estimate="<<timing.maximum_temporary_bytes
        <<" checksum="<<timing.checksum<<'\n';
}
void Run(int target_count)
{
    const auto fixture_started=Clock::now();
    rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",target_count));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    std::size_t accepted_state_bytes{};
    if(!joint::CheckedMultiplyBytes(static_cast<std::size_t>(data.domain.rows),2*sizeof(double),
        accepted_state_bytes)) throw std::overflow_error("Replay accepted-state memory estimate overflow.");
    Vector eta=Vector::Constant(static_cast<Eigen::Index>(data.layout.full_atoms.size()),std::log(.55));
    Vector beta(2*eta.size());
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {beta(2*atom)=2.0; beta(2*atom+1)=.2;}
    const auto initial=joint::FullReplayForTesting(data.domain,data.layout,data.y,eta,beta,data.context.scale);
    joint::BlockCoordinateState accepted{eta,beta,initial.prediction,initial.residual,initial.objective};
    joint::StructuralTopology topology(data.domain,data.context.atom_ids,data.layout);
    auto partition=topology.BuildCores(12); topology.PopulateMetadata(partition);
    std::vector<ReplayCandidate> candidates;
    candidates.reserve(partition.cores.size());
    for(std::size_t block=0;block<partition.cores.size();++block)
    {
        const auto & core=partition.cores[block];
        ReplayCandidate candidate;
        candidate.core_atoms=core.atoms; candidate.neighbor_atoms=core.neighbor_atoms;
        candidate.affected_rows=core.affected_rows; candidate.eta=accepted.eta; candidate.beta=accepted.beta;
        const double direction=block%2 ? -1.0 : 1.0;
        for(const auto atom:core.atoms)
        {
            candidate.eta(atom)+=direction*1e-4;
            candidate.beta(2*atom)+=direction*1e-3;
            candidate.beta(2*atom+1)-=direction*2e-4;
        }
        candidates.push_back(std::move(candidate));
    }
    const double fixture_seconds=Seconds(fixture_started);
    std::vector<std::size_t> affected_counts; affected_counts.reserve(candidates.size());
    for(const auto & candidate:candidates) affected_counts.push_back(candidate.affected_rows.size());
    std::sort(affected_counts.begin(),affected_counts.end());
    std::size_t support_entries{};
    std::vector<std::size_t> row_contributors(static_cast<std::size_t>(data.domain.rows));
    for(const auto atom:data.layout.full_atoms) for(const auto & support:data.domain.atoms[atom])
    {++support_entries; ++row_contributors[static_cast<std::size_t>(support.row)];}
    const auto contributor_total=std::accumulate(row_contributors.begin(),row_contributors.end(),std::size_t{});
    const auto sequence_bytes=CandidateSequenceBytes(candidates);
    const auto csr_started=Clock::now();
    const auto global_incidence=BuildGlobalReplayCsr(data.domain,data.layout);
    const double csr_build_seconds=Seconds(csr_started);
    const auto csr_build_temporary_bytes=static_cast<std::size_t>(data.domain.rows)*
        (sizeof(std::uint8_t)+2*sizeof(std::size_t));
    const auto csr_reference_started=Clock::now();
    const auto csr_reference=BuildGlobalReplayCsrReference(data.domain,data.layout);
    const double csr_reference_build_seconds=Seconds(csr_reference_started);
    std::size_t e2_offset_bytes{},e2_entry_bytes{};
    if(!joint::CheckedMultiplyBytes(csr_reference.offsets.capacity(),sizeof(std::size_t),e2_offset_bytes) ||
        !joint::CheckedMultiplyBytes(csr_reference.entries.capacity(),sizeof(GlobalReplayCsrEntry),e2_entry_bytes))
        throw std::overflow_error("E2 CSR component memory estimate overflow.");
    std::size_t e2_csr_storage_bytes{},e2_csr_build_peak_bytes{};
    if(!joint::CheckedAddBytes(e2_offset_bytes,e2_entry_bytes,e2_csr_storage_bytes) ||
        !joint::CheckedAddBytes(e2_csr_storage_bytes,csr_reference.build_temporary_bytes,e2_csr_build_peak_bytes) ||
        !joint::CheckedAddBytes(e2_csr_build_peak_bytes,accepted_state_bytes,e2_csr_build_peak_bytes))
        throw std::overflow_error("E2 CSR build peak estimate overflow.");
    std::cout<<std::setprecision(17)
        <<"fixture_targets="<<target_count<<" observations="<<data.domain.rows
        <<" informative_rows="<<data.layout.informative_rows.size()<<" support_entries="<<support_entries
        <<" contributors_mean="<<static_cast<double>(contributor_total)/static_cast<double>(row_contributors.size())
        <<" blocks="<<candidates.size()<<" affected_rows_min_median_max="
        <<affected_counts.front()<<":"<<affected_counts[affected_counts.size()/2]<<":"<<affected_counts.back()
        <<" candidate_sequence_bytes="<<sequence_bytes<<" fixture_and_sequence_seconds="<<fixture_seconds
        <<" accepted_state_bytes="<<accepted_state_bytes
        <<" d2_index_build_seconds="<<csr_build_seconds
        <<" d2_index_bytes="<<global_incidence.StorageBytes()
        <<" d2_index_build_peak_bytes_estimate="<<global_incidence.StorageBytes()+csr_build_temporary_bytes
        <<" e2_csr_build_seconds="<<csr_reference_build_seconds
        <<" e2_csr_offsets_bytes="<<e2_offset_bytes
        <<" e2_csr_entries_bytes="<<e2_entry_bytes
        <<" e2_csr_storage_bytes="<<e2_csr_storage_bytes
        <<" e2_csr_build_temporary_bytes="<<csr_reference.build_temporary_bytes
        <<" e2_csr_build_peak_bytes_estimate="<<e2_csr_build_peak_bytes<<'\n';

    std::array<MethodAccuracy,Methods.size()> accuracies{};
    const auto verification_started=Clock::now();
    for(const auto & candidate:candidates)
    {
        const auto reference=joint::FullReplayForTesting(data.domain,data.layout,data.y,
            candidate.eta,candidate.beta,data.context.scale);
        const auto incidence=second_stage_test::PrepareReplayIncidence(data.domain,data.layout,
            candidate.core_atoms,candidate.neighbor_atoms,candidate.affected_rows);
        for(std::size_t method_index=0;method_index<Methods.size();++method_index)
        {
            const auto & method=Methods[method_index];
            const auto result=method.global_csr
                ? EvaluateGlobalCsrCandidate(data.layout,data.y,accepted,
                    candidate.eta,candidate.beta,data.context.scale,candidate.affected_rows,global_incidence)
                : method.global_csr_reference
                    ? EvaluateGlobalCsrReferenceCandidate(data.domain,data.layout,data.y,accepted,
                        candidate.eta,candidate.beta,data.context.scale,candidate.core_atoms,
                        candidate.neighbor_atoms,candidate.affected_rows,csr_reference)
                    : second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
                        accepted,candidate.eta,candidate.beta,data.context.scale,candidate.core_atoms,
                        candidate.affected_rows,method.prediction,method.objective,
                        method.exact ? &incidence : nullptr);
            CompareMethod(data.layout,accepted,candidate,reference,result,accuracies[method_index]);
        }
    }
    const double verification_seconds=Seconds(verification_started);
    std::cout<<"full_reference_verification_seconds="<<verification_seconds<<'\n';
    for(std::size_t i=0;i<Methods.size();++i) ReportAccuracy(Methods[i].name,accuracies[i]);
    if(accuracies[6].prediction.bitwise_mismatches || accuracies[6].residual.bitwise_mismatches)
        throw std::runtime_error("E2 global CSR changed exact affected-row prediction or residual values.");

    constexpr std::size_t samples=5;
    std::array<MethodTiming,Methods.size()+1> timings{};
    for(std::size_t method_index=0;method_index<timings.size();++method_index)
    {
        auto & timing=timings[method_index];
        timing.preparation.reserve(samples); timing.candidate.reserve(samples); timing.total.reserve(samples);
        for(std::size_t sample=0;sample<=samples;++sample)
        {
            double preparation{},candidate_time{},checksum{};
            std::size_t max_temporary{};
            const auto total_started=Clock::now();
            for(const auto & replay_candidate:candidates)
            {
                if(method_index==0)
                {
                    const auto started=Clock::now();
                    auto replay=joint::FullReplayForTesting(data.domain,data.layout,data.y,
                        replay_candidate.eta,replay_candidate.beta,data.context.scale);
                    candidate_time+=Seconds(started); checksum+=replay.objective;
                    max_temporary=std::max(max_temporary,accepted_state_bytes);
                    continue;
                }
                const auto & method=Methods[method_index-1];
                second_stage_test::ReplayIncidence incidence;
                if(method.exact)
                {
                    const auto prep_started=Clock::now();
                    incidence=second_stage_test::PrepareReplayIncidence(data.domain,data.layout,
                        replay_candidate.core_atoms,replay_candidate.neighbor_atoms,replay_candidate.affected_rows);
                    preparation+=Seconds(prep_started);
                }
                const auto started=Clock::now();
                const auto result=method.global_csr
                    ? EvaluateGlobalCsrCandidate(data.layout,data.y,accepted,
                        replay_candidate.eta,replay_candidate.beta,data.context.scale,
                        replay_candidate.affected_rows,global_incidence)
                    : method.global_csr_reference
                        ? EvaluateGlobalCsrReferenceCandidate(data.domain,data.layout,data.y,accepted,
                            replay_candidate.eta,replay_candidate.beta,data.context.scale,
                            replay_candidate.core_atoms,replay_candidate.neighbor_atoms,
                            replay_candidate.affected_rows,csr_reference)
                        : second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
                            accepted,replay_candidate.eta,replay_candidate.beta,data.context.scale,
                            replay_candidate.core_atoms,replay_candidate.affected_rows,method.prediction,method.objective,
                            method.exact ? &incidence : nullptr);
                candidate_time+=Seconds(started); checksum+=result.objective;
                max_temporary=std::max(max_temporary,result.scratch_bytes);
            }
            const double total=Seconds(total_started);
            if(sample)
            {
                timing.preparation.push_back(preparation);
                timing.candidate.push_back(candidate_time);
                timing.total.push_back(total);
                timing.checksum=checksum;
                timing.maximum_temporary_bytes=std::max(timing.maximum_temporary_bytes,max_temporary);
            }
        }
    }
    ReportTimes("A",timings[0],candidates.size());
    for(std::size_t i=0;i<Methods.size();++i) ReportTimes(Methods[i].name,timings[i+1],candidates.size());
    const auto median_seconds=[](std::vector<double> values) {
        std::sort(values.begin(),values.end()); return values[values.size()/2];
    };
    const double d1_per_candidate=median_seconds(timings[5].total)/
        static_cast<double>(candidates.size());
    const double d2_per_candidate=median_seconds(timings[4].total)/
        static_cast<double>(candidates.size());
    const double d2_amortized=d2_per_candidate+csr_build_seconds/static_cast<double>(candidates.size());
    std::cout<<"d2_amortized_ms_per_candidate="<<d2_amortized*1000.0
        <<" build_amortized_ms_per_candidate="<<csr_build_seconds*1000.0/static_cast<double>(candidates.size())
        <<" d2_vs_d1_savings_ms_per_candidate="<<(d1_per_candidate-d2_amortized)*1000.0;
    if(d1_per_candidate>d2_per_candidate)
        std::cout<<" break_even_candidates_vs_d1="<<std::ceil(csr_build_seconds/
            (d1_per_candidate-d2_per_candidate))
            <<" break_even_sweeps_vs_d1="<<csr_build_seconds/
                ((d1_per_candidate-d2_per_candidate)*static_cast<double>(candidates.size()));
    else
        std::cout<<" break_even_candidates_vs_d1=never break_even_sweeps_vs_d1=never";
    const auto e1_candidate=median_seconds(timings[6].total)/static_cast<double>(candidates.size());
    const auto e2_candidate=median_seconds(timings[7].total)/static_cast<double>(candidates.size());
    std::cout<<"\ne2_break_even_candidates_vs_e1=";
    if(e1_candidate>e2_candidate)
        std::cout<<std::ceil(csr_reference_build_seconds/(e1_candidate-e2_candidate))
            <<" break_even_sweeps="<<csr_reference_build_seconds/
                ((e1_candidate-e2_candidate)*static_cast<double>(candidates.size()));
    else std::cout<<"never break_even_sweeps=never";
    std::cout<<" e2_candidate_savings_ms="<<(e1_candidate-e2_candidate)*1000.0<<'\n';
    std::cout<<'\n';
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)==0)
    {
#if defined(__APPLE__)
        const auto peak=static_cast<std::uint64_t>(usage.ru_maxrss);
#else
        const auto peak=static_cast<std::uint64_t>(usage.ru_maxrss)*1024;
#endif
        std::cout<<"process_peak_rss_bytes="<<peak<<'\n';
    }
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=2) throw std::invalid_argument("usage: rhbm_fixed_neighbor_replay_benchmark <target-count>");
        const int target_count=std::stoi(argv[1]);
        Run(target_count);
    }
    catch(const std::exception & error)
    {std::cerr<<"fixed-neighbor replay benchmark failed: "<<error.what()<<'\n'; return 1;}
    return 0;
}
