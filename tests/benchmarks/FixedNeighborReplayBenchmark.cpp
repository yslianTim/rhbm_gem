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
    bool exact;
};
constexpr std::array<Method,5> Methods{{
    {"B1",second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false},
    {"B2",second_stage_test::ReplayPredictionStrategy::DeltaSubtractThenAdd,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false},
    {"C",second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayObjectiveStrategy::Incremental,false},
    {"D1",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::FullScan,true},
    {"E",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::Incremental,true}
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
    std::cout<<std::setprecision(17)
        <<"fixture_targets="<<target_count<<" observations="<<data.domain.rows
        <<" informative_rows="<<data.layout.informative_rows.size()<<" support_entries="<<support_entries
        <<" contributors_mean="<<static_cast<double>(contributor_total)/static_cast<double>(row_contributors.size())
        <<" blocks="<<candidates.size()<<" affected_rows_min_median_max="
        <<affected_counts.front()<<":"<<affected_counts[affected_counts.size()/2]<<":"<<affected_counts.back()
        <<" candidate_sequence_bytes="<<sequence_bytes<<" fixture_and_sequence_seconds="<<fixture_seconds
        <<" accepted_state_bytes="<<2*static_cast<std::size_t>(data.domain.rows)*sizeof(double)<<'\n';

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
            const auto result=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
                accepted,candidate.eta,candidate.beta,data.context.scale,candidate.core_atoms,
                candidate.affected_rows,method.prediction,method.objective,
                method.exact ? &incidence : nullptr);
            CompareMethod(data.layout,accepted,candidate,reference,result,accuracies[method_index]);
        }
    }
    const double verification_seconds=Seconds(verification_started);
    std::cout<<"full_reference_verification_seconds="<<verification_seconds<<'\n';
    for(std::size_t i=0;i<Methods.size();++i) ReportAccuracy(Methods[i].name,accuracies[i]);

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
                    max_temporary=std::max(max_temporary,
                        2*static_cast<std::size_t>(data.domain.rows)*sizeof(double));
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
                const auto result=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
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
