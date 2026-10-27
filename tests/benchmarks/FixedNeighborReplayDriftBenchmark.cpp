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
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;
using Vector=Eigen::VectorXd;

struct Method
{
    const char * name;
    second_stage_test::ReplayPredictionStrategy prediction;
    second_stage_test::ReplayObjectiveStrategy objective;
    bool full_replay{},exact{},tiled{};
    std::size_t tile_size{};
    bool hierarchical{};
};
const std::vector<Method> Methods{
    {"A",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::FullScan,true,true},
    {"B1",second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false,false},
    {"B2",second_stage_test::ReplayPredictionStrategy::DeltaSubtractThenAdd,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false,false},
    {"C",second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayObjectiveStrategy::Incremental,false,false},
    {"D1",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::FullScan,false,true},
    {"E",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::Incremental,false,true},
    {"E_scaled",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::IncrementalScaled,false,true},
    {"E_kahan",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::IncrementalKahan,false,true},
    {"E_neumaier",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::IncrementalNeumaier,false,true},
    {"E_kahan_scaled",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::IncrementalKahanScaled,false,true},
    {"E_neumaier_scaled",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::IncrementalNeumaierScaled,false,true},
    {"tile_256",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,256,false},
    {"tile_512",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,512,false},
    {"tile_1024",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,1024,false},
    {"tile_2048",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,2048,false},
    {"tree_256",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,256,true},
    {"tree_512",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,512,true},
    {"tree_1024",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,1024,true},
    {"tree_2048",second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::CachedExternally,false,true,true,2048,true}
};
struct Candidate {Vector eta,beta; bool accept{};};
struct Metrics
{
    std::size_t prediction_mismatches{},residual_mismatches{},objective_mismatches{},acceptance_disagreements{},
        rollback_checks{},refresh_count{},refresh_mismatches{};
    Eigen::Index first_prediction_mismatch=-1,first_residual_mismatch=-1;
    std::uint64_t maximum_prediction_ulp{},maximum_residual_ulp{},maximum_objective_ulp{};
    double maximum_prediction_absolute{},maximum_prediction_relative{},maximum_residual_absolute{},
        maximum_residual_relative{},maximum_objective_absolute{},maximum_objective_relative{},
        maximum_refresh_pre_prediction{},maximum_refresh_pre_residual{},maximum_refresh_pre_objective{};
    double first_half_max_prediction{},second_half_max_prediction{},first_half_max_objective{},second_half_max_objective{};
    std::size_t temporary_bytes{};
    std::size_t objective_cache_bytes{};
    std::size_t objective_cache_update_bytes{};
    std::size_t objective_cache_tiles_recomputed{};
    double objective_cache_evaluate_seconds{},objective_cache_commit_seconds{};
    double maximum_delta_cancellation{};
};

double Seconds(Clock::time_point started)
{return std::chrono::duration<double>(Clock::now()-started).count();}
bool BitwiseEqual(double left,double right)
{return std::bit_cast<std::uint64_t>(left)==std::bit_cast<std::uint64_t>(right);}
bool BitwiseEqual(const Vector & left,const Vector & right)
{
    if(left.size()!=right.size()) return false;
    for(Eigen::Index i=0;i<left.size();++i) if(!BitwiseEqual(left(i),right(i))) return false;
    return true;
}
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
Eigen::Index AffectedPosition(const joint::Indices & rows,Eigen::Index row)
{
    const auto found=std::lower_bound(rows.begin(),rows.end(),row);
    return found==rows.end() || *found!=row ? -1 : static_cast<Eigen::Index>(found-rows.begin());
}
std::vector<Candidate> CandidateSequence(const Vector & initial_eta,const Vector & initial_beta,
    const joint::Indices & core_atoms,double eta_upper_bound)
{
    constexpr std::size_t attempts=503;
    std::vector<Candidate> output; output.reserve(attempts);
    Vector eta=initial_eta,beta=initial_beta;
    for(std::size_t attempt=1;attempt<=attempts;++attempt)
    {
        Vector candidate_eta=eta,candidate_beta=beta;
        const double direction=((attempt/50)%2) ? -1.0 : 1.0;
        for(const auto atom_value:core_atoms)
        {
            const auto atom=static_cast<Eigen::Index>(atom_value);
            if(attempt<=200)
            {
                candidate_beta(2*atom)+=std::clamp(2.0-candidate_beta(2*atom),-1e-3,1e-3);
                candidate_beta(2*atom+1)+=std::clamp(.2-candidate_beta(2*atom+1),-1e-3,1e-3);
            }
            else
            {
                candidate_beta(2*atom)+=direction*1e-12;
                candidate_beta(2*atom+1)-=direction*1e-12;
            }
            candidate_eta(atom)+=attempt%50==0 ? direction*.2 : direction*1e-12;
            if(attempt%97==0) candidate_eta(atom)=eta_upper_bound-1e-12;
        }
        const bool accept=attempt%17!=0;
        output.push_back({candidate_eta,candidate_beta,accept});
        if(accept) {eta=std::move(candidate_eta); beta=std::move(candidate_beta);}
    }
    return output;
}
void CompareCandidate(const rhbm_gem::JointParameterLayout & layout,const joint::BlockCoordinateState & accepted,
    const joint::Indices & affected_rows,const joint::FixedNeighborReplayResult & reference,
    const second_stage_test::ReplayCandidateResult & result,bool full_replay,Metrics & metrics,std::size_t attempt)
{
    double current_prediction_max{},current_objective{};
    for(const auto row_value:layout.informative_rows)
    {
        const auto row=static_cast<Eigen::Index>(row_value);
        const auto local=AffectedPosition(affected_rows,row);
        const double prediction=full_replay ? result.prediction(row) :
            (local>=0 ? result.prediction(local) : accepted.prediction(row));
        const double residual=full_replay ? result.residual(row) :
            (local>=0 ? result.residual(local) : accepted.residual(row));
        const double prediction_error=std::abs(prediction-reference.prediction(row));
        const double residual_error=std::abs(residual-reference.residual(row));
        metrics.maximum_prediction_absolute=std::max(metrics.maximum_prediction_absolute,prediction_error);
        metrics.maximum_prediction_relative=std::max(metrics.maximum_prediction_relative,
            prediction_error/std::max(std::abs(reference.prediction(row)),std::numeric_limits<double>::min()));
        metrics.maximum_residual_absolute=std::max(metrics.maximum_residual_absolute,residual_error);
        metrics.maximum_residual_relative=std::max(metrics.maximum_residual_relative,
            residual_error/std::max(std::abs(reference.residual(row)),std::numeric_limits<double>::min()));
        current_prediction_max=std::max(current_prediction_max,prediction_error);
        if(!BitwiseEqual(prediction,reference.prediction(row)))
        {
            ++metrics.prediction_mismatches;
            if(metrics.first_prediction_mismatch<0) metrics.first_prediction_mismatch=row;
        }
        if(!BitwiseEqual(residual,reference.residual(row)))
        {
            ++metrics.residual_mismatches;
            if(metrics.first_residual_mismatch<0) metrics.first_residual_mismatch=row;
        }
        if(std::isnormal(prediction) && std::isnormal(reference.prediction(row)))
            metrics.maximum_prediction_ulp=std::max(metrics.maximum_prediction_ulp,
                UlpDistance(prediction,reference.prediction(row)));
        if(std::isnormal(residual) && std::isnormal(reference.residual(row)))
            metrics.maximum_residual_ulp=std::max(metrics.maximum_residual_ulp,
                UlpDistance(residual,reference.residual(row)));
    }
    double delta_sum{},delta_absolute_sum{};
    for(std::size_t k=0;k<affected_rows.size();++k)
    {
        const double old_residual=accepted.residual(affected_rows[k]);
        const double new_residual=full_replay ? reference.residual(affected_rows[k]) :
            result.residual(static_cast<Eigen::Index>(k));
        const double delta=new_residual*new_residual-old_residual*old_residual;
        delta_sum+=delta; delta_absolute_sum+=std::abs(delta);
    }
    if(delta_absolute_sum>0)
        metrics.maximum_delta_cancellation=std::max(metrics.maximum_delta_cancellation,
            delta_sum==0 ? std::numeric_limits<double>::infinity() : delta_absolute_sum/std::abs(delta_sum));
    const double objective_error=std::abs(result.objective-reference.objective);
    metrics.maximum_objective_absolute=std::max(metrics.maximum_objective_absolute,objective_error);
    metrics.maximum_objective_relative=std::max(metrics.maximum_objective_relative,
        objective_error/std::max(std::abs(reference.objective),std::numeric_limits<double>::min()));
    current_objective=objective_error;
    if(!BitwiseEqual(result.objective,reference.objective)) ++metrics.objective_mismatches;
    if(std::isnormal(result.objective) && std::isnormal(reference.objective))
        metrics.maximum_objective_ulp=std::max(metrics.maximum_objective_ulp,
            UlpDistance(result.objective,reference.objective));
    if(attempt<=251)
    {
        metrics.first_half_max_prediction=std::max(metrics.first_half_max_prediction,current_prediction_max);
        metrics.first_half_max_objective=std::max(metrics.first_half_max_objective,current_objective);
    }
    else
    {
        metrics.second_half_max_prediction=std::max(metrics.second_half_max_prediction,current_prediction_max);
        metrics.second_half_max_objective=std::max(metrics.second_half_max_objective,current_objective);
    }
}
void Run(int target_count)
{
    constexpr std::array<std::size_t,5> refresh_intervals{1,10,50,100,11};
    constexpr std::size_t attempts=503;
    rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",target_count));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    Vector eta=Vector::Constant(static_cast<Eigen::Index>(data.layout.full_atoms.size()),std::log(.5));
    Vector beta=Vector::Zero(2*eta.size());
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {beta(2*atom)=2.0; beta(2*atom+1)=atom%2 ? -.2 : .2;}
    joint::StructuralTopology topology(data.domain,data.context.atom_ids,data.layout);
    auto partition=topology.BuildCores(12); topology.PopulateMetadata(partition);
    if(partition.cores.empty()) throw std::runtime_error("Synthetic workload produced no blocks.");
    const auto & block=partition.cores.front();
    const double eta_upper_bound=std::log(2.2);
    const auto sequence=CandidateSequence(eta,beta,block.atoms,eta_upper_bound);
    const auto initial=joint::FullReplayForTesting(data.domain,data.layout,data.y,eta,beta,data.context.scale);
    std::size_t support_entries{};
    for(const auto atom:data.layout.full_atoms) support_entries+=data.domain.atoms[static_cast<std::size_t>(atom)].size();
    std::size_t accepted_attempts{},rejected_attempts{};
    for(const auto & candidate:sequence) candidate.accept ? ++accepted_attempts : ++rejected_attempts;
    std::cout<<std::setprecision(17)<<"fixture_targets="<<target_count
        <<" observations="<<data.domain.rows<<" informative_rows="<<data.layout.informative_rows.size()
        <<" support_entries="<<support_entries<<" blocks="<<partition.cores.size()
        <<" repeated_block=1 affected_rows="<<block.affected_rows.size()
        <<" attempts="<<attempts<<" scripted_accepted_updates="<<accepted_attempts
        <<" scripted_rejections="<<rejected_attempts<<" refresh_interval_sweep=11\n";

    for(std::size_t method_index=0;method_index<Methods.size();++method_index)
    {
        const auto & method=Methods[method_index];
        for(const auto interval:refresh_intervals)
        {
            if(method.full_replay && interval!=refresh_intervals.front()) continue;
            joint::BlockCoordinateState state{eta,beta,initial.prediction,initial.residual,initial.objective};
            double accepted_reference_objective=initial.objective;
            std::size_t accepted_since_refresh{};
            std::optional<second_stage_test::TiledObjectiveCache> objective_cache;
            if(method.tiled)
                objective_cache.emplace(data.layout.informative_rows,state.residual,
                    method.tile_size,method.hierarchical);
            Metrics metrics;
            double measured_seconds{};
            for(std::size_t attempt=0;attempt<sequence.size();++attempt)
            {
                const auto & candidate=sequence[attempt];
                const auto reference=joint::FullReplayForTesting(data.domain,data.layout,data.y,
                    candidate.eta,candidate.beta,data.context.scale);
                const double reference_delta=reference.objective-accepted_reference_objective;
                const bool full_objective_acceptance=reference_delta<=0;
                const auto before_rejection=state;
                second_stage_test::ReplayCandidateResult result;
                second_stage_test::ObjectiveCacheCandidate objective_cache_candidate;
                const auto candidate_started=Clock::now();
                if(method.full_replay)
                {
                    const auto full=joint::FullReplayForTesting(data.domain,data.layout,data.y,
                        candidate.eta,candidate.beta,data.context.scale);
                    result.prediction=full.prediction; result.residual=full.residual; result.objective=full.objective;
                    result.scratch_bytes=2*static_cast<std::size_t>(data.domain.rows)*sizeof(double);
                }
                else
                {
                    second_stage_test::ReplayIncidence incidence;
                    if(method.exact)
                    {
                        incidence=second_stage_test::PrepareReplayIncidence(data.domain,data.layout,
                            block.atoms,block.neighbor_atoms,block.affected_rows);
                    }
                    result=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
                        state,candidate.eta,candidate.beta,data.context.scale,block.atoms,block.affected_rows,
                        method.prediction,method.tiled ? second_stage_test::ReplayObjectiveStrategy::CachedExternally :
                            method.objective,method.exact ? &incidence : nullptr);
                    if(method.tiled)
                    {
                        const auto cache_started=Clock::now();
                        objective_cache_candidate=objective_cache->Evaluate(state.residual,block.affected_rows,
                            result.residual,data.context.scale);
                        metrics.objective_cache_evaluate_seconds+=Seconds(cache_started);
                        metrics.objective_cache_tiles_recomputed+=objective_cache_candidate.tiles_recomputed;
                        result.objective=objective_cache_candidate.objective;
                        metrics.objective_cache_update_bytes=std::max(metrics.objective_cache_update_bytes,
                            objective_cache_candidate.updates.capacity()*sizeof(
                                std::pair<std::size_t,double>));
                    }
                }
                const double method_delta=result.objective-state.objective;
                const double method_scale=std::max(std::abs(result.objective),std::abs(state.objective));
                const double method_enclosure=joint::BlockObjectiveReplayEnclosure(method_scale);
                const double local_before=1.0,local_after=local_before+reference_delta;
                const double consistency_error=std::abs(method_delta+local_before-local_after);
                const bool method_objective_acceptance=joint::WithinBlockObjectiveReplay(consistency_error,method_scale) &&
                    (method_delta<=0 || (local_after<=local_before && method_delta<=method_enclosure));
                metrics.acceptance_disagreements+=method_objective_acceptance!=full_objective_acceptance;
                measured_seconds+=Seconds(candidate_started);
                metrics.temporary_bytes=std::max(metrics.temporary_bytes,result.scratch_bytes);
                CompareCandidate(data.layout,state,block.affected_rows,reference,result,method.full_replay,metrics,attempt+1);

                if(candidate.accept)
                {
                    const auto state_update_started=Clock::now();
                    state.eta=candidate.eta; state.beta=candidate.beta;
                    if(method.full_replay)
                    {
                        state.prediction=result.prediction; state.residual=result.residual; state.objective=result.objective;
                    }
                    else
                    {
                        for(std::size_t k=0;k<block.affected_rows.size();++k)
                        {
                            const auto row=block.affected_rows[k];
                            state.prediction(row)=result.prediction(static_cast<Eigen::Index>(k));
                            state.residual(row)=result.residual(static_cast<Eigen::Index>(k));
                        }
                        state.objective=result.objective;
                    }
                    if(method.tiled)
                    {
                        const auto cache_started=Clock::now();
                        objective_cache->Commit(objective_cache_candidate);
                        metrics.objective_cache_commit_seconds+=Seconds(cache_started);
                    }
                    accepted_reference_objective=reference.objective;
                    ++accepted_since_refresh;
                    if(!method.full_replay && accepted_since_refresh==interval)
                    {
                        const auto refreshed=joint::FullReplayForTesting(data.domain,data.layout,data.y,
                            state.eta,state.beta,data.context.scale);
                        const double refresh_prediction_error=(state.prediction-refreshed.prediction).cwiseAbs().maxCoeff();
                        const double refresh_residual_error=(state.residual-refreshed.residual).cwiseAbs().maxCoeff();
                        const double refresh_objective_error=std::abs(state.objective-refreshed.objective);
                        metrics.maximum_refresh_pre_prediction=std::max(metrics.maximum_refresh_pre_prediction,
                            refresh_prediction_error);
                        metrics.maximum_refresh_pre_residual=std::max(metrics.maximum_refresh_pre_residual,
                            refresh_residual_error);
                        metrics.maximum_refresh_pre_objective=std::max(metrics.maximum_refresh_pre_objective,
                            refresh_objective_error);
                        state.prediction=refreshed.prediction; state.residual=refreshed.residual;
                        state.objective=refreshed.objective; accepted_since_refresh=0; ++metrics.refresh_count;
                        if(method.tiled)
                            objective_cache.emplace(data.layout.informative_rows,state.residual,
                                method.tile_size,method.hierarchical);
                        metrics.refresh_mismatches+=!BitwiseEqual(state.prediction,refreshed.prediction) ||
                            !BitwiseEqual(state.residual,refreshed.residual) ||
                            !BitwiseEqual(state.objective,refreshed.objective);
                    }
                    measured_seconds+=Seconds(state_update_started);
                }
                else
                {
                    ++metrics.rollback_checks;
                    if(!BitwiseEqual(state.eta,before_rejection.eta) || !BitwiseEqual(state.beta,before_rejection.beta) ||
                        !BitwiseEqual(state.prediction,before_rejection.prediction) ||
                        !BitwiseEqual(state.residual,before_rejection.residual) ||
                        !BitwiseEqual(state.objective,before_rejection.objective))
                        throw std::runtime_error("Rejected synthetic candidate changed accepted shadow state.");
                }
            }
            const auto resident_state_bytes=2*static_cast<std::size_t>(data.domain.rows)*sizeof(double);
            const auto full_replay_scratch_bytes=resident_state_bytes;
            if(objective_cache) metrics.objective_cache_bytes=objective_cache->StorageBytes();
            const auto estimated_peak_bytes=resident_state_bytes+metrics.objective_cache_bytes+
                std::max(metrics.temporary_bytes+metrics.objective_cache_update_bytes,full_replay_scratch_bytes);
            std::cout<<"method="<<method.name<<" refresh_interval="<<interval
                <<" refresh_count="<<metrics.refresh_count<<" elapsed_seconds="<<measured_seconds
                <<" milliseconds_per_attempt="<<measured_seconds*1000.0/static_cast<double>(attempts)
                <<" candidate_temporary_bytes_estimate="<<metrics.temporary_bytes
                <<" objective_cache_bytes="<<metrics.objective_cache_bytes
                <<" objective_cache_update_bytes="<<metrics.objective_cache_update_bytes
                <<" objective_cache_tiles_recomputed="<<metrics.objective_cache_tiles_recomputed
                <<" objective_cache_evaluate_seconds="<<metrics.objective_cache_evaluate_seconds
                <<" objective_cache_commit_seconds="<<metrics.objective_cache_commit_seconds
                <<" estimated_peak_bytes="<<estimated_peak_bytes
                <<" prediction_mismatches="<<metrics.prediction_mismatches
                <<" residual_mismatches="<<metrics.residual_mismatches
                <<" objective_mismatches="<<metrics.objective_mismatches
                <<" max_prediction_abs="<<metrics.maximum_prediction_absolute
                <<" max_prediction_rel="<<metrics.maximum_prediction_relative
                <<" max_prediction_ulp="<<metrics.maximum_prediction_ulp
                <<" max_residual_abs="<<metrics.maximum_residual_absolute
                <<" max_residual_rel="<<metrics.maximum_residual_relative
                <<" max_residual_ulp="<<metrics.maximum_residual_ulp
                <<" max_objective_abs="<<metrics.maximum_objective_absolute
                <<" max_objective_rel="<<metrics.maximum_objective_relative
                <<" max_objective_ulp="<<metrics.maximum_objective_ulp
                <<" max_delta_cancellation_ratio="<<metrics.maximum_delta_cancellation
                <<" first_half_max_prediction="<<metrics.first_half_max_prediction
                <<" second_half_max_prediction="<<metrics.second_half_max_prediction
                <<" first_half_max_objective="<<metrics.first_half_max_objective
                <<" second_half_max_objective="<<metrics.second_half_max_objective
                <<" max_pre_refresh_prediction_error="<<metrics.maximum_refresh_pre_prediction
                <<" max_pre_refresh_residual_error="<<metrics.maximum_refresh_pre_residual
                <<" max_pre_refresh_objective_error="<<metrics.maximum_refresh_pre_objective
                <<" objective_acceptance_disagreements="<<metrics.acceptance_disagreements
                <<" rollback_checks="<<metrics.rollback_checks
                <<" refresh_mismatches="<<metrics.refresh_mismatches
                <<" first_prediction_mismatch_row="<<metrics.first_prediction_mismatch
                <<" first_residual_mismatch_row="<<metrics.first_residual_mismatch<<'\n';
        }
    }
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
        if(argc!=2) throw std::invalid_argument("usage: rhbm_fixed_neighbor_replay_drift_benchmark <target-count>");
        Run(std::stoi(argv[1]));
    }
    catch(const std::exception & error)
    {std::cerr<<"fixed-neighbor replay drift benchmark failed: "<<error.what()<<'\n'; return 1;}
    return 0;
}
