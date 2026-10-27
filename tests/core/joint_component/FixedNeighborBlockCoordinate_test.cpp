#include <gtest/gtest.h>
#include <boost/hash2/sha2.hpp>
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/JointComponentChecks.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include "support/JointReplayExperiment.hpp"
#include <algorithm>
#include <array>
#include "support/JointSyntheticWorkload.hpp"
#include <bit>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using rhbm_gem::core::JointProblem;
void CheckSingleBlock(const std::string & topology,int atoms)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,atoms));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(atoms,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=30;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    std::vector<n::FixedNeighborBlockSweep> observed;
    n::BlockCoordinateState observed_state; std::size_t state_callbacks{};
    policy.sweep_observer=[&](const auto & sweep){observed.push_back(sweep);};
    policy.state_observer=[&](std::size_t,const auto & state,const auto &,const auto &) {
        observed_state=state; ++state_callbacks;
    };
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    auto unobserved_policy=policy; unobserved_policy.sweep_observer={}; unobserved_policy.state_observer={};
    const auto unobserved=n::RunFixedNeighborDiagnostic(problem,eta,unobserved_policy);
    ASSERT_TRUE(result.search_converged)<<result.reason;
    ASSERT_TRUE(result.endpoint_certified)<<result.endpoint_trust.reason;
    EXPECT_GE(result.sweeps.size(),2u);
    EXPECT_GT(result.first_order_stationarity_sweep,0u);
    EXPECT_GE(result.confirmed_stationarity_sweep,result.first_order_stationarity_sweep);
    ASSERT_EQ(result.fit.components.size(),1u);
    EXPECT_EQ(result.fit.components.front().search_sweeps,result.sweeps.size());
    EXPECT_EQ(result.fit.components.front().first_order_stationarity_sweep,result.first_order_stationarity_sweep);
    EXPECT_EQ(result.fit.components.front().confirmed_stationarity_sweep,result.confirmed_stationarity_sweep);
    EXPECT_EQ(result.fit.components.front().endpoint_certification_sweep,result.sweeps.size());
    EXPECT_EQ(result.fit.components.front().endpoint_certification_sweep,result.confirmed_stationarity_sweep);
    ASSERT_EQ(observed.size(),result.sweeps.size());
    EXPECT_EQ(state_callbacks,result.sweeps.size());
    EXPECT_TRUE((observed_state.eta.array()==result.state.eta.array()).all());
    EXPECT_TRUE((observed_state.beta.array()==result.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(observed.back().global_ac_kkt,result.sweeps.back().global_ac_kkt);
    EXPECT_LE(result.sweeps.back().global_ac_kkt,policy.convergence_criteria.global_ac_kkt);
    EXPECT_LE(result.sweeps.back().global_width_gradient_inf_norm,policy.convergence_criteria.global_width_gradient);
    EXPECT_LE(result.sweeps.back().cache_replay_error,2e-12);
    EXPECT_TRUE(result.assessment.primary.valid);
    EXPECT_TRUE(result.assessment.independent_replay.passed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(result.fit.search_completed);
    const auto bounded=second_stage_test::matched::joint_abc::FitBoundedJointComponentsForTesting(
        problem,std::vector<double>(static_cast<std::size_t>(atoms),.55),2.2);
    ASSERT_EQ(bounded.components.size(),1u);
    ASSERT_TRUE(bounded.components.front().candidate_diagnostic);
    EXPECT_EQ(bounded.components.front().candidate_diagnostic->convergence_criteria_id,"practical-v1");
    ASSERT_TRUE(bounded.components.front().state);
    ASSERT_TRUE(result.fit.components.front().state);
    ASSERT_TRUE(bounded.objective && result.fit.objective);
    EXPECT_EQ(bounded.components.front().search_sweeps,result.fit.components.front().search_sweeps);
    EXPECT_EQ(bounded.components.front().stop_reason,"block-stationary");
    EXPECT_LE(std::abs(bounded.objective.value()-result.fit.objective.value()),1e-12);
    for(std::size_t k=0;k<bounded.components.front().state->ac.size();++k)
        EXPECT_LE(std::abs(bounded.components.front().state->ac[k]-result.fit.components.front().state->ac[k]),
            1e-10*(1+std::abs(result.fit.components.front().state->ac[k])));
    for(std::size_t k=0;k<bounded.components.front().state->log_b.size();++k)
        EXPECT_LE(std::abs(bounded.components.front().state->log_b[k]-result.fit.components.front().state->log_b[k]),1e-10);
    EXPECT_EQ(result.search_converged,unobserved.search_converged);
    EXPECT_EQ(result.endpoint_certified,unobserved.endpoint_certified);
    EXPECT_EQ(result.reason,unobserved.reason);
    EXPECT_EQ(result.sweeps.size(),unobserved.sweeps.size());
    EXPECT_EQ(result.state.eta.size(),unobserved.state.eta.size());
    EXPECT_TRUE((result.state.eta.array()==unobserved.state.eta.array()).all());
    EXPECT_TRUE((result.state.beta.array()==unobserved.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(result.state.objective,unobserved.state.objective);
    auto stopped_search=result.fit; stopped_search.search_completed=false;
    EXPECT_EQ(stopped_search.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    for(const auto & sweep:result.sweeps)
        EXPECT_LE(sweep.objective_after,sweep.objective_before+n::BlockObjectiveReplayEnclosure(sweep.objective_before));
    EXPECT_EQ(data.partition.components.size(),1u);
}

JointProblem PureExclusionChain(std::size_t atoms)
{
    auto input=second_stage_test::SyntheticJointWorkload("chain",static_cast<int>(atoms));
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(rhbm_gem::PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices.resize(atoms);
    std::iota(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),0);
    input.selection_domain->pure_halo_support_exclusion=rhbm_gem::JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    return JointProblem(std::move(input));
}

bool BitwiseEqual(double left,double right)
{return std::bit_cast<std::uint64_t>(left)==std::bit_cast<std::uint64_t>(right);}
bool BitwiseEqual(const n::Vector & left,const n::Vector & right)
{
    if(left.size()!=right.size()) return false;
    for(Eigen::Index i=0;i<left.size();++i)
        if(!BitwiseEqual(left(i),right(i))) return false;
    return true;
}
std::string Precise(double value)
{std::ostringstream out; out<<std::setprecision(17)<<value; return out.str();}
std::uint64_t UlpDistance(double left,double right)
{
    const auto ordered=[](double value) {
        const auto bits=std::bit_cast<std::uint64_t>(value);
        return bits&(std::uint64_t{1}<<63) ? ~bits : bits|(std::uint64_t{1}<<63);
    };
    const auto a=ordered(left),b=ordered(right);
    return a>b ? a-b : b-a;
}
bool ObjectiveDecision(double before,double after,double local_before,double local_after,bool unchanged)
{
    const double delta=after-before;
    const double enclosure=n::BlockObjectiveReplayEnclosure(std::max(std::abs(before),std::abs(after)));
    const double consistency=std::abs(delta+local_before-local_after);
    return n::WithinBlockObjectiveReplay(consistency,std::max(std::abs(before),std::abs(after))) &&
        (delta<=0 || (local_after<=local_before && delta<=enclosure)) && !unchanged;
}
struct IncrementalObjectiveAudit
{
    std::size_t candidates{},accepted{},rejected{},fallbacks{},estimated_reference_disagreements{},
        actual_reference_disagreements{},reference_positive_accepts{};
    double maximum_absolute{},maximum_relative{},maximum_cache_drift{},minimum_decision_margin{
        std::numeric_limits<double>::infinity()};
    std::uint64_t maximum_ulp{};
};
void AuditIncrementalCandidate(const rhbm_gem::JointParameterLayout & layout,double scale,
    const n::FixedNeighborCandidateEvent & event,IncrementalObjectiveAudit & audit)
{
    if(!event.candidate_replay_ran || !std::isfinite(event.incremental_objective)) return;
    ++audit.candidates; audit.accepted+=event.accepted; audit.rejected+=!event.accepted;
    audit.fallbacks+=event.objective_full_scan_fallback;
    double old_reference_squared{},candidate_reference_squared{};
    for(const auto row_value:layout.informative_rows)
    {
        const auto row=static_cast<Eigen::Index>(row_value);
        const double old_residual=event.accepted_before.residual(row);
        const double candidate_residual=event.full_residual(row);
        old_reference_squared+=old_residual*old_residual;
        candidate_reference_squared+=candidate_residual*candidate_residual;
    }
    const double old_reference=.5*old_reference_squared/(scale*scale);
    const double candidate_reference=.5*candidate_reference_squared/(scale*scale);
    const double error=std::abs(event.incremental_objective-candidate_reference);
    audit.maximum_absolute=std::max(audit.maximum_absolute,error);
    audit.maximum_relative=std::max(audit.maximum_relative,
        error/std::max(std::abs(candidate_reference),std::numeric_limits<double>::min()));
    audit.maximum_cache_drift=std::max(audit.maximum_cache_drift,
        std::abs(event.accepted_before.objective-old_reference));
    if(std::isnormal(event.incremental_objective) && std::isnormal(candidate_reference))
        audit.maximum_ulp=std::max(audit.maximum_ulp,UlpDistance(event.incremental_objective,candidate_reference));
    const bool estimated_decision=ObjectiveDecision(event.accepted_before.objective,event.incremental_objective,
        event.local_objective_before,event.local_objective_after,event.unchanged);
    const bool reference_decision=ObjectiveDecision(old_reference,candidate_reference,
        event.local_objective_before,event.local_objective_after,event.unchanged);
    audit.estimated_reference_disagreements+=estimated_decision!=reference_decision;
    audit.actual_reference_disagreements+=event.accepted!=reference_decision;
    audit.reference_positive_accepts+=reference_decision && candidate_reference>old_reference;
    const double delta=candidate_reference-old_reference;
    const double enclosure=n::BlockObjectiveReplayEnclosure(std::max(std::abs(old_reference),
        std::abs(candidate_reference)));
    const double consistency=std::abs(delta+event.local_objective_before-event.local_objective_after);
    audit.minimum_decision_margin=std::min({audit.minimum_decision_margin,std::abs(delta),
        std::abs(delta-enclosure),std::abs(event.local_objective_after-event.local_objective_before),
        std::abs(consistency-enclosure)});
}
std::string VectorProperty(n::VectorRef values)
{
    std::ostringstream out; out<<std::setprecision(17)<<'[';
    for(Eigen::Index k=0;k<values.size();++k) {if(k) out<<','; out<<values(k);}
    out<<']'; return out.str();
}
std::string ReplayDigest(const n::Vector & values)
{
    boost::hash2::sha2_256 hash;
    for(Eigen::Index i=0;i<values.size();++i)
    {
        const auto bits=std::bit_cast<std::uint64_t>(values(i));
        hash.update(&bits,sizeof(bits));
    }
    return boost::hash2::to_string(hash.result());
}

TEST(JointFixedNeighborBlockCoordinateTest, ProfileEndpointMatchesFixedStateEvaluation)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    const auto profile=n::EvaluateProfile(data.domain,data.y,eta,&data.context);
    ASSERT_TRUE(profile.valid)<<profile.reason;
    const auto state=n::EvaluateState(data.domain,data.y,profile.eta,profile.beta,data.context);
    ASSERT_TRUE(state.valid)<<state.reason;

    EXPECT_TRUE(BitwiseEqual(profile.eta,state.eta));
    EXPECT_TRUE(BitwiseEqual(profile.beta,state.beta));
    EXPECT_TRUE(BitwiseEqual(profile.gradient,state.gradient));
    EXPECT_TRUE(BitwiseEqual(profile.certificate.objective,state.certificate.objective));
    EXPECT_TRUE(BitwiseEqual(profile.certificate.projected_kkt,state.certificate.projected_kkt));
    EXPECT_EQ(profile.certificate.evaluated,state.certificate.evaluated);
    EXPECT_EQ(profile.certificate.available,state.certificate.available);
    EXPECT_EQ(profile.certificate.feasible,state.certificate.feasible);
    EXPECT_EQ(profile.certificate.kkt_passed,state.certificate.kkt_passed);
    EXPECT_EQ(profile.certificate.width_bound_applied,state.certificate.width_bound_applied);
    EXPECT_EQ(profile.certificate.width_bound_feasible,state.certificate.width_bound_feasible);
    EXPECT_TRUE(BitwiseEqual(profile.certificate.width_bound_projected_kkt,
        state.certificate.width_bound_projected_kkt));
}

TEST(JointFixedNeighborBlockCoordinateTest, ReplayWorkspacePreservesExactSearchResultsAcrossRuns)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=4;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.capture_local_trajectory=true; policy.eta_upper_bound=std::log(2.2);
    const auto first=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    const auto second=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    EXPECT_EQ(first.reason,second.reason);
    EXPECT_EQ(first.total_block_solves,second.total_block_solves);
    EXPECT_EQ(first.total_profile_evaluations,second.total_profile_evaluations);
    ASSERT_EQ(first.blocks.size(),second.blocks.size());
    ASSERT_EQ(first.sweeps.size(),second.sweeps.size());
    EXPECT_TRUE(BitwiseEqual(first.state.eta,second.state.eta));
    EXPECT_TRUE(BitwiseEqual(first.state.beta,second.state.beta));
    EXPECT_TRUE(BitwiseEqual(first.state.prediction,second.state.prediction));
    EXPECT_TRUE(BitwiseEqual(first.state.residual,second.state.residual));
    EXPECT_TRUE(BitwiseEqual(first.state.objective,second.state.objective));

    for(std::size_t i=0;i<first.blocks.size();++i)
    {
        const auto & left=first.blocks[i]; const auto & right=second.blocks[i];
        EXPECT_EQ(left.sweep,right.sweep); EXPECT_EQ(left.block,right.block);
        EXPECT_EQ(left.status,right.status); EXPECT_EQ(left.reason,right.reason);
        EXPECT_EQ(left.accepted,right.accepted); EXPECT_EQ(left.recovery_retry,right.recovery_retry);
        EXPECT_EQ(left.recovery_success,right.recovery_success);
        EXPECT_TRUE(BitwiseEqual(left.local_objective_before,right.local_objective_before));
        EXPECT_TRUE(BitwiseEqual(left.local_objective_after,right.local_objective_after));
        EXPECT_TRUE(BitwiseEqual(left.global_replay_delta,right.global_replay_delta));
        ASSERT_EQ(left.profile_trials.size(),right.profile_trials.size());
        for(std::size_t j=0;j<left.profile_trials.size();++j)
        {
            EXPECT_EQ(left.profile_trials[j].profile_evaluation,right.profile_trials[j].profile_evaluation);
            EXPECT_EQ(left.profile_trials[j].accepted,right.profile_trials[j].accepted);
            EXPECT_EQ(left.profile_trials[j].accepted_update,right.profile_trials[j].accepted_update);
            EXPECT_EQ(left.profile_trials[j].acceptance_reason,right.profile_trials[j].acceptance_reason);
            EXPECT_TRUE(BitwiseEqual(left.profile_trials[j].local_objective_after,
                right.profile_trials[j].local_objective_after));
        }
    }
    for(std::size_t i=0;i<first.sweeps.size();++i)
    {
        EXPECT_EQ(first.sweeps[i].accepted_blocks,second.sweeps[i].accepted_blocks);
        EXPECT_EQ(first.sweeps[i].unchanged_blocks,second.sweeps[i].unchanged_blocks);
        EXPECT_EQ(first.sweeps[i].block_solves,second.sweeps[i].block_solves);
        EXPECT_TRUE(BitwiseEqual(first.sweeps[i].objective_after,second.sweeps[i].objective_after));
        EXPECT_TRUE(BitwiseEqual(first.sweeps[i].cache_replay_error,second.sweeps[i].cache_replay_error));
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, CandidateObserverDoesNotChangeFullReplaySearch)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=4;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.eta_upper_bound=std::log(2.2);
    policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
        const n::Indices &,n::VectorRef candidate_eta,n::VectorRef candidate_beta) {
        return n::FullReplayForTesting(data.domain,data.layout,data.y,candidate_eta,candidate_beta,data.context.scale);
    };
    const auto baseline=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    std::size_t candidate_count{},full_replay_count{},accepted_count{};
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        ++candidate_count; full_replay_count+=event.candidate_replay_ran; accepted_count+=event.accepted;
        EXPECT_EQ(event.candidate_eta.size(),event.accepted_before.eta.size());
        EXPECT_EQ(event.candidate_beta.size(),event.accepted_before.beta.size());
        EXPECT_EQ(event.full_prediction.size(),event.accepted_before.prediction.size());
        EXPECT_EQ(event.full_residual.size(),event.accepted_before.residual.size());
        EXPECT_GT(event.core_atoms.size(),0u);
        EXPECT_GT(event.affected_rows.size(),0u);
        EXPECT_TRUE(std::isfinite(event.full_objective));
    };
    const auto observed=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    EXPECT_GT(candidate_count,0u);
    EXPECT_GT(full_replay_count,0u);
    EXPECT_GT(accepted_count,0u);
    EXPECT_EQ(baseline.reason,observed.reason);
    EXPECT_EQ(baseline.total_block_solves,observed.total_block_solves);
    EXPECT_EQ(baseline.sweep_count,observed.sweep_count);
    EXPECT_TRUE(BitwiseEqual(baseline.state.eta,observed.state.eta));
    EXPECT_TRUE(BitwiseEqual(baseline.state.beta,observed.state.beta));
    EXPECT_TRUE(BitwiseEqual(baseline.state.prediction,observed.state.prediction));
    EXPECT_TRUE(BitwiseEqual(baseline.state.residual,observed.state.residual));
    EXPECT_TRUE(BitwiseEqual(baseline.state.objective,observed.state.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, AffectedRowReplayPreservesFullContributorAndObjectiveOrder)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::Vector beta(6);
    const auto cancel_basis=n::EvaluateKernel(0.0,.55,2.5);
    beta<<1.0,-cancel_basis.gaussian/cancel_basis.charge*(1.0+1e-15),0.0,-.2,1.4,.2;
    n::BlockCoordinateState accepted; accepted.eta=eta; accepted.beta=beta;
    auto initial=n::FullReplayForTesting(data.domain,data.layout,data.y,eta,beta,data.context.scale);
    accepted.prediction=std::move(initial.prediction); accepted.residual=std::move(initial.residual);
    accepted.objective=initial.objective;
    n::StructuralTopology topology(data.domain,data.context.atom_ids,data.layout);
    auto partition=topology.BuildCores(2); topology.PopulateMetadata(partition);
    ASSERT_FALSE(partition.cores.empty());
    ASSERT_EQ(partition.cores.front().atoms.size(),2u);
    std::size_t delta_method_prediction_mismatches{};
    for(std::size_t mode=0;mode<4;++mode)
    {
        const auto & core=partition.cores[mode%partition.cores.size()];
        n::Vector candidate_eta=accepted.eta,candidate_beta=accepted.beta;
        for(const auto atom:core.atoms)
        {
            if(mode==0 || mode==3) candidate_beta(2*atom)+=.01;
            if(mode==1 || mode==3) candidate_beta(2*atom+1)-=.007;
            if(mode==2 || mode==3) candidate_eta(atom)+=.002;
        }
        const auto incidence=second_stage_test::PrepareReplayIncidence(data.domain,data.layout,
            core.atoms,core.neighbor_atoms,core.affected_rows);
        ASSERT_TRUE(second_stage_test::VerifyReplayIncidenceCoverage(data.domain,data.layout,
            core.affected_rows,incidence));
        const auto full=n::FullReplayForTesting(data.domain,data.layout,data.y,
            candidate_eta,candidate_beta,data.context.scale);
        const auto exact=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
            accepted,candidate_eta,candidate_beta,data.context.scale,core.atoms,core.affected_rows,
            second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
            second_stage_test::ReplayObjectiveStrategy::FullScan,&incidence);
        const auto b1=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
            accepted,candidate_eta,candidate_beta,data.context.scale,core.atoms,core.affected_rows,
            second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
            second_stage_test::ReplayObjectiveStrategy::FullScan);
        const auto b2=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
            accepted,candidate_eta,candidate_beta,data.context.scale,core.atoms,core.affected_rows,
            second_stage_test::ReplayPredictionStrategy::DeltaSubtractThenAdd,
            second_stage_test::ReplayObjectiveStrategy::FullScan);
        const auto incremental_exact=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
            accepted,candidate_eta,candidate_beta,data.context.scale,core.atoms,core.affected_rows,
            second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
            second_stage_test::ReplayObjectiveStrategy::Incremental,&incidence);
        for(const auto objective_strategy:std::array{
            second_stage_test::ReplayObjectiveStrategy::IncrementalScaled,
            second_stage_test::ReplayObjectiveStrategy::IncrementalKahan,
            second_stage_test::ReplayObjectiveStrategy::IncrementalNeumaier,
            second_stage_test::ReplayObjectiveStrategy::IncrementalKahanScaled,
            second_stage_test::ReplayObjectiveStrategy::IncrementalNeumaierScaled})
        {
            const auto compensated=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
                accepted,candidate_eta,candidate_beta,data.context.scale,core.atoms,core.affected_rows,
                second_stage_test::ReplayPredictionStrategy::AffectedRowExact,objective_strategy,&incidence);
            EXPECT_TRUE(std::isfinite(compensated.objective));
        }
        std::size_t local{};
        for(const auto row_value:data.layout.informative_rows)
        {
            const auto row=static_cast<Eigen::Index>(row_value);
            const auto found=std::lower_bound(core.affected_rows.begin(),core.affected_rows.end(),row);
            const bool affected=found!=core.affected_rows.end() && *found==row;
            const auto index=static_cast<Eigen::Index>(found-core.affected_rows.begin());
            const double expected_prediction=affected ? exact.prediction(index) : accepted.prediction(row);
            const double expected_residual=affected ? exact.residual(index) : accepted.residual(row);
            EXPECT_TRUE(BitwiseEqual(full.prediction(row),expected_prediction));
            EXPECT_TRUE(BitwiseEqual(full.residual(row),expected_residual));
            if(affected) ++local;
        }
        EXPECT_EQ(local,core.affected_rows.size());
        EXPECT_TRUE(BitwiseEqual(full.objective,exact.objective));
        EXPECT_TRUE(std::isfinite(b1.objective)); EXPECT_TRUE(std::isfinite(b2.objective));
        EXPECT_TRUE(std::isfinite(incremental_exact.objective));
        for(Eigen::Index k=0;k<exact.prediction.size();++k)
            delta_method_prediction_mismatches+=!BitwiseEqual(b1.prediction(k),b2.prediction(k));
    }
    RecordProperty("B1_B2_prediction_bitwise_mismatch_count",std::to_string(delta_method_prediction_mismatches));
}

TEST(JointFixedNeighborBlockCoordinateTest, TiledObjectiveCacheUpdatesAffectedTilesAndRollsBackByOmission)
{
    const std::vector<std::size_t> informative{0,1,2,3,4,5,6};
    n::Vector accepted(7); accepted<<1.0,-.5,1e-14,2.0,-1.25,.75,1e-8;
    n::Indices affected{1,2,6};
    n::Vector candidate_residuals(3); candidate_residuals<<-.5000001,-2e-14,3e-8;
    constexpr double scale=1.7;
    const auto full_objective=[&](const n::Vector & residuals) {
        double squared{};
        for(const auto row:informative)
        {
            const auto index=static_cast<Eigen::Index>(row);
            squared+=residuals(index)*residuals(index);
        }
        return .5*squared/(scale*scale);
    };
    second_stage_test::TiledObjectiveCache tiled(informative,accepted,2,false);
    second_stage_test::TiledObjectiveCache tree(informative,accepted,2,true);
    const auto tiled_rejected=tiled.Evaluate(accepted,affected,candidate_residuals,scale);
    const auto tree_rejected=tree.Evaluate(accepted,affected,candidate_residuals,scale);
    EXPECT_EQ(tiled_rejected.tiles_recomputed,3u);
    EXPECT_EQ(tree_rejected.tiles_recomputed,3u);
    EXPECT_GT(tree.StorageBytes(),tiled.StorageBytes());
    EXPECT_NEAR(tiled_rejected.objective,tree_rejected.objective,1e-14);
    EXPECT_NEAR(tiled_rejected.objective,[&] {
        n::Vector candidate=accepted;
        for(std::size_t k=0;k<affected.size();++k)
            candidate(affected[k])=candidate_residuals(static_cast<Eigen::Index>(k));
        return full_objective(candidate);
    }(),1e-14);

    const auto repeated_tiled=tiled.Evaluate(accepted,affected,candidate_residuals,scale);
    const auto repeated_tree=tree.Evaluate(accepted,affected,candidate_residuals,scale);
    EXPECT_TRUE(BitwiseEqual(tiled_rejected.objective,repeated_tiled.objective));
    EXPECT_TRUE(BitwiseEqual(tree_rejected.objective,repeated_tree.objective));
    tiled.Commit(tiled_rejected); tree.Commit(tree_rejected);
    for(std::size_t k=0;k<affected.size();++k)
        accepted(affected[k])=candidate_residuals(static_cast<Eigen::Index>(k));

    const n::Indices second_affected{0,4};
    n::Vector second_residuals(2); second_residuals<<.9,-1.1;
    const auto next_tiled=tiled.Evaluate(accepted,second_affected,second_residuals,scale);
    const auto next_tree=tree.Evaluate(accepted,second_affected,second_residuals,scale);
    n::Vector second_candidate=accepted;
    for(std::size_t k=0;k<second_affected.size();++k)
        second_candidate(second_affected[k])=second_residuals(static_cast<Eigen::Index>(k));
    EXPECT_NEAR(next_tiled.objective,full_objective(second_candidate),1e-14);
    EXPECT_NEAR(next_tree.objective,full_objective(second_candidate),1e-14);
}

TEST(JointFixedNeighborBlockCoordinateTest, LocalIncidencePreservesDuplicateSupportEntries)
{
    n::Domain domain(2,{{{0,0.0},{0,0.0},{1,.25}},{{0,0.0},{1,.25}}});
    rhbm_gem::JointParameterLayout layout;
    layout.full_atoms={0,1}; layout.informative_rows={0,1};
    n::Vector observations(2); observations<<0.1,-0.2;
    n::BlockCoordinateState accepted;
    accepted.eta.resize(2); accepted.eta<<std::log(.7),std::log(.9);
    accepted.beta.resize(4); accepted.beta<<1.0,-.2,1.1,-0.0;
    auto replay=n::FullReplayForTesting(domain,layout,observations,accepted.eta,accepted.beta,1.0);
    accepted.prediction=std::move(replay.prediction); accepted.residual=std::move(replay.residual);
    accepted.objective=replay.objective;
    const n::Indices core{0},neighbors{1},affected{0,1};
    const auto incidence=second_stage_test::PrepareReplayIncidence(domain,layout,core,neighbors,affected);
    ASSERT_TRUE(second_stage_test::VerifyReplayIncidenceCoverage(domain,layout,affected,incidence));
    ASSERT_EQ(incidence.rows[0].size(),3u);
    n::Vector candidate_eta=accepted.eta,candidate_beta=accepted.beta;
    candidate_eta(0)+=.01; candidate_beta(0)+=.05;
    const auto reference=n::FullReplayForTesting(domain,layout,observations,
        candidate_eta,candidate_beta,1.0);
    const auto exact=second_stage_test::EvaluateReplayCandidate(domain,layout,observations,accepted,
        candidate_eta,candidate_beta,1.0,core,affected,
        second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayObjectiveStrategy::FullScan,&incidence);
    for(Eigen::Index row=0;row<2;++row)
    {
        EXPECT_TRUE(BitwiseEqual(reference.prediction(row),exact.prediction(row)));
        EXPECT_TRUE(BitwiseEqual(reference.residual(row),exact.residual(row)));
    }
    EXPECT_TRUE(BitwiseEqual(reference.objective,exact.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, ExactReplayShadowMatchesProductionCandidateSearch)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",30));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(30,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=5;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.eta_upper_bound=std::log(2.2);
    policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
        const n::Indices &,n::VectorRef candidate_eta,n::VectorRef candidate_beta) {
        return n::FullReplayForTesting(data.domain,data.layout,data.y,candidate_eta,candidate_beta,data.context.scale);
    };
    const auto baseline=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    n::StructuralTopology topology(data.domain,data.context.atom_ids,data.layout);
    auto partition=topology.BuildCores(policy.core_atoms); topology.PopulateMetadata(partition);
    n::BlockCoordinateState shadow; shadow.eta=eta; shadow.beta=n::Vector::Zero(2*eta.size());
    auto initial=n::FullReplayForTesting(data.domain,data.layout,data.y,shadow.eta,shadow.beta,data.context.scale);
    shadow.prediction=std::move(initial.prediction); shadow.residual=std::move(initial.residual);
    shadow.objective=initial.objective;
    constexpr std::array prediction_strategies{
        second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayPredictionStrategy::DeltaSubtractThenAdd,
        second_stage_test::ReplayPredictionStrategy::DeltaAddDifference,
        second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
        second_stage_test::ReplayPredictionStrategy::AffectedRowExact};
    constexpr std::array objective_strategies{
        second_stage_test::ReplayObjectiveStrategy::FullScan,
        second_stage_test::ReplayObjectiveStrategy::FullScan,
        second_stage_test::ReplayObjectiveStrategy::Incremental,
        second_stage_test::ReplayObjectiveStrategy::FullScan,
        second_stage_test::ReplayObjectiveStrategy::Incremental};
    constexpr std::array<const char *,5> method_names{"B1","B2","C","D1","E"};
    struct ShadowMetrics
    {
        std::size_t prediction_mismatches{},residual_mismatches{},objective_mismatches{},acceptance_disagreements{};
        Eigen::Index first_prediction_mismatch=-1,first_residual_mismatch=-1,first_objective_mismatch=-1;
        std::uint64_t maximum_prediction_ulp{},maximum_residual_ulp{},maximum_objective_ulp{};
        double maximum_prediction_absolute{},maximum_prediction_relative{},maximum_residual_absolute{},
            maximum_residual_relative{},maximum_objective_absolute{},maximum_objective_relative{};
    };
    std::array<n::BlockCoordinateState,5> experiment_states;
    for(auto & state:experiment_states) state=shadow;
    std::array<ShadowMetrics,5> experiment_metrics{};
    std::size_t candidates{},full_replays{},accepted_updates{},acceptance_disagreements{},
        prediction_mismatches{},residual_mismatches{},objective_mismatches{},accepted_since_refresh{},
        maximum_accepted_since_refresh{};
    Eigen::Index first_prediction_mismatch=-1,first_residual_mismatch=-1;
    std::uint64_t maximum_prediction_ulp{},maximum_residual_ulp{},maximum_objective_ulp{};
    double maximum_prediction_absolute{},maximum_prediction_relative{},maximum_residual_absolute{},
        maximum_residual_relative{},maximum_objective_absolute{};
    double minimum_acceptance_margin=std::numeric_limits<double>::infinity();
    double minimum_acceptance_margin_ratio=std::numeric_limits<double>::infinity();
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        ++candidates; full_replays+=event.candidate_replay_ran;
        EXPECT_TRUE(BitwiseEqual(shadow.eta,event.accepted_before.eta));
        EXPECT_TRUE(BitwiseEqual(shadow.beta,event.accepted_before.beta));
        EXPECT_TRUE(BitwiseEqual(shadow.prediction,event.accepted_before.prediction));
        EXPECT_TRUE(BitwiseEqual(shadow.residual,event.accepted_before.residual));
        EXPECT_TRUE(BitwiseEqual(shadow.objective,event.accepted_before.objective));
        ASSERT_LT(event.block,partition.cores.size()+1);
        const auto & structural=partition.cores[event.block-1];
        EXPECT_EQ(event.core_atoms,structural.atoms);
        const auto incidence=second_stage_test::PrepareReplayIncidence(data.domain,data.layout,
            event.core_atoms,structural.neighbor_atoms,event.affected_rows);
        const auto exact=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
            shadow,event.candidate_eta,event.candidate_beta,data.context.scale,event.core_atoms,
            event.affected_rows,second_stage_test::ReplayPredictionStrategy::AffectedRowExact,
            second_stage_test::ReplayObjectiveStrategy::FullScan,&incidence);
        for(const auto row_value:data.layout.informative_rows)
        {
            const auto row=static_cast<Eigen::Index>(row_value);
            const auto found=std::lower_bound(event.affected_rows.begin(),event.affected_rows.end(),row);
            const bool affected=found!=event.affected_rows.end() && *found==row;
            const auto local=static_cast<Eigen::Index>(found-event.affected_rows.begin());
            const double prediction=affected ? exact.prediction(local) : shadow.prediction(row);
            const double residual=affected ? exact.residual(local) : shadow.residual(row);
            const double prediction_error=std::abs(prediction-event.full_prediction(row));
            const double residual_error=std::abs(residual-event.full_residual(row));
            maximum_prediction_absolute=std::max(maximum_prediction_absolute,prediction_error);
            maximum_residual_absolute=std::max(maximum_residual_absolute,residual_error);
            maximum_prediction_relative=std::max(maximum_prediction_relative,
                prediction_error/std::max(std::abs(event.full_prediction(row)),std::numeric_limits<double>::min()));
            maximum_residual_relative=std::max(maximum_residual_relative,
                residual_error/std::max(std::abs(event.full_residual(row)),std::numeric_limits<double>::min()));
            if(std::isnormal(prediction) && std::isnormal(event.full_prediction(row)))
                maximum_prediction_ulp=std::max(maximum_prediction_ulp,UlpDistance(prediction,event.full_prediction(row)));
            if(std::isnormal(residual) && std::isnormal(event.full_residual(row)))
                maximum_residual_ulp=std::max(maximum_residual_ulp,UlpDistance(residual,event.full_residual(row)));
            if(!BitwiseEqual(prediction,event.full_prediction(row)))
            {++prediction_mismatches; if(first_prediction_mismatch<0) first_prediction_mismatch=row;}
            if(!BitwiseEqual(residual,event.full_residual(row)))
            {++residual_mismatches; if(first_residual_mismatch<0) first_residual_mismatch=row;}
        }
        for(std::size_t method=0;method<experiment_states.size();++method)
        {
            auto & method_state=experiment_states[method];
            EXPECT_TRUE(BitwiseEqual(method_state.eta,event.accepted_before.eta));
            EXPECT_TRUE(BitwiseEqual(method_state.beta,event.accepted_before.beta));
            const auto result=second_stage_test::EvaluateReplayCandidate(data.domain,data.layout,data.y,
                method_state,event.candidate_eta,event.candidate_beta,data.context.scale,event.core_atoms,
                event.affected_rows,prediction_strategies[method],objective_strategies[method],
                method>=3 ? &incidence : nullptr);
            auto & metrics=experiment_metrics[method];
            for(const auto row_value:data.layout.informative_rows)
            {
                const auto row=static_cast<Eigen::Index>(row_value);
                const auto found=std::lower_bound(event.affected_rows.begin(),event.affected_rows.end(),row);
                const bool affected=found!=event.affected_rows.end() && *found==row;
                const auto local=static_cast<Eigen::Index>(found-event.affected_rows.begin());
                const double prediction=affected ? result.prediction(local) : method_state.prediction(row);
                const double residual=affected ? result.residual(local) : method_state.residual(row);
                const double prediction_error=std::abs(prediction-event.full_prediction(row));
                const double residual_error=std::abs(residual-event.full_residual(row));
                EXPECT_TRUE(std::isfinite(prediction)); EXPECT_TRUE(std::isfinite(residual));
                metrics.maximum_prediction_absolute=std::max(metrics.maximum_prediction_absolute,prediction_error);
                metrics.maximum_prediction_relative=std::max(metrics.maximum_prediction_relative,
                    prediction_error/std::max(std::abs(event.full_prediction(row)),std::numeric_limits<double>::min()));
                metrics.maximum_residual_absolute=std::max(metrics.maximum_residual_absolute,residual_error);
                metrics.maximum_residual_relative=std::max(metrics.maximum_residual_relative,
                    residual_error/std::max(std::abs(event.full_residual(row)),std::numeric_limits<double>::min()));
                if(std::isnormal(prediction) && std::isnormal(event.full_prediction(row)))
                    metrics.maximum_prediction_ulp=std::max(metrics.maximum_prediction_ulp,
                        UlpDistance(prediction,event.full_prediction(row)));
                if(std::isnormal(residual) && std::isnormal(event.full_residual(row)))
                    metrics.maximum_residual_ulp=std::max(metrics.maximum_residual_ulp,
                        UlpDistance(residual,event.full_residual(row)));
                if(!BitwiseEqual(prediction,event.full_prediction(row)))
                {
                    ++metrics.prediction_mismatches;
                    if(metrics.first_prediction_mismatch<0) metrics.first_prediction_mismatch=row;
                }
                if(!BitwiseEqual(residual,event.full_residual(row)))
                {
                    ++metrics.residual_mismatches;
                    if(metrics.first_residual_mismatch<0) metrics.first_residual_mismatch=row;
                }
            }
            const double objective_error=std::abs(result.objective-event.full_objective);
            EXPECT_TRUE(std::isfinite(result.objective));
            metrics.maximum_objective_absolute=std::max(metrics.maximum_objective_absolute,objective_error);
            metrics.maximum_objective_relative=std::max(metrics.maximum_objective_relative,
                objective_error/std::max(std::abs(event.full_objective),std::numeric_limits<double>::min()));
            if(std::isnormal(result.objective) && std::isnormal(event.full_objective))
                metrics.maximum_objective_ulp=std::max(metrics.maximum_objective_ulp,
                    UlpDistance(result.objective,event.full_objective));
            if(!BitwiseEqual(result.objective,event.full_objective))
            {
                ++metrics.objective_mismatches;
                if(metrics.first_objective_mismatch<0) metrics.first_objective_mismatch=static_cast<Eigen::Index>(candidates-1);
            }
            const double reference=std::max(std::abs(result.objective),std::abs(method_state.objective));
            const double enclosure=n::BlockObjectiveReplayEnclosure(reference);
            const double delta=result.objective-method_state.objective;
            const double consistency_error=std::abs(delta+event.local_objective_before-event.local_objective_after);
            const bool method_acceptance=!event.unchanged && n::WithinBlockObjectiveReplay(consistency_error,reference) &&
                (delta<=0 || (event.local_objective_after<=event.local_objective_before && delta<=enclosure));
            metrics.acceptance_disagreements+=method_acceptance!=event.accepted;
            if(event.accepted)
            {
                method_state.eta=event.candidate_eta; method_state.beta=event.candidate_beta;
                for(std::size_t k=0;k<event.affected_rows.size();++k)
                {
                    const auto row=event.affected_rows[k];
                    method_state.prediction(row)=result.prediction(static_cast<Eigen::Index>(k));
                    method_state.residual(row)=result.residual(static_cast<Eigen::Index>(k));
                }
                method_state.objective=result.objective;
            }
        }
        const double objective_error=std::abs(exact.objective-event.full_objective);
        maximum_objective_absolute=std::max(maximum_objective_absolute,objective_error);
        if(std::isnormal(exact.objective) && std::isnormal(event.full_objective))
            maximum_objective_ulp=std::max(maximum_objective_ulp,UlpDistance(exact.objective,event.full_objective));
        objective_mismatches+=!BitwiseEqual(exact.objective,event.full_objective);
        const double reference=std::max(std::abs(exact.objective),std::abs(shadow.objective));
        const double enclosure=n::BlockObjectiveReplayEnclosure(reference);
        const double delta=exact.objective-shadow.objective;
        const double consistency_error=std::abs(delta+event.local_objective_before-event.local_objective_after);
        const bool objective_rule=delta<=0 || (event.local_objective_after<=event.local_objective_before &&
            delta<=enclosure);
        const bool exact_acceptance=!event.unchanged && n::WithinBlockObjectiveReplay(consistency_error,reference) &&
            objective_rule;
        acceptance_disagreements+=exact_acceptance!=event.accepted;
        const double margin=std::min({enclosure-delta,
            event.local_objective_before-event.local_objective_after,enclosure-consistency_error});
        minimum_acceptance_margin=std::min(minimum_acceptance_margin,margin);
        minimum_acceptance_margin_ratio=std::min(minimum_acceptance_margin_ratio,margin/enclosure);
        if(event.accepted)
        {
            ++accepted_updates; ++accepted_since_refresh;
            shadow.eta=event.candidate_eta; shadow.beta=event.candidate_beta;
            for(std::size_t k=0;k<event.affected_rows.size();++k)
            {
                const auto row=event.affected_rows[k];
                shadow.prediction(row)=exact.prediction(static_cast<Eigen::Index>(k));
                shadow.residual(row)=exact.residual(static_cast<Eigen::Index>(k));
            }
            shadow.objective=exact.objective;
            maximum_accepted_since_refresh=std::max(maximum_accepted_since_refresh,accepted_since_refresh);
        }
    };
    policy.state_observer=[&](std::size_t,const n::BlockCoordinateState & state,const auto &,const auto &) {
        EXPECT_TRUE(BitwiseEqual(shadow.eta,state.eta));
        EXPECT_TRUE(BitwiseEqual(shadow.beta,state.beta));
        EXPECT_TRUE(BitwiseEqual(shadow.prediction,state.prediction));
        EXPECT_TRUE(BitwiseEqual(shadow.residual,state.residual));
        EXPECT_TRUE(BitwiseEqual(shadow.objective,state.objective));
        for(std::size_t method=0;method<experiment_states.size();++method)
        {
            EXPECT_TRUE(BitwiseEqual(experiment_states[method].eta,state.eta));
            EXPECT_TRUE(BitwiseEqual(experiment_states[method].beta,state.beta));
            if(method==3)
            {
                EXPECT_TRUE(BitwiseEqual(experiment_states[method].prediction,state.prediction));
                EXPECT_TRUE(BitwiseEqual(experiment_states[method].residual,state.residual));
                EXPECT_TRUE(BitwiseEqual(experiment_states[method].objective,state.objective));
            }
        }
    };
    const auto shadow_result=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    EXPECT_GT(candidates,0u); EXPECT_GT(full_replays,0u); EXPECT_GT(accepted_updates,0u);
    EXPECT_EQ(acceptance_disagreements,0u);
    EXPECT_EQ(prediction_mismatches,0u); EXPECT_EQ(residual_mismatches,0u);
    EXPECT_EQ(objective_mismatches,0u);
    EXPECT_EQ(baseline.reason,shadow_result.reason);
    ASSERT_EQ(baseline.sweeps.size(),shadow_result.sweeps.size());
    for(std::size_t k=0;k<baseline.sweeps.size();++k)
    {
        EXPECT_EQ(baseline.sweeps[k].accepted_blocks,shadow_result.sweeps[k].accepted_blocks);
        EXPECT_EQ(baseline.sweeps[k].unchanged_blocks,shadow_result.sweeps[k].unchanged_blocks);
        EXPECT_TRUE(BitwiseEqual(baseline.sweeps[k].objective_after,shadow_result.sweeps[k].objective_after));
        EXPECT_TRUE(BitwiseEqual(baseline.sweeps[k].cache_replay_error,shadow_result.sweeps[k].cache_replay_error));
    }
    EXPECT_TRUE(BitwiseEqual(baseline.state.eta,shadow_result.state.eta));
    EXPECT_TRUE(BitwiseEqual(baseline.state.beta,shadow_result.state.beta));
    EXPECT_TRUE(BitwiseEqual(baseline.state.prediction,shadow_result.state.prediction));
    EXPECT_TRUE(BitwiseEqual(baseline.state.residual,shadow_result.state.residual));
    EXPECT_TRUE(BitwiseEqual(baseline.state.objective,shadow_result.state.objective));
    RecordProperty("D1_shadow_candidates",std::to_string(candidates));
    RecordProperty("D1_shadow_full_replays",std::to_string(full_replays));
    RecordProperty("D1_shadow_accepted_updates",std::to_string(accepted_updates));
    RecordProperty("D1_shadow_mismatches",std::to_string(prediction_mismatches+residual_mismatches+objective_mismatches));
    RecordProperty("D1_shadow_max_prediction_absolute",Precise(maximum_prediction_absolute));
    RecordProperty("D1_shadow_max_prediction_relative",Precise(maximum_prediction_relative));
    RecordProperty("D1_shadow_max_residual_absolute",Precise(maximum_residual_absolute));
    RecordProperty("D1_shadow_max_residual_relative",Precise(maximum_residual_relative));
    RecordProperty("D1_shadow_max_prediction_ulp",std::to_string(maximum_prediction_ulp));
    RecordProperty("D1_shadow_max_residual_ulp",std::to_string(maximum_residual_ulp));
    RecordProperty("D1_shadow_max_objective_absolute",Precise(maximum_objective_absolute));
    RecordProperty("D1_shadow_max_objective_ulp",std::to_string(maximum_objective_ulp));
    RecordProperty("D1_shadow_min_acceptance_margin",Precise(minimum_acceptance_margin));
    RecordProperty("D1_shadow_min_acceptance_margin_ratio",Precise(minimum_acceptance_margin_ratio));
    RecordProperty("D1_shadow_max_accepted_updates_without_refresh",std::to_string(maximum_accepted_since_refresh));
    RecordProperty("D1_shadow_first_prediction_mismatch",std::to_string(first_prediction_mismatch));
    RecordProperty("D1_shadow_first_residual_mismatch",std::to_string(first_residual_mismatch));
    for(std::size_t method=0;method<experiment_metrics.size();++method)
    {
        const auto prefix=std::string(method_names[method])+"_shadow_";
        const auto & metrics=experiment_metrics[method];
        RecordProperty(prefix+"prediction_mismatches",std::to_string(metrics.prediction_mismatches));
        RecordProperty(prefix+"residual_mismatches",std::to_string(metrics.residual_mismatches));
        RecordProperty(prefix+"objective_mismatches",std::to_string(metrics.objective_mismatches));
        RecordProperty(prefix+"acceptance_disagreements",std::to_string(metrics.acceptance_disagreements));
        RecordProperty(prefix+"max_prediction_abs",Precise(metrics.maximum_prediction_absolute));
        RecordProperty(prefix+"max_prediction_rel",Precise(metrics.maximum_prediction_relative));
        RecordProperty(prefix+"max_residual_abs",Precise(metrics.maximum_residual_absolute));
        RecordProperty(prefix+"max_residual_rel",Precise(metrics.maximum_residual_relative));
        RecordProperty(prefix+"max_objective_abs",Precise(metrics.maximum_objective_absolute));
        RecordProperty(prefix+"max_objective_rel",Precise(metrics.maximum_objective_relative));
        RecordProperty(prefix+"max_prediction_ulp",std::to_string(metrics.maximum_prediction_ulp));
        RecordProperty(prefix+"max_residual_ulp",std::to_string(metrics.maximum_residual_ulp));
        RecordProperty(prefix+"max_objective_ulp",std::to_string(metrics.maximum_objective_ulp));
        RecordProperty(prefix+"first_prediction_mismatch",std::to_string(metrics.first_prediction_mismatch));
        RecordProperty(prefix+"first_residual_mismatch",std::to_string(metrics.first_residual_mismatch));
        RecordProperty(prefix+"first_objective_mismatch_candidate",std::to_string(metrics.first_objective_mismatch));
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, AffectedRowExactReplayMatchesClosedLoopFullReplay)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy full_policy; full_policy.core_atoms=2; full_policy.maximum_sweeps=40;
    full_policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    full_policy.eta_upper_bound=std::log(2.2); full_policy.capture_local_trajectory=true;
    full_policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
        const n::Indices &,n::VectorRef candidate_eta,n::VectorRef candidate_beta) {
        return n::FullReplayForTesting(data.domain,data.layout,data.y,candidate_eta,candidate_beta,data.context.scale);
    };
    struct CandidateTrace
    {
        std::size_t sweep{},block{};
        n::Indices core_atoms;
        n::Vector eta_before,beta_before,candidate_eta,candidate_beta;
        std::string before_prediction_digest,before_residual_digest,prediction_digest,residual_digest;
        double objective_before{},objective_after{},local_objective_before{},local_objective_after{},
            global_delta{},consistency_error{},replay_enclosure{};
        bool accepted{},unchanged{},replay_ran{};
    };
    std::vector<CandidateTrace> full_trace,d1_trace;
    std::vector<n::BlockCoordinateState> full_sweep_states,d1_sweep_states;
    const auto trace_candidate=[&](std::vector<CandidateTrace> & trace) {
        auto * const trace_ptr=&trace;
        return [trace_ptr](const n::FixedNeighborCandidateEvent & event) {
            trace_ptr->push_back({event.sweep,event.block,event.core_atoms,event.accepted_before.eta,
                event.accepted_before.beta,event.candidate_eta,event.candidate_beta,
                ReplayDigest(event.accepted_before.prediction),ReplayDigest(event.accepted_before.residual),
                ReplayDigest(event.full_prediction),ReplayDigest(event.full_residual),
                event.accepted_before.objective,event.full_objective,event.local_objective_before,
                event.local_objective_after,event.global_replay_delta,event.local_global_delta_error,
                event.objective_replay_enclosure,event.accepted,event.unchanged,event.candidate_replay_ran});
        };
    };
    full_policy.candidate_observer=trace_candidate(full_trace);
    full_policy.state_observer=[&](std::size_t,const n::BlockCoordinateState & state,const auto &,const auto &) {
        full_sweep_states.push_back(state);
    };
    const auto full=n::RunFixedNeighborDiagnostic(problem,eta,full_policy);

    auto d1_policy=full_policy;
    d1_policy.candidate_replay_for_testing={};
    d1_policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::FullScan;
    d1_policy.candidate_observer=trace_candidate(d1_trace);
    d1_policy.state_observer=[&](std::size_t,const n::BlockCoordinateState & state,const auto &,const auto &) {
        d1_sweep_states.push_back(state);
    };
    const auto d1=n::RunFixedNeighborDiagnostic(problem,eta,d1_policy);
    const auto d1_candidate_replays=static_cast<std::size_t>(std::count_if(d1_trace.begin(),d1_trace.end(),
        [](const auto & candidate){return candidate.replay_ran;}));

    ASSERT_EQ(full.reason,d1.reason);
    EXPECT_EQ(full.search_converged,d1.search_converged);
    EXPECT_EQ(full.endpoint_certified,d1.endpoint_certified);
    EXPECT_EQ(full.sweep_count,d1.sweep_count);
    ASSERT_EQ(full.blocks.size(),d1.blocks.size());
    ASSERT_EQ(full.sweeps.size(),d1.sweeps.size());
    ASSERT_EQ(full_trace.size(),d1_trace.size());
    ASSERT_EQ(full_sweep_states.size(),d1_sweep_states.size());
    for(std::size_t k=0;k<full_trace.size();++k)
    {
        const auto & a=full_trace[k]; const auto & b=d1_trace[k];
        EXPECT_EQ(a.sweep,b.sweep); EXPECT_EQ(a.block,b.block); EXPECT_EQ(a.core_atoms,b.core_atoms);
        EXPECT_TRUE(BitwiseEqual(a.eta_before,b.eta_before)); EXPECT_TRUE(BitwiseEqual(a.beta_before,b.beta_before));
        EXPECT_TRUE(BitwiseEqual(a.candidate_eta,b.candidate_eta)); EXPECT_TRUE(BitwiseEqual(a.candidate_beta,b.candidate_beta));
        EXPECT_EQ(a.before_prediction_digest,b.before_prediction_digest);
        EXPECT_EQ(a.before_residual_digest,b.before_residual_digest);
        EXPECT_EQ(a.prediction_digest,b.prediction_digest); EXPECT_EQ(a.residual_digest,b.residual_digest);
        EXPECT_TRUE(BitwiseEqual(a.objective_before,b.objective_before));
        EXPECT_TRUE(BitwiseEqual(a.objective_after,b.objective_after));
        EXPECT_TRUE(BitwiseEqual(a.local_objective_before,b.local_objective_before));
        EXPECT_TRUE(BitwiseEqual(a.local_objective_after,b.local_objective_after));
        EXPECT_TRUE(BitwiseEqual(a.global_delta,b.global_delta));
        EXPECT_TRUE(BitwiseEqual(a.consistency_error,b.consistency_error));
        EXPECT_TRUE(BitwiseEqual(a.replay_enclosure,b.replay_enclosure));
        EXPECT_EQ(a.accepted,b.accepted); EXPECT_EQ(a.unchanged,b.unchanged); EXPECT_EQ(a.replay_ran,b.replay_ran);
    }
    for(std::size_t k=0;k<full.blocks.size();++k)
    {
        const auto & a=full.blocks[k]; const auto & b=d1.blocks[k];
        EXPECT_EQ(a.sweep,b.sweep); EXPECT_EQ(a.block,b.block); EXPECT_EQ(a.core_atoms,b.core_atoms);
        EXPECT_EQ(a.status,b.status); EXPECT_EQ(a.reason,b.reason); EXPECT_EQ(a.accepted,b.accepted);
        EXPECT_EQ(a.recovery_retry,b.recovery_retry); EXPECT_EQ(a.recovery_success,b.recovery_success);
        EXPECT_EQ(a.profile_evaluations,b.profile_evaluations); EXPECT_EQ(a.accepted_updates,b.accepted_updates);
        EXPECT_TRUE(BitwiseEqual(a.local_objective_before,b.local_objective_before));
        EXPECT_TRUE(BitwiseEqual(a.local_objective_after,b.local_objective_after));
        EXPECT_TRUE(BitwiseEqual(a.global_replay_delta,b.global_replay_delta));
        EXPECT_TRUE(BitwiseEqual(a.local_global_delta_error,b.local_global_delta_error));
        EXPECT_TRUE(BitwiseEqual(a.objective_after,b.objective_after));
        EXPECT_TRUE(BitwiseEqual(a.local_final_gradient_inf_norm,b.local_final_gradient_inf_norm));
        EXPECT_TRUE(BitwiseEqual(a.local_final_ac_kkt,b.local_final_ac_kkt));
        ASSERT_EQ(a.profile_trials.size(),b.profile_trials.size());
        for(std::size_t trial=0;trial<a.profile_trials.size();++trial)
        {
            const auto & x=a.profile_trials[trial]; const auto & y=b.profile_trials[trial];
            EXPECT_EQ(x.trial_index,y.trial_index); EXPECT_EQ(x.profile_evaluation,y.profile_evaluation);
            EXPECT_EQ(x.accepted,y.accepted); EXPECT_EQ(x.accepted_update,y.accepted_update);
            EXPECT_TRUE(BitwiseEqual(x.local_objective_before,y.local_objective_before));
            EXPECT_TRUE(BitwiseEqual(x.local_objective_after,y.local_objective_after));
            EXPECT_TRUE(BitwiseEqual(x.eta_change_inf,y.eta_change_inf));
            EXPECT_TRUE(BitwiseEqual(x.gradient_inf_norm,y.gradient_inf_norm));
            EXPECT_TRUE(BitwiseEqual(x.lm_actual_decrease,y.lm_actual_decrease));
            EXPECT_TRUE(BitwiseEqual(x.lm_predicted_decrease,y.lm_predicted_decrease));
            EXPECT_TRUE(BitwiseEqual(x.lm_ratio,y.lm_ratio));
        }
    }
    for(std::size_t k=0;k<full.sweeps.size();++k)
    {
        const auto & a=full.sweeps[k]; const auto & b=d1.sweeps[k];
        EXPECT_EQ(a.sweep,b.sweep); EXPECT_EQ(a.accepted_blocks,b.accepted_blocks);
        EXPECT_EQ(a.unchanged_blocks,b.unchanged_blocks); EXPECT_EQ(a.block_solves,b.block_solves);
        EXPECT_TRUE(BitwiseEqual(a.objective_before,b.objective_before));
        EXPECT_TRUE(BitwiseEqual(a.objective_after,b.objective_after));
        EXPECT_TRUE(BitwiseEqual(a.global_ac_kkt,b.global_ac_kkt));
        EXPECT_TRUE(BitwiseEqual(a.projected_width_gradient_inf_norm,b.projected_width_gradient_inf_norm));
        EXPECT_TRUE(BitwiseEqual(a.cache_replay_error,b.cache_replay_error));
        EXPECT_TRUE(BitwiseEqual(a.objective_replay_error,b.objective_replay_error));
        EXPECT_TRUE(BitwiseEqual(full_sweep_states[k].eta,d1_sweep_states[k].eta));
        EXPECT_TRUE(BitwiseEqual(full_sweep_states[k].beta,d1_sweep_states[k].beta));
        EXPECT_TRUE(BitwiseEqual(full_sweep_states[k].prediction,d1_sweep_states[k].prediction));
        EXPECT_TRUE(BitwiseEqual(full_sweep_states[k].residual,d1_sweep_states[k].residual));
        EXPECT_TRUE(BitwiseEqual(full_sweep_states[k].objective,d1_sweep_states[k].objective));
    }
    EXPECT_TRUE(BitwiseEqual(full.state.eta,d1.state.eta));
    EXPECT_TRUE(BitwiseEqual(full.state.beta,d1.state.beta));
    EXPECT_TRUE(BitwiseEqual(full.state.prediction,d1.state.prediction));
    EXPECT_TRUE(BitwiseEqual(full.state.residual,d1.state.residual));
    EXPECT_TRUE(BitwiseEqual(full.state.objective,d1.state.objective));
    auto context=data.context; context.eta_upper_bound=full_policy.eta_upper_bound;
    const auto full_primary=n::EvaluateStateStreaming(data.domain,data.y,full.state.eta,full.state.beta,context,
        n::ProfileEvaluationRole::AcceptedEndpoint);
    const auto d1_primary=n::EvaluateStateStreaming(data.domain,data.y,d1.state.eta,d1.state.beta,context,
        n::ProfileEvaluationRole::AcceptedEndpoint);
    ASSERT_TRUE(full_primary.valid)<<full_primary.reason; ASSERT_TRUE(d1_primary.valid)<<d1_primary.reason;
    const auto full_independent=n::CheckIndependentReplay(data.domain,data.y,full_primary,context,
        nullptr,nullptr,full_policy.convergence_criteria);
    const auto d1_independent=n::CheckIndependentReplay(data.domain,data.y,d1_primary,context,
        nullptr,nullptr,d1_policy.convergence_criteria);
    ASSERT_TRUE(full_independent.available && full_independent.valid)<<full_independent.reason;
    ASSERT_TRUE(d1_independent.available && d1_independent.valid)<<d1_independent.reason;
    EXPECT_EQ(full_independent.passed,d1_independent.passed);
    EXPECT_TRUE(BitwiseEqual(full_independent.primary_objective,d1_independent.primary_objective));
    EXPECT_TRUE(BitwiseEqual(full_independent.independent_objective,d1_independent.independent_objective));
    EXPECT_TRUE(BitwiseEqual(full_independent.primary_ac_kkt,d1_independent.primary_ac_kkt));
    EXPECT_TRUE(BitwiseEqual(full_independent.independent_ac_kkt,d1_independent.independent_ac_kkt));
    EXPECT_EQ(full.fit.RuntimeConvergence(),d1.fit.RuntimeConvergence());
    ASSERT_TRUE(full.endpoint_certified)<<full.endpoint_trust.reason;
    ASSERT_TRUE(d1.endpoint_certified)<<d1.endpoint_trust.reason;
    ASSERT_TRUE(full.endpoint_trust.independent.has_value());
    ASSERT_TRUE(d1.endpoint_trust.independent.has_value());
    EXPECT_EQ(full_independent.passed,full.endpoint_trust.independent->passed);
    EXPECT_EQ(d1_independent.passed,d1.endpoint_trust.independent->passed);
    const auto full_accepted=static_cast<std::size_t>(std::count_if(full_trace.begin(),full_trace.end(),
        [](const auto & candidate){return candidate.accepted;}));
    const auto profile_trials=std::accumulate(full.blocks.begin(),full.blocks.end(),std::size_t{},
        [](std::size_t count,const auto & block){return count+block.profile_trials.size();});
    EXPECT_GT(d1_candidate_replays,0u);
    EXPECT_EQ(d1_candidate_replays,static_cast<std::size_t>(std::count_if(d1_trace.begin(),d1_trace.end(),
        [](const auto & candidate){return candidate.replay_ran;})));
    RecordProperty("N4_candidates",std::to_string(full_trace.size()));
    RecordProperty("N4_accepted_candidates",std::to_string(full_accepted));
    RecordProperty("N4_replay_calls",std::to_string(d1_candidate_replays));
    RecordProperty("N4_sweeps",std::to_string(full.sweeps.size()));
    RecordProperty("N4_blocks",std::to_string(full.blocks.size()));
    RecordProperty("N4_profile_trials",std::to_string(profile_trials));
    RecordProperty("N4_converged",std::to_string(full.search_converged));
    RecordProperty("N4_stop_reason",full.reason);
    RecordProperty("N4_endpoint_certified",std::to_string(full.endpoint_certified));
    RecordProperty("N4_independent_replay_passed",std::to_string(full_independent.passed));
    RecordProperty("N4_global_ac_kkt",Precise(full.sweeps.back().global_ac_kkt));
    RecordProperty("N4_projected_width_kkt",Precise(full.sweeps.back().projected_width_gradient_inf_norm));
    RecordProperty("N4_final_objective",Precise(full.state.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveControlsClosedLoopAndGuardedFallback)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy base; base.core_atoms=2; base.maximum_sweeps=40;
    base.convergence_criteria=n::PracticalJointConvergenceCriteria;
    base.eta_upper_bound=std::log(2.2); base.collect_telemetry=true;

    auto full_replay_policy=base;
    full_replay_policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
        const n::Indices &,n::VectorRef candidate_eta,n::VectorRef candidate_beta) {
        return n::FullReplayForTesting(data.domain,data.layout,data.y,candidate_eta,candidate_beta,data.context.scale);
    };
    const auto full_replay=n::RunFixedNeighborDiagnostic(problem,initial_eta,full_replay_policy);
    auto d1_policy=base;
    d1_policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::FullScan;
    const auto d1=n::RunFixedNeighborDiagnostic(problem,initial_eta,d1_policy);
    ASSERT_TRUE(full_replay.search_converged)<<full_replay.reason;
    ASSERT_TRUE(d1.search_converged)<<d1.reason;
    ASSERT_TRUE(full_replay.endpoint_certified)<<full_replay.endpoint_trust.reason;
    ASSERT_TRUE(d1.endpoint_certified)<<d1.endpoint_trust.reason;

    const auto run_incremental=[&](n::FixedNeighborObjectiveStrategyForTesting strategy,double guard_fraction,
        bool refresh,const std::string & prefix) {
        auto policy=base; policy.objective_strategy_for_testing=strategy;
        policy.objective_guard_band_fraction_for_testing=guard_fraction;
        policy.refresh_objective_after_sweep_for_testing=refresh;
        IncrementalObjectiveAudit audit;
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            AuditIncrementalCandidate(data.layout,data.context.scale,event,audit);
        };
        const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        RecordProperty(prefix+"_candidates",std::to_string(audit.candidates));
        RecordProperty(prefix+"_accepted",std::to_string(audit.accepted));
        RecordProperty(prefix+"_rejected",std::to_string(audit.rejected));
        RecordProperty(prefix+"_fallbacks",std::to_string(audit.fallbacks));
        RecordProperty(prefix+"_estimated_reference_disagreements",std::to_string(audit.estimated_reference_disagreements));
        RecordProperty(prefix+"_actual_reference_disagreements",std::to_string(audit.actual_reference_disagreements));
        RecordProperty(prefix+"_max_objective_abs_error",Precise(audit.maximum_absolute));
        RecordProperty(prefix+"_max_objective_relative_error",Precise(audit.maximum_relative));
        RecordProperty(prefix+"_max_accepted_cache_drift",Precise(audit.maximum_cache_drift));
        RecordProperty(prefix+"_max_objective_ulp",std::to_string(audit.maximum_ulp));
        RecordProperty(prefix+"_minimum_reference_decision_margin",Precise(audit.minimum_decision_margin));
        RecordProperty(prefix+"_search_converged",std::to_string(result.search_converged));
        RecordProperty(prefix+"_stop_reason",result.reason);
        RecordProperty(prefix+"_sweeps",std::to_string(result.sweep_count));
        RecordProperty(prefix+"_objective_full_scan_fallbacks",std::to_string(result.work.objective_full_scan_fallbacks));
        RecordProperty(prefix+"_objective_cache_refreshes",std::to_string(std::count_if(result.sweeps.begin(),
            result.sweeps.end(),[](const auto & sweep){return sweep.objective_cache_refreshed;})));
        if(result.final_sweep)
        {
            RecordProperty(prefix+"_global_ac_kkt",Precise(result.final_sweep->global_ac_kkt));
            RecordProperty(prefix+"_projected_width_kkt",Precise(result.final_sweep->projected_width_gradient_inf_norm));
            RecordProperty(prefix+"_objective_replay_error",Precise(result.final_sweep->objective_replay_error));
        }
        if(result.endpoint_certified)
        {
            EXPECT_TRUE(result.assessment.primary.valid);
            EXPECT_TRUE(result.assessment.independent_replay.passed);
            EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        }
        EXPECT_EQ(audit.candidates,result.work.incremental_objective_evaluations);
        EXPECT_EQ(audit.fallbacks,result.work.objective_full_scan_fallbacks);
        EXPECT_GT(audit.candidates,0u);
        RecordProperty(prefix+"_final_objective",Precise(result.state.objective));
        RecordProperty(prefix+"_final_eta",VectorProperty(result.state.eta));
        RecordProperty(prefix+"_final_width",VectorProperty(result.state.eta.array().exp()));
        RecordProperty(prefix+"_final_beta",VectorProperty(result.state.beta));
        double amplitude_squared{},charge_squared{},width_squared{};
        for(Eigen::Index atom=0;atom<result.state.eta.size();++atom)
        {
            amplitude_squared+=std::pow(result.state.beta(2*atom)-2.0,2);
            charge_squared+=std::pow(result.state.beta(2*atom+1)-0.2,2);
            width_squared+=std::pow(std::exp(result.state.eta(atom))-0.5,2);
        }
        const double atom_count=static_cast<double>(result.state.eta.size());
        RecordProperty(prefix+"_truth_rmse_A",Precise(std::sqrt(amplitude_squared/atom_count)));
        RecordProperty(prefix+"_truth_rmse_C",Precise(std::sqrt(charge_squared/atom_count)));
        RecordProperty(prefix+"_truth_rmse_width",Precise(std::sqrt(width_squared/atom_count)));
        if(result.fit.objective) RecordProperty(prefix+"_independent_final_objective",Precise(*result.fit.objective));
        return result;
    };

    const auto direct=run_incremental(n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect,
        0.0,false,"E1_direct");
    const auto narrow_guard=run_incremental(n::FixedNeighborObjectiveStrategyForTesting::IncrementalGuarded,
        0.01,true,"E1_guarded_001");
    const auto guarded=run_incremental(n::FixedNeighborObjectiveStrategyForTesting::IncrementalGuarded,
        0.1,true,"E1_guarded");
    const auto wide_guard=run_incremental(n::FixedNeighborObjectiveStrategyForTesting::IncrementalGuarded,
        1.0,true,"E1_guarded_1");
    EXPECT_TRUE(direct.search_converged)<<direct.reason;
    EXPECT_TRUE(guarded.search_converged)<<guarded.reason;
    EXPECT_TRUE(direct.endpoint_certified)<<direct.endpoint_trust.reason;
    EXPECT_TRUE(guarded.endpoint_certified)<<guarded.endpoint_trust.reason;
    EXPECT_TRUE(direct.assessment.independent_replay.passed);
    EXPECT_TRUE(guarded.assessment.independent_replay.passed);
    EXPECT_TRUE(narrow_guard.search_converged)<<narrow_guard.reason;
    EXPECT_TRUE(wide_guard.search_converged)<<wide_guard.reason;
    EXPECT_TRUE(narrow_guard.endpoint_certified)<<narrow_guard.endpoint_trust.reason;
    EXPECT_TRUE(wide_guard.endpoint_certified)<<wide_guard.endpoint_trust.reason;
    const auto parameter_differences=[&](const n::BlockCoordinateState & candidate,const std::string & prefix) {
        double max_amplitude{},max_charge{},max_width{};
        for(Eigen::Index atom=0;atom<d1.state.eta.size();++atom)
        {
            max_amplitude=std::max(max_amplitude,std::abs(candidate.beta(2*atom)-d1.state.beta(2*atom)));
            max_charge=std::max(max_charge,std::abs(candidate.beta(2*atom+1)-d1.state.beta(2*atom+1)));
            max_width=std::max(max_width,std::abs(std::exp(candidate.eta(atom))-std::exp(d1.state.eta(atom))));
        }
        RecordProperty(prefix+"_max_abs_A_difference",Precise(max_amplitude));
        RecordProperty(prefix+"_max_abs_C_difference",Precise(max_charge));
        RecordProperty(prefix+"_max_abs_width_difference",Precise(max_width));
    };
    parameter_differences(direct.state,"E1_direct");
    parameter_differences(guarded.state,"E1_guarded");
    parameter_differences(narrow_guard.state,"E1_guarded_001");
    parameter_differences(wide_guard.state,"E1_guarded_1");
    EXPECT_LE(std::abs(direct.state.objective-direct.fit.objective.value()),
        n::BlockObjectiveReplayEnclosure(direct.fit.objective.value()));
    EXPECT_LE(std::abs(guarded.state.objective-guarded.fit.objective.value()),
        n::BlockObjectiveReplayEnclosure(guarded.fit.objective.value()));
    EXPECT_LE(std::abs(narrow_guard.state.objective-narrow_guard.fit.objective.value()),
        n::BlockObjectiveReplayEnclosure(narrow_guard.fit.objective.value()));
    EXPECT_LE(std::abs(wide_guard.state.objective-wide_guard.fit.objective.value()),
        n::BlockObjectiveReplayEnclosure(wide_guard.fit.objective.value()));
    RecordProperty("E0_D1_sweeps",std::to_string(d1.sweep_count));
    RecordProperty("E0_D1_candidates",std::to_string(d1.blocks.size()));
    RecordProperty("E0_D1_final_objective",Precise(d1.state.objective));
    RecordProperty("E0_D1_global_ac_kkt",Precise(d1.final_sweep->global_ac_kkt));
    RecordProperty("E0_D1_projected_width_kkt",Precise(d1.final_sweep->projected_width_gradient_inf_norm));
    EXPECT_TRUE(BitwiseEqual(full_replay.state.eta,d1.state.eta));
    EXPECT_TRUE(BitwiseEqual(full_replay.state.beta,d1.state.beta));
    EXPECT_TRUE(BitwiseEqual(full_replay.state.prediction,d1.state.prediction));
    EXPECT_TRUE(BitwiseEqual(full_replay.state.residual,d1.state.residual));
    EXPECT_TRUE(BitwiseEqual(full_replay.state.objective,d1.state.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveQualifiesAtActiveWidthUpperBound)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(3,std::log(.5));
    const auto run=[&](n::FixedNeighborObjectiveStrategyForTesting strategy,bool full_replay_reference) {
        n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=40;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(.5); policy.collect_telemetry=true;
        policy.objective_strategy_for_testing=strategy;
        policy.refresh_objective_after_sweep_for_testing=strategy==
            n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
        if(full_replay_reference)
            policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
                const n::Indices &,n::VectorRef eta,n::VectorRef beta) {
                return n::FullReplayForTesting(data.domain,data.layout,data.y,eta,beta,data.context.scale);
            };
        return n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    };
    const auto f1=run(n::FixedNeighborObjectiveStrategyForTesting::FullScan,true);
    const auto d1=run(n::FixedNeighborObjectiveStrategyForTesting::FullScan,false);
    const auto e1=run(n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect,false);
    for(const auto & [name,result]:std::array<std::pair<const char *,const n::FixedNeighborDiagnosticResult *>,3>{{
        {"F1",&f1},{"D1",&d1},{"E1",&e1}}})
    {
        const auto prefix=std::string("active_width_")+name+"_";
        RecordProperty((prefix+"objective").c_str(),Precise(result->fit.objective.value()));
        RecordProperty((prefix+"global_ac_kkt").c_str(),Precise(result->final_sweep->global_ac_kkt));
        RecordProperty((prefix+"projected_width_kkt").c_str(),
            Precise(result->final_sweep->projected_width_gradient_inf_norm));
        RecordProperty((prefix+"active_bound_count").c_str(),
            std::to_string(result->final_sweep->active_width_bound_count));
    }
    for(const auto * result:{&f1,&d1,&e1})
    {
        ASSERT_TRUE(result->search_converged)<<result->reason;
        ASSERT_TRUE(result->endpoint_certified)<<result->endpoint_trust.reason;
        ASSERT_TRUE(result->assessment.independent_replay.passed);
        ASSERT_TRUE(result->final_sweep);
        EXPECT_TRUE(result->final_sweep->width_bound_applied);
        EXPECT_GT(result->final_sweep->active_width_bound_count,0u);
        EXPECT_LE(result->final_sweep->maximum_width_bound_violation,0.0);
        EXPECT_LE(result->final_sweep->projected_width_gradient_inf_norm,
            n::PracticalJointConvergenceCriteria.global_width_gradient);
    }
    EXPECT_TRUE(BitwiseEqual(f1.state.eta,d1.state.eta));
    EXPECT_TRUE(BitwiseEqual(f1.state.beta,d1.state.beta));
    EXPECT_TRUE(BitwiseEqual(d1.state.eta,e1.state.eta));
    EXPECT_TRUE(BitwiseEqual(d1.state.beta,e1.state.beta));
    EXPECT_EQ(e1.work.objective_full_scan_fallbacks,0u);
    EXPECT_GT(e1.work.incremental_objective_evaluations,0u);
    EXPECT_EQ(std::count_if(e1.sweeps.begin(),e1.sweeps.end(),
        [](const auto & sweep){return sweep.objective_cache_refreshed;}),
        static_cast<std::ptrdiff_t>(e1.sweeps.size()));
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveRefreshIntervalsReconcileCache)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",300));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(300,std::log(.55));
    n::FixedNeighborPolicy base; base.core_atoms=12; base.maximum_sweeps=30;
    base.convergence_criteria=n::PracticalJointConvergenceCriteria;
    base.stop_after_stationarity=true; base.certify_endpoint=true;
    base.assess_local_conditioning=true; base.assess_parameter_stability=true;
    base.eta_upper_bound=std::log(2.2); base.collect_telemetry=true;
    const auto run_interval=[&](std::size_t interval,const std::string & prefix) {
        auto policy=base;
        policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
        policy.objective_refresh_interval_for_testing=interval;
        IncrementalObjectiveAudit audit;
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            AuditIncrementalCandidate(data.layout,data.context.scale,event,audit);
        };
        const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        RecordProperty(prefix+"_candidates",std::to_string(audit.candidates));
        RecordProperty(prefix+"_accepted",std::to_string(audit.accepted));
        RecordProperty(prefix+"_rejected",std::to_string(audit.rejected));
        RecordProperty(prefix+"_decision_disagreements",std::to_string(audit.actual_reference_disagreements));
        RecordProperty(prefix+"_maximum_absolute_objective_error",Precise(audit.maximum_absolute));
        RecordProperty(prefix+"_maximum_relative_objective_error",Precise(audit.maximum_relative));
        RecordProperty(prefix+"_maximum_cache_drift",Precise(audit.maximum_cache_drift));
        RecordProperty(prefix+"_maximum_objective_ulp",std::to_string(audit.maximum_ulp));
        RecordProperty(prefix+"_refreshes",std::to_string(result.work.objective_cache_refreshes));
        RecordProperty(prefix+"_refresh_max_error",Precise(result.work.maximum_objective_refresh_error));
        RecordProperty(prefix+"_refresh_max_accepted_updates",std::to_string(
            result.work.maximum_accepted_updates_between_objective_refreshes));
        RecordProperty(prefix+"_refresh_positive_deltas",std::to_string(result.work.objective_refresh_positive_deltas));
        RecordProperty(prefix+"_refresh_failures",std::to_string(result.work.objective_refresh_failures));
        RecordProperty(prefix+"_refresh_seconds",Precise(result.work.objective_cache_refresh_seconds));
        RecordProperty(prefix+"_sweeps",std::to_string(result.sweep_count));
        RecordProperty(prefix+"_stop_reason",result.reason);
        ASSERT_TRUE(result.search_converged)<<result.reason;
        ASSERT_TRUE(result.endpoint_certified)<<result.endpoint_trust.reason;
        EXPECT_TRUE(result.assessment.independent_replay.passed);
        EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        EXPECT_EQ(audit.candidates,result.work.incremental_objective_evaluations);
        EXPECT_EQ(audit.actual_reference_disagreements,0u);
        EXPECT_EQ(result.work.objective_refresh_failures,0u);
        EXPECT_GE(result.work.objective_cache_refreshes,result.sweep_count);
        EXPECT_EQ(result.work.objective_reconciliations,result.sweep_count);
        EXPECT_LE(result.work.maximum_accepted_updates_between_objective_refreshes,interval);
        EXPECT_TRUE(result.final_sweep.has_value());
        EXPECT_LE(std::abs(result.state.objective-result.fit.objective.value()),
            n::BlockObjectiveReplayEnclosure(result.fit.objective.value()));
    };
    run_interval(10,"E2_refresh_10");
    run_interval(50,"E2_refresh_50");
    run_interval(100,"E2_refresh_100");

    auto sweep_policy=base;
    sweep_policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
    sweep_policy.refresh_objective_after_sweep_for_testing=true;
    const auto every_sweep=n::RunFixedNeighborDiagnostic(problem,initial_eta,sweep_policy);
    RecordProperty("E2_refresh_every_sweep_count",std::to_string(every_sweep.work.objective_cache_refreshes));
    RecordProperty("E2_refresh_every_sweep_max_error",Precise(every_sweep.work.maximum_objective_refresh_error));
    RecordProperty("E2_refresh_every_sweep_max_updates",std::to_string(
        every_sweep.work.maximum_accepted_updates_between_objective_refreshes));
    ASSERT_TRUE(every_sweep.search_converged)<<every_sweep.reason;
    ASSERT_TRUE(every_sweep.endpoint_certified)<<every_sweep.endpoint_trust.reason;
    EXPECT_TRUE(every_sweep.assessment.independent_replay.passed);
    EXPECT_EQ(every_sweep.work.objective_cache_refreshes,every_sweep.sweep_count);
    EXPECT_EQ(every_sweep.work.objective_reconciliations,every_sweep.sweep_count);
    EXPECT_GT(every_sweep.work.sweep_replay_seconds,0.0);
    EXPECT_GT(every_sweep.work.sweep_replay_audit_seconds,0.0);
    EXPECT_GT(every_sweep.work.sweep_global_state_seconds,0.0);
    EXPECT_GT(every_sweep.work.objective_cache_refresh_seconds,0.0);
    EXPECT_LE(every_sweep.work.sweep_replay_seconds+every_sweep.work.sweep_replay_audit_seconds+
        every_sweep.work.sweep_global_state_seconds+every_sweep.work.objective_cache_refresh_seconds,
        every_sweep.search_seconds);
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveRemainsQualifiedAcrossInitialWidths)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    for(const double width:std::array<double,3>{.35,.55,.9})
    {
        const n::Vector initial_eta=n::Vector::Constant(3,std::log(width));
        n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=40;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
        policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
        auto f1_policy=policy;
        auto d1_policy=policy;
        d1_policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::FullScan;
        f1_policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
            const n::Indices &,n::VectorRef candidate_eta,n::VectorRef candidate_beta) {
            return n::FullReplayForTesting(data.domain,data.layout,data.y,
                candidate_eta,candidate_beta,data.context.scale);
        };
        const auto f1=n::RunFixedNeighborDiagnostic(problem,initial_eta,f1_policy);
        const auto d1=n::RunFixedNeighborDiagnostic(problem,initial_eta,d1_policy);
        policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
        const auto e1=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        ASSERT_TRUE(f1.search_converged)<<f1.reason;
        ASSERT_TRUE(d1.search_converged)<<d1.reason;
        ASSERT_TRUE(e1.search_converged)<<e1.reason;
        ASSERT_TRUE(f1.endpoint_certified)<<f1.endpoint_trust.reason;
        ASSERT_TRUE(d1.endpoint_certified)<<d1.endpoint_trust.reason;
        ASSERT_TRUE(e1.endpoint_certified)<<e1.endpoint_trust.reason;
        const auto prefix="initial_width_"+std::to_string(width);
        RecordProperty(prefix+"_F1_objective",Precise(f1.fit.objective.value()));
        RecordProperty(prefix+"_D1_objective",Precise(d1.fit.objective.value()));
        RecordProperty(prefix+"_E1_objective",Precise(e1.fit.objective.value()));
        RecordProperty(prefix+"_D1_sweeps",std::to_string(d1.sweep_count));
        RecordProperty(prefix+"_E1_sweeps",std::to_string(e1.sweep_count));
        EXPECT_TRUE(f1.assessment.independent_replay.passed);
        EXPECT_TRUE(d1.assessment.independent_replay.passed);
        EXPECT_TRUE(e1.assessment.independent_replay.passed);
        EXPECT_EQ(e1.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        EXPECT_GT(e1.work.incremental_objective_evaluations,0u);
        EXPECT_EQ(e1.final_sweep->global_ac_kkt,d1.final_sweep->global_ac_kkt);
        EXPECT_EQ(e1.final_sweep->projected_width_gradient_inf_norm,
            d1.final_sweep->projected_width_gradient_inf_norm);
        EXPECT_TRUE((e1.state.eta.array()==d1.state.eta.array()).all());
        EXPECT_TRUE((e1.state.beta.array()==d1.state.beta.array()).all());
        EXPECT_TRUE(BitwiseEqual(f1.state.eta,d1.state.eta));
        EXPECT_TRUE(BitwiseEqual(f1.state.beta,d1.state.beta));
        (void)data;
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveHandlesSignedChargeCancellationAndLowAmplitude)
{
    auto input=second_stage_test::SyntheticJointWorkload("chain",3);
    for(std::size_t atom=0;atom<3;++atom)
        for(const auto & support:input.support[atom])
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,.5,2.5);
            if(atom==0) input.observations[support.row]-=2.0*basis.gaussian;
            if(atom!=1) input.observations[support.row]-=.4*basis.charge;
        }
    JointProblem problem(std::move(input));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=500;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
    policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
    policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::FullScan;
    const auto d1=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
    IncrementalObjectiveAudit signed_audit;
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        AuditIncrementalCandidate(data.layout,data.context.scale,event,signed_audit);
    };
    const auto e1=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    RecordProperty("signed_low_amplitude_D1_converged",std::to_string(d1.search_converged));
    RecordProperty("signed_low_amplitude_E1_converged",std::to_string(e1.search_converged));
    RecordProperty("signed_low_amplitude_D1_stop_reason",d1.reason);
    RecordProperty("signed_low_amplitude_E1_stop_reason",e1.reason);
    RecordProperty("signed_low_amplitude_D1_sweeps",std::to_string(d1.sweep_count));
    RecordProperty("signed_low_amplitude_E1_sweeps",std::to_string(e1.sweep_count));
    RecordProperty("signed_low_amplitude_E1_incremental_candidates",
        std::to_string(e1.work.incremental_objective_evaluations));
    RecordProperty("signed_low_amplitude_D1_accepted_updates",std::to_string(d1.total_accepted_local_updates));
    RecordProperty("signed_low_amplitude_E1_accepted_updates",std::to_string(e1.total_accepted_local_updates));
    RecordProperty("signed_low_amplitude_E1_candidates",std::to_string(signed_audit.candidates));
    RecordProperty("signed_low_amplitude_E1_candidate_accepted",std::to_string(signed_audit.accepted));
    RecordProperty("signed_low_amplitude_E1_candidate_rejected",std::to_string(signed_audit.rejected));
    RecordProperty("signed_low_amplitude_E1_reference_decision_disagreements",
        std::to_string(signed_audit.actual_reference_disagreements));
    RecordProperty("signed_low_amplitude_E1_maximum_objective_absolute_error",Precise(signed_audit.maximum_absolute));
    RecordProperty("signed_low_amplitude_E1_maximum_cache_drift",Precise(signed_audit.maximum_cache_drift));
    RecordProperty("signed_low_amplitude_E1_maximum_objective_ulp",std::to_string(signed_audit.maximum_ulp));
    RecordProperty("signed_low_amplitude_E1_recovery_retries",std::to_string(std::count_if(
        e1.blocks.begin(),e1.blocks.end(),[](const auto & block){return block.recovery_retry;})));
    RecordProperty("signed_low_amplitude_E1_recovery_successes",std::to_string(std::count_if(
        e1.blocks.begin(),e1.blocks.end(),[](const auto & block){return block.recovery_success;})));
    RecordProperty("signed_low_amplitude_E1_final_objective_cache",Precise(e1.state.objective));
    RecordProperty("signed_low_amplitude_E1_final_eta",VectorProperty(e1.state.eta));
    RecordProperty("signed_low_amplitude_E1_final_width",VectorProperty(e1.state.eta.array().exp()));
    RecordProperty("signed_low_amplitude_E1_final_beta",VectorProperty(e1.state.beta));
    if(d1.final_sweep)
    {
        RecordProperty("signed_low_amplitude_D1_global_ac_kkt",Precise(d1.final_sweep->global_ac_kkt));
        RecordProperty("signed_low_amplitude_D1_projected_width_kkt",
            Precise(d1.final_sweep->projected_width_gradient_inf_norm));
    }
    if(e1.final_sweep)
    {
        RecordProperty("signed_low_amplitude_E1_global_ac_kkt",Precise(e1.final_sweep->global_ac_kkt));
        RecordProperty("signed_low_amplitude_E1_projected_width_kkt",
            Precise(e1.final_sweep->projected_width_gradient_inf_norm));
        RecordProperty("signed_low_amplitude_E1_objective_replay_error",
            Precise(e1.final_sweep->objective_replay_error));
        RecordProperty("signed_low_amplitude_E1_prediction_replay_error",
            Precise(e1.final_sweep->cache_replay_error));
    }
    if(!d1.search_converged && !e1.search_converged)
    {
        EXPECT_EQ(d1.reason,"block-sweep-budget");
        EXPECT_EQ(e1.reason,"block-sweep-budget");
        ASSERT_TRUE(d1.final_sweep); ASSERT_TRUE(e1.final_sweep);
        EXPECT_GT(d1.final_sweep->projected_width_gradient_inf_norm,
            policy.convergence_criteria.global_width_gradient);
        EXPECT_GT(e1.final_sweep->projected_width_gradient_inf_norm,
            policy.convergence_criteria.global_width_gradient);
        EXPECT_FALSE(d1.endpoint_certified); EXPECT_FALSE(e1.endpoint_certified);
        EXPECT_NE(d1.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        EXPECT_NE(e1.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        EXPECT_NEAR(d1.final_sweep->projected_width_gradient_inf_norm,
            e1.final_sweep->projected_width_gradient_inf_norm,1e-12);
        EXPECT_EQ(signed_audit.actual_reference_disagreements,0u);
        RecordProperty("signed_low_amplitude_classification","legacy-unqualified-non-regression");
        return;
    }
    if(!d1.search_converged)
    {
        ASSERT_TRUE(e1.search_converged)<<e1.reason;
        ASSERT_TRUE(e1.endpoint_certified)<<e1.endpoint_trust.reason;
        EXPECT_TRUE(e1.assessment.independent_replay.passed);
        EXPECT_EQ(e1.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        RecordProperty("signed_low_amplitude_classification","E1-improvement");
        return;
    }
    ASSERT_TRUE(e1.search_converged)<<"E1 regression: "<<e1.reason;
    ASSERT_TRUE(d1.endpoint_certified)<<d1.endpoint_trust.reason;
    ASSERT_TRUE(e1.endpoint_certified)<<e1.endpoint_trust.reason;
    EXPECT_TRUE(d1.assessment.independent_replay.passed);
    EXPECT_TRUE(e1.assessment.independent_replay.passed);
    EXPECT_EQ(e1.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_GT(e1.work.incremental_objective_evaluations,0u);
    EXPECT_LE(std::abs(e1.state.beta(0)),1e-5);
    EXPECT_LT(std::abs(e1.state.beta(2*2+1)),.01);
    RecordProperty("signed_low_amplitude_D1_objective",Precise(*d1.fit.objective));
    RecordProperty("signed_low_amplitude_E1_objective",Precise(*e1.fit.objective));
    RecordProperty("signed_low_amplitude_E1_direct_cache_error",
        Precise(std::abs(e1.state.objective-*e1.fit.objective)));
    RecordProperty("signed_low_amplitude_local_conditioning_blocks",std::to_string(e1.local_conditioning.size()));
    (void)data;
}

TEST(JointFixedNeighborBlockCoordinateTest, SignedCancellationFixtureMatchesHistoricalInputDigest)
{
    auto input=second_stage_test::SyntheticJointWorkload("chain",3);
    for(std::size_t atom=0;atom<3;++atom)
        for(const auto & support:input.support[atom])
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,.5,2.5);
            if(atom==0) input.observations[support.row]-=2.0*basis.gaussian;
            if(atom!=1) input.observations[support.row]-=.4*basis.charge;
        }
    std::size_t support_entries{};
    for(const auto & support:input.support) support_entries+=support.size();
    EXPECT_EQ(input.atom_ids,(std::vector<std::string>{"1","2","3"}));
    EXPECT_EQ(input.row_ids.size(),1425u);
    EXPECT_EQ(input.observations.size(),1425u);
    EXPECT_EQ(support_entries,1545u);
    EXPECT_EQ(second_stage_test::SyntheticJointWorkloadHash(input),
        "0b5c3798029192eb18c137134db3214e9647fa87868433366bd744210cce34eb");
}

TEST(JointFixedNeighborBlockCoordinateTest, StandardQualificationFixturesMatchHistoricalInputDigests)
{
    const std::array<std::pair<int,std::string>,4> fixtures{{
        {3,"4b3381463b3e7d97099342a9de86df01bf435cf2d64b8ea96af9ac0ee2eab176"},
        {30,"13cb46a3575b92368d8393f168cf45cc7db9615e0e2cce06b7f97326a8580c57"},
        {300,"bfb80c5d252cee06f7e7f57906c39123adc007d37bccbe34fe546390268d8544"},
        {1000,"3d81218fbf38f8c661551631fbe4038694c67254b26283f51d10d49299ee33a3"},
    }};
    for(const auto & [count,expected_digest]:fixtures)
    {
        const auto input=second_stage_test::SyntheticJointWorkload("chain",count);
        EXPECT_EQ(second_stage_test::SyntheticJointWorkloadHash(input),expected_digest)<<count;
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, HardCaseD1E1ObjectiveNonRegression)
{
    using Input=rhbm_gem::core::JointProblemInput;
    struct Fixture {std::string name,digest; Input input; std::size_t core_atoms;};
    std::vector<Fixture> fixtures;
    fixtures.push_back({"high-overlap-cube-8",
        "2a1894dc5dc5f7789fbc802477e77e08f520a6c25f0c7eac6b42a4caf3f717b2",
        second_stage_test::SyntheticJointWorkload("cube",8),4});
    auto rank_deficient=second_stage_test::SyntheticJointWorkload("chain",3);
    for(auto & support:rank_deficient.support)
        for(auto & entry:support) entry.squared_distance=0.0;
    fixtures.push_back({"rank-deficient-constant-support",
        "fd39aa13434e60a2193641f97ff9bdd05669fb6bc1b57a4d52f454bdacaa7151",
        std::move(rank_deficient),2});
    Input near_singular;
    near_singular.atom_ids={"1","2"};
    near_singular.row_ids={"0","1","2","3","4"};
    near_singular.observations.assign(5,0.0);
    near_singular.support={{{0,.1},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,.1+1e-8},{1,.5-1e-8},{2,1.0+1e-8},{3,1.5-1e-8},{4,2.0+1e-8}}};
    near_singular.selection_domain.emplace(); near_singular.selection_domain->target_indices={0,1};
    for(std::size_t row=0;row<near_singular.observations.size();++row)
    {
        const auto first=n::EvaluateKernel(near_singular.support[0][row].squared_distance,.7,2.5);
        const auto second=n::EvaluateKernel(near_singular.support[1][row].squared_distance,.7,2.5);
        near_singular.observations[row]=1.0*first.gaussian+.2*first.charge+
            1.4*second.gaussian-.1*second.charge;
    }
    fixtures.push_back({"near-singular-overlapping-support",
        "fc9e6c0ec7f746d9de6221d1b680334bd870d7797012a666c32cdcdb467454ef",
        std::move(near_singular),2});
    auto low_signal=second_stage_test::SyntheticJointWorkload("chain",3);
    for(std::size_t row=0;row<low_signal.observations.size();++row)
        low_signal.observations[row]=1e-5*std::sin(.17*static_cast<double>(row));
    fixtures.push_back({"low-signal-noise-only",
        "0d2ecf601cb21d6ce3d4e23f1aa7f37b78f17fe7ec58a19c49b397374ad99d2f",
        std::move(low_signal),2});
    Input weak_width;
    weak_width.atom_ids={"1","2"};
    for(Eigen::Index row=0;row<12;++row)
    {
        weak_width.row_ids.push_back(std::to_string(row));
        const double first_square=.22*static_cast<double>(row);
        const double second_radius=.12*static_cast<double>(row)+.31;
        weak_width.support.resize(2);
        weak_width.support[0].push_back({static_cast<std::size_t>(row),first_square});
        weak_width.support[1].push_back({static_cast<std::size_t>(row),second_radius*second_radius});
    }
    const std::array<double,2> truth_widths{.72,1.08},truth_coefficients{1.25,.19};
    const std::array<double,2> second_coefficients{.83,.11};
    weak_width.observations.resize(12);
    for(Eigen::Index row=0;row<12;++row)
    {
        const double first_square=weak_width.support[0][static_cast<std::size_t>(row)].squared_distance;
        const double second_square=weak_width.support[1][static_cast<std::size_t>(row)].squared_distance;
        const auto first=n::EvaluateKernel(first_square,truth_widths[0],2.5);
        const auto second=n::EvaluateKernel(second_square,truth_widths[1],2.5);
        weak_width.observations[static_cast<std::size_t>(row)]=truth_coefficients[0]*first.gaussian+
            truth_coefficients[1]*first.charge+second_coefficients[0]*second.gaussian+
            second_coefficients[1]*second.charge;
    }
    fixtures.push_back({"weak-width-direction",
        "aa1d50f8e6dbcdaa22e401be857e849573ca0c6457d714583029003f23f9ef9c",
        std::move(weak_width),2});

    for(auto & fixture:fixtures)
    {
        EXPECT_EQ(second_stage_test::SyntheticJointWorkloadHash(fixture.input),fixture.digest)<<fixture.name;
        JointProblem problem(std::move(fixture.input));
        const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
        n::Vector initial_eta;
        if(fixture.name=="weak-width-direction")
        {initial_eta.resize(2); initial_eta<<std::log(.72),std::log(1.08);}
        else initial_eta=n::Vector::Constant(static_cast<Eigen::Index>(data.context.atom_ids.size()),std::log(.55));
        const auto run=[&](n::FixedNeighborObjectiveStrategyForTesting strategy,bool full_replay_reference,
            IncrementalObjectiveAudit * audit) {
            n::FixedNeighborPolicy policy; policy.core_atoms=fixture.core_atoms; policy.maximum_sweeps=100;
            policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
            policy.stop_after_stationarity=true; policy.certify_endpoint=true;
            policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
            policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
            policy.objective_strategy_for_testing=strategy;
            policy.refresh_objective_after_sweep_for_testing=strategy==
                n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
            if(full_replay_reference)
                policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
                    const n::Indices &,n::VectorRef eta,n::VectorRef beta) {
                    return n::FullReplayForTesting(data.domain,data.layout,data.y,eta,beta,data.context.scale);
                };
            if(audit)
                policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
                    AuditIncrementalCandidate(data.layout,data.context.scale,event,*audit);
                };
            return n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        };
        IncrementalObjectiveAudit e1_audit;
        const auto f1=run(n::FixedNeighborObjectiveStrategyForTesting::FullScan,true,nullptr);
        const auto d1=run(n::FixedNeighborObjectiveStrategyForTesting::FullScan,false,nullptr);
        const auto e1=run(n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect,false,&e1_audit);
        if(f1.state.eta.size()!=initial_eta.size() || d1.state.eta.size()!=initial_eta.size() ||
            e1.state.eta.size()!=initial_eta.size())
        {
            ADD_FAILURE()<<fixture.name<<" returned no complete accepted state: F1="<<f1.reason
                <<", D1="<<d1.reason<<", E1="<<e1.reason;
            continue;
        }
        const auto f1_replay=n::FullReplayForTesting(data.domain,data.layout,data.y,
            f1.state.eta,f1.state.beta,data.context.scale);
        const auto d1_replay=n::FullReplayForTesting(data.domain,data.layout,data.y,
            d1.state.eta,d1.state.beta,data.context.scale);
        const auto e1_replay=n::FullReplayForTesting(data.domain,data.layout,data.y,
            e1.state.eta,e1.state.beta,data.context.scale);
        const auto prefix="hard_case_"+fixture.name+"_";
        RecordProperty((prefix+"F1_reason").c_str(),f1.reason);
        RecordProperty((prefix+"D1_reason").c_str(),d1.reason);
        RecordProperty((prefix+"E1_reason").c_str(),e1.reason);
        RecordProperty((prefix+"digest").c_str(),fixture.digest);
        RecordProperty((prefix+"F1_objective").c_str(),Precise(f1_replay.objective));
        RecordProperty((prefix+"D1_objective").c_str(),Precise(d1_replay.objective));
        RecordProperty((prefix+"E1_objective").c_str(),Precise(e1_replay.objective));
        RecordProperty((prefix+"F1_converged").c_str(),std::to_string(f1.search_converged));
        RecordProperty((prefix+"D1_converged").c_str(),std::to_string(d1.search_converged));
        RecordProperty((prefix+"E1_converged").c_str(),std::to_string(e1.search_converged));
        RecordProperty((prefix+"F1_endpoint").c_str(),std::to_string(f1.endpoint_certified));
        RecordProperty((prefix+"D1_endpoint").c_str(),std::to_string(d1.endpoint_certified));
        RecordProperty((prefix+"E1_endpoint").c_str(),std::to_string(e1.endpoint_certified));
        RecordProperty((prefix+"F1_independent_replay").c_str(),
            std::to_string(f1.assessment.independent_replay.passed));
        RecordProperty((prefix+"D1_independent_replay").c_str(),
            std::to_string(d1.assessment.independent_replay.passed));
        RecordProperty((prefix+"E1_independent_replay").c_str(),
            std::to_string(e1.assessment.independent_replay.passed));
        RecordProperty((prefix+"E1_candidates").c_str(),std::to_string(e1_audit.candidates));
        RecordProperty((prefix+"E1_accepted").c_str(),std::to_string(e1_audit.accepted));
        RecordProperty((prefix+"E1_rejected").c_str(),std::to_string(e1_audit.rejected));
        RecordProperty((prefix+"E1_decision_disagreements").c_str(),
            std::to_string(e1_audit.actual_reference_disagreements));
        RecordProperty((prefix+"E1_candidate_max_abs_error").c_str(),Precise(e1_audit.maximum_absolute));
        RecordProperty((prefix+"F1_D1_max_beta_difference").c_str(),
            Precise((f1.state.beta-d1.state.beta).cwiseAbs().maxCoeff()));
        RecordProperty((prefix+"E1_D1_max_beta_difference").c_str(),
            Precise((e1.state.beta-d1.state.beta).cwiseAbs().maxCoeff()));
        RecordProperty((prefix+"F1_D1_max_width_difference").c_str(),
            Precise((f1.state.eta.array().exp()-d1.state.eta.array().exp()).abs().maxCoeff()));
        RecordProperty((prefix+"E1_D1_max_width_difference").c_str(),
            Precise((e1.state.eta.array().exp()-d1.state.eta.array().exp()).abs().maxCoeff()));
        if(f1.final_sweep && d1.final_sweep && e1.final_sweep)
        {
            RecordProperty((prefix+"sweeps").c_str(),std::to_string(e1.sweep_count));
            RecordProperty((prefix+"E1_eta_change_inf").c_str(),Precise(e1.final_sweep->eta_change_inf));
            RecordProperty((prefix+"E1_beta_scaled_change").c_str(),Precise(e1.final_sweep->beta_scaled_change));
            RecordProperty((prefix+"E1_max_log_width_step").c_str(),
                Precise(e1.final_sweep->maximum_accepted_log_width_step));
            RecordProperty((prefix+"E1_active_width_bounds").c_str(),
                std::to_string(e1.final_sweep->active_width_bound_count));
            RecordProperty((prefix+"F1_ac_kkt").c_str(),Precise(f1.final_sweep->global_ac_kkt));
            RecordProperty((prefix+"D1_ac_kkt").c_str(),Precise(d1.final_sweep->global_ac_kkt));
            RecordProperty((prefix+"E1_ac_kkt").c_str(),Precise(e1.final_sweep->global_ac_kkt));
            RecordProperty((prefix+"F1_width_kkt").c_str(),
                Precise(f1.final_sweep->projected_width_gradient_inf_norm));
            RecordProperty((prefix+"D1_width_kkt").c_str(),
                Precise(d1.final_sweep->projected_width_gradient_inf_norm));
            RecordProperty((prefix+"E1_width_kkt").c_str(),
                Precise(e1.final_sweep->projected_width_gradient_inf_norm));
        }
        for(const auto * result:{&f1,&d1,&e1})
        {
            EXPECT_TRUE(result->state.eta.allFinite())<<fixture.name<<':'<<result->reason;
            EXPECT_TRUE(result->state.beta.allFinite())<<fixture.name<<':'<<result->reason;
            EXPECT_TRUE(std::isfinite(result->state.objective))<<fixture.name<<':'<<result->reason;
            EXPECT_EQ(result->global_state_memory_estimate.status,n::MemoryStatus::Pass)<<fixture.name;
            if(result->search_converged)
            {
                EXPECT_TRUE(result->endpoint_certified)<<fixture.name<<':'<<result->endpoint_trust.reason;
                EXPECT_TRUE(result->assessment.independent_replay.passed)<<fixture.name;
            }
        }
        EXPECT_LE(std::abs(f1.state.objective-f1_replay.objective),
            n::BlockObjectiveReplayEnclosure(f1_replay.objective))<<fixture.name;
        EXPECT_LE(std::abs(d1.state.objective-d1_replay.objective),
            n::BlockObjectiveReplayEnclosure(d1_replay.objective))<<fixture.name;
        EXPECT_LE(std::abs(e1.state.objective-e1_replay.objective),
            n::BlockObjectiveReplayEnclosure(e1_replay.objective))<<fixture.name;
        EXPECT_EQ(e1_audit.actual_reference_disagreements,0u)<<fixture.name;
        if(d1.search_converged)
        {
            EXPECT_TRUE(e1.search_converged)<<fixture.name<<':'<<e1.reason;
            EXPECT_TRUE(e1.endpoint_certified)<<fixture.name<<':'<<e1.endpoint_trust.reason;
            EXPECT_TRUE(e1.assessment.independent_replay.passed)<<fixture.name;
        }
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, SignedCancellationHistoricalFailureIsNonRegressionAcrossReplayModes)
{
    auto input=second_stage_test::SyntheticJointWorkload("chain",3);
    for(std::size_t atom=0;atom<3;++atom)
        for(const auto & support:input.support[atom])
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,.5,2.5);
            if(atom==0) input.observations[support.row]-=2.0*basis.gaussian;
            if(atom!=1) input.observations[support.row]-=.4*basis.charge;
        }
    const auto input_hash=second_stage_test::SyntheticJointWorkloadHash(input);
    ASSERT_EQ(input_hash,"0b5c3798029192eb18c137134db3214e9647fa87868433366bd744210cce34eb");
    JointProblem problem(std::move(input));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(3,std::log(.55));
    const std::vector<double> truth_beta{0,-.2,2,.2,2,-.2};
    const char * trace_path=std::getenv("RHBM_SIGNED_TRACE_PATH");
    std::ofstream trace;
    if(trace_path) trace.open(trace_path);
    if(trace_path) ASSERT_TRUE(trace.good());
    trace << std::setprecision(17);
    trace << "record|method|budget|sweep|block|full_objective|cache_objective|global_ac_kkt|width_kkt|"
             "eta_change|ac_change|max_log_width_change|accepted|rejected|unchanged|recovery_attempts|"
             "recovery_successes|objective_replay_error|prediction_replay_error|active_width_bounds|"
             "local_stop|block_reason|candidate_before|candidate_objective|incremental_objective|"
             "global_delta|local_global_error|candidate_accepted|candidate_unchanged|profile_evaluation|"
             "profile_accepted|lm_ratio|acceptance_reason|stop_reason|details\n";
    const auto run=[&](const std::string & method,std::size_t budget,
        n::FixedNeighborObjectiveStrategyForTesting strategy,bool full_replay_reference,
        IncrementalObjectiveAudit * audit_out) {
        n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=budget;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
        policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
        policy.capture_local_trajectory=true;
        policy.objective_strategy_for_testing=strategy;
        IncrementalObjectiveAudit audit;
        std::vector<std::size_t> candidate_accepts(budget+1),candidate_rejects(budget+1),candidate_unchanged(budget+1);
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            AuditIncrementalCandidate(data.layout,data.context.scale,event,audit);
            if(event.unchanged) ++candidate_unchanged.at(event.sweep);
            else if(event.accepted) ++candidate_accepts.at(event.sweep);
            else ++candidate_rejects.at(event.sweep);
            trace << "candidate|" << method << '|' << budget << '|' << event.sweep << '|' << event.block
                << "|||||||||||||||||" << event.accepted_before.objective << '|' << event.full_objective
                << '|' << event.incremental_objective << '|' << event.global_replay_delta << '|'
                << event.local_global_delta_error << '|' << event.accepted << '|' << event.unchanged
                << "||||||" << (event.objective_full_scan_fallback ? "full-scan-fallback" : "direct") << '\n';
        };
        policy.state_observer=[&](std::size_t sweep_index,const n::BlockCoordinateState & state,
            const n::FixedNeighborBlockSweep & sweep,const std::vector<n::FixedNeighborBlockRecord> & blocks) {
            const auto replay=n::FullReplayForTesting(data.domain,data.layout,data.y,
                state.eta,state.beta,data.context.scale);
            std::size_t block_accepted{},block_unchanged{},block_failed{};
            for(const auto & block:blocks) if(block.sweep==sweep_index)
            {
                block_accepted+=block.accepted;
                block_unchanged+=block.status=="unchanged";
                block_failed+=block.status=="failed";
                trace << "block|" << method << '|' << budget << '|' << block.sweep << '|' << block.block
                    << "||||||||" << block.accepted << '|' << (block.status=="failed") << '|'
                    << (block.status=="unchanged") << '|' << block.recovery_attempt << '|'
                    << block.recovery_success << "||||||||||||" << block.local_search_stop_reason << '|'
                    << block.reason << "||||||||" << block.profile_trials.size() << '\n';
                for(const auto & trial:block.profile_trials)
                    trace << "trial|" << method << '|' << budget << '|' << block.sweep << '|' << block.block
                        << "|||||||||||||||||||||||||" << trial.profile_evaluation << '|' << trial.accepted
                        << '|' << trial.lm_ratio << '|' << trial.acceptance_reason << "||\n";
            }
            trace << "sweep|" << method << '|' << budget << '|' << sweep_index << "|"
                << "|" << replay.objective << '|' << state.objective << '|' << sweep.global_ac_kkt << '|'
                << sweep.projected_width_gradient_inf_norm << '|' << sweep.eta_change_inf << '|'
                << sweep.beta_scaled_change << '|' << sweep.maximum_accepted_log_width_step << '|'
                << candidate_accepts.at(sweep_index) << '|' << candidate_rejects.at(sweep_index) << '|'
                << candidate_unchanged.at(sweep_index) << '|' << sweep.recovery_attempts << '|'
                << sweep.recovery_successes << '|' << sweep.objective_replay_error << '|'
                << sweep.cache_replay_error << '|' << sweep.active_width_bound_count
                << "|||||||||||||full-replay-audit,block_accepts=" << block_accepted
                << ",block_unchanged=" << block_unchanged << ",block_failures=" << block_failed << '\n';
        };
        if(full_replay_reference)
            policy.candidate_replay_for_testing=[&](const n::BlockCoordinateState &,const n::Indices &,
                const n::Indices &,n::VectorRef eta,n::VectorRef beta) {
                return n::FullReplayForTesting(data.domain,data.layout,data.y,eta,beta,data.context.scale);
            };
        const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        const auto full=n::FullReplayForTesting(data.domain,data.layout,data.y,
            result.state.eta,result.state.beta,data.context.scale);
        const double beta_error=(result.state.beta-n::Vector::Map(truth_beta.data(),6)).norm()/std::sqrt(6.0);
        const double width_error=(result.state.eta.array().exp()-.5).matrix().norm()/std::sqrt(3.0);
        if(!result.final_sweep) {ADD_FAILURE() << "Missing final sweep for " << method; return result;}
        trace << "summary|" << method << '|' << budget << '|' << result.sweep_count << "||"
            << full.objective << '|' << result.state.objective << '|' << result.final_sweep->global_ac_kkt << '|'
            << result.final_sweep->projected_width_gradient_inf_norm << '|' << result.final_sweep->eta_change_inf << '|'
            << result.final_sweep->beta_scaled_change << '|' << result.final_sweep->maximum_accepted_log_width_step
            << "||||" << result.final_sweep->recovery_attempts << '|' << result.final_sweep->recovery_successes << "||||"
            << result.final_sweep->active_width_bound_count << "||||||||||||" << result.reason << '|'
            << "converged=" << result.search_converged << ",endpoint=" << result.endpoint_certified
            << ",independent_replay=" << result.assessment.independent_replay.passed
            << ",objective_error=" << std::abs(result.state.objective-full.objective)
            << ",beta_truth_rmse=" << beta_error << ",width_truth_rmse=" << width_error
            << ",candidate_count=" << audit.candidates << ",candidate_accepts=" << audit.accepted
            << ",candidate_rejects=" << audit.rejected
            << ",actual_decision_disagreements=" << audit.actual_reference_disagreements
            << ",max_candidate_objective_error=" << audit.maximum_absolute
            << ",max_candidate_cache_drift=" << audit.maximum_cache_drift
            << ",max_candidate_objective_ulp=" << audit.maximum_ulp
            << ",input_sha256=" << input_hash << ",scale=" << data.context.scale
            << ",memory_status=" << static_cast<int>(result.global_state_memory_estimate.status)
            << ",memory_peak=" << result.global_state_memory_estimate.estimated_peak_bytes << '\n';
        EXPECT_TRUE(result.state.eta.allFinite()); EXPECT_TRUE(result.state.beta.allFinite());
        EXPECT_TRUE(std::isfinite(result.state.objective)); EXPECT_TRUE(std::isfinite(full.objective));
        EXPECT_EQ(result.global_state_memory_estimate.status,n::MemoryStatus::Pass);
        EXPECT_FALSE(result.work.full_replay_memory_fallback);
        if(result.search_converged)
        {
            EXPECT_TRUE(result.endpoint_certified);
            EXPECT_TRUE(result.assessment.independent_replay.passed);
        }
        else
        {
            EXPECT_EQ(result.reason,"block-sweep-budget");
            EXPECT_GT(result.final_sweep->projected_width_gradient_inf_norm,
                policy.convergence_criteria.global_width_gradient);
        }
        if(audit_out) *audit_out=audit;
        return result;
    };
    for(const std::size_t budget:{200u,500u,1000u})
    {
        const auto f1=run("F1",budget,n::FixedNeighborObjectiveStrategyForTesting::FullScan,true,nullptr);
        const auto d1=run("D1",budget,n::FixedNeighborObjectiveStrategyForTesting::FullScan,false,nullptr);
        IncrementalObjectiveAudit e1_audit;
        const auto e1=run("E1-Direct",budget,n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect,
            false,&e1_audit);
        const auto f1_full=n::FullReplayForTesting(data.domain,data.layout,data.y,f1.state.eta,f1.state.beta,data.context.scale);
        const auto d1_full=n::FullReplayForTesting(data.domain,data.layout,data.y,d1.state.eta,d1.state.beta,data.context.scale);
        const auto e1_full=n::FullReplayForTesting(data.domain,data.layout,data.y,e1.state.eta,e1.state.beta,data.context.scale);
        ASSERT_EQ(f1.search_converged,d1.search_converged);
        ASSERT_EQ(f1.reason,d1.reason);
        EXPECT_LE(std::abs(d1_full.objective-f1_full.objective),n::BlockObjectiveReplayEnclosure(f1_full.objective));
        EXPECT_LE(std::abs(e1_full.objective-d1_full.objective),n::BlockObjectiveReplayEnclosure(d1_full.objective));
        EXPECT_LE(std::abs(e1.state.objective-e1_full.objective),n::BlockObjectiveReplayEnclosure(e1_full.objective));
        EXPECT_EQ(e1_audit.candidates,e1.work.incremental_objective_evaluations);
        EXPECT_EQ(e1_audit.actual_reference_disagreements,0u);
        if(!e1.search_converged)
        {
            EXPECT_EQ(e1.reason,d1.reason);
            EXPECT_NEAR(e1.final_sweep->global_ac_kkt,d1.final_sweep->global_ac_kkt,1e-12);
            EXPECT_NEAR(e1.final_sweep->projected_width_gradient_inf_norm,
                d1.final_sweep->projected_width_gradient_inf_norm,1e-12);
        }
        const double maximum_beta_difference=(e1.state.beta-d1.state.beta).cwiseAbs().maxCoeff();
        const double maximum_width_difference=(e1.state.eta.array().exp()-d1.state.eta.array().exp()).abs().maxCoeff();
        const double f1_d1_beta_difference=(f1.state.beta-d1.state.beta).cwiseAbs().maxCoeff();
        const double f1_d1_width_difference=(f1.state.eta.array().exp()-d1.state.eta.array().exp()).abs().maxCoeff();
        RecordProperty(("signed_"+std::to_string(budget)+"_f1_outcome").c_str(),f1.reason);
        RecordProperty(("signed_"+std::to_string(budget)+"_d1_outcome").c_str(),d1.reason);
        RecordProperty(("signed_"+std::to_string(budget)+"_e1_outcome").c_str(),e1.reason);
        RecordProperty(("signed_"+std::to_string(budget)+"_F1_full_objective").c_str(),Precise(f1_full.objective));
        RecordProperty(("signed_"+std::to_string(budget)+"_D1_full_objective").c_str(),Precise(d1_full.objective));
        RecordProperty(("signed_"+std::to_string(budget)+"_E1_full_objective").c_str(),Precise(e1_full.objective));
        RecordProperty(("signed_"+std::to_string(budget)+"_E1_cache_error").c_str(),
            Precise(std::abs(e1.state.objective-e1_full.objective)));
        RecordProperty(("signed_"+std::to_string(budget)+"_E1_max_beta_difference").c_str(),Precise(maximum_beta_difference));
        RecordProperty(("signed_"+std::to_string(budget)+"_E1_max_width_difference").c_str(),Precise(maximum_width_difference));
        RecordProperty(("signed_"+std::to_string(budget)+"_F1_D1_max_beta_difference").c_str(),Precise(f1_d1_beta_difference));
        RecordProperty(("signed_"+std::to_string(budget)+"_F1_D1_max_width_difference").c_str(),Precise(f1_d1_width_difference));
        RecordProperty(("signed_"+std::to_string(budget)+"_E1_candidate_objective_error").c_str(),
            Precise(e1_audit.maximum_absolute));
        RecordProperty(("signed_"+std::to_string(budget)+"_E1_candidate_decision_disagreements").c_str(),
            std::to_string(e1_audit.actual_reference_disagreements));
        EXPECT_LE(f1_d1_beta_difference,1e-12);
        EXPECT_LE(f1_d1_width_difference,1e-12);
        EXPECT_LE(maximum_beta_difference,1e-12);
        EXPECT_LE(maximum_width_difference,1e-12);
        EXPECT_LE(e1_audit.maximum_absolute,n::BlockObjectiveReplayEnclosure(d1_full.objective));
    }
    trace.close();
    if(trace_path) ASSERT_TRUE(trace.good());
}

TEST(JointFixedNeighborBlockCoordinateTest, AffectedReplayFallsBackBeforeScratchExceedsMemoryBudget)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=5;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.collect_telemetry=true; policy.eta_upper_bound=std::log(2.2);
    const auto affected=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    ASSERT_FALSE(affected.work.full_replay_memory_fallback);
    ASSERT_EQ(affected.global_state_memory_estimate.status,n::MemoryStatus::Pass);

    auto full_replay_shape=affected.global_state_memory_estimate.shape;
    full_replay_shape.affected_replay_scratch=false;
    const auto full_replay_plan=n::PlanMemory(full_replay_shape,n::MemoryStage::GlobalStateEvaluation);
    ASSERT_EQ(full_replay_plan.status,n::MemoryStatus::Pass);
    const auto affected_at_full_budget=n::PlanMemory(affected.global_state_memory_estimate.shape,
        n::MemoryStage::GlobalStateEvaluation,full_replay_plan.estimated_peak_bytes);
    EXPECT_EQ(affected_at_full_budget.status,n::MemoryStatus::Reject);

    auto fallback_policy=policy;
    fallback_policy.global_state_memory_budget_for_testing=full_replay_plan.estimated_peak_bytes;
    fallback_policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
    std::size_t fallback_candidate_count{};
    fallback_policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        if(!event.candidate_replay_ran) return;
        ++fallback_candidate_count;
        const auto before=n::FullReplayForTesting(data.domain,data.layout,data.y,
            event.accepted_before.eta,event.accepted_before.beta,data.context.scale);
        const auto candidate=n::FullReplayForTesting(data.domain,data.layout,data.y,
            event.candidate_eta,event.candidate_beta,data.context.scale);
        EXPECT_LE(std::abs(event.accepted_before.objective-before.objective),
            n::BlockObjectiveReplayEnclosure(before.objective));
        EXPECT_DOUBLE_EQ(event.full_objective,candidate.objective);
        EXPECT_FALSE(std::isfinite(event.incremental_objective));
    };
    const auto fallback=n::RunFixedNeighborDiagnostic(problem,eta,fallback_policy);
    ASSERT_TRUE(fallback.work.full_replay_memory_fallback);
    EXPECT_GT(fallback.work.full_candidate_replays,0u);
    EXPECT_EQ(fallback.work.affected_candidate_replays,0u);
    EXPECT_EQ(fallback.work.incremental_objective_evaluations,0u);
    EXPECT_GT(fallback_candidate_count,0u);
    EXPECT_EQ(affected.reason,fallback.reason);
    EXPECT_EQ(affected.sweep_count,fallback.sweep_count);
    EXPECT_TRUE(BitwiseEqual(affected.state.eta,fallback.state.eta));
    EXPECT_TRUE(BitwiseEqual(affected.state.beta,fallback.state.beta));
    EXPECT_TRUE(BitwiseEqual(affected.state.prediction,fallback.state.prediction));
    EXPECT_TRUE(BitwiseEqual(affected.state.residual,fallback.state.residual));
    EXPECT_TRUE(BitwiseEqual(affected.state.objective,fallback.state.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, E1MemoryPreflightRejectsBeforeCreatingAcceptedState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=2;
    policy.eta_upper_bound=std::log(2.2);
    policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
    policy.global_state_memory_budget_for_testing=1;
    std::size_t candidates{};
    policy.candidate_observer=[&](const auto &) {++candidates;};
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        eta,data.context,policy);

    EXPECT_NE(result.reason.find("resource-limit"),std::string::npos);
    EXPECT_EQ(result.memory_estimate.status,n::MemoryStatus::Reject);
    EXPECT_EQ(result.state.eta.size(),0);
    EXPECT_EQ(result.state.beta.size(),0);
    EXPECT_EQ(candidates,0u);
}

TEST(JointFixedNeighborBlockCoordinateTest, E1NonfiniteObservationsFailBeforeAcceptedStateCreation)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    n::Vector observations=data.y;
    observations(0)=std::numeric_limits<double>::quiet_NaN();
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=2;
    policy.eta_upper_bound=std::log(2.2);
    policy.objective_strategy_for_testing=n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
    std::size_t candidates{};
    policy.candidate_observer=[&](const auto &) {++candidates;};
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,observations,
        eta,data.context,policy);

    EXPECT_EQ(result.reason,"block-invalid-partition");
    EXPECT_EQ(result.state.eta.size(),0);
    EXPECT_EQ(result.state.beta.size(),0);
    EXPECT_EQ(candidates,0u);
}

TEST(JointFixedNeighborBlockCoordinateTest, E1NonfiniteObjectiveDoesNotAcceptOrPartiallyCommitCandidate)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector observations=n::Vector::Constant(data.y.size(),1e200);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    struct Run
    {
        n::FixedNeighborSearchResult result;
        std::size_t candidates{},accepted{};
    };
    const auto run=[&](n::FixedNeighborObjectiveStrategyForTesting strategy) {
        Run output;
        n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=2;
        policy.eta_upper_bound=std::log(2.2);
        policy.objective_strategy_for_testing=strategy;
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            if(event.candidate_replay_ran) ++output.candidates;
            output.accepted+=event.accepted;
        };
        output.result=n::SearchFixedNeighborComponent(data.domain,data.layout,observations,
            eta,data.context,policy);
        return output;
    };
    const auto d1=run(n::FixedNeighborObjectiveStrategyForTesting::FullScan);
    const auto e1=run(n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect);
    const auto initial_replay=n::FullReplayForTesting(data.domain,data.layout,observations,
        eta,n::Vector::Zero(6),data.context.scale);

    ASSERT_TRUE(std::isinf(initial_replay.objective));
    EXPECT_FALSE(d1.result.search_converged);
    EXPECT_FALSE(e1.result.search_converged);
    EXPECT_EQ(d1.result.reason,"block-cache-replay-failed");
    EXPECT_EQ(d1.result.reason,e1.result.reason);
    EXPECT_TRUE(std::isinf(d1.result.state.objective));
    EXPECT_TRUE(std::isinf(e1.result.state.objective));
    EXPECT_EQ(d1.candidates,0u);
    EXPECT_EQ(e1.candidates,0u);
    EXPECT_EQ(d1.accepted,0u);
    EXPECT_EQ(e1.accepted,0u);
    EXPECT_TRUE(BitwiseEqual(d1.result.state.eta,eta));
    EXPECT_TRUE(BitwiseEqual(e1.result.state.eta,eta));
    EXPECT_TRUE(BitwiseEqual(d1.result.state.beta,n::Vector::Zero(6)));
    EXPECT_TRUE(BitwiseEqual(e1.result.state.beta,n::Vector::Zero(6)));
}

TEST(JointFixedNeighborBlockCoordinateTest, UnchangedCandidateSkipsOnlyCandidateReplay)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=30;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.collect_telemetry=true; policy.eta_upper_bound=std::log(2.2);
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    EXPECT_GT(result.work.candidate_replay_skipped_unchanged,0u);
    EXPECT_LT(result.work.full_candidate_replays+result.work.affected_candidate_replays,result.total_block_solves);
    EXPECT_LE(result.work.full_candidate_replays+result.work.affected_candidate_replays+
        result.work.candidate_replay_skipped_unchanged,
        result.total_block_solves);
    EXPECT_EQ(result.work.sweep_replays,result.sweeps.size());
    std::size_t unchanged_records{};
    for(const auto & block:result.blocks)
        if(block.reason=="block-unchanged")
        {
            ++unchanged_records;
            EXPECT_TRUE(BitwiseEqual(block.global_replay_delta,0.0));
            EXPECT_LE(block.local_global_delta_error,block.objective_replay_enclosure);
        }
    EXPECT_EQ(unchanged_records,result.work.candidate_replay_skipped_unchanged);
}


TEST(JointFixedNeighborBlockCoordinateTest, ProductionBoundedSearchReleasesTheStartingUpperFace)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const double eta_max=std::log(.55);
    const n::Vector initial=n::Vector::Constant(3,eta_max);
    n::FixedNeighborPolicy policy;
    policy.core_atoms=4;
    policy.maximum_sweeps=1;
    policy.certify_endpoint=false;
    policy.stop_after_stationarity=false;
    policy.capture_local_trajectory=true;
    policy.eta_upper_bound=eta_max;
    const auto result=n::RunFixedNeighborDiagnostic(problem,initial,policy);
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto & bounded_lm=result.work.local_profile_work.bounded_lm;
    EXPECT_GT(bounded_lm.subproblem_solves,0u);
    EXPECT_EQ(bounded_lm.failed_subproblems,0u);
    EXPECT_GT(bounded_lm.active_set_iterations,0u);
    EXPECT_TRUE(std::isfinite(bounded_lm.maximum_subproblem_kkt));
    ASSERT_EQ(result.blocks.size(),1u)<<result.reason;
    ASSERT_EQ(result.state.eta.size(),3);
    EXPECT_LE(result.state.eta.maxCoeff(),eta_max+1e-14);
    bool bounded_subproblem_seen=false,released=false;
    for(const auto & trial:result.blocks.front().profile_trials)
        if(trial.bounded_subproblem)
        {
            bounded_subproblem_seen=true;
            EXPECT_EQ(trial.bounded_subproblem_status,"solved");
            EXPECT_LE(trial.lm_subproblem_kkt_inf_norm,1e-9);
            released|=trial.active_set_releases>0;
        }
    EXPECT_TRUE(bounded_subproblem_seen);
    EXPECT_TRUE(released);
    ASSERT_TRUE(result.final_sweep);
    EXPECT_TRUE(result.final_sweep->width_bound_applied);
    auto context=data.context; context.eta_upper_bound=eta_max;
    const auto primary=n::EvaluateStateStreaming(data.domain,data.y,result.state.eta,result.state.beta,context);
    ASSERT_TRUE(primary.valid)<<primary.reason;
    const auto replay=n::CheckIndependentReplay(data.domain,data.y,primary,context);
    ASSERT_TRUE(replay.available && replay.valid)<<replay.reason;
    EXPECT_TRUE(replay.width_bound_applied);
    EXPECT_TRUE(replay.width_bound_feasible);
    EXPECT_DOUBLE_EQ(result.final_sweep->global_width_gradient_inf_norm,
        primary.gradient.lpNorm<Eigen::Infinity>());
    EXPECT_DOUBLE_EQ(result.final_sweep->projected_width_gradient_inf_norm,
        replay.independent_width_projected_kkt);
    EXPECT_DOUBLE_EQ(result.final_sweep->stationarity_width_gradient_inf_norm,
        result.final_sweep->projected_width_gradient_inf_norm);
    EXPECT_EQ(result.final_sweep->active_width_bound_count,replay.active_width_bounds);
}

TEST(JointFixedNeighborBlockCoordinateTest, AssembledBoundedWidthStationarityUsesIndependentProjectedKkt)
{
    auto problem=PureExclusionChain(3);
    const double eta_max=std::log(.4);
    n::FixedNeighborSearchPolicy policy;
    policy.eta_upper_bound=eta_max;
    const auto fit=n::FitFixedNeighborComponentsImpl(problem,std::vector<double>(3,.4),policy);
    const auto stationarity=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="width-stationarity" &&
            check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    ASSERT_NE(stationarity,fit.evidence.end());
    EXPECT_EQ(stationarity->reason,"independent-bound-aware-projected-width-kkt");
    ASSERT_TRUE(stationarity->value.has_value());
    const auto bound_feasibility=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="width-bound-feasibility" &&
            check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    ASSERT_NE(bound_feasibility,fit.evidence.end());
    EXPECT_EQ(bound_feasibility->status,rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(fit.assembled_state.has_value());
    double projected_inf{},raw_inf{};
    bool active_bound{};
    for(std::size_t atom=0;atom<fit.assembled_state->log_b.size();++atom)
    {
        const double eta=fit.assembled_state->log_b[atom];
        const double gradient=fit.assembled_state->width_gradient[atom];
        projected_inf=std::max(projected_inf,std::abs(eta-std::min(eta-gradient,eta_max)));
        raw_inf=std::max(raw_inf,std::abs(gradient));
        active_bound|=eta>=eta_max-1e-12;
    }
    EXPECT_TRUE(active_bound);
    EXPECT_GT(raw_inf,n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_NEAR(*stationarity->value,projected_inf,1e-12);
    EXPECT_LE(*stationarity->value,n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_EQ(stationarity->status,rhbm_gem::JointCheckStatus::Passed);
}

void CheckMultiBlockBudgetFailure(const std::string & topology)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    n::StructuralTopology structural_topology(data.domain,data.context.atom_ids,data.layout);
    const auto partition=structural_topology.BuildCores(16);
    ASSERT_EQ(partition.cores.size(),2u);
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=0;
    policy.eta_upper_bound=std::log(2.2);
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-sweep-budget");
    EXPECT_FALSE(result.endpoint_certified);
    EXPECT_FALSE(result.fit.search_completed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_TRUE(result.fit.components.empty());
}
}
TEST(JointFixedNeighborBlockCoordinateTest, SmallSingleBlockFixturesUseExistingEndpointCertification)
{
    CheckSingleBlock("chain",8); CheckSingleBlock("cube",8);
    CheckSingleBlock("chain",32); CheckSingleBlock("cube",32);
}
TEST(JointFixedNeighborBlockCoordinateTest, MultiBlockBudgetFailureRemainsDistinctFromRuntimeConvergence)
{
    CheckMultiBlockBudgetFailure("chain"); CheckMultiBlockBudgetFailure("cube");
}
TEST(JointFixedNeighborBlockCoordinateTest, CheapStationarityNeedsCompleteSweepConfirmation)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=1;
    policy.eta_upper_bound=std::log(2.2);
    const auto candidate=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    EXPECT_FALSE(candidate.search_converged);
    EXPECT_EQ(candidate.reason,"block-sweep-budget");
    EXPECT_EQ(candidate.sweeps.size(),1u);
    EXPECT_EQ(candidate.first_order_stationarity_sweep,0u);
    EXPECT_EQ(candidate.confirmed_stationarity_sweep,0u);

    policy.maximum_sweeps=30;
    const auto confirmed=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    ASSERT_TRUE(confirmed.search_converged)<<confirmed.reason;
    EXPECT_GT(confirmed.first_order_stationarity_sweep,0u);
    EXPECT_GT(confirmed.confirmed_stationarity_sweep,confirmed.first_order_stationarity_sweep);
    EXPECT_EQ(confirmed.confirmed_stationarity_sweep,confirmed.sweeps.size());
}

TEST(JointFixedNeighborBlockCoordinateTest, StabilityHistoryKeepsOnlyTheLastThreeCompletedSweeps)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    for(std::size_t sweep_budget=0;sweep_budget<=5;++sweep_budget)
    {
        n::FixedNeighborPolicy policy;
        policy.maximum_sweeps=sweep_budget;
        policy.eta_upper_bound=std::log(2.2);
        policy.stop_after_stationarity=false;
        policy.collect_diagnostics=true;
        policy.assess_parameter_stability=true;
        const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
        ASSERT_EQ(result.sweep_count,sweep_budget);
        ASSERT_EQ(result.sweeps.size(),sweep_budget);
        const auto recent_count=std::min<std::size_t>(sweep_budget,3);
        EXPECT_EQ(result.recent_sweep_count,recent_count);
        double maximum_eta_change{},maximum_scaled_ac_change{};
        for(std::size_t k=sweep_budget-recent_count;k<sweep_budget;++k)
        {
            maximum_eta_change=std::max(maximum_eta_change,result.sweeps[k].eta_change_inf);
            maximum_scaled_ac_change=std::max(maximum_scaled_ac_change,result.sweeps[k].beta_scaled_change);
        }
        EXPECT_DOUBLE_EQ(result.maximum_eta_change,maximum_eta_change);
        EXPECT_DOUBLE_EQ(result.maximum_scaled_ac_change,maximum_scaled_ac_change);
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, FailedBlockRetriesAfterLaterBlocksUpdateTheNeighborState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(2,std::log(.55));
    struct RecoveryRun
    {
        n::FixedNeighborSearchResult result;
        std::optional<n::BlockCoordinateState> first_sweep_state,retry_candidate_before;
        double retry_candidate_objective{},retry_candidate_reference{};
        bool retry_candidate_seen{};
    };
    const auto run=[&](n::FixedNeighborObjectiveStrategyForTesting strategy) {
        RecoveryRun output;
        n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
        policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
        policy.objective_strategy_for_testing=strategy;
        policy.refresh_objective_after_sweep_for_testing=
            strategy==n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
        policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t sweep) {
            return block==0 && sweep==0;
        };
        policy.state_observer=[&](std::size_t sweep,const auto & state,const auto &,const auto &) {
            if(sweep==1) output.first_sweep_state=state;
        };
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            if(event.sweep==2 && event.block==1)
            {
                output.retry_candidate_seen=true;
                output.retry_candidate_before=event.accepted_before;
                output.retry_candidate_objective=event.full_objective;
                const auto replay=n::FullReplayForTesting(data.domain,data.layout,data.y,
                    event.candidate_eta,event.candidate_beta,data.context.scale);
                output.retry_candidate_reference=replay.objective;
            }
        };
        output.result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
            initial_eta,data.context,policy);
        return output;
    };
    const auto d1=run(n::FixedNeighborObjectiveStrategyForTesting::FullScan);
    const auto e1=run(n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect);
    for(const auto * current:{&d1,&e1})
    {
        ASSERT_TRUE(current->first_sweep_state.has_value());
        ASSERT_TRUE(current->retry_candidate_seen);
        ASSERT_TRUE(current->retry_candidate_before.has_value());
        EXPECT_GT(current->first_sweep_state->beta.norm(),0.0);
        EXPECT_EQ(current->result.recovery_attempts,1u);
        EXPECT_EQ(current->result.recovery_successes,1u);
        const auto first_failure=std::find_if(current->result.blocks.begin(),current->result.blocks.end(),[](const auto & block) {
            return block.block==1 && block.sweep==1;
        });
        ASSERT_NE(first_failure,current->result.blocks.end());
        EXPECT_EQ(first_failure->status,"retry-pending");
        const auto retried=std::find_if(current->result.blocks.begin(),current->result.blocks.end(),[](const auto & block) {
            return block.block==1 && block.sweep==2;
        });
        ASSERT_NE(retried,current->result.blocks.end());
        EXPECT_TRUE(retried->recovery_retry);
        EXPECT_EQ(retried->recovery_attempt,1u);
        EXPECT_LE(std::abs(retried->objective_before-current->first_sweep_state->objective),
            n::BlockObjectiveReplayEnclosure(current->first_sweep_state->objective));
        const auto trusted=n::FullReplayForTesting(data.domain,data.layout,data.y,
            current->first_sweep_state->eta,current->first_sweep_state->beta,data.context.scale);
        EXPECT_LE(std::abs(current->first_sweep_state->objective-trusted.objective),
            n::BlockObjectiveReplayEnclosure(trusted.objective));
        EXPECT_LE(std::abs(current->retry_candidate_before->objective-trusted.objective),
            n::BlockObjectiveReplayEnclosure(trusted.objective));
        EXPECT_LE(std::abs(current->retry_candidate_objective-current->retry_candidate_reference),
            n::BlockObjectiveReplayEnclosure(current->retry_candidate_reference));
    }
    EXPECT_GT(e1.result.work.incremental_objective_evaluations,0u);
    EXPECT_TRUE(BitwiseEqual(d1.result.state.eta,e1.result.state.eta));
    EXPECT_TRUE(BitwiseEqual(d1.result.state.beta,e1.result.state.beta));
    EXPECT_TRUE(BitwiseEqual(d1.result.state.objective,e1.result.state.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, ExhaustedBlockRecoveryLeavesTheLastAcceptedComponentStateUntouched)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(2,std::log(.55));
    std::array<n::BlockCoordinateState,2> accepted_states;
    std::array<n::FixedNeighborSearchResult,2> results;
    for(std::size_t i=0;i<results.size();++i)
    {
        n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
        policy.eta_upper_bound=std::log(2.2);
        policy.objective_strategy_for_testing=i==0 ? n::FixedNeighborObjectiveStrategyForTesting::FullScan :
            n::FixedNeighborObjectiveStrategyForTesting::IncrementalDirect;
        policy.refresh_objective_after_sweep_for_testing=i==1;
        policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t) {return block==0;};
        policy.state_observer=[&,i](std::size_t sweep,const auto & state,const auto &,const auto &) {
            if(sweep==1) accepted_states[i]=state;
        };
        results[i]=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
            initial_eta,data.context,policy);
    }
    for(std::size_t i=0;i<results.size();++i)
    {
        const auto & result=results[i];
        EXPECT_FALSE(result.search_converged);
        EXPECT_EQ(result.reason,"block-recovery-exhausted");
        EXPECT_EQ(result.recovery_attempts,1u);
        EXPECT_EQ(result.recovery_successes,0u);
        EXPECT_EQ(result.sweep_count,1u);
        EXPECT_EQ(result.first_order_stationarity_sweep,0u);
        EXPECT_TRUE((result.state.eta.array()==accepted_states[i].eta.array()).all());
        EXPECT_TRUE((result.state.beta.array()==accepted_states[i].beta.array()).all());
        const auto replay=n::FullReplayForTesting(data.domain,data.layout,data.y,
            result.state.eta,result.state.beta,data.context.scale);
        EXPECT_LE(std::abs(result.state.objective-replay.objective),
            n::BlockObjectiveReplayEnclosure(replay.objective));
        EXPECT_FALSE(result.certification.has_value());
    }
    EXPECT_TRUE(BitwiseEqual(results[0].state.eta,results[1].state.eta));
    EXPECT_TRUE(BitwiseEqual(results[0].state.beta,results[1].state.beta));
    EXPECT_TRUE(BitwiseEqual(results[0].state.objective,results[1].state.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, SearchResultAndIndependentCertificationStaySeparate)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(32,std::log(.55));
    const n::Domain domain=data.domain;
    const auto context=data.context;
    n::FixedNeighborPolicy quiet_policy; quiet_policy.core_atoms=16; quiet_policy.maximum_sweeps=3;
    quiet_policy.eta_upper_bound=std::log(2.2);
    const auto quiet=n::SearchFixedNeighborComponent(domain,data.layout,data.y,initial_eta,context,quiet_policy);
    auto diagnostic_policy=quiet_policy; diagnostic_policy.collect_diagnostics=true;
    const auto diagnostic=n::SearchFixedNeighborComponent(domain,data.layout,data.y,initial_eta,context,diagnostic_policy);
    const auto quiet_state=quiet.state;
    const auto quiet_certification=n::CertifyFixedNeighborEndpoint(domain,data.y,quiet.state,context);
    const auto diagnostic_certification=n::CertifyFixedNeighborEndpoint(domain,data.y,diagnostic.state,context);
    const auto resource_limited=n::CertifyFixedNeighborEndpoint(domain,data.y,quiet.state,context,nullptr,
        n::JointMemoryPlanningBudgetBytes+1);

    EXPECT_TRUE(quiet.blocks.empty());
    EXPECT_TRUE(quiet.sweeps.empty());
    ASSERT_FALSE(diagnostic.sweeps.empty());
    ASSERT_TRUE(quiet.final_sweep.has_value());
    EXPECT_EQ(quiet.sweep_count,diagnostic.sweep_count);
    EXPECT_EQ(quiet.total_block_solves,diagnostic.total_block_solves);
    EXPECT_EQ(quiet.total_profile_evaluations,diagnostic.total_profile_evaluations);
    EXPECT_EQ(quiet.block_preparations,2u);
    EXPECT_EQ(quiet.domain_preparations,quiet.block_preparations);
    EXPECT_EQ(quiet.mapping_preparations,quiet.block_preparations);
    EXPECT_EQ(diagnostic.block_preparations,quiet.block_preparations);
    EXPECT_GT(quiet.work.local_profile_work.total.evaluations,0u);
    EXPECT_EQ(quiet.work.local_profile_work.total.evaluations,
        diagnostic.work.local_profile_work.total.evaluations);
    EXPECT_EQ(quiet.reason,diagnostic.reason);
    EXPECT_EQ(quiet.search_converged,diagnostic.search_converged);
    EXPECT_EQ(quiet_certification.verified(),diagnostic_certification.verified());
    EXPECT_TRUE((quiet.state.eta.array()==quiet_state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==quiet_state.beta.array()).all());
    EXPECT_TRUE((quiet.state.eta.array()==diagnostic.state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==diagnostic.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(quiet.state.objective,diagnostic.state.objective);
    EXPECT_DOUBLE_EQ(quiet.final_sweep->global_ac_kkt,diagnostic.sweeps.back().global_ac_kkt);
    EXPECT_DOUBLE_EQ(quiet.final_sweep->global_width_gradient_inf_norm,
        diagnostic.sweeps.back().global_width_gradient_inf_norm);
    EXPECT_TRUE(resource_limited.resource_limited());
    EXPECT_EQ(resource_limited.verification_memory_estimate.status,n::MemoryStatus::Reject);
    EXPECT_EQ(resource_limited.resource_limit_reason(),
        "resource-limit:endpoint-verification:estimated-peak-exceeds-budget");
    EXPECT_TRUE((quiet_certification.assessment.primary.eta.array()==
        diagnostic_certification.assessment.primary.eta.array()).all());
    EXPECT_TRUE((quiet_certification.assessment.primary.eta.array()==quiet.state.eta.array()).all());
    EXPECT_TRUE((quiet_certification.assessment.primary.beta.array()==quiet.state.beta.array()).all());
}

TEST(JointFixedNeighborBlockCoordinateTest, ReplayStructureFollowsComponentLocalRows)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"a","b"};
    input.row_ids={"r0","r1","r2","r3","r4","r5","r6","r7","r8","r9","constant"};
    input.observations.resize(11);
    input.support.resize(2);
    const std::vector<double> squares{0.0,.25,1.0,2.25,4.0};
    for(std::size_t k=0;k<squares.size();++k)
        for(std::size_t atom=0;atom<2;++atom)
        {
            const auto row=2*k+atom;
            input.support[atom].push_back({row,squares[k]});
            const auto basis=n::EvaluateKernel(squares[k],.6,2.5);
            input.observations[row]=(atom==0 ? 1.3 : 2.1)*basis.gaussian+
                (atom==0 ? .1 : -.2)*basis.charge;
        }
    input.observations.back()=.4;

    const JointProblem problem(std::move(input));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.partition.components.size(),2u);
    EXPECT_EQ(data.partition.components[0].rows,(n::Indices{0,2,4,6,8}));
    EXPECT_EQ(data.partition.components[1].rows,(n::Indices{1,3,5,7,9}));
    EXPECT_EQ(data.partition.constant_rows,(n::Indices{10}));

    const auto fit=second_stage_test::matched::joint_abc::FitBoundedJointComponentsForTesting(
        problem,{.55,.55},2.2);
    ASSERT_EQ(fit.components.size(),2u);
    for(const auto & component:fit.components)
    {
        ASSERT_TRUE(component.candidate_diagnostic) << component.stop_reason;
        EXPECT_EQ(component.candidate_diagnostic->endpoint_verification_status,
            rhbm_gem::JointCheckStatus::Passed);
        ASSERT_TRUE(component.state) << component.stop_reason;
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, EmptyOrInvalidReplayRowMappingCannotAcceptState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(2,std::log(.55));
    auto check_invalid=[&](const rhbm_gem::JointParameterLayout & layout) {
        n::FixedNeighborPolicy policy; policy.eta_upper_bound=std::log(2.2);
        const auto result=n::SearchFixedNeighborComponent(data.domain,layout,data.y,eta,data.context,policy);
        EXPECT_EQ(result.reason,"block-invalid-partition");
        EXPECT_EQ(result.state.eta.size(),0);
        EXPECT_FALSE(result.search_converged);
    };
    auto empty_layout=data.layout;
    empty_layout.informative_rows.clear();
    check_invalid(empty_layout);
    auto invalid_layout=data.layout;
    invalid_layout.informative_rows.front()=static_cast<std::size_t>(data.domain.rows);
    check_invalid(invalid_layout);
}

TEST(JointFixedNeighborBlockCoordinateTest, BlockSinkReceivesRecordsWithoutRetainingThemInSearchResult)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(2,std::log(.55));
    n::FixedNeighborDiagnosticSink sink;
    std::vector<n::FixedNeighborBlockRecord> records;
    sink.block=[&](const auto & record){records.push_back(record);};
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=1;
    policy.eta_upper_bound=std::log(2.2);
    policy.diagnostic_sink=&sink;
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    EXPECT_TRUE(result.blocks.empty());
    ASSERT_EQ(records.size(),2u);
    EXPECT_EQ(records[0].core_atoms.size(),1u);
    EXPECT_EQ(records[1].core_atoms.size(),1u);
}

TEST(JointFixedNeighborBlockCoordinateTest, FixedNeighborAcceptsAtMostOneTrustedLocalUpdatePerBlockVisit)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=2;
    policy.eta_upper_bound=std::log(2.2);
    policy.certify_endpoint=false; policy.capture_local_trajectory=true;
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    ASSERT_FALSE(result.blocks.empty());
    for(const auto & block:result.blocks)
    {
        std::size_t accepted_updates{},rejected_trials{};
        for(const auto & trial:block.profile_trials)
        {
            if(trial.accepted && trial.accepted_update && *trial.accepted_update>0) ++accepted_updates;
            if(!trial.accepted) ++rejected_trials;
        }
        EXPECT_LE(accepted_updates,1u);
        EXPECT_EQ(accepted_updates,static_cast<std::size_t>(block.accepted_updates));
        EXPECT_LE(block.accepted_updates,1);
        if(rejected_trials>0) EXPECT_LE(accepted_updates,1u);
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, StrictAndPracticalStationarityCriteriaHaveInclusiveFiniteBoundaries)
{
    const auto & strict=n::StrictJointConvergenceCriteria;
    const auto & practical=n::PracticalJointConvergenceCriteria;
    EXPECT_TRUE(n::IsValidJointConvergenceCriteria(strict));
    EXPECT_TRUE(n::IsValidJointConvergenceCriteria(practical));

    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(5e-8,0,0,true,strict),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(5e-8,0,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,5e-10,0,true,strict),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,5e-10,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,0,5e-7,true,strict),
        n::FixedNeighborStationarityState::CandidateStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,0,5e-7,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);

    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(practical.global_ac_kkt,0,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(
        std::nextafter(practical.global_ac_kkt,0.0),0,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(
        std::nextafter(practical.global_ac_kkt,std::numeric_limits<double>::infinity()),0,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,practical.global_width_gradient,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,
        std::nextafter(practical.global_width_gradient,0.0),0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,
        std::nextafter(practical.global_width_gradient,std::numeric_limits<double>::infinity()),0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,0,practical.confirmed_log_width_change,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_TRUE(n::IsFixedNeighborEtaChangeConfirmed(
        std::nextafter(practical.confirmed_log_width_change,0.0),true,practical));
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(
        std::nextafter(practical.confirmed_log_width_change,std::numeric_limits<double>::infinity()),true,practical));

    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double infinity=std::numeric_limits<double>::infinity();
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(nan,0,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,infinity,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(-1,0,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,-1,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(-1,true,practical));
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(infinity,true,practical));
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(nan,true,practical));
    EXPECT_NE(n::ClassifyFixedNeighborStationarity(0,0,-1,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_NE(n::ClassifyFixedNeighborStationarity(0,0,nan,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);

    auto invalid=practical; invalid.global_ac_kkt=nan;
    EXPECT_FALSE(n::IsValidJointConvergenceCriteria(invalid));
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,0,0,true,invalid),
        n::FixedNeighborStationarityState::NotStationary);
    invalid=practical; invalid.global_width_gradient=infinity;
    EXPECT_FALSE(n::IsValidJointConvergenceCriteria(invalid));
    invalid=practical; invalid.confirmed_log_width_change=-1;
    EXPECT_FALSE(n::IsValidJointConvergenceCriteria(invalid));
}
TEST(JointFixedNeighborBlockCoordinateTest, PracticalCriteriaDoNotQualifySweepsWithFailedBlocks)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(2,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t) {return block==0;};
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-recovery-exhausted");
    EXPECT_EQ(result.first_order_stationarity_sweep,0u);
    EXPECT_EQ(result.confirmed_stationarity_sweep,0u);
    EXPECT_FALSE(result.certification.has_value());
}
TEST(JointFixedNeighborBlockCoordinateTest, ProductionSearchDefaultsToValidatedPracticalCriteria)
{
    const auto production=n::PracticalJointConvergenceCriteria;
    EXPECT_DOUBLE_EQ(production.global_ac_kkt,n::PracticalJointConvergenceCriteria.global_ac_kkt);
    EXPECT_DOUBLE_EQ(production.global_width_gradient,n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_DOUBLE_EQ(production.confirmed_log_width_change,
        n::PracticalJointConvergenceCriteria.confirmed_log_width_change);
    const auto local=n::FixedNeighborPolicy{}.convergence_criteria;
    EXPECT_DOUBLE_EQ(local.global_ac_kkt,n::StrictJointConvergenceCriteria.global_ac_kkt);
    EXPECT_DOUBLE_EQ(local.global_width_gradient,n::StrictJointConvergenceCriteria.global_width_gradient);
    EXPECT_DOUBLE_EQ(local.confirmed_log_width_change,n::StrictJointConvergenceCriteria.confirmed_log_width_change);
}

TEST(JointFixedNeighborBlockCoordinateTest, BudgetExhaustionPreservesCandidateWithoutInventingTrustedState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    n::FixedNeighborSearchPolicy policy; policy.eta_upper_bound=std::log(2.2);
    policy.configure_search_for_testing=[](n::FixedNeighborPolicy & search) {search.maximum_sweeps=1;};
    const auto fit=n::FitFixedNeighborComponents(problem,std::vector<double>(8,.55),policy);
    ASSERT_EQ(fit.components.size(),1u);
    const auto & component=fit.components.front();
    EXPECT_FALSE(fit.search_completed);
    EXPECT_FALSE(component.search_completed);
    EXPECT_EQ(component.stop_reason,"block-sweep-budget");
    EXPECT_EQ(component.search_sweeps,1u);
    EXPECT_EQ(component.endpoint_certification_sweep,1u);
    ASSERT_TRUE(component.candidate_diagnostic);
    EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
    const bool qualified=component.candidate_diagnostic->numerical_qualification_status==rhbm_gem::JointCheckStatus::Passed;
    EXPECT_EQ(component.state.has_value(),qualified);
    EXPECT_EQ(fit.assembled_state.has_value(),qualified);
    EXPECT_EQ(fit.objective.has_value(),qualified);
    if(!qualified) EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    const auto trust=std::find_if(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
        return check.name=="independent-replay";
    });
    ASSERT_NE(trust,component.evidence.end());
    if(qualified)
        EXPECT_EQ(trust->status,rhbm_gem::JointCheckStatus::Passed)<<trust->reason;
    else
    {
        EXPECT_NE(trust->status,rhbm_gem::JointCheckStatus::Passed)<<trust->reason;
        EXPECT_EQ(trust->reason,"primary-bounded-width-kkt-failed");
    }

    const auto copy_result=rhbm_gem::core::CaptureJointAnalysisResult(fit);
    auto movable_fit=fit;
    const auto move_result=rhbm_gem::core::CaptureJointAnalysisResult(std::move(movable_fit));
    EXPECT_EQ(move_result.search_completed,copy_result.search_completed);
    EXPECT_FALSE(copy_result.search_completed);
    ASSERT_EQ(copy_result.components.size(),1u);
    ASSERT_EQ(move_result.components.size(),copy_result.components.size());
    EXPECT_EQ(move_result.components.front().stop_reason,copy_result.components.front().stop_reason);
    ASSERT_TRUE(copy_result.components.front().candidate_diagnostic);
    ASSERT_TRUE(move_result.components.front().candidate_diagnostic);
    EXPECT_EQ(move_result.components.front().candidate_diagnostic->candidate_available,
        copy_result.components.front().candidate_diagnostic->candidate_available);
    EXPECT_EQ(move_result.components.front().candidate_diagnostic->numerical_qualification_status,
        copy_result.components.front().candidate_diagnostic->numerical_qualification_status);
    if(!qualified) EXPECT_FALSE(copy_result.components.front().state);
}
