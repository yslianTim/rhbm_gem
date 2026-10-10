#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using rhbm_gem::core::JointProblem;
void CheckSingleBlock(const std::string & topology,int atoms)
{
    JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(atoms,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=30;
    std::vector<n::FixedNeighborBlockSweep> observed;
    n::BlockCoordinateState observed_state; std::size_t state_callbacks{};
    policy.sweep_observer=[&](const auto & sweep){observed.push_back(sweep);};
    policy.state_observer=[&](std::size_t,const auto & state,const auto &,const auto &) {
        observed_state=state; ++state_callbacks;
    };
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    auto unobserved_policy=policy; unobserved_policy.sweep_observer={}; unobserved_policy.state_observer={};
    const auto unobserved=n::SearchFixedNeighbor(problem,eta,unobserved_policy);
    ASSERT_TRUE(result.search_converged)<<result.reason;
    ASSERT_TRUE(result.endpoint_certified)<<result.endpoint_trust.reason;
    EXPECT_GE(result.sweeps.size(),2u);
    EXPECT_GT(result.first_order_stationarity_sweep,0u);
    EXPECT_GT(result.confirmed_stationarity_sweep,result.first_order_stationarity_sweep);
    EXPECT_EQ(result.confirmed_stationarity_sweep,result.sweeps.size());
    ASSERT_EQ(observed.size(),result.sweeps.size());
    EXPECT_EQ(state_callbacks,result.sweeps.size());
    EXPECT_TRUE((observed_state.eta.array()==result.state.eta.array()).all());
    EXPECT_TRUE((observed_state.beta.array()==result.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(observed.back().global_ac_kkt,result.sweeps.back().global_ac_kkt);
    EXPECT_LE(result.sweeps.back().global_ac_kkt,1e-10);
    EXPECT_LE(result.sweeps.back().global_width_gradient_inf_norm,1e-12);
    EXPECT_LE(result.sweeps.back().cache_replay_error,2e-12);
    EXPECT_TRUE(result.assessment.primary.valid);
    EXPECT_TRUE(result.assessment.reference.valid);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(result.fit.search_completed);
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
void CheckMultiBlockBudgetFailure(const std::string & topology)
{
    JointProblem problem(second_stage_test::OperatorWorkload(topology,32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto partition=n::BuildStructuralBlockPartition(*data.input,data.layout,16);
    ASSERT_EQ(partition.cores.size(),2u);
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=0;
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-sweep-budget");
    EXPECT_FALSE(result.endpoint_certified);
    EXPECT_FALSE(result.fit.search_completed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(result.fit.components.size(),1u);
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
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=1;
    const auto candidate=n::SearchFixedNeighbor(problem,eta,policy);
    EXPECT_FALSE(candidate.search_converged);
    EXPECT_EQ(candidate.reason,"block-sweep-budget");
    EXPECT_EQ(candidate.sweeps.size(),1u);
    EXPECT_EQ(candidate.first_order_stationarity_sweep,0u);
    EXPECT_EQ(candidate.confirmed_stationarity_sweep,0u);

    policy.maximum_sweeps=30;
    const auto confirmed=n::SearchFixedNeighbor(problem,eta,policy);
    ASSERT_TRUE(confirmed.search_converged)<<confirmed.reason;
    EXPECT_GT(confirmed.first_order_stationarity_sweep,0u);
    EXPECT_GT(confirmed.confirmed_stationarity_sweep,confirmed.first_order_stationarity_sweep);
    EXPECT_EQ(confirmed.confirmed_stationarity_sweep,confirmed.sweeps.size());
}
TEST(JointFixedNeighborBlockCoordinateTest, SearchOnlySkipsFinalAssessmentWithoutChangingTrajectory)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    const auto assessed=n::SearchFixedNeighbor(problem,eta);
    n::FixedNeighborPolicy policy; policy.assess_final_endpoint=false;
    const auto search_only=n::SearchFixedNeighbor(problem,eta,policy);
    ASSERT_TRUE(assessed.endpoint_certified);
    ASSERT_TRUE(search_only.search_converged)<<search_only.reason;
    EXPECT_FALSE(search_only.endpoint_certified);
    EXPECT_TRUE(search_only.fit.components.empty());
    EXPECT_EQ(search_only.reason,assessed.reason);
    ASSERT_EQ(search_only.sweeps.size(),assessed.sweeps.size());
    EXPECT_TRUE((search_only.state.eta.array()==assessed.state.eta.array()).all());
    EXPECT_TRUE((search_only.state.beta.array()==assessed.state.beta.array()).all());
    EXPECT_TRUE((search_only.state.prediction.array()==assessed.state.prediction.array()).all());
    EXPECT_TRUE((search_only.state.residual.array()==assessed.state.residual.array()).all());
    EXPECT_DOUBLE_EQ(search_only.state.objective,assessed.state.objective);
    for(std::size_t k=0;k<search_only.sweeps.size();++k)
    {
        EXPECT_EQ(search_only.sweeps[k].sweep,assessed.sweeps[k].sweep);
        EXPECT_DOUBLE_EQ(search_only.sweeps[k].objective_after,assessed.sweeps[k].objective_after);
        EXPECT_DOUBLE_EQ(search_only.sweeps[k].global_ac_kkt,assessed.sweeps[k].global_ac_kkt);
        EXPECT_DOUBLE_EQ(search_only.sweeps[k].global_width_gradient_inf_norm,
            assessed.sweeps[k].global_width_gradient_inf_norm);
    }
}
TEST(JointFixedNeighborBlockCoordinateTest, LocalTrajectoryTelemetryIsOptInAndDoesNotChangeSearchState)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy policy; policy.assess_final_endpoint=false; policy.maximum_sweeps=1;
    policy.capture_local_trajectory=true;
    const auto captured=n::SearchFixedNeighbor(problem,eta,policy);
    policy.capture_local_trajectory=false;
    const auto quiet=n::SearchFixedNeighbor(problem,eta,policy);
    ASSERT_FALSE(captured.blocks.empty());
    ASSERT_EQ(captured.blocks.size(),quiet.blocks.size());
    ASSERT_FALSE(captured.blocks.front().profile_trials.empty());
    for(std::size_t k=0;k<captured.blocks.size();++k)
    {
        const auto & block=captured.blocks[k];
        double factor_seconds{};
        std::size_t accepted_updates{};
        for(const auto & trial:block.profile_trials)
        {
            factor_seconds+=trial.factor_seconds;
            if(trial.accepted && trial.accepted_update && *trial.accepted_update>0) ++accepted_updates;
            EXPECT_GE(trial.cumulative_factor_seconds,factor_seconds);
        }
        EXPECT_DOUBLE_EQ(block.profile_factor_seconds,factor_seconds);
        EXPECT_EQ(accepted_updates,static_cast<std::size_t>(block.accepted_updates));
        EXPECT_TRUE(block.local_objective_before==block.local_objective_before);
        EXPECT_EQ(block.profile_evaluations,static_cast<int>(block.profile_trials.size()));
    }
    EXPECT_EQ(captured.reason,quiet.reason);
    EXPECT_EQ(captured.sweeps.size(),quiet.sweeps.size());
    EXPECT_TRUE((captured.state.eta.array()==quiet.state.eta.array()).all());
    EXPECT_TRUE((captured.state.beta.array()==quiet.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(captured.state.objective,quiet.state.objective);
}
TEST(JointFixedNeighborBlockCoordinateTest, AggregateWorkTelemetryIsOptInAndDoesNotChangeSearchState)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=2;
    policy.assess_final_endpoint=false; policy.collect_telemetry=true;
    const auto measured=n::SearchFixedNeighbor(problem,eta,policy);
    policy.collect_telemetry=false;
    const auto quiet=n::SearchFixedNeighbor(problem,eta,policy);

    ASSERT_GT(measured.work.block_solves,0u);
    std::size_t sweep_block_solves{};
    for(const auto & sweep:measured.sweeps) sweep_block_solves+=sweep.block_solves;
    EXPECT_EQ(measured.work.block_solves,sweep_block_solves);
    EXPECT_EQ(measured.work.old_core_basis_builds,measured.work.block_solves);
    EXPECT_EQ(measured.work.full_candidate_replays,measured.work.candidate_state_full_copies);
    EXPECT_GT(measured.work.affected_row_updates,0u);
    for(const double seconds:{measured.work.old_core_seconds,measured.work.effective_response_seconds,
        measured.work.local_search_seconds,measured.work.local_state_seconds,measured.work.candidate_copy_seconds,
        measured.work.candidate_replay_seconds,measured.work.cache_update_seconds,
        measured.work.sweep_replay_seconds,measured.work.sweep_global_state_seconds})
        EXPECT_GE(seconds,0.0);

    EXPECT_EQ(quiet.work.block_solves,0u);
    EXPECT_EQ(quiet.work.full_candidate_replays,0u);
    EXPECT_EQ(quiet.work.candidate_state_full_copies,0u);
    EXPECT_EQ(quiet.work.affected_row_updates,0u);
    EXPECT_EQ(quiet.work.old_core_basis_builds,0u);
    EXPECT_DOUBLE_EQ(quiet.work.old_core_seconds,0.0);
    EXPECT_DOUBLE_EQ(quiet.work.candidate_replay_seconds,0.0);
    EXPECT_DOUBLE_EQ(quiet.work.sweep_replay_seconds,0.0);
    for(const double seconds:{quiet.work.effective_response_seconds,quiet.work.local_search_seconds,
        quiet.work.local_state_seconds,quiet.work.candidate_copy_seconds,quiet.work.cache_update_seconds,
        quiet.work.sweep_global_state_seconds})
        EXPECT_DOUBLE_EQ(seconds,0.0);
    EXPECT_EQ(measured.reason,quiet.reason);
    EXPECT_EQ(measured.search_converged,quiet.search_converged);
    EXPECT_EQ(measured.sweeps.size(),quiet.sweeps.size());
    EXPECT_TRUE((measured.state.eta.array()==quiet.state.eta.array()).all());
    EXPECT_TRUE((measured.state.beta.array()==quiet.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(measured.state.objective,quiet.state.objective);
}
TEST(JointFixedNeighborBlockCoordinateTest, LocalProfileAttributionIsAggregateAndOptIn)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",16));
    const n::Vector eta=n::Vector::Constant(16,std::log(.55));
    n::FixedNeighborPolicy measured_policy; measured_policy.core_atoms=8; measured_policy.maximum_sweeps=1;
    measured_policy.assess_final_endpoint=false; measured_policy.collect_telemetry=true;
    const auto measured=n::SearchFixedNeighbor(problem,eta,measured_policy);
    auto quiet_policy=measured_policy; quiet_policy.collect_telemetry=false;
    const auto quiet=n::SearchFixedNeighbor(problem,eta,quiet_policy);
    const auto & profile=measured.work.local_profile_work;
    ASSERT_GT(profile.total.evaluations,0u);
    EXPECT_EQ(profile.initial_profile.evaluations+profile.trial_profile.evaluations,
        std::accumulate(measured.sweeps.begin(),measured.sweeps.end(),std::size_t{},
            [](std::size_t total,const auto & sweep){return total+sweep.profile_evaluations;}));
    EXPECT_EQ(profile.initial_profile.evaluations,static_cast<std::size_t>(measured.sweeps.front().block_solves));
    EXPECT_EQ(profile.total.evaluations,profile.initial_profile.evaluations+profile.trial_profile.evaluations+
        profile.accepted_endpoint.evaluations+profile.reference_evaluation.evaluations);
    for(const double seconds:{profile.total.evaluation_seconds,profile.total.basis_seconds,
        profile.total.linear_matrix_preparation_seconds,profile.total.linear_rhs_solve_seconds,
        profile.total.linear_certificate_seconds,profile.total.derivative_prepare_seconds,
        profile.total.derivative_raw_assembly_seconds,profile.total.derivative_free_design_assembly_seconds,
        profile.total.derivative_factor_match_seconds,profile.total.derivative_factor_build_seconds,
        profile.total.derivative_factor_compact_seconds,profile.total.derivative_rank_seconds,
        profile.total.derivative_least_squares_seconds,profile.total.derivative_normal_solve_seconds,
        profile.total.derivative_cancellation_check_seconds,profile.total.derivative_cancellation_fallback_seconds,
        profile.total.derivative_reduce_seconds,profile.total.derivative_rows_seconds,
        profile.total.derivative_jacobian_qr_seconds,profile.total.derivative_norms_seconds,
        profile.total.derivative_outer_overhead_seconds,profile.total.tiled_qr_assembly_copy_seconds,
        profile.total.tiled_qr_householder_seconds,profile.total.tiled_qr_rhs_transform_seconds,
        profile.total.replay_trust_seconds,
        profile.lm_overhead_seconds}) EXPECT_GE(seconds,0.0);
    EXPECT_EQ(quiet.work.local_profile_work.total.evaluations,0u);
    EXPECT_DOUBLE_EQ(quiet.work.local_profile_work.total.basis_seconds,0.0);
    EXPECT_DOUBLE_EQ(quiet.work.local_profile_work.lm_overhead_seconds,0.0);
    EXPECT_EQ(measured.reason,quiet.reason);
    EXPECT_TRUE((measured.state.eta.array()==quiet.state.eta.array()).all());
    EXPECT_TRUE((measured.state.beta.array()==quiet.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(measured.state.objective,quiet.state.objective);
}
TEST(JointFixedNeighborBlockCoordinateTest, ProductionSearchKeepsOnlyMinimalOutputWithoutDiagnostics)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(32,std::log(.55));
    const n::Domain domain=n::ProfileDomain(data.domain,data.layout);
    const auto context=n::ProfileContext(data.context,data.layout,data.domain.rows);
    n::FixedNeighborPolicy quiet_policy; quiet_policy.core_atoms=16; quiet_policy.maximum_sweeps=3;
    const auto quiet=n::SearchFixedNeighborComponent(*data.input,data.layout,domain,data.y,data.y,
        initial_eta,context,quiet_policy);
    auto diagnostic_policy=quiet_policy; diagnostic_policy.collect_diagnostics=true;
    const auto diagnostic=n::SearchFixedNeighborComponent(*data.input,data.layout,domain,data.y,data.y,
        initial_eta,context,diagnostic_policy);

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
    EXPECT_EQ(quiet.work.local_profile_work.total.evaluations,0u);
    EXPECT_EQ(quiet.reason,diagnostic.reason);
    EXPECT_EQ(quiet.search_converged,diagnostic.search_converged);
    EXPECT_EQ(quiet.endpoint_certified,diagnostic.endpoint_certified);
    EXPECT_TRUE((quiet.state.eta.array()==diagnostic.state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==diagnostic.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(quiet.state.objective,diagnostic.state.objective);
    EXPECT_DOUBLE_EQ(quiet.final_sweep->global_ac_kkt,diagnostic.sweeps.back().global_ac_kkt);
    EXPECT_DOUBLE_EQ(quiet.final_sweep->global_width_gradient_inf_norm,
        diagnostic.sweeps.back().global_width_gradient_inf_norm);
    EXPECT_TRUE((quiet.assessment.primary.eta.array()==diagnostic.assessment.primary.eta.array()).all());
}
TEST(JointFixedNeighborBlockCoordinateTest, FixedNeighborAcceptsAtMostOneTrustedLocalUpdatePerBlockVisit)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=2;
    policy.assess_final_endpoint=false; policy.capture_local_trajectory=true;
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
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
TEST(JointFixedNeighborBlockCoordinateTest, EtaConfirmationUsesInclusiveExistingThreshold)
{
    constexpr double threshold=1e-10;
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(threshold,false));
    EXPECT_TRUE(n::IsFixedNeighborEtaChangeConfirmed(threshold,true));
    EXPECT_TRUE(n::IsFixedNeighborEtaChangeConfirmed(std::nextafter(threshold,0.0),true));
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(std::nextafter(threshold,
        std::numeric_limits<double>::infinity()),true));
    EXPECT_FALSE(n::IsFixedNeighborEtaChangeConfirmed(std::numeric_limits<double>::quiet_NaN(),true));
}
TEST(JointFixedNeighborBlockCoordinateTest, SweepCoordinateTelemetryUsesCompleteSweepStates)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",32));
    const n::Vector initial_eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=5;
    policy.stop_after_stationarity=false;
    n::Vector previous_eta=initial_eta,previous_beta=n::Vector::Zero(64);
    std::size_t observed_sweeps{};
    policy.state_observer=[&](std::size_t sweep,const n::BlockCoordinateState & state,
        const n::FixedNeighborBlockSweep & telemetry,const std::vector<n::FixedNeighborBlockRecord> & blocks) {
        ++observed_sweeps;
        EXPECT_EQ(telemetry.coordinate_confirmation_available,sweep>1);
        EXPECT_DOUBLE_EQ(telemetry.eta_change_inf,(state.eta-previous_eta).lpNorm<Eigen::Infinity>());
        const double beta_change=((state.beta-previous_beta).array().abs() /
            (1.0+state.beta.array().abs().max(previous_beta.array().abs()))).maxCoeff();
        EXPECT_DOUBLE_EQ(telemetry.beta_scaled_change,beta_change);
        std::size_t accepted{},unchanged{};
        for(const auto & block:blocks) if(block.sweep==sweep)
        {
            accepted+=block.accepted ? 1u : 0u;
            unchanged+=block.status=="unchanged" ? 1u : 0u;
        }
        EXPECT_EQ(telemetry.accepted_blocks,accepted);
        EXPECT_EQ(telemetry.unchanged_blocks,unchanged);
        previous_eta=state.eta; previous_beta=state.beta;
    };
    const auto observed=n::SearchFixedNeighbor(problem,initial_eta,policy);
    auto quiet_policy=policy; quiet_policy.state_observer={};
    const auto quiet=n::SearchFixedNeighbor(problem,initial_eta,quiet_policy);
    ASSERT_EQ(observed_sweeps,observed.sweeps.size());
    ASSERT_EQ(observed.sweeps.size(),quiet.sweeps.size());
    EXPECT_EQ(observed.search_converged,quiet.search_converged);
    EXPECT_EQ(observed.reason,quiet.reason);
    EXPECT_TRUE((observed.state.eta.array()==quiet.state.eta.array()).all());
    EXPECT_TRUE((observed.state.beta.array()==quiet.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(observed.state.objective,quiet.state.objective);
    for(std::size_t k=0;k<observed.sweeps.size();++k)
    {
        EXPECT_DOUBLE_EQ(observed.sweeps[k].eta_change_inf,quiet.sweeps[k].eta_change_inf);
        EXPECT_DOUBLE_EQ(observed.sweeps[k].beta_scaled_change,quiet.sweeps[k].beta_scaled_change);
        EXPECT_EQ(observed.sweeps[k].accepted_blocks,quiet.sweeps[k].accepted_blocks);
        EXPECT_EQ(observed.sweeps[k].unchanged_blocks,quiet.sweeps[k].unchanged_blocks);
    }
}
TEST(JointFixedNeighborBlockCoordinateTest, InvalidCorePolicyKeepsFailureSemantics)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55)); n::FixedNeighborPolicy policy; policy.core_atoms=0;
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    EXPECT_FALSE(result.search_converged); EXPECT_EQ(result.reason,"block-invalid-partition");
    EXPECT_FALSE(result.endpoint_certified); EXPECT_EQ(result.sweeps.size(),0u);
}
TEST(JointFixedNeighborBlockCoordinateTest, SameEtaRawPrimaryAndReferenceDiagnosticsPreserveSearchState)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(8,std::log(.55));
    const auto search=n::SearchFixedNeighbor(problem,initial_eta);
    ASSERT_TRUE(search.search_converged);
    const n::Vector eta=search.state.eta,beta=search.state.beta;
    const n::Domain domain=n::ProfileDomain(data.domain,data.layout);
    const auto context=n::ProfileContext(data.context,data.layout,data.domain.rows);
    const n::Vector y=data.y;
    const auto raw=n::EvaluateState(domain,y,eta,beta,context);
    const auto primary=n::EvaluateProfile(domain,y,eta,false,&context);
    const auto reference=n::EvaluateProfile(domain,y,eta,true,&context);
    const auto raw_assessment=n::AssessEvaluated(domain,y,raw,reference,context,true);
    const auto primary_assessment=n::AssessEvaluated(domain,y,primary,reference,context,false);
    const auto reference_assessment=n::AssessEvaluated(domain,y,reference,primary,context,false);
    const auto raw_trust=n::CheckTrust(domain,y,raw,context,reference);
    const auto primary_trust=n::CheckTrust(domain,y,primary,context,reference);
    const auto reference_trust=n::CheckTrust(domain,y,reference,context,primary);

    EXPECT_TRUE((search.state.eta.array()==eta.array()).all());
    EXPECT_TRUE((raw.eta.array()==eta.array()).all());
    EXPECT_TRUE((primary.eta.array()==eta.array()).all());
    EXPECT_TRUE((reference.eta.array()==eta.array()).all());
    EXPECT_TRUE(raw.valid); EXPECT_TRUE(primary.valid); EXPECT_TRUE(reference.valid);
    EXPECT_TRUE(std::isfinite(raw.gradient.lpNorm<Eigen::Infinity>()));
    EXPECT_TRUE(std::isfinite(primary.gradient.lpNorm<Eigen::Infinity>()));
    EXPECT_TRUE(std::isfinite(reference.gradient.lpNorm<Eigen::Infinity>()));
    const auto difference=[](const n::Vector & lhs,const n::Vector & rhs) {
        return ((lhs-rhs).array().abs()/(1+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
    };
    EXPECT_TRUE(std::isfinite(difference(beta,primary.beta)));
    EXPECT_TRUE(std::isfinite(difference(primary.beta,reference.beta)));
    EXPECT_TRUE(std::isfinite(difference(beta,reference.beta)));
    for(const auto * assessment:{&raw_assessment,&primary_assessment,&reference_assessment})
    {
        EXPECT_TRUE(assessment->primary.valid);
        EXPECT_TRUE(assessment->reference.valid);
        EXPECT_TRUE(assessment->design.has_value());
        EXPECT_TRUE(assessment->widths.has_value());
        EXPECT_TRUE(assessment->jacobian.has_value());
        EXPECT_TRUE(assessment->normalized_widths.has_value());
        EXPECT_EQ(assessment->gradient,assessment->primary.gradient.lpNorm<Eigen::Infinity>()<=1e-12 &&
            assessment->reference.gradient.lpNorm<Eigen::Infinity>()<=1e-12);
    }
    EXPECT_TRUE(raw_trust.primary_valid); EXPECT_TRUE(primary_trust.primary_valid); EXPECT_TRUE(reference_trust.primary_valid);
}
