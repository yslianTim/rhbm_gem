#include <gtest/gtest.h>
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <fstream>
#include "support/JointSyntheticWorkload.hpp"
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using rhbm_gem::core::JointProblem;
std::pair<JointProblem,std::vector<double>> FixedNeighborRetryFixture()
{
    const std::string root="joint_component/fixed-neighbor-retry/";
    auto model=rhbm_gem::ReadModel(command_test::TestDataPath(root+"model.cif"));
    auto map=rhbm_gem::ReadMap(command_test::TestDataPath(root+"map.map"));
    model->ApplySimulationMetadata(.5); model->SelectAllAtoms();
    model->ApplyElementSelection(Element::HYDROGEN,true);
    JointProblem problem=rhbm_gem::core::BuildJointProblem(*map,*model);
    std::ifstream file(command_test::TestDataPath(root+"initial_widths.txt"));
    std::vector<double> widths; double width{};
    while(file>>width) widths.push_back(width);
    return {std::move(problem),std::move(widths)};
}
void CheckSingleBlock(const std::string & topology,int atoms)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,atoms));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(atoms,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=30;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
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
    const auto production=rhbm_gem::core::FitJointComponents(problem,std::vector<double>(static_cast<std::size_t>(atoms),.55));
    ASSERT_EQ(production.components.size(),1u);
    ASSERT_TRUE(production.components.front().candidate_diagnostic);
    EXPECT_EQ(production.components.front().candidate_diagnostic->convergence_criteria_id,"practical-v1");
    ASSERT_TRUE(production.components.front().state);
    ASSERT_TRUE(result.fit.components.front().state);
    ASSERT_TRUE(production.objective && result.fit.objective);
    EXPECT_EQ(production.components.front().search_sweeps,result.fit.components.front().search_sweeps);
    EXPECT_EQ(production.components.front().stop_reason,"block-stationary");
    EXPECT_LE(std::abs(production.objective.value()-result.fit.objective.value()),1e-12);
    for(std::size_t k=0;k<production.components.front().state->ac.size();++k)
        EXPECT_LE(std::abs(production.components.front().state->ac[k]-result.fit.components.front().state->ac[k]),
            1e-10*(1+std::abs(result.fit.components.front().state->ac[k])));
    for(std::size_t k=0;k<production.components.front().state->log_b.size();++k)
        EXPECT_LE(std::abs(production.components.front().state->log_b[k]-result.fit.components.front().state->log_b[k]),1e-10);
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
TEST(JointFixedNeighborBlockCoordinateTest, IncrementalCandidateDeltaMatchesAffectedRowFullReplay)
{
    n::Vector residual(4); residual<<1.25,-.75,2e-8,-2.0;
    const n::Indices affected_rows{0,2};
    n::Vector old_core(2),new_core(2); old_core<<.5,-.3; new_core<<.25,-.30000001;
    const double scale=.75;
    const double current=.5*residual.squaredNorm()/(scale*scale);
    const auto incremental=n::EvaluateCandidateDelta(residual,affected_rows,old_core,new_core,current,scale);
    n::Vector full_residual=residual;
    for(std::size_t k=0;k<affected_rows.size();++k)
        full_residual(static_cast<Eigen::Index>(affected_rows[k]))+=new_core(static_cast<Eigen::Index>(k))-
            old_core(static_cast<Eigen::Index>(k));
    const double full_objective=.5*full_residual.squaredNorm()/(scale*scale);
    ASSERT_EQ(incremental.residual_delta.size(),2);
    EXPECT_DOUBLE_EQ(incremental.residual_delta(0),-.25);
    EXPECT_NEAR(incremental.residual_delta(1),-1e-8,1e-16);
    EXPECT_NEAR(incremental.objective,full_objective,2e-15);
    EXPECT_NEAR(incremental.objective_delta,full_objective-current,2e-15);
    EXPECT_TRUE(std::isfinite(incremental.delta_condition));
}
TEST(JointFixedNeighborBlockCoordinateTest, IncrementalCandidateExperimentMatchesFullReplayOnChainAndCube)
{
    for(const auto & topology:{"chain","cube"})
    {
        JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,8));
        const n::Vector eta=n::Vector::Constant(8,std::log(.55));
        n::FixedNeighborPolicy policy; policy.maximum_sweeps=1; policy.certify_endpoint=false;
        policy.collect_telemetry=true; policy.retain_diagnostic_records=false;
        policy.compare_incremental_replay=true;
        const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
        EXPECT_EQ(result.sweep_count,1u)<<topology;
        EXPECT_EQ(result.work.incremental_candidate_replays,result.work.full_candidate_replays)<<topology;
        EXPECT_EQ(result.work.incremental_acceptance_mismatches,0u)<<topology;
        EXPECT_LE(result.work.maximum_incremental_objective_difference,1e-10)<<topology;
        EXPECT_LE(result.work.maximum_incremental_residual_difference,1e-12)<<topology;
        ASSERT_TRUE(result.final_sweep.has_value())<<topology;
        EXPECT_LE(result.final_sweep->cache_replay_error,2e-12)<<topology;
        EXPECT_TRUE(result.blocks.empty())<<topology;
        EXPECT_TRUE(result.sweeps.empty())<<topology;
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, ResearchBoundedSearchReleasesTheStartingUpperFace)
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
    policy.eta_upper_bound_for_research=eta_max;
    const auto result=n::RunFixedNeighborDiagnostic(problem,initial,policy);
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
}

void CheckMultiBlockBudgetFailure(const std::string & topology)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto partition=n::BuildStructuralBlockPartition(*data.input,data.layout,16);
    ASSERT_EQ(partition.cores.size(),2u);
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=0;
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
TEST(JointFixedNeighborBlockCoordinateTest, DiagnosticFacadeCanSkipCertificationWithoutChangingCandidate)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    const auto assessed=n::RunFixedNeighborDiagnostic(problem,eta);
    n::FixedNeighborPolicy policy; policy.certify_endpoint=false;
    const auto search_only=n::RunFixedNeighborDiagnostic(problem,eta,policy);
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
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy policy; policy.certify_endpoint=false; policy.maximum_sweeps=1;
    policy.capture_local_trajectory=true;
    const auto captured=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    policy.capture_local_trajectory=false;
    const auto quiet=n::RunFixedNeighborDiagnostic(problem,eta,policy);
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

TEST(JointFixedNeighborBlockCoordinateTest, WidthTrajectoryCaptureHonorsAtomEventAndRowLimits)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    ASSERT_FALSE(data.input->atom_ids.empty());
    n::FixedNeighborWidthTrajectoryCapture capture;
    capture.atom_initial_widths.emplace_back(data.input->atom_ids.front(),.55);
    capture.maximum_events=1; capture.maximum_atoms=1; capture.maximum_rows=0;
    n::FixedNeighborDiagnosticSink sink; std::size_t observed_events{};
    std::size_t retained_blocks{},retained_sweeps{};
    sink.block=[&](const n::FixedNeighborBlockRecord & block) {
        for(const auto & trial:block.profile_trials) if(!trial.tracked_atoms.empty())
        {
            ++observed_events;
            EXPECT_EQ(trial.tracked_atoms.size(),1u);
            EXPECT_EQ(trial.tracked_atoms.front().atom_id,data.input->atom_ids.front());
            EXPECT_EQ(trial.tracked_atoms.front().local_design_rank,-1);
            EXPECT_EQ(trial.width_trust_evidence_unavailable_reason,"row-limit");
        }
    };
    n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=1;
    policy.configure_search_for_testing=[&](n::FixedNeighborPolicy & local) {
        local.width_trajectory=&capture; local.diagnostic_sink=&sink;
    };
    policy.search_result_observer=[&](const n::FixedNeighborSearchResult & result) {
        retained_blocks=result.blocks.size(); retained_sweeps=result.sweeps.size();
    };
    const auto traced=n::FitFixedNeighborComponents(problem,std::vector<double>(8,.55),policy);
    n::FixedNeighborSearchPolicy quiet_policy; quiet_policy.maximum_sweeps=1;
    const auto quiet=n::FitFixedNeighborComponents(problem,std::vector<double>(8,.55),quiet_policy);
    EXPECT_EQ(capture.events_emitted,1u);
    EXPECT_EQ(observed_events,1u);
    EXPECT_EQ(retained_blocks,0u);
    EXPECT_EQ(retained_sweeps,0u);
    ASSERT_EQ(traced.components.size(),1u);
    ASSERT_EQ(quiet.components.size(),1u);
    EXPECT_EQ(traced.components.front().stop_reason,quiet.components.front().stop_reason);
    EXPECT_EQ(traced.components.front().search_sweeps,quiet.components.front().search_sweeps);
    EXPECT_EQ(traced.components.front().state.has_value(),quiet.components.front().state.has_value());
    EXPECT_EQ(traced.assembled_state.has_value(),quiet.assembled_state.has_value());
    if(traced.components.front().state && quiet.components.front().state)
    {
        EXPECT_EQ(traced.components.front().state->log_b,quiet.components.front().state->log_b);
        EXPECT_EQ(traced.components.front().state->ac,quiet.components.front().state->ac);
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, WidthTrajectoryRecordsLocalIdentifiabilityAndCancellation)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    n::FixedNeighborWidthTrajectoryCapture capture;
    capture.atom_initial_widths.emplace_back(data.input->atom_ids.front(),.55);
    capture.maximum_events=8; capture.maximum_atoms=1; capture.maximum_rows=4096;
    n::FixedNeighborDiagnosticSink sink; std::size_t observed{};
    sink.block=[&](const n::FixedNeighborBlockRecord & block) {
        for(const auto & trial:block.profile_trials) if(!trial.tracked_atoms.empty())
        {
            ++observed;
            ASSERT_EQ(trial.tracked_atoms.size(),1u);
            const auto & atom=trial.tracked_atoms.front();
            EXPECT_EQ(atom.atom_id,data.input->atom_ids.front());
            EXPECT_TRUE(std::isfinite(atom.proposed_eta_step));
            EXPECT_TRUE(std::isfinite(atom.accepted_eta_step));
            EXPECT_DOUBLE_EQ(atom.accepted_eta_step,trial.accepted ? atom.proposed_eta_step : 0.0);
            EXPECT_TRUE(std::isfinite(atom.ac_cancellation_ratio));
            EXPECT_GE(atom.ac_cancellation_ratio,1.0-1e-12);
            EXPECT_FALSE(trial.projected_width_singular_values.empty());
            for(const auto singular:trial.projected_width_singular_values)
                EXPECT_TRUE(std::isfinite(singular));
        }
    };
    n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=1;
    policy.configure_search_for_testing=[&](n::FixedNeighborPolicy & local) {
        local.width_trajectory=&capture; local.diagnostic_sink=&sink;
    };
    const auto traced=n::FitFixedNeighborComponents(problem,std::vector<double>(2,.55),policy);
    EXPECT_EQ(observed,capture.events_emitted);
    EXPECT_GT(observed,0u);
}

TEST(JointFixedNeighborBlockCoordinateTest, FixedWidthOracleOptimizesAcWithoutChangingTheRestrictedModel)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",1));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(1,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=4; policy.certify_endpoint=false;
    policy.fixed_width_atom_ids=data.input->atom_ids;
    const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    EXPECT_TRUE(result.search_converged);
    EXPECT_FALSE(result.endpoint_certified);
    EXPECT_TRUE((result.state.eta.array()==initial_eta.array()).all());
    ASSERT_TRUE(result.final_sweep.has_value());
    EXPECT_LE(result.final_sweep->global_ac_kkt,1e-10);
    EXPECT_GT(result.final_sweep->global_width_gradient_inf_norm,1e-12);
}

TEST(JointFixedNeighborBlockCoordinateTest, WeakWidthSafeguardRequiresLargeWeakAndReplayScaleStep)
{
    n::Matrix jacobian=n::Matrix::Zero(2,2); jacobian(0,0)=1; jacobian(1,1)=1e-16;
    n::Vector step(2); step<<0,12;
    const n::RankPolicy rank{100,2,2};
    const auto assessment=n::AssessWeakWidthStep(jacobian,step,rank);
    ASSERT_TRUE(assessment.available);
    EXPECT_NEAR(assessment.weak_step_fraction,1.0,1e-14);
    EXPECT_TRUE(n::RejectWeakWidthTrial(12,assessment.weak_step_fraction,1.0,1.0));
    EXPECT_FALSE(n::RejectWeakWidthTrial(12,assessment.weak_step_fraction,1.0,0.99));
    EXPECT_FALSE(n::RejectWeakWidthTrial(2,assessment.weak_step_fraction,1.0,1.0));
    EXPECT_FALSE(n::RejectWeakWidthTrial(12,0.0,1.0,1.0));
    EXPECT_FALSE(n::RejectWeakWidthTrial(12,assessment.weak_step_fraction,1.0,1.01));
}

TEST(JointFixedNeighborBlockCoordinateTest, FailedBlockRetriesAfterLaterBlocksUpdateTheNeighborState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(2,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
    policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t sweep) {
        return block==0 && sweep==0;
    };
    std::optional<n::BlockCoordinateState> first_sweep_state;
    policy.state_observer=[&](std::size_t sweep,const auto & state,const auto &,const auto &) {
        if(sweep==1) first_sweep_state=state;
    };
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    ASSERT_TRUE(first_sweep_state.has_value());
    EXPECT_GT(first_sweep_state->beta.norm(),0.0);
    EXPECT_EQ(result.recovery_attempts,1u);
    EXPECT_EQ(result.recovery_successes,1u);
    const auto first_failure=std::find_if(result.blocks.begin(),result.blocks.end(),[](const auto & block) {
        return block.block==1 && block.sweep==1;
    });
    ASSERT_NE(first_failure,result.blocks.end());
    EXPECT_EQ(first_failure->status,"retry-pending");
    const auto retried=std::find_if(result.blocks.begin(),result.blocks.end(),[](const auto & block) {
        return block.block==1 && block.sweep==2;
    });
    ASSERT_NE(retried,result.blocks.end());
    EXPECT_TRUE(retried->recovery_retry);
    EXPECT_EQ(retried->recovery_attempt,1u);
    EXPECT_DOUBLE_EQ(retried->objective_before,first_sweep_state->objective);
}

TEST(JointFixedNeighborBlockCoordinateTest, ExhaustedBlockRecoveryLeavesTheLastAcceptedComponentStateUntouched)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(2,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
    policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t) {return block==0;};
    std::optional<n::BlockCoordinateState> last_complete_sweep_state;
    policy.state_observer=[&](std::size_t sweep,const auto & state,const auto &,const auto &) {
        if(sweep==1) last_complete_sweep_state=state;
    };
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    ASSERT_TRUE(last_complete_sweep_state.has_value());
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-recovery-exhausted");
    EXPECT_EQ(result.recovery_attempts,1u);
    EXPECT_EQ(result.recovery_successes,0u);
    EXPECT_EQ(result.sweep_count,1u);
    EXPECT_EQ(result.first_order_stationarity_sweep,0u);
    EXPECT_TRUE((result.state.eta.array()==last_complete_sweep_state->eta.array()).all());
    EXPECT_TRUE((result.state.beta.array()==last_complete_sweep_state->beta.array()).all());
    EXPECT_FALSE(result.certification.has_value());
}

TEST(JointFixedNeighborBlockCoordinateTest, FailedBlockReportsLocalDesignAndHaloSupportDiagnostics)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"z-target","a-weak-halo"}; input.row_ids={"r0","r1"};
    input.observations={1.0,2.0}; input.support.resize(2);
    input.support[0]={{0,.1},{1,.2}}; input.support[1]={{0,.8},{1,.8}};
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};
    JointProblem problem(std::move(input));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto eta=n::Vector::Constant(2,std::log(.6));
    const auto context=data.context;
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=1; policy.collect_diagnostics=true;
    policy.component_index=4; policy.component_id="weak-component";
    policy.halo_census=[](std::size_t & count,std::size_t & minimum) {
        count=1; minimum=2;
    };
    const auto diagnostic=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,context,policy);
    auto quiet_policy=policy; quiet_policy.collect_diagnostics=false;
    const auto quiet=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,context,quiet_policy);

    EXPECT_EQ(diagnostic.reason,"block-recovery-sweep-budget");
    EXPECT_EQ(quiet.reason,diagnostic.reason);
    EXPECT_TRUE((quiet.state.eta.array()==diagnostic.state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==diagnostic.state.beta.array()).all());
    ASSERT_FALSE(diagnostic.blocks.empty());
    ASSERT_TRUE(diagnostic.blocks.front().failure_diagnostics);
    EXPECT_EQ(diagnostic.blocks.front().status,"retry-pending");
    EXPECT_EQ(diagnostic.blocks.front().reason,"block-retry-pending");
    const auto & failure=*diagnostic.blocks.front().failure_diagnostics;
    EXPECT_EQ(failure.component_index,4u);
    EXPECT_EQ(failure.component_id,"weak-component");
    EXPECT_EQ(failure.sweep_index,1u);
    EXPECT_EQ(failure.core_atom_ids,(std::vector<std::string>{"a-weak-halo"}));
    EXPECT_EQ(failure.core_atom_count,1u);
    EXPECT_EQ(failure.affected_row_count,2u);
    EXPECT_EQ(failure.halo_count,1u);
    EXPECT_EQ(failure.minimum_halo_support_rows,2u);
    EXPECT_FALSE(failure.initial_profile_valid);
    EXPECT_FALSE(failure.initial_profile_trusted);
    EXPECT_EQ(failure.local_design_rows,2);
    EXPECT_EQ(failure.local_design_columns,2);
    ASSERT_TRUE(failure.local_design_rank);
    EXPECT_EQ(*failure.local_design_rank,1);
    EXPECT_EQ(failure.local_search_stop_reason,"inner-invalid-input");
}

TEST(JointFixedNeighborBlockCoordinateTest, PreparedRowMappingsScaleWithMembershipsRatherThanRowsTimesBlocks)
{
    constexpr int atoms=192;
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",atoms));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    n::FixedNeighborPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=0;
    const n::Vector eta=n::Vector::Constant(atoms,std::log(.55));
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    const auto blocks=n::BuildStructuralBlockPartition(*data.input,data.layout,policy.core_atoms);
    ASSERT_GT(blocks.cores.size(),1u);
    EXPECT_EQ(result.mapping_preparations,blocks.cores.size());
    EXPECT_EQ(result.memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_EQ(result.global_state_memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_EQ(result.global_state_memory_estimate.stage,n::MemoryStage::GlobalStateEvaluation);
    ASSERT_GT(result.prepared_mapping_bytes,0u);
    std::size_t legacy_bytes{};
    ASSERT_TRUE(n::CheckedMultiplyBytes(static_cast<std::size_t>(data.domain.rows),blocks.cores.size(),legacy_bytes));
    ASSERT_TRUE(n::CheckedMultiplyBytes(legacy_bytes,sizeof(Eigen::Index),legacy_bytes));
    EXPECT_LT(result.prepared_mapping_bytes,legacy_bytes);
    EXPECT_LT(result.prepared_mapping_entries,
        static_cast<std::size_t>(data.domain.rows)*blocks.cores.size());
    EXPECT_EQ(result.reason,"block-sweep-budget");
}
TEST(JointFixedNeighborBlockCoordinateTest, ProductionScaleSingleComponentPreparationStaysBelowBudget)
{
    constexpr Eigen::Index rows=5'938'087;
    constexpr std::size_t atoms=37'406;
    constexpr std::size_t target_atoms=27'406;
    constexpr std::size_t halo_atoms=atoms-target_atoms;
    auto input=std::make_shared<rhbm_gem::core::JointProblemInput>();
    input->observations.assign(static_cast<std::size_t>(rows),0.0);
    input->row_ids.resize(static_cast<std::size_t>(rows));
    input->atom_ids.reserve(atoms); input->support.resize(atoms);
    for(std::size_t atom=0;atom<target_atoms;++atom)
    {
        input->atom_ids.push_back("target-"+std::to_string(atom));
        input->support[atom].push_back({atom,0.0});
        input->support[atom].push_back({atom+1,0.0});
    }
    for(std::size_t halo=0;halo<halo_atoms;++halo)
    {
        const auto atom=target_atoms+halo;
        input->atom_ids.push_back("halo-"+std::to_string(halo));
        input->support[atom].push_back({0,0.0});
        input->support[atom].push_back({target_atoms,0.0});
    }
    input->selection_domain.emplace();
    input->selection_domain->contributor_policy=std::string(rhbm_gem::FixedHaloContributorPolicy);
    input->selection_domain->fixed_halo_selection=rhbm_gem::JointFixedHaloSelectionProvenance{
        1,halo_atoms,halo_atoms,0};
    input->selection_domain->target_indices.resize(target_atoms);
    std::iota(input->selection_domain->target_indices.begin(),input->selection_domain->target_indices.end(),std::size_t{0});
    const std::shared_ptr<const rhbm_gem::core::JointProblemInput> snapshot=input;
    n::Domain all_contributors(snapshot); all_contributors.rows=rows;
    n::Indices target_indices(target_atoms);
    std::iota(target_indices.begin(),target_indices.end(),Eigen::Index{0});
    const auto domain=all_contributors.Select(target_indices,rows,n::CompactRowMapping::Identity(rows));
    rhbm_gem::JointParameterLayout layout; layout.full_atoms.resize(target_atoms);
    std::iota(layout.full_atoms.begin(),layout.full_atoms.end(),std::size_t{0});
    layout.informative_rows.resize(static_cast<std::size_t>(rows));
    std::iota(layout.informative_rows.begin(),layout.informative_rows.end(),Eigen::Index{0});
    n::EvaluationContext context; context.scale=1;
    context.rank={rows,2*static_cast<Eigen::Index>(target_atoms),static_cast<Eigen::Index>(target_atoms)};
    context.linear.rank_relative=context.rank.Relative(context.rank.design_columns);
    auto effective=std::make_shared<std::vector<double>>(input->observations);
    const std::shared_ptr<const std::vector<double>> effective_snapshot=effective;
    context.observations=n::Observe(effective_snapshot);
    context.atom_ids=n::Identities(std::shared_ptr<const std::vector<std::string>>(snapshot,&snapshot->atom_ids))
        .Select(target_indices);
    context.row_ids=n::Identities(std::shared_ptr<const std::vector<std::string>>(snapshot,&snapshot->row_ids));
    n::VectorMap y(effective->data(),rows);
    const n::Vector eta=n::Vector::Constant(static_cast<Eigen::Index>(target_atoms),std::log(.55));
    const Eigen::VectorXd fixed_prediction=Eigen::VectorXd::Zero(rows);
    const auto target_partition=n::Partition(domain,context.atom_ids);
    ASSERT_EQ(target_partition.components.size(),1u);
    EXPECT_EQ(target_partition.components.front().atoms.size(),target_atoms);
    n::FixedNeighborPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=0; policy.certify_endpoint=false;
    const auto preparation=n::SearchFixedNeighborComponent(domain,layout,y,eta,context,policy);
    EXPECT_EQ(preparation.reason,"block-sweep-budget");
    EXPECT_EQ(preparation.memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_EQ(preparation.global_state_memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_LE(preparation.global_state_memory_estimate.estimated_peak_bytes,n::JointMemoryPlanningBudgetBytes);
    EXPECT_EQ(preparation.sweep_count,0u);
    EXPECT_EQ(preparation.memory_estimate.shape.component_free_columns,2*target_atoms);
    EXPECT_EQ(preparation.memory_estimate.shape.component_blocks,
        (target_atoms+policy.core_atoms-1)/policy.core_atoms);

    const auto support_entries=n::CountSupportEntries(*snapshot);
    ASSERT_TRUE(support_entries);
    EXPECT_EQ(*support_entries,2*atoms);
    std::size_t target_support_entries{};
    for(std::size_t atom=0;atom<target_atoms;++atom) target_support_entries+=input->support[atom].size();
    EXPECT_EQ(target_support_entries,2*target_atoms);
    n::BlockCoordinateState accepted; accepted.eta=eta; accepted.beta=n::Vector::Zero(2*target_atoms);
    EXPECT_EQ(static_cast<std::size_t>(accepted.eta.size()+accepted.beta.size()),3*target_atoms);
    n::AssessmentWorkForTesting()={}; n::NumericsWorkForTesting()={};
    const auto certification=n::CertifyFixedNeighborEndpoint(domain,y,accepted,context);
    ASSERT_TRUE(certification.qualified())<<certification.assessment.failure;
    EXPECT_TRUE(certification.assessment.numerical_qualified);
    EXPECT_TRUE(certification.assessment.endpoint_verified);
    const auto & stages=n::AssessmentWorkForTesting().stages;
    EXPECT_EQ(std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="global-rank-assessment" || stage.name=="projected-width-rank-assessment";
    }),stages.end());
    const auto evidence=n::AssessmentEvidence(certification.assessment,rhbm_gem::JointEvidenceScope::ComponentLocal);
    EXPECT_EQ(n::ConvergenceStatus(evidence,rhbm_gem::JointEvidenceScope::ComponentLocal),
        rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,0);
    const auto rss=n::CurrentRssBytes();
    const auto peak=n::ProcessPeakRssBytes();
    RecordProperty("synthetic_support_entries",std::to_string(*support_entries));
    RecordProperty("target_count",std::to_string(target_atoms));
    RecordProperty("fixed_halo_count",std::to_string(halo_atoms));
    RecordProperty("target_only_blocks",std::to_string(preparation.memory_estimate.shape.component_blocks));
    RecordProperty("estimated_endpoint_verification_peak_bytes",
        std::to_string(certification.verification_memory_estimate.estimated_peak_bytes));
    std::cout<<"joint-large-synthetic observations="<<rows<<" contributors="<<atoms
        <<" targets="<<target_atoms<<" fixed_halos="<<halo_atoms
        <<" target_components="<<target_partition.components.size()
        <<" target_only_blocks="<<preparation.memory_estimate.shape.component_blocks
        <<" ac_design_columns="<<2*target_atoms
        <<" target_parameters="<<3*target_atoms
        <<" fixed_halo_prediction_bytes="<<static_cast<std::size_t>(fixed_prediction.size())*sizeof(double)
        <<" effective_observation_bytes="<<effective->size()*sizeof(double)
        <<" support_entries="<<*support_entries<<" target_support_entries="<<target_support_entries
        <<" preparation_peak_estimate_bytes="<<preparation.global_state_memory_estimate.estimated_peak_bytes
        <<" preflight_status="<<(preparation.global_state_memory_estimate.status==n::MemoryStatus::Pass ? "pass" : "fail")
        <<" preflight_budget_bytes="<<n::JointMemoryPlanningBudgetBytes
        <<" endpoint_peak_estimate_bytes="<<certification.verification_memory_estimate.estimated_peak_bytes
        <<" observed_current_rss_bytes="<<(rss ? std::to_string(*rss) : "unavailable")
        <<" observed_process_peak_rss_bytes="<<(peak ? std::to_string(*peak) : "unavailable")
        <<" global_dense_matrix=none"<<std::endl;
}
TEST(JointFixedNeighborBlockCoordinateTest, EndpointReplayUsesStreamingVerificationForLargeComponents)
{
    constexpr Eigen::Index atoms=37'406;
    std::vector<std::vector<n::Support>> support(static_cast<std::size_t>(atoms));
    for(auto & atom:support) atom.push_back({0,.5});
    const n::Domain domain(8,std::move(support));
    const n::Vector y=n::Vector::Zero(8),eta=n::Vector::Constant(atoms,std::log(.55));
    const auto context=n::CreateContext(y,atoms);
    n::AssessmentWorkForTesting()={};
    n::ProfileSearchWork work;
    const auto reference=n::EvaluateProfile(domain,y,eta,true,&context,nullptr,
        n::ProfileEvaluationRole::Reference,&work);
    EXPECT_FALSE(reference.valid);
    EXPECT_EQ(reference.reason,"resource-limit:reference-solve");
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,1);
    EXPECT_DOUBLE_EQ(work.total.basis_seconds,0.0);

    n::BlockCoordinateState state; state.eta=eta; state.beta=n::Vector::Zero(2*atoms);
    n::AssessmentWorkForTesting()={};
    const auto certification=n::CertifyFixedNeighborEndpoint(domain,y,state,context);
    EXPECT_EQ(certification.verification_memory_estimate.stage,n::MemoryStage::EndpointVerification);
    EXPECT_EQ(certification.verification_memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_TRUE(certification.trust.primary_valid);
    EXPECT_TRUE(certification.verified())<<certification.trust.reason;
    EXPECT_TRUE(certification.assessment.endpoint_verified);
    EXPECT_TRUE(certification.assessment.numerical_qualified);
    EXPECT_TRUE(certification.qualified());
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,0);
}
TEST(JointFixedNeighborBlockCoordinateTest, AggregateWorkTelemetryIsOptInAndDoesNotChangeSearchState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=2;
    policy.certify_endpoint=false; policy.collect_telemetry=true;
    const auto measured=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    policy.collect_telemetry=false;
    const auto quiet=n::RunFixedNeighborDiagnostic(problem,eta,policy);

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
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",16));
    const n::Vector eta=n::Vector::Constant(16,std::log(.55));
    n::FixedNeighborPolicy measured_policy; measured_policy.core_atoms=8; measured_policy.maximum_sweeps=1;
    measured_policy.certify_endpoint=false; measured_policy.collect_telemetry=true;
    const auto measured=n::RunFixedNeighborDiagnostic(problem,eta,measured_policy);
    auto quiet_policy=measured_policy; quiet_policy.collect_telemetry=false;
    const auto quiet=n::RunFixedNeighborDiagnostic(problem,eta,quiet_policy);
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
TEST(JointFixedNeighborBlockCoordinateTest, SearchResultAndIndependentCertificationStaySeparate)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(32,std::log(.55));
    const n::Domain domain=data.domain;
    const auto context=data.context;
    n::FixedNeighborPolicy quiet_policy; quiet_policy.core_atoms=16; quiet_policy.maximum_sweeps=3;
    const auto quiet=n::SearchFixedNeighborComponent(domain,data.layout,data.y,initial_eta,context,quiet_policy);
    auto diagnostic_policy=quiet_policy; diagnostic_policy.collect_diagnostics=true;
    const auto diagnostic=n::SearchFixedNeighborComponent(domain,data.layout,data.y,initial_eta,context,diagnostic_policy);
    const auto quiet_state=quiet.state;
    const auto quiet_certification=n::CertifyFixedNeighborEndpoint(domain,data.y,quiet.state,context);
    const auto diagnostic_certification=n::CertifyFixedNeighborEndpoint(domain,data.y,diagnostic.state,context);

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
    EXPECT_EQ(quiet_certification.verified(),diagnostic_certification.verified());
    EXPECT_TRUE((quiet.state.eta.array()==quiet_state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==quiet_state.beta.array()).all());
    EXPECT_TRUE((quiet.state.eta.array()==diagnostic.state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==diagnostic.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(quiet.state.objective,diagnostic.state.objective);
    EXPECT_DOUBLE_EQ(quiet.final_sweep->global_ac_kkt,diagnostic.sweeps.back().global_ac_kkt);
    EXPECT_DOUBLE_EQ(quiet.final_sweep->global_width_gradient_inf_norm,
        diagnostic.sweeps.back().global_width_gradient_inf_norm);
    EXPECT_TRUE((quiet_certification.assessment.primary.eta.array()==
        diagnostic_certification.assessment.primary.eta.array()).all());
    EXPECT_TRUE((quiet_certification.assessment.primary.eta.array()==quiet.state.eta.array()).all());
    EXPECT_TRUE((quiet_certification.assessment.primary.beta.array()==quiet.state.beta.array()).all());
}
TEST(JointFixedNeighborBlockCoordinateTest, FixedNeighborAcceptsAtMostOneTrustedLocalUpdatePerBlockVisit)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=2;
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
    const auto production=n::FixedNeighborSearchPolicy{}.convergence_criteria;
    EXPECT_DOUBLE_EQ(production.global_ac_kkt,n::PracticalJointConvergenceCriteria.global_ac_kkt);
    EXPECT_DOUBLE_EQ(production.global_width_gradient,n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_DOUBLE_EQ(production.confirmed_log_width_change,
        n::PracticalJointConvergenceCriteria.confirmed_log_width_change);
    const auto local=n::FixedNeighborPolicy{}.convergence_criteria;
    EXPECT_DOUBLE_EQ(local.global_ac_kkt,n::StrictJointConvergenceCriteria.global_ac_kkt);
    EXPECT_DOUBLE_EQ(local.global_width_gradient,n::StrictJointConvergenceCriteria.global_width_gradient);
    EXPECT_DOUBLE_EQ(local.confirmed_log_width_change,n::StrictJointConvergenceCriteria.confirmed_log_width_change);
}
TEST(JointFixedNeighborBlockCoordinateTest, CertificationRetryWaitsTenSweepsWhileStationarityPersists)
{
    std::size_t next_retry_sweep{};
    const auto stationary=n::ClassifyFixedNeighborStationarity(1e-10,1e-12,0.0,true);
    ASSERT_EQ(stationary,n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(stationary,5,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,15u);
    EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(stationary,5,next_retry_sweep));
    for(std::size_t sweep=6;sweep<15;++sweep)
        EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(stationary,sweep,next_retry_sweep));
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(stationary,15,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,25u);
}
TEST(JointFixedNeighborBlockCoordinateTest, CertificationRetryResetsAfterKktLoss)
{
    std::size_t next_retry_sweep{};
    ASSERT_TRUE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::ConfirmedStationary,5,next_retry_sweep));
    const auto lost=n::ClassifyFixedNeighborStationarity(1.1e-10,5e-13,0.0,true);
    ASSERT_EQ(lost,n::FixedNeighborStationarityState::NotStationary);
    EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(lost,10,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,0u);
    const auto restored=n::ClassifyFixedNeighborStationarity(5e-11,5e-13,0.0,true);
    EXPECT_EQ(restored,n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(restored,11,next_retry_sweep));
}
TEST(JointFixedNeighborBlockCoordinateTest, CertificationRetryResetsAfterWidthGradientLoss)
{
    std::size_t next_retry_sweep{};
    ASSERT_TRUE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::ConfirmedStationary,5,next_retry_sweep));
    const auto lost=n::ClassifyFixedNeighborStationarity(5e-11,1.1e-12,0.0,true);
    ASSERT_EQ(lost,n::FixedNeighborStationarityState::NotStationary);
    EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(lost,10,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,0u);
    const auto restored=n::ClassifyFixedNeighborStationarity(5e-11,5e-13,0.0,true);
    EXPECT_EQ(restored,n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(restored,11,next_retry_sweep));
}
TEST(JointFixedNeighborBlockCoordinateTest, EtaChangeLossResetsCertificationRetryEligibility)
{
    std::size_t next_retry_sweep{};
    ASSERT_TRUE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::ConfirmedStationary,5,next_retry_sweep));
    const auto candidate=n::ClassifyFixedNeighborStationarity(1e-10,1e-12,1.1e-10,true);
    ASSERT_EQ(candidate,n::FixedNeighborStationarityState::CandidateStationary);
    EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(candidate,10,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,0u);
    const auto reconfirmed=n::ClassifyFixedNeighborStationarity(1e-10,1e-12,1e-10,true);
    EXPECT_EQ(reconfirmed,n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(reconfirmed,11,next_retry_sweep));
}
TEST(JointFixedNeighborBlockCoordinateTest, SweepCoordinateTelemetryUsesCompleteSweepStates)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
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
    const auto observed=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    auto quiet_policy=policy; quiet_policy.state_observer={};
    const auto quiet=n::RunFixedNeighborDiagnostic(problem,initial_eta,quiet_policy);
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
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55)); n::FixedNeighborPolicy policy; policy.core_atoms=0;
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    EXPECT_FALSE(result.search_converged); EXPECT_EQ(result.reason,"block-invalid-partition");
    EXPECT_FALSE(result.endpoint_certified); EXPECT_EQ(result.sweeps.size(),0u);
}
TEST(JointFixedNeighborBlockCoordinateTest, BudgetExhaustionPreservesCandidateWithoutInventingTrustedState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=1;
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
    EXPECT_EQ(trust->status,rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointFixedNeighborBlockCoordinateTest, InvalidEndpointCertificationStopsWithoutRetry)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto & view=data.partition.components.front();
    const auto context=n::ChildContext(data.context,view,true);
    n::Vector y(static_cast<Eigen::Index>(view.rows.size()));
    for(std::size_t k=0;k<view.rows.size();++k) y(static_cast<Eigen::Index>(k))=data.y(view.rows[k]);
    n::BlockCoordinateState invalid;
    invalid.eta=n::Vector::Constant(8,std::log(.55)); invalid.beta=n::Vector::Zero(16);
    invalid.beta(0)=std::numeric_limits<double>::quiet_NaN();
    const auto certification=n::CertifyFixedNeighborEndpoint(view.domain,y,invalid,context);
    EXPECT_FALSE(certification.verified());
    EXPECT_FALSE(certification.trust.primary_valid);
    EXPECT_EQ(certification.trust.reason,"invalid-primary");

    n::FixedNeighborPolicy policy; policy.maximum_sweeps=30;
    std::size_t certification_attempts{};
    policy.certification_for_testing=[&](const auto & endpoint_domain,n::VectorRef endpoint_y,
        const n::BlockCoordinateState & state,const n::EvaluationContext & endpoint_context,
        n::ProfileSearchWork * work) {
        ++certification_attempts;
        auto result=n::CertifyFixedNeighborEndpoint(endpoint_domain,endpoint_y,state,endpoint_context,work);
        result.trust.passed=false; result.trust.primary_valid=false; result.trust.reason="invalid-primary";
        result.assessment.independent_replay.available=true;
        result.assessment.independent_replay.valid=false;
        result.assessment.independent_replay.passed=false;
        result.assessment.independent_replay.reason="invalid-primary";
        result.assessment.endpoint_verified=false;
        result.assessment.numerical_qualified=false;
        return result;
    };
    const auto search=n::RunFixedNeighborDiagnostic(problem,n::Vector::Constant(8,std::log(.55)),policy);
    EXPECT_FALSE(search.search_converged);
    EXPECT_EQ(search.reason,"invalid-primary");
    EXPECT_EQ(certification_attempts,1u);
    EXPECT_GT(search.confirmed_stationarity_sweep,0u);
    EXPECT_EQ(search.sweeps.size(),search.confirmed_stationarity_sweep);
    EXPECT_FALSE(search.fit.search_completed);
    ASSERT_EQ(search.fit.components.size(),1u);
    EXPECT_EQ(search.fit.components.front().stop_reason,"invalid-primary");
    const auto endpoint_trust=std::find_if(search.fit.components.front().evidence.begin(),
        search.fit.components.front().evidence.end(),[](const auto & check){return check.name=="independent-replay";});
    ASSERT_NE(endpoint_trust,search.fit.components.front().evidence.end());
    EXPECT_EQ(endpoint_trust->status,rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(endpoint_trust->reason,"invalid-primary");
    EXPECT_FALSE(search.fit.components.front().state);
    EXPECT_FALSE(search.fit.assembled_state);
    EXPECT_FALSE(search.fit.prediction);
    EXPECT_FALSE(search.fit.objective);
    EXPECT_EQ(search.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
}
TEST(JointFixedNeighborBlockCoordinateTest, RecoverableWidthStationarityFailureRetriesAtBoundedInterval)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=40;
    std::size_t current_sweep{},attempts{}; std::vector<std::size_t> certification_sweeps;
    policy.sweep_observer=[&](const auto & sweep){current_sweep=sweep.sweep;};
    policy.certification_for_testing=[&](const auto & domain,n::VectorRef y,const n::BlockCoordinateState & state,
        const n::EvaluationContext & context,n::ProfileSearchWork * work) {
        ++attempts; certification_sweeps.push_back(current_sweep);
        auto result=n::CertifyFixedNeighborEndpoint(domain,y,state,context,work);
        if(attempts==1)
        {
            result.assessment.gradient=false;
            result.assessment.numerical_qualified=false;
            result.assessment.failure="global-width-gradient-failed";
        }
        return result;
    };
    const auto search=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    ASSERT_TRUE(search.search_converged)<<search.reason;
    ASSERT_TRUE(search.endpoint_certified);
    ASSERT_EQ(attempts,2u);
    ASSERT_EQ(certification_sweeps.size(),2u);
    EXPECT_EQ(certification_sweeps.front(),search.confirmed_stationarity_sweep);
    EXPECT_GE(certification_sweeps.back(),certification_sweeps.front()+10);
    EXPECT_EQ(search.reason,"block-stationary");
    EXPECT_TRUE(search.assessment.numerical_qualified);
    ASSERT_EQ(search.fit.components.size(),1u);
    EXPECT_TRUE(search.fit.components.front().search_completed);
    EXPECT_TRUE(search.fit.components.front().state);
    EXPECT_EQ(search.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointFixedNeighborBlockCoordinateTest, CertificationAtLastBudgetSweepIsNotRepeated)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy search_only; search_only.certify_endpoint=false;
    const auto confirmed=n::RunFixedNeighborDiagnostic(problem,eta,search_only);
    ASSERT_TRUE(confirmed.search_converged)<<confirmed.reason;
    ASSERT_GT(confirmed.confirmed_stationarity_sweep,0u);

    n::FixedNeighborPolicy policy; policy.maximum_sweeps=confirmed.confirmed_stationarity_sweep;
    std::size_t current_sweep{}; std::vector<std::size_t> certification_sweeps;
    policy.sweep_observer=[&](const auto & sweep){current_sweep=sweep.sweep;};
    policy.certification_for_testing=[&](const auto & domain,n::VectorRef y,const n::BlockCoordinateState & state,
        const n::EvaluationContext & context,n::ProfileSearchWork * work) {
        certification_sweeps.push_back(current_sweep);
        auto result=n::CertifyFixedNeighborEndpoint(domain,y,state,context,work);
        result.trust.passed=false; result.trust.primary_valid=true; result.trust.reason="independent-replay-disagreement";
        result.assessment.independent_replay.passed=false;
        result.assessment.independent_replay.reason=result.trust.reason;
        result.assessment.endpoint_verified=false;
        result.assessment.numerical_qualified=false;
        return result;
    };
    const auto budget=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    EXPECT_FALSE(budget.search_converged);
    EXPECT_EQ(budget.reason,"block-sweep-budget");
    EXPECT_EQ(budget.sweeps.size(),confirmed.confirmed_stationarity_sweep);
    ASSERT_EQ(certification_sweeps.size(),1u);
    EXPECT_EQ(certification_sweeps.front(),budget.sweeps.size());
    ASSERT_EQ(budget.fit.components.size(),1u);
    EXPECT_EQ(budget.fit.components.front().endpoint_certification_sweep,budget.sweeps.size());
    EXPECT_TRUE((budget.state.eta.array()==budget.assessment.primary.eta.array()).all());
    EXPECT_TRUE((budget.state.beta.array()==budget.assessment.primary.beta.array()).all());
    EXPECT_FALSE(budget.fit.components.front().state);
    EXPECT_FALSE(budget.fit.assembled_state);
}
TEST(JointFixedNeighborBlockCoordinateTest, ReferenceFree168StopsAtFirstCorrectionFreeQualification)
{
    auto [problem,widths]=FixedNeighborRetryFixture();
    ASSERT_EQ(widths.size(),168u);
    ASSERT_EQ(problem.Input().observations.size(),30108u);
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.partition.components.size(),1u);
    const auto blocks=n::BuildStructuralBlockPartition(*data.input,data.layout,12);
    ASSERT_EQ(blocks.cores.size(),23u);
    struct SweepPoint { n::FixedNeighborBlockSweep sweep; n::BlockCoordinateState state; };
    struct SearchRun
    {
        rhbm_gem::core::JointFitResult fit;
        std::vector<std::size_t> certification_sweeps;
        std::vector<n::FixedNeighborEndpointCertification> certifications;
        std::vector<SweepPoint> trajectory;
        int production_reference_evaluations{};
    };
    const auto run=[&](std::size_t maximum_sweeps,bool stop_after_stationarity) {
        SearchRun result;
        n::AssessmentWorkForTesting()={};
        n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=maximum_sweeps;
        policy.convergence_criteria=n::StrictJointConvergenceCriteria;
        std::size_t current_sweep{};
        policy.configure_search_for_testing=[&](n::FixedNeighborPolicy & search_policy) {
            search_policy.stop_after_stationarity=stop_after_stationarity;
            search_policy.sweep_observer=[&](const auto & sweep){current_sweep=sweep.sweep;};
            search_policy.state_observer=[&](std::size_t,const auto & state,const auto & sweep,const auto &) {
                result.trajectory.push_back({sweep,state});
            };
            search_policy.certification_for_testing=[&](const auto & domain,n::VectorRef y,
                const n::BlockCoordinateState & state,const n::EvaluationContext & context,
                n::ProfileSearchWork * work) {
                auto certification=n::CertifyFixedNeighborEndpoint(domain,y,state,context,work);
                result.certification_sweeps.push_back(current_sweep);
                result.certifications.push_back(certification);
                return certification;
            };
        };
        result.fit=n::FitFixedNeighborComponents(problem,widths,policy);
        result.production_reference_evaluations=n::AssessmentWorkForTesting().reference_evaluations;
        return result;
    };
    const auto qualified_run=run(200,true);
    const auto continued_run=run(127,false);
    EXPECT_EQ(qualified_run.production_reference_evaluations,0);
    EXPECT_EQ(continued_run.production_reference_evaluations,0);
    ASSERT_EQ(qualified_run.fit.components.size(),1u);
    const auto & component=qualified_run.fit.components.front();
    ASSERT_TRUE(component.candidate_diagnostic);
    EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
    ASSERT_TRUE(component.state);
    ASSERT_TRUE(qualified_run.fit.assembled_state);
    EXPECT_TRUE(qualified_run.fit.search_completed);
    EXPECT_EQ(component.search_sweeps,107u);
    EXPECT_EQ(component.first_order_stationarity_sweep,106u);
    EXPECT_EQ(component.confirmed_stationarity_sweep,107u);
    EXPECT_EQ(component.endpoint_certification_sweep,107u);
    EXPECT_EQ(qualified_run.certification_sweeps,(std::vector<std::size_t>{107u}));
    ASSERT_EQ(qualified_run.certifications.size(),1u);
    EXPECT_EQ(qualified_run.trajectory.size(),107u);
    EXPECT_EQ(component.candidate_diagnostic->numerical_qualification_status,
        rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(component.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(qualified_run.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_FALSE(std::any_of(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
        return check.name=="local-correction";
    }));

    const auto evidence_status=[&](const char * name) {
        const auto found=std::find_if(component.evidence.begin(),component.evidence.end(),[&](const auto & check) {
            return check.name==name;
        });
        EXPECT_NE(found,component.evidence.end())<<name;
        return found==component.evidence.end() ? rhbm_gem::JointCheckStatus::Unavailable : found->status;
    };
    for(const auto * name:{"ac-feasibility","conditional-ac-optimality","independent-ac-optimality",
            "independent-kkt-consistency","independent-replay","width-gradient-consistency",
            "width-stationarity","endpoint-verification","numerical-qualification",
            "complete-sweep-stationarity"})
        EXPECT_EQ(evidence_status(name),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(component.target_evidence);
    const auto target_validation=std::find_if(component.target_evidence->checks.begin(),
        component.target_evidence->checks.end(),[](const auto & check){return check.name=="target-point-validation";});
    ASSERT_NE(target_validation,component.target_evidence->checks.end());
    EXPECT_EQ(target_validation->status,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(component.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);

    ASSERT_EQ(continued_run.trajectory.size(),127u);
    EXPECT_TRUE(continued_run.certification_sweeps.empty());
    for(std::size_t k=0;k<qualified_run.trajectory.size();++k)
    {
        const auto & a=qualified_run.trajectory[k]; const auto & b=continued_run.trajectory[k];
        EXPECT_EQ(a.sweep.sweep,b.sweep.sweep);
        EXPECT_EQ(a.sweep.accepted_blocks,b.sweep.accepted_blocks);
        EXPECT_EQ(a.sweep.accepted_local_updates,b.sweep.accepted_local_updates);
        EXPECT_EQ(a.sweep.profile_evaluations,b.sweep.profile_evaluations);
        EXPECT_DOUBLE_EQ(a.sweep.objective_after,b.sweep.objective_after);
        EXPECT_DOUBLE_EQ(a.sweep.global_ac_kkt,b.sweep.global_ac_kkt);
        EXPECT_DOUBLE_EQ(a.sweep.global_width_gradient_inf_norm,b.sweep.global_width_gradient_inf_norm);
        EXPECT_DOUBLE_EQ(a.sweep.cache_replay_error,b.sweep.cache_replay_error);
        EXPECT_DOUBLE_EQ(a.state.objective,b.state.objective);
        EXPECT_TRUE((a.state.eta.array()==b.state.eta.array()).all());
        EXPECT_TRUE((a.state.beta.array()==b.state.beta.array()).all());
        EXPECT_TRUE((a.state.prediction.array()==b.state.prediction.array()).all());
        EXPECT_TRUE((a.state.residual.array()==b.state.residual.array()).all());
    }
    const auto & sweep117=continued_run.trajectory[116].state;
    const auto sweep117_certification=n::CertifyFixedNeighborEndpoint(data.domain,data.y,sweep117,data.context);
    EXPECT_TRUE(sweep117_certification.verified());
    EXPECT_TRUE(sweep117_certification.qualified());
    const auto & old_qualified=continued_run.trajectory[126].state;
    const auto sweep127_certification=n::CertifyFixedNeighborEndpoint(data.domain,data.y,old_qualified,data.context);
    EXPECT_TRUE(sweep127_certification.verified());
    EXPECT_TRUE(sweep127_certification.qualified());

    const auto & early=qualified_run.trajectory.back().state;
    const auto & early_assessment=qualified_run.certifications.front().assessment;
    EXPECT_TRUE(early_assessment.numerical_qualified);
    const auto early_primary=n::EvaluateState(data.domain,data.y,early.eta,early.beta,data.context);
    ASSERT_TRUE(early_primary.valid)<<early_primary.reason;
    const auto early_replay=n::CheckIndependentReplay(data.domain,data.y,early_primary,data.context);
    ASSERT_TRUE(early_replay.available && early_replay.valid);
    EXPECT_TRUE(early_replay.passed);
    EXPECT_EQ(component.state->ac.size(),static_cast<std::size_t>(early.beta.size()));
    EXPECT_EQ(component.state->log_b.size(),static_cast<std::size_t>(early.eta.size()));
    for(std::size_t k=0;k<component.state->ac.size();++k)
        EXPECT_DOUBLE_EQ(component.state->ac[k],early.beta(static_cast<Eigen::Index>(k)));
    for(std::size_t k=0;k<component.state->log_b.size();++k)
        EXPECT_DOUBLE_EQ(component.state->log_b[k],early.eta(static_cast<Eigen::Index>(k)));

    double amplitude_difference{},charge_difference{},width_difference{};
    for(Eigen::Index k=0;k<early.beta.size();k+=2)
    {
        amplitude_difference=std::max(amplitude_difference,std::abs(early.beta(k)-old_qualified.beta(k)));
        charge_difference=std::max(charge_difference,std::abs(early.beta(k+1)-old_qualified.beta(k+1)));
        width_difference=std::max(width_difference,
            std::abs(std::exp(early.eta(k/2))-std::exp(old_qualified.eta(k/2))));
    }
    const double objective_difference=std::abs(early.objective-old_qualified.objective);
    const double residual_difference=(early_primary.residual-old_qualified.residual).norm();
    RecordProperty("correction_free_168_new_qualified_sweep",std::to_string(component.search_sweeps));
    RecordProperty("correction_free_168_old_contract_qualified_sweep","127");
    RecordProperty("correction_baseline_168_sweep107_local_correction","9.0237250078785534e-10");
    RecordProperty("correction_baseline_168_sweep117_mapped_correction","1.5377372185134297e-10");
    RecordProperty("correction_baseline_168_old_objective","2.6038831254099719e-16");
    RecordProperty("correction_baseline_168_old_ac_kkt","1.300543006621524e-12");
    RecordProperty("correction_baseline_168_old_width_gradient","1.9630561012973708e-14");
    RecordProperty("correction_baseline_168_old_replay_prediction_difference","3.5527136788005009e-15");
    RecordProperty("correction_baseline_168_old_replay_residual_difference","3.5527136788005009e-15");
    RecordProperty("correction_baseline_168_old_replay_objective_difference","1.3676001633950818e-21");
    RecordProperty("correction_free_168_max_amplitude_difference",std::to_string(amplitude_difference));
    RecordProperty("correction_free_168_max_charge_difference",std::to_string(charge_difference));
    RecordProperty("correction_free_168_max_width_difference",std::to_string(width_difference));
    RecordProperty("correction_free_168_objective_difference",std::to_string(objective_difference));
    RecordProperty("correction_free_168_residual_difference",std::to_string(residual_difference));
    RecordProperty("correction_free_168_ac_kkt",std::to_string(early_primary.certificate.projected_kkt));
    RecordProperty("correction_free_168_width_gradient",std::to_string(early_primary.gradient.lpNorm<Eigen::Infinity>()));
    RecordProperty("correction_free_168_replay_prediction_difference",std::to_string(early_replay.prediction_difference));
    RecordProperty("correction_free_168_replay_objective_difference",std::to_string(early_replay.objective_difference));
    EXPECT_EQ(component.candidate_diagnostic->endpoint_verification_status,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(component.state);
    EXPECT_TRUE(qualified_run.fit.assembled_state);
    const auto saved=rhbm_gem::core::CaptureJointAnalysisResult(qualified_run.fit);
    EXPECT_EQ(saved.components.front().runtime_convergence,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(saved.runtime_convergence,rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointFixedNeighborBlockCoordinateTest, SameEtaNumericalAssessmentsPreserveSearchState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(8,std::log(.55));
    const auto search=n::RunFixedNeighborDiagnostic(problem,initial_eta);
    ASSERT_TRUE(search.search_converged);
    const n::Vector eta=search.state.eta,beta=search.state.beta;
    const auto & domain=data.domain; const auto & context=data.context; const auto & y=data.y;
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
    ASSERT_TRUE(raw.certificate.available && reference.certificate.available);
    const auto difference=[](const n::Vector & lhs,const n::Vector & rhs) {
        return ((lhs-rhs).array().abs()/(1+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
    };
    EXPECT_LE(difference(beta,primary.beta),1e-10);
    EXPECT_LE(difference(primary.beta,reference.beta),1e-10);
    EXPECT_LE(difference(beta,reference.beta),1e-10);
    EXPECT_TRUE(raw_assessment.inner && primary_assessment.inner && reference_assessment.inner);
    EXPECT_TRUE(raw_assessment.gradient && primary_assessment.gradient && reference_assessment.gradient);
    EXPECT_TRUE(raw_assessment.numerical_qualified && primary_assessment.numerical_qualified &&
        reference_assessment.numerical_qualified);
    EXPECT_TRUE(raw_trust.passed && primary_trust.passed && reference_trust.passed);
    EXPECT_LE(raw.certificate.projected_kkt,1e-10);
    EXPECT_LE(reference.certificate.projected_kkt,1e-10);
    EXPECT_LE(raw.gradient.lpNorm<Eigen::Infinity>(),1e-12);
    EXPECT_LE(reference.gradient.lpNorm<Eigen::Infinity>(),1e-12);
    RecordProperty("numerical_qualification",raw_assessment.numerical_qualified ? "passed" : "failed");
    RecordProperty("endpoint_replay",raw_trust.passed ? "passed" : "failed");
}
