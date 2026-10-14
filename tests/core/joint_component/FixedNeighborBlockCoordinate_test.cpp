#include <gtest/gtest.h>
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "support/BoundedRankResearch.hpp"
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
    EXPECT_GT(result.confirmed_stationarity_sweep,result.first_order_stationarity_sweep);
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
    EXPECT_LE(result.sweeps.back().global_ac_kkt,1e-10);
    EXPECT_LE(result.sweeps.back().global_width_gradient_inf_norm,1e-12);
    EXPECT_LE(result.sweeps.back().cache_replay_error,2e-12);
    EXPECT_TRUE(result.assessment.primary.valid);
    EXPECT_TRUE(result.assessment.independent_replay.passed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(result.fit.search_completed);
    const auto production=rhbm_gem::core::FitJointComponents(problem,std::vector<double>(static_cast<std::size_t>(atoms),.55));
    ASSERT_EQ(production.components.size(),1u);
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

TEST(JointFixedNeighborBlockCoordinateTest, RankFailureReportsLocalDesignAndHaloDiagnostics)
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
    policy.halo_census=[](std::size_t & count,std::size_t & minimum,std::size_t & deficient) {
        count=1; minimum=2; deficient=1;
    };
    const auto diagnostic=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,context,policy);
    auto quiet_policy=policy; quiet_policy.collect_diagnostics=false;
    const auto quiet=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,context,quiet_policy);

    EXPECT_EQ(diagnostic.reason,"block-search-failed");
    EXPECT_EQ(quiet.reason,diagnostic.reason);
    EXPECT_TRUE((quiet.state.eta.array()==diagnostic.state.eta.array()).all());
    EXPECT_TRUE((quiet.state.beta.array()==diagnostic.state.beta.array()).all());
    ASSERT_FALSE(diagnostic.blocks.empty());
    ASSERT_TRUE(diagnostic.blocks.front().failure_diagnostics);
    const auto & failure=*diagnostic.blocks.front().failure_diagnostics;
    EXPECT_EQ(failure.component_index,4u);
    EXPECT_EQ(failure.component_id,"weak-component");
    EXPECT_EQ(failure.sweep_index,1u);
    EXPECT_EQ(failure.core_atom_ids,(std::vector<std::string>{"a-weak-halo"}));
    EXPECT_EQ(failure.core_atom_count,1u);
    EXPECT_EQ(failure.affected_row_count,2u);
    EXPECT_EQ(failure.halo_count,1u);
    EXPECT_EQ(failure.minimum_halo_support_rows,2u);
    EXPECT_EQ(failure.rank_deficient_halo_count,1u);
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
    auto input=std::make_shared<rhbm_gem::core::JointProblemInput>();
    input->observations.assign(static_cast<std::size_t>(rows),0.0);
    input->row_ids.resize(static_cast<std::size_t>(rows));
    input->atom_ids.reserve(atoms); input->support.resize(atoms);
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        input->atom_ids.push_back("atom-"+std::to_string(atom));
        input->support[atom].push_back({atom,0.0});
        if(atom+1<atoms) input->support[atom].push_back({atom+1,0.0});
    }
    input->selection_domain.emplace(); input->selection_domain->target_indices.resize(atoms);
    std::iota(input->selection_domain->target_indices.begin(),input->selection_domain->target_indices.end(),std::size_t{0});
    const std::shared_ptr<const rhbm_gem::core::JointProblemInput> snapshot=input;
    n::Domain domain(snapshot); domain.rows=rows;
    rhbm_gem::JointParameterLayout layout; layout.full_atoms.resize(atoms);
    std::iota(layout.full_atoms.begin(),layout.full_atoms.end(),std::size_t{0});
    layout.informative_rows.resize(static_cast<std::size_t>(rows));
    std::iota(layout.informative_rows.begin(),layout.informative_rows.end(),Eigen::Index{0});
    n::EvaluationContext context; context.scale=1; context.target_columns=3*atoms;
    context.rank={rows,2*static_cast<Eigen::Index>(atoms),static_cast<Eigen::Index>(atoms)};
    context.linear.rank_relative=context.rank.Relative(context.rank.design_columns);
    context.observations=n::Observe(snapshot);
    context.atom_ids=n::Identities(std::shared_ptr<const std::vector<std::string>>(snapshot,&snapshot->atom_ids));
    context.row_ids=n::Identities(std::shared_ptr<const std::vector<std::string>>(snapshot,&snapshot->row_ids));
    const n::Vector y=n::Vector::Zero(rows),eta=n::Vector::Constant(static_cast<Eigen::Index>(atoms),std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=0; policy.certify_endpoint=false;
    const auto result=n::SearchFixedNeighborComponent(domain,layout,y,eta,context,policy);

    const auto support_entries=n::CountSupportEntries(*snapshot);
    ASSERT_TRUE(support_entries);
    EXPECT_EQ(*support_entries,2*atoms-1);
    EXPECT_EQ(result.memory_estimate.shape.component_atoms,atoms);
    EXPECT_EQ(result.memory_estimate.shape.component_observations,static_cast<std::size_t>(rows));
    EXPECT_EQ(result.memory_estimate.shape.component_support_entries,*support_entries);
    EXPECT_EQ(result.memory_estimate.shape.component_blocks,3'118u);
    EXPECT_EQ(result.memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_EQ(result.global_state_memory_estimate.status,n::MemoryStatus::Pass);
    EXPECT_LE(result.global_state_memory_estimate.estimated_peak_bytes,n::JointMemoryPlanningBudgetBytes);
    EXPECT_EQ(result.sweep_count,0u);
    EXPECT_EQ(result.reason,"block-sweep-budget");

    std::size_t legacy_bytes{};
    ASSERT_TRUE(n::CheckedMultiplyBytes(static_cast<std::size_t>(rows),result.memory_estimate.shape.component_blocks,legacy_bytes));
    ASSERT_TRUE(n::CheckedMultiplyBytes(legacy_bytes,sizeof(Eigen::Index),legacy_bytes));
    EXPECT_EQ(legacy_bytes,148'119'642'128u);
    ASSERT_GT(result.prepared_mapping_bytes,0u);
    EXPECT_LT(result.prepared_mapping_bytes,legacy_bytes);
    EXPECT_LT(result.prepared_mapping_bytes,n::JointMemoryPlanningBudgetBytes);
    n::Domain endpoint_domain(snapshot); endpoint_domain.rows=rows;
    n::EvaluationContext endpoint_context; endpoint_context.scale=1;
    endpoint_context.target_columns=3*atoms;
    endpoint_context.rank={rows,2*static_cast<Eigen::Index>(atoms),static_cast<Eigen::Index>(atoms)};
    endpoint_context.linear.rank_relative=endpoint_context.rank.Relative(endpoint_context.rank.design_columns);
    endpoint_context.observations=n::Observe(snapshot);
    endpoint_context.atom_ids=n::Identities(std::shared_ptr<const std::vector<std::string>>(snapshot,&snapshot->atom_ids));
    endpoint_context.row_ids=n::Identities(std::shared_ptr<const std::vector<std::string>>(snapshot,&snapshot->row_ids));
    const auto endpoint_memory=n::EndpointVerificationMemoryEstimate(endpoint_domain,y,endpoint_context);
    EXPECT_EQ(endpoint_memory.stage,n::MemoryStage::EndpointVerification);
    EXPECT_EQ(endpoint_memory.status,n::MemoryStatus::Pass);
    EXPECT_LT(endpoint_memory.estimated_peak_bytes,n::JointMemoryPlanningBudgetBytes);
    n::BlockCoordinateState accepted; accepted.eta=eta; accepted.beta=n::Vector::Zero(2*atoms);
    n::AssessmentWorkForTesting()={};
    n::NumericsWorkForTesting()={};
    const auto certification=n::CertifyFixedNeighborEndpoint(endpoint_domain,y,accepted,endpoint_context);
    EXPECT_TRUE(certification.trust.primary_valid)<<certification.trust.reason;
    EXPECT_TRUE(certification.verified())<<certification.trust.reason;
    EXPECT_TRUE(certification.assessment.endpoint_verified);
    EXPECT_FALSE(certification.assessment.qualified);
    EXPECT_EQ(certification.assessment.rank_memory_estimate.status,n::MemoryStatus::Reject);
    EXPECT_FALSE(certification.assessment.design.has_value());
    ASSERT_TRUE(certification.assessment.design_structure.has_value());
    const auto & structure=*certification.assessment.design_structure;
    EXPECT_TRUE(structure.available);
    EXPECT_FALSE(structure.matching_available);
    EXPECT_TRUE(structure.proven_deficient);
    EXPECT_EQ(structure.rows,rows);
    EXPECT_EQ(structure.columns,74'812);
    EXPECT_EQ(structure.nonzero_entries,149'622);
    EXPECT_EQ(structure.supported_rows,37'406);
    EXPECT_EQ(structure.rank_upper_bound,37'406);
    EXPECT_EQ(certification.assessment.rank_memory_estimate.shape.component_free_columns,37'406u);
    std::vector<Eigen::Triplet<double>> free_face_entries; free_face_entries.reserve(*support_entries);
    for(Eigen::Index atom=0;atom<static_cast<Eigen::Index>(atoms);++atom)
        for(const auto & support:endpoint_domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto basis=n::EvaluateKernel(support.square,.55,2.5);
            if(basis.charge!=0) free_face_entries.emplace_back(support.row,atom,basis.charge);
        }
    n::Sparse free_face(rows,static_cast<Eigen::Index>(atoms));
    free_face.setFromTriplets(free_face_entries.begin(),free_face_entries.end()); free_face.makeCompressed();
    std::size_t free_face_bytes{},free_face_outer{};
    ASSERT_TRUE(n::CheckedMultiplyBytes(static_cast<std::size_t>(free_face.nonZeros()),
        sizeof(double)+sizeof(Eigen::Index),free_face_bytes));
    ASSERT_TRUE(n::CheckedAddBytes(static_cast<std::size_t>(free_face.outerSize()),1,free_face_outer));
    ASSERT_TRUE(n::CheckedMultiplyBytes(free_face_outer,sizeof(Eigen::Index),free_face_outer));
    ASSERT_TRUE(n::CheckedAddBytes(free_face_bytes,free_face_outer,free_face_bytes));
    ASSERT_TRUE(n::CheckedAddBytes(free_face_bytes,certification.assessment.rank_memory_estimate.resident_bytes,
        free_face_bytes));
    const auto free_structure=n::AssessStructuralRank(free_face,free_face_bytes);
    ASSERT_TRUE(free_structure.matching_available);
    EXPECT_TRUE(free_structure.structurally_full_rank);
    EXPECT_FALSE(free_structure.proven_deficient);
    EXPECT_EQ(free_structure.columns,37'406);
    EXPECT_EQ(free_structure.supported_rows,37'406);
    EXPECT_EQ(free_structure.matching_size,37'406);
    std::vector<Eigen::Triplet<double>> full_design_entries; full_design_entries.reserve(2*(*support_entries));
    for(Eigen::Index atom=0;atom<static_cast<Eigen::Index>(atoms);++atom)
        for(const auto & support:endpoint_domain.atoms[static_cast<std::size_t>(atom)])
        {
            const auto basis=n::EvaluateKernel(support.square,.55,2.5);
            if(basis.gaussian!=0) full_design_entries.emplace_back(support.row,2*atom,basis.gaussian);
            if(basis.charge!=0) full_design_entries.emplace_back(support.row,2*atom+1,basis.charge);
        }
    n::Sparse full_design(rows,2*static_cast<Eigen::Index>(atoms));
    full_design.setFromTriplets(full_design_entries.begin(),full_design_entries.end()); full_design.makeCompressed();
    std::size_t full_design_bytes{},full_design_outer{};
    ASSERT_TRUE(n::CheckedMultiplyBytes(static_cast<std::size_t>(full_design.nonZeros()),
        sizeof(double)+sizeof(Eigen::Index),full_design_bytes));
    ASSERT_TRUE(n::CheckedAddBytes(static_cast<std::size_t>(full_design.outerSize()),1,full_design_outer));
    ASSERT_TRUE(n::CheckedMultiplyBytes(full_design_outer,sizeof(Eigen::Index),full_design_outer));
    ASSERT_TRUE(n::CheckedAddBytes(full_design_bytes,full_design_outer,full_design_bytes));
    ASSERT_TRUE(n::CheckedAddBytes(full_design_bytes,free_face_bytes,full_design_bytes));
    const auto giant_streaming=second_stage_test::bounded_rank::ScalarGershgorin(full_design,
        endpoint_context.rank,full_design_bytes);
    EXPECT_TRUE(giant_streaming.available);
    EXPECT_FALSE(giant_streaming.certified_full_rank);
    EXPECT_LE(giant_streaming.estimated_peak_bytes,n::JointMemoryPlanningBudgetBytes);
    const auto giant_sparseqr=second_stage_test::bounded_rank::SparseQRRank(full_design,
        endpoint_context.rank,1.0,second_stage_test::bounded_rank::SparseQROrdering::COLAMD,
        n::CurrentRssBytes().value_or(full_design_bytes));
    EXPECT_FALSE(giant_sparseqr.preflight_passed);
    EXPECT_FALSE(giant_sparseqr.valid);
    EXPECT_EQ(giant_sparseqr.status,"ResourceLimit");
    EXPECT_EQ(giant_sparseqr.reason,"estimated-peak-exceeds-budget");
    const auto rank_evidence=n::AssessmentEvidence(certification.assessment,rhbm_gem::JointEvidenceScope::ComponentLocal);
    const auto design_check=std::find_if(rank_evidence.begin(),rank_evidence.end(),[](const auto & check) {
        return check.name=="design-identifiability";
    });
    ASSERT_NE(design_check,rank_evidence.end());
    EXPECT_EQ(design_check->status,rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(n::ConvergenceStatus(rank_evidence,rhbm_gem::JointEvidenceScope::ComponentLocal),
        rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,0);
    const auto structural_rss_before=n::CurrentRssBytes();
    const auto structural_started=std::chrono::steady_clock::now();
    const auto structural_repeat=n::AssessDesignStructuralRank(endpoint_domain,eta,
        certification.assessment.rank_memory_estimate.resident_bytes,
        certification.assessment.rank_memory_estimate.budget_bytes,false);
    const double structural_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-structural_started).count();
    const auto structural_rss_after=n::CurrentRssBytes();
    const auto structural_process_peak=n::ProcessPeakRssBytes();
    ASSERT_TRUE(structural_repeat.proven_deficient);
    EXPECT_EQ(structural_repeat.rank_upper_bound,structure.rank_upper_bound);
    const auto & rank_stages=certification.assessment.rank_memory_estimate.attributed_stages;
    ASSERT_EQ(rank_stages.size(),1u);
    auto interior_rank_shape=certification.assessment.rank_memory_estimate.shape;
    interior_rank_shape.component_free_columns=interior_rank_shape.component_design_columns;
    const auto interior_rank_memory=n::PlanMemory(interior_rank_shape,n::MemoryStage::ProjectedWidthAssessment,
        certification.assessment.rank_memory_estimate.budget_bytes);
    EXPECT_EQ(interior_rank_memory.status,n::MemoryStatus::Reject);
    auto mixed_rank_shape=certification.assessment.rank_memory_estimate.shape;
    mixed_rank_shape.component_free_columns=atoms+atoms/2;
    const auto mixed_rank_memory=n::PlanMemory(mixed_rank_shape,n::MemoryStage::ProjectedWidthAssessment,
        certification.assessment.rank_memory_estimate.budget_bytes);
    EXPECT_EQ(mixed_rank_memory.status,n::MemoryStatus::Reject);
    const auto target_rss_before=n::CurrentRssBytes();
    const auto target_started=std::chrono::steady_clock::now();
    const auto target_memory=n::TargetIdentifiabilityMemoryEstimate(*snapshot,layout,rows);
    const double target_preflight_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-target_started).count();
    const auto target_rss_after=n::CurrentRssBytes();
    const auto target_process_peak=n::ProcessPeakRssBytes();
    ASSERT_EQ(target_memory.shape.component_target_atoms,atoms);
    ASSERT_EQ(target_memory.shape.component_halo_atoms,0u);
    EXPECT_EQ(target_memory.status,n::MemoryStatus::Pass);
    const auto & assessment_stages=n::AssessmentWorkForTesting().stages;
    const auto find_assessment_stage=[&](const std::string & name)->const n::AssessmentStageWork * {
        const auto found=std::find_if(assessment_stages.begin(),assessment_stages.end(),[&](const auto & stage) {
            return stage.name==name;
        });
        return found==assessment_stages.end() ? nullptr : &*found;
    };
    const auto rank_aggregate=std::find_if(assessment_stages.begin(),assessment_stages.end(),[](const auto & stage) {
        return stage.name=="numerical-identifiability";
    });
    ASSERT_NE(rank_aggregate,assessment_stages.end());
    EXPECT_EQ(rank_aggregate->memory_status,n::MemoryStatus::Reject);
    EXPECT_EQ(rank_aggregate->attributed_memory_stages.size(),1u);
    EXPECT_EQ(std::find_if(assessment_stages.begin(),assessment_stages.end(),[](const auto & stage) {
        return stage.name=="design-spectrum";
    }),assessment_stages.end());
    EXPECT_EQ(std::find_if(assessment_stages.begin(),assessment_stages.end(),[](const auto & stage) {
        return stage.name=="projected-width-assessment";
    }),assessment_stages.end());
    EXPECT_EQ(std::find_if(assessment_stages.begin(),assessment_stages.end(),[](const auto & stage) {
        return stage.name=="local-correction";
    }),assessment_stages.end());
    const auto current_rss=n::CurrentRssBytes();
    const auto process_peak_rss=n::ProcessPeakRssBytes();
    EXPECT_TRUE(target_memory.attributed_stages.empty());
    EXPECT_EQ(target_memory.dense_factor_bytes,0u);
    EXPECT_EQ(target_memory.right_singular_vector_bytes,0u);
    const auto snapshot_bytes=n::EstimateSnapshotBytes(*snapshot);
    ASSERT_TRUE(snapshot_bytes);
    const auto * primary_stage=find_assessment_stage("endpoint-primary-evaluation");
    const auto * replay_stage=find_assessment_stage("independent-streaming-replay");
    ASSERT_NE(primary_stage,nullptr); ASSERT_NE(replay_stage,nullptr);
    RecordProperty("synthetic_support_entries",std::to_string(*support_entries));
    RecordProperty("prepared_mapping_bytes",std::to_string(result.prepared_mapping_bytes));
    RecordProperty("estimated_fixed_neighbor_peak_bytes",
        std::to_string(result.global_state_memory_estimate.estimated_peak_bytes));
    RecordProperty("estimated_endpoint_verification_peak_bytes",std::to_string(endpoint_memory.estimated_peak_bytes));
    std::cout<<"joint-large-synthetic rows="<<rows<<" atoms="<<atoms
        <<" support_entries="<<*support_entries
        <<" structural_nonzeros="<<structure.nonzero_entries
        <<" unique_supported_rows="<<structure.supported_rows
        <<" full_design_columns="<<structure.columns
        <<" structural_rank_upper_bound="<<structure.rank_upper_bound
        <<" structural_memory_estimate_bytes="<<structure.memory_estimate_bytes
        <<" structural_repeat_seconds="<<structural_seconds
        <<" structural_repeat_rss_before_bytes="<<(structural_rss_before ? std::to_string(*structural_rss_before) : "unavailable")
        <<" structural_repeat_rss_after_bytes="<<(structural_rss_after ? std::to_string(*structural_rss_after) : "unavailable")
        <<" structural_repeat_process_peak_rss_bytes="<<(structural_process_peak ? std::to_string(*structural_process_peak) : "unavailable")
        <<" structural_repeat_status="<<structural_repeat.reason
        <<" free_face_columns="<<free_structure.columns
        <<" free_face_nonzeros="<<free_structure.nonzero_entries
        <<" free_face_supported_rows="<<free_structure.supported_rows
        <<" free_face_matching="<<free_structure.matching_size
        <<" free_face_structural_status="<<free_structure.reason
        <<" free_face_memory_estimate_bytes="<<free_structure.memory_estimate_bytes
        <<" streaming_lambda_min_lower="<<giant_streaming.lambda_min_lower
        <<" streaming_lambda_max_upper="<<giant_streaming.lambda_max_upper
        <<" streaming_bound_gap="<<giant_streaming.bound_gap
        <<" streaming_status="<<giant_streaming.reason
        <<" streaming_estimated_peak_bytes="<<giant_streaming.estimated_peak_bytes
        <<" streaming_rss_before_bytes="<<(giant_streaming.current_rss_before_bytes ? std::to_string(*giant_streaming.current_rss_before_bytes) : "unavailable")
        <<" streaming_rss_after_bytes="<<(giant_streaming.current_rss_after_bytes ? std::to_string(*giant_streaming.current_rss_after_bytes) : "unavailable")
        <<" streaming_process_peak_rss_bytes="<<(giant_streaming.process_peak_rss_bytes ? std::to_string(*giant_streaming.process_peak_rss_bytes) : "unavailable")
        <<" streaming_seconds="<<giant_streaming.seconds
        <<" sparseqr_estimated_temporary_bytes="<<giant_sparseqr.estimated_temporary_bytes
        <<" sparseqr_estimated_peak_bytes="<<giant_sparseqr.estimated_peak_bytes
        <<" sparseqr_rss_before_bytes="<<(giant_sparseqr.current_rss_before_bytes ?
            std::to_string(*giant_sparseqr.current_rss_before_bytes) : "unavailable")
        <<" sparseqr_rss_after_bytes="<<(giant_sparseqr.current_rss_after_bytes ?
            std::to_string(*giant_sparseqr.current_rss_after_bytes) : "unavailable")
        <<" sparseqr_process_peak_rss_bytes="<<(giant_sparseqr.process_peak_rss_bytes ?
            std::to_string(*giant_sparseqr.process_peak_rss_bytes) : "unavailable")
        <<" sparseqr_status="<<giant_sparseqr.status
        <<" sparseqr_reason="<<giant_sparseqr.reason
        <<" sparseqr_factorization=NotRun"
        <<" endpoint_estimated_peak_bytes="<<endpoint_memory.estimated_peak_bytes
        <<" rank_estimated_peak_bytes="<<certification.assessment.rank_memory_estimate.estimated_peak_bytes
        <<" projected_width_assessment_peak_bytes="<<rank_stages[0].estimated_peak_bytes
        <<" all_interior_projected_width_peak_bytes="<<interior_rank_memory.estimated_peak_bytes
        <<" all_interior_resource_status="<<n::MemoryStatusName(interior_rank_memory.status)
        <<" all_interior_temporary_bytes="<<interior_rank_memory.temporary_bytes
        <<" all_interior_sparse_design_bytes="<<interior_rank_memory.sparse_design_bytes
        <<" all_interior_qr_tile_workspace_bytes="<<interior_rank_memory.qr_tile_workspace_bytes
        <<" all_interior_free_design_factor_bytes="<<interior_rank_memory.free_design_factor_bytes
        <<" all_interior_free_design_svd_workspace_bytes="<<interior_rank_memory.free_design_svd_workspace_bytes
        <<" all_interior_projection_coefficient_bytes="<<interior_rank_memory.projection_coefficient_bytes
        <<" all_interior_projected_width_factor_bytes="<<interior_rank_memory.projected_width_factor_bytes
        <<" all_interior_projected_width_svd_workspace_bytes="<<interior_rank_memory.projected_width_svd_workspace_bytes
        <<" mixed_active_free_design_columns="<<mixed_rank_memory.shape.component_free_columns
        <<" mixed_active_projected_width_peak_bytes="<<mixed_rank_memory.estimated_peak_bytes
        <<" mixed_active_resource_status="<<n::MemoryStatusName(mixed_rank_memory.status)
        <<" mixed_active_temporary_bytes="<<mixed_rank_memory.temporary_bytes
        <<" mixed_active_sparse_design_bytes="<<mixed_rank_memory.sparse_design_bytes
        <<" mixed_active_qr_tile_workspace_bytes="<<mixed_rank_memory.qr_tile_workspace_bytes
        <<" mixed_active_free_design_factor_bytes="<<mixed_rank_memory.free_design_factor_bytes
        <<" mixed_active_free_design_svd_workspace_bytes="<<mixed_rank_memory.free_design_svd_workspace_bytes
        <<" mixed_active_projection_coefficient_bytes="<<mixed_rank_memory.projection_coefficient_bytes
        <<" mixed_active_projected_width_factor_bytes="<<mixed_rank_memory.projected_width_factor_bytes
        <<" mixed_active_projected_width_svd_workspace_bytes="<<mixed_rank_memory.projected_width_svd_workspace_bytes
        <<" free_design_factor_bytes="<<rank_stages[0].free_design_factor_bytes
        <<" free_design_svd_workspace_bytes="<<rank_stages[0].free_design_svd_workspace_bytes
        <<" projection_coefficient_bytes="<<rank_stages[0].projection_coefficient_bytes
        <<" projected_width_factor_bytes="<<rank_stages[0].projected_width_factor_bytes
        <<" projected_width_svd_workspace_bytes="<<rank_stages[0].projected_width_svd_workspace_bytes
        <<" rank_assessment_temporary_bytes="<<certification.assessment.rank_memory_estimate.temporary_bytes
        <<" rank_assessment_sparse_design_bytes="<<rank_stages[0].sparse_design_bytes
        <<" rank_assessment_qr_tile_workspace_bytes="<<rank_stages[0].qr_tile_workspace_bytes
        <<" target_identifiability_peak_bytes="<<target_memory.estimated_peak_bytes
        <<" overall_estimated_peak_bytes="<<std::max({result.memory_estimate.estimated_peak_bytes,
            result.global_state_memory_estimate.estimated_peak_bytes,endpoint_memory.estimated_peak_bytes,
            certification.assessment.rank_memory_estimate.estimated_peak_bytes,target_memory.estimated_peak_bytes})
        <<" problem_resident_bytes="<<*snapshot_bytes
        <<" preparation_estimated_peak_bytes="<<result.memory_estimate.estimated_peak_bytes
        <<" global_state_estimated_peak_bytes="<<result.global_state_memory_estimate.estimated_peak_bytes
        <<" prepared_block_bytes="<<result.prepared_mapping_bytes
        <<" affected_rows="<<result.memory_estimate.shape.affected_row_entries
        <<" total_block_support_entries="<<result.memory_estimate.shape.total_block_support_entries
        <<" maximum_block_rows="<<result.memory_estimate.shape.maximum_block_rows
        <<" maximum_block_support_entries="<<result.memory_estimate.shape.maximum_block_support_entries
        <<" blocks="<<result.memory_estimate.shape.component_blocks
        <<" components=1 target_atoms="<<target_memory.shape.component_target_atoms
        <<" halo_atoms="<<target_memory.shape.component_halo_atoms
        <<" free_ac_columns="<<certification.assessment.rank_memory_estimate.shape.component_free_columns
        <<" design_columns="<<certification.assessment.rank_memory_estimate.shape.component_design_columns
        <<" width_columns="<<certification.assessment.rank_memory_estimate.shape.component_width_columns
        <<" target_columns="<<target_memory.shape.component_target_columns
        <<" target_evidence_atoms="<<target_memory.shape.component_target_atoms
        <<" rank_factor_dimensions="<<rank_stages[0].factor_rows<<"x"<<rank_stages[0].factor_columns
        <<" rank_factor_workspace_bytes="<<rank_stages[0].factor_workspace_bytes
        <<" sparse_design_bytes="<<rank_stages[0].sparse_design_bytes
        <<" projected_width_right_vector_bytes="<<rank_stages[0].right_singular_vector_bytes
        <<" primary_endpoint_seconds="<<certification.assessment.primary_endpoint_seconds
        <<" independent_replay_seconds="<<certification.assessment.independent_replay_seconds
        <<" rank_assessment_seconds="<<certification.assessment.rank_assessment_seconds
        <<" primary_stage_rss_before_bytes="<<(primary_stage->current_rss_before_bytes ?
            std::to_string(*primary_stage->current_rss_before_bytes) : "unavailable")
        <<" primary_stage_rss_after_bytes="<<(primary_stage->current_rss_after_bytes ?
            std::to_string(*primary_stage->current_rss_after_bytes) : "unavailable")
        <<" replay_stage_rss_before_bytes="<<(replay_stage->current_rss_before_bytes ?
            std::to_string(*replay_stage->current_rss_before_bytes) : "unavailable")
        <<" replay_stage_rss_after_bytes="<<(replay_stage->current_rss_after_bytes ?
            std::to_string(*replay_stage->current_rss_after_bytes) : "unavailable")
        <<" rank_preflight_seconds="<<rank_aggregate->seconds
        <<" rank_preflight_rss_before_bytes="<<(rank_aggregate->current_rss_before_bytes ?
            std::to_string(*rank_aggregate->current_rss_before_bytes) : "unavailable")
        <<" rank_preflight_rss_after_bytes="<<(rank_aggregate->current_rss_after_bytes ?
            std::to_string(*rank_aggregate->current_rss_after_bytes) : "unavailable")
        <<" target_preflight_seconds="<<target_preflight_seconds
        <<" target_preflight_rss_before_bytes="<<(target_rss_before ?
            std::to_string(*target_rss_before) : "unavailable")
        <<" target_preflight_rss_after_bytes="<<(target_rss_after ?
            std::to_string(*target_rss_after) : "unavailable")
        <<" target_preflight_process_peak_rss_bytes="<<(target_process_peak ?
            std::to_string(*target_process_peak) : "unavailable")
        <<" current_rss_bytes="<<(current_rss ? std::to_string(*current_rss) : "unavailable")
        <<" process_peak_rss_bytes="<<(process_peak_rss ? std::to_string(*process_peak_rss) : "unavailable")
        <<" reference_qr_factorization_calls="<<n::NumericsWorkForTesting().reference
        <<" endpoint_verification_status=Passed design_rank_status=ResourceLimit"
        <<" projected_width_status=NotRun"
        <<" numerical_identifiability_status=Unavailable target_identifiability_status=ResourceLimit"
        <<" runtime_convergence_status=Unavailable physical_input=not-measured"
        <<" endpoint_verification="<<(certification.verified()?"passed":"failed")
        <<" numerical_qualification="<<(certification.assessment.qualified?"passed":"unavailable")
        <<" reference_solve_calls="<<n::AssessmentWorkForTesting().reference_evaluations<<std::endl;
}
TEST(JointFixedNeighborBlockCoordinateTest, EndpointReplayIsReferenceFreeAndRankPreflightRemainsUnavailable)
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
    EXPECT_FALSE(certification.assessment.qualified);
    EXPECT_EQ(certification.assessment.rank_memory_estimate.stage,n::MemoryStage::RankAssessment);
    EXPECT_EQ(certification.assessment.rank_memory_estimate.status,n::MemoryStatus::Reject);
    EXPECT_FALSE(certification.assessment.design.has_value());
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
        result.assessment.qualified=false;
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
TEST(JointFixedNeighborBlockCoordinateTest, VerifiedWidthIdentifiabilityFailureRetriesAtBoundedInterval)
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
            result.assessment.identified=false;
            result.assessment.qualified=false;
            result.assessment.failure="width-unidentified";
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
    EXPECT_TRUE(search.assessment.qualified);
    ASSERT_EQ(search.fit.components.size(),1u);
    EXPECT_TRUE(search.fit.components.front().search_completed);
    EXPECT_TRUE(search.fit.components.front().state);
    EXPECT_EQ(search.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointFixedNeighborBlockCoordinateTest, RankResourceLimitStopsQualificationRetryAndPreservesCandidate)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    n::FixedNeighborSearchPolicy production_policy; production_policy.maximum_sweeps=40;
    std::size_t current_sweep{},attempts{};
    production_policy.configure_search_for_testing=[&](n::FixedNeighborPolicy & policy) {
        policy.sweep_observer=[&](const auto & sweep){current_sweep=sweep.sweep;};
        policy.certification_for_testing=[&](const auto & domain,n::VectorRef y,
            const n::BlockCoordinateState & state,const n::EvaluationContext & context,
            n::ProfileSearchWork * work) {
            ++attempts;
            auto result=n::CertifyFixedNeighborEndpoint(domain,y,state,context,work);
            result.assessment.rank_memory_estimate.status=n::MemoryStatus::Reject;
            result.assessment.rank_memory_estimate.reason="estimated-peak-exceeds-budget";
            result.assessment.qualified=false;
            result.assessment.failure="resource-limit:rank-assessment";
            return result;
        };
    };
    const auto fit=n::FitFixedNeighborComponents(problem,std::vector<double>(8,.55),production_policy);
    ASSERT_EQ(fit.components.size(),1u);
    const auto & component=fit.components.front();
    EXPECT_FALSE(fit.search_completed);
    EXPECT_FALSE(component.search_completed);
    EXPECT_EQ(component.stop_reason,"resource-limit:rank-assessment:estimated-peak-exceeds-budget");
    EXPECT_EQ(attempts,1u);
    EXPECT_EQ(current_sweep,component.confirmed_stationarity_sweep);
    EXPECT_EQ(component.endpoint_certification_sweep,component.confirmed_stationarity_sweep);
    ASSERT_TRUE(component.candidate_diagnostic);
    EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
    EXPECT_EQ(component.candidate_diagnostic->resource_limit_reason,
        "resource-limit:rank-assessment:estimated-peak-exceeds-budget");
    EXPECT_EQ(component.candidate_diagnostic->endpoint_verification_status,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(component.candidate_diagnostic->numerical_qualification_status,rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_FALSE(component.state);
    EXPECT_FALSE(fit.assembled_state);
    EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
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
        result.assessment.qualified=false;
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
        std::size_t production_design_spectrum_calls{};
    };
    const auto run=[&](std::size_t maximum_sweeps,bool stop_after_stationarity) {
        SearchRun result;
        n::AssessmentWorkForTesting()={};
        n::NumericsWorkForTesting().design_spectrum_calls=0;
        n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=maximum_sweeps;
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
        result.production_design_spectrum_calls=n::NumericsWorkForTesting().design_spectrum_calls;
        return result;
    };
    const auto qualified_run=run(200,true);
    const auto continued_run=run(127,false);
    EXPECT_EQ(qualified_run.production_reference_evaluations,0);
    EXPECT_EQ(continued_run.production_reference_evaluations,0);
    EXPECT_EQ(qualified_run.production_design_spectrum_calls,0u);
    EXPECT_EQ(continued_run.production_design_spectrum_calls,0u);
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
    for(const auto * name:{"conditional-ac-optimality","independent-replay","width-stationarity",
            "design-identifiability","numerical-identifiability"})
        EXPECT_EQ(evidence_status(name),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(component.target_evidence);
    const auto target_identifiability=std::find_if(component.target_evidence->checks.begin(),
        component.target_evidence->checks.end(),[](const auto & check){return check.name=="target-identifiability";});
    ASSERT_NE(target_identifiability,component.target_evidence->checks.end());
    EXPECT_EQ(target_identifiability->status,rhbm_gem::JointCheckStatus::Passed);

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
    ASSERT_TRUE(sweep127_certification.assessment.design && sweep127_certification.assessment.widths);
    EXPECT_EQ(sweep127_certification.assessment.design->rank,336);
    EXPECT_EQ(sweep127_certification.assessment.widths->rank,168);

    const auto & early=qualified_run.trajectory.back().state;
    const auto & early_assessment=qualified_run.certifications.front().assessment;
    ASSERT_TRUE(early_assessment.design && early_assessment.widths);
    EXPECT_EQ(early_assessment.design->rank,336);
    EXPECT_EQ(early_assessment.widths->rank,168);
    EXPECT_TRUE(early_assessment.identified);
    EXPECT_TRUE(early_assessment.qualified);
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
    RecordProperty("correction_baseline_168_old_design_rank","336/336");
    RecordProperty("correction_baseline_168_old_projected_width_rank","168/168");
    RecordProperty("correction_baseline_168_old_target_identifiability","passed");
    RecordProperty("correction_baseline_168_old_target_local_correction","3.1549147319454434e-11");
    RecordProperty("correction_free_168_max_amplitude_difference",std::to_string(amplitude_difference));
    RecordProperty("correction_free_168_max_charge_difference",std::to_string(charge_difference));
    RecordProperty("correction_free_168_max_width_difference",std::to_string(width_difference));
    RecordProperty("correction_free_168_objective_difference",std::to_string(objective_difference));
    RecordProperty("correction_free_168_residual_difference",std::to_string(residual_difference));
    RecordProperty("correction_free_168_ac_kkt",std::to_string(early_primary.certificate.projected_kkt));
    RecordProperty("correction_free_168_width_gradient",std::to_string(early_primary.gradient.lpNorm<Eigen::Infinity>()));
    RecordProperty("correction_free_168_replay_prediction_difference",std::to_string(early_replay.prediction_difference));
    RecordProperty("correction_free_168_replay_objective_difference",std::to_string(early_replay.objective_difference));
    RecordProperty("correction_free_168_design_rank",std::to_string(early_assessment.design->rank));
    RecordProperty("correction_free_168_projected_width_rank",std::to_string(early_assessment.widths->rank));
    EXPECT_EQ(component.candidate_diagnostic->endpoint_verification_status,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(component.state);
    EXPECT_TRUE(qualified_run.fit.assembled_state);
    const auto saved=rhbm_gem::core::CaptureJointAnalysisResult(qualified_run.fit);
    EXPECT_EQ(saved.components.front().runtime_convergence,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(saved.runtime_convergence,rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointFixedNeighborBlockCoordinateTest, SameEtaRawPrimaryAndReferenceDiagnosticsPreserveSearchState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(8,std::log(.55));
    const auto search=n::RunFixedNeighborDiagnostic(problem,initial_eta);
    ASSERT_TRUE(search.search_converged);
    const n::Vector eta=search.state.eta,beta=search.state.beta;
    const n::Domain domain=data.domain;
    const auto context=data.context;
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
    ASSERT_TRUE(raw.certificate.available && reference.certificate.available);
    ASSERT_TRUE(raw_trust.reference.has_value());
    ASSERT_TRUE(raw_assessment.design && raw_assessment.widths);
    EXPECT_EQ(raw_assessment.design->rank,raw_assessment.design->columns);
    EXPECT_TRUE(std::isfinite(raw_assessment.design->condition));
    EXPECT_LE(raw_assessment.design->condition,1e8);
    EXPECT_EQ(raw.certificate.feasible,reference.certificate.feasible);
    EXPECT_EQ(raw.certificate.active_atoms,reference.certificate.active_atoms);
    EXPECT_LE((raw.beta-reference.beta).lpNorm<Eigen::Infinity>(),1e-10*(1.+reference.beta.lpNorm<Eigen::Infinity>()));
    EXPECT_LE(std::abs(raw.certificate.objective-reference.certificate.objective),
        1e-12*(1.+std::abs(reference.certificate.objective)));
    EXPECT_LE((raw.residual-reference.residual).norm(),1e-10*(1.+reference.residual.norm()));
    EXPECT_LE(raw.certificate.projected_kkt,1e-10);
    EXPECT_LE(reference.certificate.projected_kkt,1e-10);
    EXPECT_LE(raw.gradient.lpNorm<Eigen::Infinity>(),1e-12);
    EXPECT_LE(reference.gradient.lpNorm<Eigen::Infinity>(),1e-12);
    EXPECT_LE(raw_trust.coefficient_difference,1e-10);
    EXPECT_LE(raw_trust.kkt_difference,1e-13);
    EXPECT_LE(raw_trust.prediction_difference,2e-12+2e-13*(raw.x*raw.beta).cwiseAbs().maxCoeff());
    EXPECT_LE(raw_trust.gradient_difference,1e-13+2e-9*raw.gradient.lpNorm<Eigen::Infinity>());
    EXPECT_EQ(raw_assessment.design->rank,raw.x.cols());
    EXPECT_EQ(raw_assessment.widths->rank,eta.size());
    EXPECT_TRUE(raw_assessment.inner);
    EXPECT_TRUE(raw_assessment.identified);
    EXPECT_TRUE(raw_assessment.gradient);
    EXPECT_TRUE(raw_trust.passed);
    const auto number=[](double value) {
        std::ostringstream out; out<<std::setprecision(17)<<value; return out.str();
    };
    RecordProperty("legacy_primary_ac_kkt",number(raw.certificate.projected_kkt));
    RecordProperty("legacy_primary_width_gradient",number(raw.gradient.lpNorm<Eigen::Infinity>()));
    RecordProperty("legacy_primary_feasibility",raw.certificate.feasible ? "true" : "false");
    RecordProperty("legacy_reference_ac_kkt",number(reference.certificate.projected_kkt));
    RecordProperty("legacy_reference_width_gradient",number(reference.gradient.lpNorm<Eigen::Infinity>()));
    RecordProperty("legacy_reference_coefficient_difference",number(raw_trust.coefficient_difference));
    RecordProperty("legacy_independent_prediction_difference",number(raw_trust.prediction_difference));
    RecordProperty("legacy_independent_kkt_difference",number(raw_trust.kkt_difference));
    RecordProperty("legacy_independent_gradient_difference",number(raw_trust.gradient_difference));
    RecordProperty("legacy_design_rank",std::to_string(raw_assessment.design->rank));
    RecordProperty("legacy_projected_width_rank",std::to_string(raw_assessment.widths->rank));
    RecordProperty("legacy_numerical_identifiability",raw_assessment.identified ? "passed" : "failed");
    RecordProperty("legacy_endpoint_trust",raw_trust.passed ? "passed" : "failed");
    RecordProperty("legacy_runtime_convergence",raw_assessment.identified && raw_assessment.inner &&
        raw_assessment.gradient && raw_trust.passed ? "passed" : "failed");
    for(const auto * assessment:{&raw_assessment,&primary_assessment,&reference_assessment})
    {
        EXPECT_TRUE(assessment->primary.valid);
        EXPECT_TRUE(assessment->reference.valid);
        EXPECT_TRUE(assessment->design.has_value());
        EXPECT_TRUE(assessment->widths.has_value());
        EXPECT_EQ(assessment->gradient,assessment->primary.gradient.lpNorm<Eigen::Infinity>()<=1e-12 &&
            assessment->reference.gradient.lpNorm<Eigen::Infinity>()<=1e-12);
    }
    EXPECT_TRUE(raw_trust.primary_valid); EXPECT_TRUE(primary_trust.primary_valid); EXPECT_TRUE(reference_trust.primary_valid);
}
