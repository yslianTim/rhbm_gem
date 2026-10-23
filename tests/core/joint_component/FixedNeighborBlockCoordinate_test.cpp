#include <gtest/gtest.h>
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <algorithm>
#include "support/JointSyntheticWorkload.hpp"
#include <cmath>
#include <limits>
#include <numeric>
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
    policy.maximum_sweeps=200;
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
    EXPECT_GT(raw_inf,policy.convergence_criteria.global_width_gradient);
    EXPECT_NEAR(*stationarity->value,projected_inf,1e-12);
    EXPECT_LE(*stationarity->value,policy.convergence_criteria.global_width_gradient);
    EXPECT_EQ(stationarity->status,rhbm_gem::JointCheckStatus::Passed);
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

TEST(JointFixedNeighborBlockCoordinateTest, StabilityHistoryKeepsOnlyTheLastThreeCompletedSweeps)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    for(std::size_t sweep_budget=0;sweep_budget<=5;++sweep_budget)
    {
        n::FixedNeighborPolicy policy;
        policy.maximum_sweeps=sweep_budget;
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

    const auto fit=rhbm_gem::core::FitJointComponents(problem,{.55,.55});
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
        const auto result=n::SearchFixedNeighborComponent(data.domain,layout,data.y,eta,data.context);
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
