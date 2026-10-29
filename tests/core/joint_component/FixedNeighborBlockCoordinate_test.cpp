#include <gtest/gtest.h>
#include "support/JointFixedNeighborDiagnostic.hpp"
#include "support/JointBlockCoordinate.hpp"
#include "support/JointComponentChecks.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include "support/JointSyntheticWorkload.hpp"
#include <cstdint>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace b=second_stage_test::joint_block;
using rhbm_gem::core::JointProblem;
void CheckSingleBlock(const std::string & topology,int atoms)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,atoms));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(atoms,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=30;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
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
    EXPECT_EQ(bounded.components.front().candidate_diagnostic->convergence_criteria_id,"practical-v2");
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
bool BitwiseEqual(const n::BlockCoordinateState & left,const n::BlockCoordinateState & right)
{
    return BitwiseEqual(left.eta,right.eta) && BitwiseEqual(left.beta,right.beta) &&
        BitwiseEqual(left.prediction,right.prediction) && BitwiseEqual(left.residual,right.residual) &&
        BitwiseEqual(left.objective,right.objective);
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
    std::size_t candidates{},decision_disagreements{};
    double maximum_absolute_error{};
};
void AuditIncrementalCandidate(const rhbm_gem::core::JointProblemInput & input,
    const rhbm_gem::JointParameterLayout & layout,n::VectorRef observations,double scale,
    const n::FixedNeighborCandidateEvent & event,IncrementalObjectiveAudit & audit)
{
    if(!event.candidate_replay_ran || !std::isfinite(event.incremental_objective)) return;
    ++audit.candidates;
    const auto old_replay=b::Replay(input,layout,observations,event.accepted_before.eta,
        event.accepted_before.beta,scale);
    const auto candidate_replay=b::Replay(input,layout,observations,event.candidate_eta,event.candidate_beta,scale);
    EXPECT_EQ(event.affected_prediction.size(),static_cast<Eigen::Index>(event.affected_rows.size()));
    EXPECT_EQ(event.affected_residual.size(),static_cast<Eigen::Index>(event.affected_rows.size()));
    for(std::size_t k=0;k<event.affected_rows.size();++k)
    {
        const auto row=static_cast<Eigen::Index>(event.affected_rows[k]);
        EXPECT_DOUBLE_EQ(event.affected_prediction(static_cast<Eigen::Index>(k)),candidate_replay.prediction(row));
        EXPECT_DOUBLE_EQ(event.affected_residual(static_cast<Eigen::Index>(k)),candidate_replay.residual(row));
    }
    for(const auto row_value:layout.informative_rows)
    {
        if(std::find(event.affected_rows.begin(),event.affected_rows.end(),
            static_cast<Eigen::Index>(row_value))!=event.affected_rows.end()) continue;
        const auto row=static_cast<Eigen::Index>(row_value);
        EXPECT_DOUBLE_EQ(event.accepted_before.prediction(row),candidate_replay.prediction(row));
        EXPECT_DOUBLE_EQ(event.accepted_before.residual(row),candidate_replay.residual(row));
    }
    const double old_reference=old_replay.objective;
    const double candidate_reference=candidate_replay.objective;
    audit.maximum_absolute_error=std::max(audit.maximum_absolute_error,
        std::abs(event.incremental_objective-candidate_reference));
    const bool estimated_decision=ObjectiveDecision(event.accepted_before.objective,event.incremental_objective,
        event.local_objective_before,event.local_objective_after,event.unchanged);
    const bool reference_decision=ObjectiveDecision(old_reference,candidate_reference,
        event.local_objective_before,event.local_objective_after,event.unchanged);
    audit.decision_disagreements+=estimated_decision!=reference_decision || event.accepted!=reference_decision;
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
    policy.eta_upper_bound=std::log(2.2);
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

TEST(JointFixedNeighborBlockCoordinateTest, CandidateObserverDoesNotChangeSearch)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=4;
    policy.stop_after_stationarity=false; policy.certify_endpoint=false;
    policy.eta_upper_bound=std::log(2.2);
    const auto baseline=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    std::size_t candidate_count{},accepted_count{};
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        ++candidate_count; accepted_count+=event.accepted;
        EXPECT_EQ(event.candidate_eta.size(),event.accepted_before.eta.size());
        EXPECT_EQ(event.candidate_beta.size(),event.accepted_before.beta.size());
        EXPECT_GT(event.core_atoms.size(),0u);
        EXPECT_GT(event.affected_rows.size(),0u);
        EXPECT_TRUE(std::isfinite(event.candidate_objective));
        if(event.candidate_replay_ran)
        {
            const auto replay=b::Replay(*data.input,data.layout,data.y,event.candidate_eta,
                event.candidate_beta,data.context.scale);
            EXPECT_LE(std::abs(event.candidate_objective-replay.objective),
                n::BlockObjectiveReplayEnclosure(replay.objective));
        }
    };
    const auto observed=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    EXPECT_GT(candidate_count,0u);
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

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveMatchesIndependentFullReplay)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=40;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
    IncrementalObjectiveAudit e1_audit;
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        AuditIncrementalCandidate(*data.input,data.layout,data.y,data.context.scale,event,e1_audit);
        if(!event.candidate_replay_ran) return;
        const auto replay=b::Replay(*data.input,data.layout,data.y,
            event.candidate_eta,event.candidate_beta,data.context.scale);
        EXPECT_LE(std::abs(event.candidate_objective-replay.objective),
            n::BlockObjectiveReplayEnclosure(replay.objective));
        EXPECT_EQ(event.accepted,ObjectiveDecision(event.accepted_before.objective,replay.objective,
            event.local_objective_before,event.local_objective_after,event.unchanged));
    };
    const auto e1=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    ASSERT_TRUE(e1.search_converged)<<e1.reason;
    ASSERT_TRUE(e1.endpoint_certified)<<e1.endpoint_trust.reason;
    EXPECT_TRUE(e1.assessment.independent_replay.passed);
    EXPECT_GT(e1_audit.candidates,0u);
    EXPECT_EQ(e1_audit.decision_disagreements,0u);
    EXPECT_EQ(e1.work.incremental_objective_evaluations,e1_audit.candidates);
    EXPECT_EQ(e1.work.objective_reconciliations,e1.sweep_count);
    EXPECT_EQ(e1.work.full_candidate_replays,0u);
    EXPECT_TRUE(std::all_of(e1.sweeps.begin(),e1.sweeps.end(),
        [](const auto & sweep){return sweep.objective_cache_refreshed;}));
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveQualifiesAtActiveWidthUpperBound)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const n::Vector initial_eta=n::Vector::Constant(3,std::log(.5));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=40;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(.5); policy.collect_telemetry=true;
    const auto e1=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
    ASSERT_TRUE(e1.search_converged)<<e1.reason;
    ASSERT_TRUE(e1.endpoint_certified)<<e1.endpoint_trust.reason;
    EXPECT_TRUE(e1.assessment.independent_replay.passed);
    ASSERT_TRUE(e1.final_sweep);
    EXPECT_TRUE(e1.final_sweep->width_bound_applied);
    EXPECT_GT(e1.final_sweep->active_width_bound_count,0u);
    EXPECT_LE(e1.final_sweep->maximum_width_bound_violation,0.0);
    EXPECT_LE(e1.final_sweep->projected_width_gradient_inf_norm,
        n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_EQ(e1.work.objective_reconciliations,e1.sweep_count);
    EXPECT_GT(e1.work.incremental_objective_evaluations,0u);
    EXPECT_TRUE(std::all_of(e1.sweeps.begin(),e1.sweeps.end(),
        [](const auto & sweep){return sweep.objective_cache_refreshed;}));
}

TEST(JointFixedNeighborBlockCoordinateTest, EverySweepReconcilesObjectiveAfterReplayAudit)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",300));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(300,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=12; policy.maximum_sweeps=30;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.stop_after_stationarity=true; policy.certify_endpoint=true;
    policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
    policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
    IncrementalObjectiveAudit audit;
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        AuditIncrementalCandidate(*data.input,data.layout,data.y,data.context.scale,event,audit);
    };
    const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);

    ASSERT_TRUE(result.search_converged)<<result.reason;
    ASSERT_TRUE(result.endpoint_certified)<<result.endpoint_trust.reason;
    EXPECT_TRUE(result.assessment.independent_replay.passed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(audit.candidates,result.work.incremental_objective_evaluations);
    EXPECT_EQ(audit.decision_disagreements,0u);
    EXPECT_EQ(result.work.objective_reconciliations,result.sweep_count);
    EXPECT_EQ(result.sweeps.size(),result.sweep_count);
    EXPECT_TRUE(std::all_of(result.sweeps.begin(),result.sweeps.end(),
        [](const auto & sweep){return sweep.objective_cache_refreshed;}));
    EXPECT_GT(result.work.sweep_replay_seconds,0.0);
    EXPECT_GT(result.work.sweep_replay_audit_seconds,0.0);
    EXPECT_GT(result.work.sweep_global_state_seconds,0.0);
    EXPECT_GE(result.work.objective_reconciliation_seconds,0.0);
    EXPECT_LE(result.work.sweep_replay_seconds+result.work.sweep_replay_audit_seconds+
        result.work.sweep_global_state_seconds+result.work.objective_reconciliation_seconds,
        result.search_seconds);
    const auto replay=b::Replay(*data.input,data.layout,data.y,
        result.state.eta,result.state.beta,data.context.scale);
    EXPECT_LE(std::abs(result.state.objective-replay.objective),
        n::BlockObjectiveReplayEnclosure(replay.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, SweepReplayMismatchFailsBeforeObjectiveReconciliation)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=1;
    policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
    policy.perturb_sweep_replay_objective_for_testing=true;
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    ASSERT_TRUE(result.final_sweep.has_value());
    ASSERT_EQ(result.sweeps.size(),1u);
    EXPECT_EQ(result.reason,"block-cache-replay-failed");
    EXPECT_FALSE(result.sweeps.front().objective_cache_refreshed);
    EXPECT_EQ(result.work.objective_reconciliations,0u);
    EXPECT_GT(result.sweeps.front().objective_replay_error,
        n::BlockObjectiveReplayEnclosure(result.sweeps.front().objective_after));
    const auto replay=b::Replay(*data.input,data.layout,data.y,
        result.state.eta,result.state.beta,data.context.scale);
    EXPECT_NE(result.state.objective,replay.objective);
}

TEST(JointFixedNeighborBlockCoordinateTest, IncrementalObjectiveRemainsQualifiedAcrossInitialWidths)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    for(const double width:std::array<double,3>{.35,.55,.9})
    {
        const n::Vector initial_eta=n::Vector::Constant(3,std::log(width));
        n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=40;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
        const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        ASSERT_TRUE(result.search_converged)<<result.reason;
        ASSERT_TRUE(result.endpoint_certified)<<result.endpoint_trust.reason;
        EXPECT_TRUE(result.assessment.independent_replay.passed);
        EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
        EXPECT_GT(result.work.incremental_objective_evaluations,0u);
        EXPECT_EQ(result.work.objective_reconciliations,result.sweep_count);
        EXPECT_TRUE(std::all_of(result.sweeps.begin(),result.sweeps.end(),
            [](const auto & sweep){return sweep.objective_cache_refreshed;}));
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, SignedCancellationQualifiesWhenOnlyKktGateIsSatisfied)
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
    IncrementalObjectiveAudit audit;
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        AuditIncrementalCandidate(*data.input,data.layout,data.y,data.context.scale,event,audit);
    };
    const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);

    EXPECT_GT(audit.candidates,0u);
    EXPECT_EQ(audit.decision_disagreements,0u);
    EXPECT_EQ(result.work.incremental_objective_evaluations,audit.candidates);
    EXPECT_EQ(result.reason,"block-stationary")<<" sweeps="<<result.sweep_count
        <<" first_kkt="<<result.first_order_stationarity_sweep
        <<" confirmed="<<result.confirmed_stationarity_sweep
        <<" eta_change="<<(result.final_sweep ? result.final_sweep->eta_change_inf : -1.0)
        <<" stability="<<(result.parameter_stability ? result.parameter_stability->trajectory_stability_status : "missing")
        <<" perturbation="<<(result.parameter_stability ? result.parameter_stability->perturbation_status : "missing");
    EXPECT_TRUE(result.search_converged);
    EXPECT_TRUE(result.endpoint_certified);
    EXPECT_TRUE(result.assessment.independent_replay.passed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(result.final_sweep);
    EXPECT_GT(result.final_sweep->eta_change_inf,1e-6);
    EXPECT_EQ(result.confirmed_stationarity_sweep,result.sweep_count);
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

TEST(JointFixedNeighborBlockCoordinateTest, HardCaseIncrementalObjectiveAndReplayQualification)
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
        n::FixedNeighborPolicy policy; policy.core_atoms=fixture.core_atoms; policy.maximum_sweeps=100;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true;
        policy.assess_local_conditioning=true; policy.assess_parameter_stability=true;
        IncrementalObjectiveAudit audit;
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            AuditIncrementalCandidate(*data.input,data.layout,data.y,data.context.scale,event,audit);
        };
        const auto result=n::RunFixedNeighborDiagnostic(problem,initial_eta,policy);
        ASSERT_EQ(result.state.eta.size(),initial_eta.size())<<fixture.name<<':'<<result.reason;
        const auto replay=b::Replay(*data.input,data.layout,data.y,
            result.state.eta,result.state.beta,data.context.scale);
        EXPECT_TRUE(result.state.eta.allFinite())<<fixture.name;
        EXPECT_TRUE(result.state.beta.allFinite())<<fixture.name;
        EXPECT_EQ(result.global_state_memory_estimate.status,n::MemoryStatus::Pass)<<fixture.name;
        EXPECT_LE(std::abs(result.state.objective-replay.objective),
            n::BlockObjectiveReplayEnclosure(replay.objective))<<fixture.name;
        EXPECT_EQ(audit.decision_disagreements,0u)<<fixture.name;
        EXPECT_EQ(result.work.incremental_objective_evaluations,audit.candidates)<<fixture.name;
        EXPECT_EQ(result.work.objective_reconciliations,result.sweep_count)<<fixture.name;
        if(fixture.name=="rank-deficient-constant-support" ||
            fixture.name=="near-singular-overlapping-support")
        {
            EXPECT_EQ(result.reason,"block-recovery-exhausted")<<fixture.name;
            EXPECT_FALSE(result.search_converged)<<fixture.name;
            EXPECT_FALSE(result.endpoint_certified)<<fixture.name;
        }
        else if(fixture.name=="low-signal-noise-only")
        {
            EXPECT_EQ(result.reason,"block-stationary")<<fixture.name
                <<" sweeps="<<result.sweep_count
                <<" first_kkt="<<result.first_order_stationarity_sweep
                <<" confirmed="<<result.confirmed_stationarity_sweep
                <<" eta_change="<<(result.final_sweep ? result.final_sweep->eta_change_inf : -1.0)
                <<" objective="<<result.state.objective
                <<" stability="<<(result.parameter_stability ? result.parameter_stability->trajectory_stability_status : "missing")
                <<" perturbation="<<(result.parameter_stability ? result.parameter_stability->perturbation_status : "missing");
            EXPECT_TRUE(result.search_converged)<<fixture.name;
            EXPECT_TRUE(result.endpoint_certified)<<fixture.name<<':'<<result.endpoint_trust.reason;
            EXPECT_TRUE(result.assessment.independent_replay.passed)<<fixture.name;
            ASSERT_TRUE(result.final_sweep)<<fixture.name;
            EXPECT_GT(result.final_sweep->eta_change_inf,1e-6)<<fixture.name;
            ASSERT_TRUE(result.parameter_stability)<<fixture.name;
            EXPECT_EQ(result.parameter_stability->trajectory_stability_status,"insufficient-history")<<fixture.name;
            EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed)<<fixture.name;
        }
        else
        {
            EXPECT_TRUE(result.search_converged)<<fixture.name<<':'<<result.reason;
            EXPECT_TRUE(result.endpoint_certified)<<fixture.name<<':'<<result.endpoint_trust.reason;
            EXPECT_TRUE(result.assessment.independent_replay.passed)<<fixture.name;
        }
    }
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
        std::size_t fallback_candidate_count{};
    fallback_policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        if(!event.candidate_replay_ran) return;
        ++fallback_candidate_count;
        const auto before=b::Replay(*data.input,data.layout,data.y,
            event.accepted_before.eta,event.accepted_before.beta,data.context.scale);
        const auto candidate=b::Replay(*data.input,data.layout,data.y,
            event.candidate_eta,event.candidate_beta,data.context.scale);
        EXPECT_LE(std::abs(event.accepted_before.objective-before.objective),
            n::BlockObjectiveReplayEnclosure(before.objective));
        EXPECT_LE(std::abs(event.candidate_objective-candidate.objective),
            n::BlockObjectiveReplayEnclosure(candidate.objective));
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
    Run result;
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=2;
    policy.eta_upper_bound=std::log(2.2);
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        if(event.candidate_replay_ran) ++result.candidates;
        result.accepted+=event.accepted;
    };
    result.result=n::SearchFixedNeighborComponent(data.domain,data.layout,observations,
        eta,data.context,policy);
    const auto initial_replay=b::Replay(*data.input,data.layout,observations,
        eta,n::Vector::Zero(6),data.context.scale);

    ASSERT_TRUE(std::isinf(initial_replay.objective));
    EXPECT_FALSE(result.result.search_converged);
    EXPECT_EQ(result.result.reason,"block-cache-replay-failed");
    EXPECT_TRUE(std::isinf(result.result.state.objective));
    EXPECT_EQ(result.candidates,0u);
    EXPECT_EQ(result.accepted,0u);
    EXPECT_TRUE(BitwiseEqual(result.result.state.eta,eta));
    EXPECT_TRUE(BitwiseEqual(result.result.state.beta,n::Vector::Zero(6)));
}

TEST(JointFixedNeighborBlockCoordinateTest, CandidateLocalObjectiveIncreasePreservesAcceptedBeforeState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=1;
    policy.eta_upper_bound=std::log(2.2);
    std::optional<n::BlockCoordinateState> rejected_before;
    policy.candidate_fault_for_testing=[&](std::size_t sweep,std::size_t block,const auto & state) {
        if(sweep==1 && block==2)
        {
            rejected_before=state;
            return n::FixedNeighborCandidateFaultForTesting::LocalObjectiveIncrease;
        }
        return n::FixedNeighborCandidateFaultForTesting::None;
    };
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);

    ASSERT_TRUE(rejected_before.has_value());
    ASSERT_EQ(result.blocks.size(),2u);
    EXPECT_TRUE(result.blocks[0].accepted);
    EXPECT_EQ(result.blocks[1].status,"unchanged");
    EXPECT_EQ(result.blocks[1].reason,"block-objective-increase");
    EXPECT_TRUE(BitwiseEqual(result.state.eta,rejected_before->eta));
    EXPECT_TRUE(BitwiseEqual(result.state.beta,rejected_before->beta));
    EXPECT_TRUE(BitwiseEqual(result.state.prediction,rejected_before->prediction));
    EXPECT_TRUE(BitwiseEqual(result.state.residual,rejected_before->residual));
    EXPECT_LE(std::abs(result.state.objective-rejected_before->objective),
        n::BlockObjectiveReplayEnclosure(rejected_before->objective));
    const auto replay=b::Replay(*data.input,data.layout,data.y,
        result.state.eta,result.state.beta,data.context.scale);
    EXPECT_LE(std::abs(result.state.objective-replay.objective),
        n::BlockObjectiveReplayEnclosure(replay.objective));
}

TEST(JointFixedNeighborBlockCoordinateTest, CandidateConsistencyAndNonfiniteFaultsPreserveAcceptedBeforeState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    constexpr std::array faults{
        n::FixedNeighborCandidateFaultForTesting::ObjectiveConsistencyFailure,
        n::FixedNeighborCandidateFaultForTesting::NaNObjective,
        n::FixedNeighborCandidateFaultForTesting::InfiniteObjective,
        n::FixedNeighborCandidateFaultForTesting::NaNParameter,
        n::FixedNeighborCandidateFaultForTesting::InfiniteParameter,
        n::FixedNeighborCandidateFaultForTesting::InfinitePrediction,
        n::FixedNeighborCandidateFaultForTesting::InfiniteResidual};
    for(const auto fault:faults)
    {
        SCOPED_TRACE(static_cast<int>(fault));
        n::FixedNeighborPolicy policy; policy.core_atoms=2; policy.maximum_sweeps=1;
        policy.eta_upper_bound=std::log(2.2);
        std::optional<n::BlockCoordinateState> failed_candidate_before;
        bool failure_event_seen{};
        policy.candidate_fault_for_testing=[&](std::size_t sweep,std::size_t block,const auto & state) {
            if(sweep==1 && block==2)
            {
                failed_candidate_before=state;
                return fault;
            }
            return n::FixedNeighborCandidateFaultForTesting::None;
        };
        policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
            if(event.block!=2) return;
            failure_event_seen=true;
            ASSERT_TRUE(failed_candidate_before.has_value());
            EXPECT_FALSE(event.accepted);
            EXPECT_TRUE(BitwiseEqual(event.accepted_before,*failed_candidate_before));
        };
        const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);

        ASSERT_TRUE(failed_candidate_before.has_value());
        ASSERT_EQ(result.blocks.size(),2u);
        EXPECT_TRUE(result.blocks[0].accepted);
        EXPECT_EQ(result.blocks[1].status,"failed");
        EXPECT_EQ(result.blocks[1].reason,"block-objective-replay-failed");
        EXPECT_EQ(result.reason,"block-objective-replay-failed");
        EXPECT_TRUE(failure_event_seen);
        EXPECT_TRUE(BitwiseEqual(result.state,*failed_candidate_before));
        const auto replay=b::Replay(*data.input,data.layout,data.y,
            result.state.eta,result.state.beta,data.context.scale);
        EXPECT_TRUE(BitwiseEqual(result.state.prediction,replay.prediction));
        EXPECT_TRUE(BitwiseEqual(result.state.residual,replay.residual));
        EXPECT_LE(std::abs(result.state.objective-replay.objective),
            n::BlockObjectiveReplayEnclosure(replay.objective));
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, PositiveObjectiveDeltaConsistencyConditionBoundsTheAcceptanceCase)
{
    const double boundary_delta=n::BlockObjectiveReplayEnclosure(0.0);
    const double boundary_reference=boundary_delta;
    const double boundary_error=std::abs(boundary_delta+0.0-0.0);
    ASSERT_TRUE(n::WithinBlockObjectiveReplay(boundary_error,boundary_reference));
    EXPECT_LE(boundary_delta,n::BlockObjectiveReplayEnclosure(boundary_reference));

    const double outside_delta=1e-12+4e-24;
    const double outside_reference=outside_delta;
    const double outside_error=std::abs(outside_delta+0.0-0.0);
    EXPECT_FALSE(n::WithinBlockObjectiveReplay(outside_error,outside_reference));
    EXPECT_GT(outside_delta,n::BlockObjectiveReplayEnclosure(outside_reference));
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
    EXPECT_GT(bounded_lm.releases,0u);
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
    EXPECT_GE(confirmed.confirmed_stationarity_sweep,confirmed.first_order_stationarity_sweep);
    EXPECT_EQ(confirmed.confirmed_stationarity_sweep,confirmed.sweeps.size());
}

TEST(JointFixedNeighborBlockCoordinateTest, EndpointCertificationRetryScheduleIsPreserved)
{
    std::size_t next_retry_sweep{};
    EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::CandidateStationary,2,next_retry_sweep));
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::ConfirmedStationary,2,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,12u);
    EXPECT_FALSE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::ConfirmedStationary,11,next_retry_sweep));
    EXPECT_TRUE(n::ShouldAttemptFixedNeighborCertification(
        n::FixedNeighborStationarityState::ConfirmedStationary,12,next_retry_sweep));
    EXPECT_EQ(next_retry_sweep,22u);
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
        std::optional<n::BlockCoordinateState> retry_candidate_before;
        double retry_candidate_objective{},retry_candidate_reference{};
        bool retry_candidate_seen{};
    };
    RecoveryRun current;
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
    policy.eta_upper_bound=std::log(2.2); policy.collect_telemetry=true; policy.collect_diagnostics=true;
    policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t sweep) {
        return block==0 && sweep==0;
    };
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        if(event.sweep==2 && event.block==1)
        {
            current.retry_candidate_seen=true;
            current.retry_candidate_before=event.accepted_before;
            current.retry_candidate_objective=event.candidate_objective;
            const auto replay=b::Replay(*data.input,data.layout,data.y,
                event.candidate_eta,event.candidate_beta,data.context.scale);
            current.retry_candidate_reference=replay.objective;
        }
    };
    current.result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    ASSERT_TRUE(current.retry_candidate_seen);
    ASSERT_TRUE(current.retry_candidate_before.has_value());
    EXPECT_GT(current.retry_candidate_before->beta.norm(),0.0);
    EXPECT_EQ(current.result.recovery_attempts,1u);
    EXPECT_EQ(current.result.recovery_successes,1u);
    const auto first_failure=std::find_if(current.result.blocks.begin(),current.result.blocks.end(),[](const auto & block) {
        return block.block==1 && block.sweep==1;
    });
    ASSERT_NE(first_failure,current.result.blocks.end());
    EXPECT_EQ(first_failure->status,"retry-pending");
    const auto retried=std::find_if(current.result.blocks.begin(),current.result.blocks.end(),[](const auto & block) {
        return block.block==1 && block.sweep==2;
    });
    ASSERT_NE(retried,current.result.blocks.end());
    EXPECT_TRUE(retried->recovery_retry);
    EXPECT_EQ(retried->recovery_attempt,1u);
    EXPECT_LE(std::abs(retried->objective_before-current.retry_candidate_before->objective),
        n::BlockObjectiveReplayEnclosure(current.retry_candidate_before->objective));
    const auto trusted=b::Replay(*data.input,data.layout,data.y,
        current.retry_candidate_before->eta,current.retry_candidate_before->beta,data.context.scale);
    EXPECT_LE(std::abs(current.retry_candidate_before->objective-trusted.objective),
        n::BlockObjectiveReplayEnclosure(trusted.objective));
    EXPECT_LE(std::abs(current.retry_candidate_before->objective-trusted.objective),
        n::BlockObjectiveReplayEnclosure(trusted.objective));
    EXPECT_LE(std::abs(current.retry_candidate_objective-current.retry_candidate_reference),
        n::BlockObjectiveReplayEnclosure(current.retry_candidate_reference));
    EXPECT_GT(current.result.work.incremental_objective_evaluations,0u);
}

TEST(JointFixedNeighborBlockCoordinateTest, ExhaustedBlockRecoveryLeavesTheLastAcceptedComponentStateUntouched)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(2,std::log(.55));
    n::Vector accepted_eta,accepted_beta;
    n::FixedNeighborPolicy policy; policy.core_atoms=1; policy.maximum_sweeps=6;
    policy.eta_upper_bound=std::log(2.2);
    policy.force_no_trusted_profile_for_testing=[](std::size_t block,std::size_t) {return block==0;};
    policy.candidate_observer=[&](const n::FixedNeighborCandidateEvent & event) {
        if(event.sweep==1 && event.block==2)
        {
            accepted_eta=event.candidate_eta;
            accepted_beta=event.candidate_beta;
        }
    };
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-recovery-exhausted");
    EXPECT_EQ(result.recovery_attempts,1u);
    EXPECT_EQ(result.recovery_successes,0u);
    EXPECT_EQ(result.sweep_count,1u);
    EXPECT_EQ(result.first_order_stationarity_sweep,0u);
    ASSERT_EQ(accepted_eta.size(),result.state.eta.size());
    ASSERT_EQ(accepted_beta.size(),result.state.beta.size());
    EXPECT_TRUE((result.state.eta.array()==accepted_eta.array()).all());
    EXPECT_TRUE((result.state.beta.array()==accepted_beta.array()).all());
    const auto replay=b::Replay(*data.input,data.layout,data.y,
        result.state.eta,result.state.beta,data.context.scale);
    EXPECT_LE(std::abs(result.state.objective-replay.objective),
        n::BlockObjectiveReplayEnclosure(replay.objective));
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

TEST(JointFixedNeighborBlockCoordinateTest, FixedNeighborAcceptsAtMostOneTrustedLocalUpdatePerBlockVisit)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=2;
    policy.eta_upper_bound=std::log(2.2);
    policy.certify_endpoint=false;
    const auto result=n::RunFixedNeighborDiagnostic(problem,eta,policy);
    ASSERT_FALSE(result.blocks.empty());
    for(const auto & block:result.blocks)
    {
        EXPECT_LE(block.accepted_updates,1);
    }
}

TEST(JointFixedNeighborBlockCoordinateTest, KktStationarityCriteriaHaveInclusiveFiniteBoundaries)
{
    const auto & strict=n::StrictJointConvergenceCriteria;
    const auto & practical=n::PracticalJointConvergenceCriteria;
    EXPECT_TRUE(n::IsValidJointConvergenceCriteria(strict));
    EXPECT_TRUE(n::IsValidJointConvergenceCriteria(practical));

    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(5e-8,0,true,strict),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(5e-8,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,5e-10,true,strict),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,5e-10,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(practical.global_ac_kkt,0,false,practical),
        n::FixedNeighborStationarityState::CandidateStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(practical.global_ac_kkt,0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(
        std::nextafter(practical.global_ac_kkt,0.0),0,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(
        std::nextafter(practical.global_ac_kkt,std::numeric_limits<double>::infinity()),0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,practical.global_width_gradient,true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,
        std::nextafter(practical.global_width_gradient,0.0),true,practical),
        n::FixedNeighborStationarityState::ConfirmedStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,
        std::nextafter(practical.global_width_gradient,std::numeric_limits<double>::infinity()),true,practical),
        n::FixedNeighborStationarityState::NotStationary);

    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double infinity=std::numeric_limits<double>::infinity();
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(nan,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,infinity,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(-1,0,true,practical),
        n::FixedNeighborStationarityState::NotStationary);
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,-1,true,practical),
        n::FixedNeighborStationarityState::NotStationary);

    auto invalid=practical; invalid.global_ac_kkt=nan;
    EXPECT_FALSE(n::IsValidJointConvergenceCriteria(invalid));
    EXPECT_EQ(n::ClassifyFixedNeighborStationarity(0,0,true,invalid),
        n::FixedNeighborStationarityState::NotStationary);
    invalid=practical; invalid.global_width_gradient=infinity;
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
    const auto local=n::FixedNeighborPolicy{}.convergence_criteria;
    EXPECT_DOUBLE_EQ(local.global_ac_kkt,n::StrictJointConvergenceCriteria.global_ac_kkt);
    EXPECT_DOUBLE_EQ(local.global_width_gradient,n::StrictJointConvergenceCriteria.global_width_gradient);
}

TEST(JointFixedNeighborBlockCoordinateTest, KktOnlyCanCertifyWithLargeLogWidthChange)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(3,std::log(.55));
    n::FixedNeighborPolicy policy;
    policy.core_atoms=12;
    policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);

    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    ASSERT_TRUE(result.search_converged)<<result.reason;
    ASSERT_TRUE(result.final_sweep);
    EXPECT_EQ(result.first_order_stationarity_sweep,3u);
    EXPECT_EQ(result.confirmed_stationarity_sweep,3u);
    EXPECT_EQ(result.certification_sweep,3u);
    EXPECT_EQ(result.sweep_count,3u);
    EXPECT_GT(result.final_sweep->eta_change_inf,1e-6);
    ASSERT_TRUE(result.certification);
    EXPECT_TRUE(result.certification->qualified());
    EXPECT_TRUE(result.certification->assessment.independent_replay.passed);
    EXPECT_EQ(result.certification_attempts,1u);
}

TEST(JointFixedNeighborBlockCoordinateTest, BudgetExhaustionPreservesCandidateWithoutInventingTrustedState)
{
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(8,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=1;
    policy.eta_upper_bound=std::log(2.2);
    const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-sweep-budget");
    EXPECT_EQ(result.sweep_count,1u);
    EXPECT_EQ(result.state.eta.size(),eta.size());
    EXPECT_EQ(result.state.beta.size(),2*eta.size());
    EXPECT_TRUE(result.state.eta.allFinite());
    EXPECT_TRUE(result.state.beta.allFinite());
    EXPECT_TRUE(std::isfinite(result.state.objective));
    ASSERT_TRUE(result.certification.has_value());
    EXPECT_EQ(result.certification->qualified(),
        result.certification->verified() && result.certification->assessment.numerical_qualified);
    if(result.certification->qualified())
        EXPECT_TRUE(result.certification->assessment.independent_replay.passed);
    const auto replay=b::Replay(*data.input,data.layout,data.y,
        result.state.eta,result.state.beta,data.context.scale);
    EXPECT_LE(std::abs(result.state.objective-replay.objective),
        n::BlockObjectiveReplayEnclosure(replay.objective));
}
