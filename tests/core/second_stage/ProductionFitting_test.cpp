#include <gtest/gtest.h>
#include "support/SecondStageTestSupport.hpp"
#include "support/SolverFailureCapture.hpp"
#include "core/detail/second_stage/FixedPointRecovery.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <cstdlib>
#include <cmath>
#include <filesystem>

namespace {
namespace d = rhbm_gem::core::detail;
namespace rg = rhbm_gem;

second_stage_test::JointPolishFixture Fixture()
{
    auto fixture{ second_stage_test::BuildJointPolishFixture(
        { {6.0, 0.55, 0.12}, {7.0, 0.52, -0.2} },
        { {6.0, 0.5, 0.1}, {7.0, 0.5, -0.1} }) };
    for (auto & atom : fixture.context.atom_list)
    {
        // Nonzero observational noise avoids a zero-variance solver boundary.
        for (std::size_t i = 0; i < atom.raw_sampling_entries.size(); ++i)
            atom.raw_sampling_entries[i].response += 0.001 * std::sin(1.0 + static_cast<double>(i));
        atom.refit_design = d::PreparedLocalGaussianDesign(atom.raw_sampling_entries, 0.0, 0.7);
    }
    return fixture;
}

d::FixedPointOperatorEvidence QualifiedEvidence(const d::FitState & state)
{
    d::FixedPointOperatorEvidence evidence;
    for (const auto & atom : state) evidence.state.push_back(atom.mdpde.GetModel());
    evidence.shape_available_atom_mask.assign(state.size(), 1);
    evidence.offset_available_atom_mask.assign(state.size(), 1);
    evidence.shape_solves.resize(state.size());
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        evidence.shape_solves[i].status = rg::RHBMEstimationStatus::SUCCESS;
        evidence.offset_solves.emplace(d::ClusterKey{i}, d::JointOffsetSolveResult{d::JointOffsetSolveStatus::Converged});
    }
    return evidence;
}
}

TEST(ProductionFittingTest, NominalQualificationBelongsToEveryEndpoint)
{
    auto fixture{ Fixture() };
    auto evidence{ QualifiedEvidence(fixture.state) };
    EXPECT_TRUE(d::AssessNominalOperator(evidence, fixture.state).certificate.StrictOperatorPassed());
    evidence.shape_solves[1].status = rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    EXPECT_FALSE(d::IsNominalOperatorSolverQualified(evidence));
    EXPECT_FALSE(d::QualifiedNominalResidualMeanSquare(evidence, fixture.state));
    evidence.shape_solves[1].status = rg::RHBMEstimationStatus::SUCCESS;
    evidence.offset_solves.at({1}).status = d::JointOffsetSolveStatus::IrlsMaximumIterationsReached;
    EXPECT_FALSE(d::IsNominalOperatorSolverQualified(evidence));
    evidence.offset_solves.at({1}).status = d::JointOffsetSolveStatus::Converged;
    evidence.shape_available_atom_mask[1] = 0;
    EXPECT_FALSE(d::AssessNominalOperator(evidence, fixture.state).certificate.StrictOperatorPassed());
}

TEST(ProductionFittingTest, NominalEvaluationIsReadOnlyAndUsesUnrestrictedSolveStatus)
{
    auto fixture{ Fixture() };
    const auto before{ fixture.state };
    const auto options{ second_stage_test::MakeSecondStageOptions() };
    const std::vector<d::ClusterKey> keys{{0}, {1}};
    const auto evidence{ d::EvaluateNominalOperator(fixture.context, keys, fixture.state, options, {1.0, 1.0}) };
    d::ClusterSolverWorkspaceMap workspaces;
    for (const auto & key : keys) workspaces.try_emplace(key);
    d::SuspiciousBlockActivity activity{{0,0}, {1,0}, {0,0}};
    const auto proposal{ d::BuildIterationProposal(fixture.context, keys, fixture.state, options,
        {1.0,1.0}, activity, workspaces) };
    for (std::size_t i = 0; i < before.size(); ++i)
    {
        second_stage_test::ExpectGaussianModelsNear(before[i].mdpde.GetModel(), fixture.state[i].mdpde.GetModel(), 0.0);
        second_stage_test::ExpectGaussianModelsNear(evidence.state[i], proposal.fixed_point_operator.state[i], 0.0);
        EXPECT_EQ(evidence.shape_solves[i].status, proposal.fixed_point_operator.shape_solves[i].status);
        EXPECT_EQ(evidence.shape_solves[i].diagnostics.iterations, proposal.fixed_point_operator.shape_solves[i].diagnostics.iterations);
    }
}

TEST(ProductionFittingTest, HistoricalPatchPreservesOtherCandidateChanges)
{
    auto fixture{ Fixture() };
    auto candidate{ d::FitStatePatch::FromState(fixture.state, {0,1}) };
    candidate.mdpde_list[1] = second_stage_test::MakeGaussianResult({8.0, 0.6, -0.4}).mdpde;
    auto best{ d::FitStatePatch::FromState(fixture.state, {0}) };
    best.mdpde_list[0] = second_stage_test::MakeGaussianResult({5.0, 0.4, 0.7}).mdpde;
    const d::FitStateView view{ fixture.state, candidate };
    const auto comparison{ d::OverlayMemberBest(view, best) };
    EXPECT_DOUBLE_EQ(comparison.Find(0)->GetModel().GetOffset(), 0.7);
    EXPECT_DOUBLE_EQ(comparison.Find(1)->GetModel().GetOffset(), -0.4);
    EXPECT_DOUBLE_EQ(fixture.state[0].mdpde.GetModel().GetOffset(), 0.12);
}

TEST(ProductionFittingTest, MemberBestRejectsPreviousOnlyProgressAndReevaluatesDomain)
{
    auto fixture{ Fixture() };
    const std::vector<d::ClusterKey> keys{{0},{1}};
    const auto domain{ d::BuildObjectiveDomain(fixture.context,
        d::BuildSecondStageModelSnapshot(fixture.context, fixture.state), keys) };
    d::MemberBestState history;
    d::UpdateMemberBestState(fixture.context, domain, fixture.state, keys, history);
    const auto original{ history.at({0}).mdpde_list[0].GetModel() };
    auto previous{ fixture.state };
    previous[0].mdpde = second_stage_test::MakeGaussianResult(original.WithOffset(original.GetOffset()+1.0)).mdpde;
    auto patch{ d::FitStatePatch::FromState(previous, {0}) };
    patch.mdpde_list[0] = second_stage_test::MakeGaussianResult(original.WithOffset(original.GetOffset()+0.9)).mdpde;
    const auto baseline{ d::BuildResidualBaseline(fixture.context, previous) };
    const auto & samples{ domain.cluster_by_key.at({0}).sample_ref_list };
    const auto score{ d::EvaluateObjectiveContribution(baseline, {0}, samples, domain) };
    ASSERT_TRUE(score);
    d::ClusterSolverWorkspaceMap solvers; d::BoundaryJointCorrectionWorkspaceMap corrections;
    d::PerformanceCounters counters(true, fixture.context, solvers, corrections);
    const d::ClusterKey key{0};
    d::LocalCandidateReference reference{d::LocalObjectivePolicy::PreviousNonRegression, key, samples,
        &*score, domain, {}, counters};
    const d::CandidateEvaluationOverlay overlay{fixture.context, baseline, previous, patch};
    EXPECT_TRUE(d::EvaluateLocalCandidate(overlay, reference).accepted);
    reference.member_best = &history.at(key);
    const auto evaluated{ d::EvaluateLocalCandidate(overlay, reference) };
    EXPECT_FALSE(evaluated.accepted);
    EXPECT_TRUE(evaluated.evidence.rejected_by_member_best);
    // Rebuilding a key creates its reference from the committed state, not a stale scalar.
    const std::vector<d::ClusterKey> merged{{0,1}};
    const auto next_domain{ d::BuildObjectiveDomain(fixture.context,
        d::BuildSecondStageModelSnapshot(fixture.context, previous), merged) };
    d::UpdateMemberBestState(fixture.context, next_domain, previous, merged, history);
    EXPECT_EQ(history.size(), 1U);
    EXPECT_EQ(history.begin()->first, merged.front());
}

TEST(ProductionFittingTest, MemberBestReevaluatesItsPatchAfterBackgroundChanges)
{
    d::SecondStageContext context;
    context.atom_list.resize(1);
    context.atom_list[0].neighbor_atom_sample_offset_list = {0,0,0,0};
    const d::FitState seed{second_stage_test::MakeGaussianResult({4.0,0.5,0.0})};
    const d::FitState historical{second_stage_test::MakeGaussianResult({6.0,0.5,0.0})};
    const d::FitState current{second_stage_test::MakeGaussianResult({5.0,0.5,0.0})};
    for (double distance : {0.15,0.35,0.60})
    {
        context.atom_list[0].raw_sampling_entries.emplace_back(LocalPotentialSample{
            2.0 * current[0].mdpde.GetModel().ResponseAtDistance(distance), SamplingPoint{distance}});
        context.atom_list[0].unselected_distance_list_by_sample.push_back({distance});
    }
    context.frozen_background = d::BuildFrozenBackground(context, seed);
    const std::vector<d::ClusterKey> keys{{0}};
    const auto domain{d::BuildObjectiveDomain(context, d::BuildSecondStageModelSnapshot(context,historical), keys)};
    d::MemberBestState history;
    d::UpdateMemberBestState(context, domain, historical, keys, history);
    d::UpdateMemberBestState(context, domain, current, keys, history);
    EXPECT_DOUBLE_EQ(history.at({0}).mdpde_list[0].GetModel().GetAmplitude(),6.0);
    context.frozen_background = d::BuildFrozenBackground(context, current);
    // Evaluating an uncommitted proposal must leave history untouched.
    const auto baseline{d::BuildResidualBaseline(context,current)};
    const d::CandidateEvaluationOverlay overlay{context,baseline,current,history.at({0})};
    EXPECT_TRUE(d::EvaluateObjectiveContribution(overlay,{0},domain.cluster_by_key.at({0}).sample_ref_list,domain));
    EXPECT_DOUBLE_EQ(history.at({0}).mdpde_list[0].GetModel().GetAmplitude(),6.0);
    // A commit in the new background can beat the same historical patch.
    d::UpdateMemberBestState(context, domain, current, keys, history);
    EXPECT_DOUBLE_EQ(history.at({0}).mdpde_list[0].GetModel().GetAmplitude(),5.0);
}

TEST(ProductionFittingTest, RecoveryRequiresResidualProgressWithinOneBestEnvelope)
{
    EXPECT_TRUE(d::IsRecoveryProgressAcceptable(100.0, 90.0, 1.0, 1.0005, 1.0));
    EXPECT_FALSE(d::IsRecoveryProgressAcceptable(100.0, 100.0, 1.0, 0.99, 1.0));
    EXPECT_FALSE(d::IsRecoveryProgressAcceptable(100.0, 90.0, 1.0, 1.0011, 1.0));
    for (int trial = 0; trial < 8; ++trial)
        EXPECT_FALSE(d::IsRecoveryProgressAcceptable(100.0, 101.0, 1.0, 1.0, std::ldexp(1.0,-trial)));
    EXPECT_FALSE(d::IsRecoveryProgressAcceptable(0.0, 0.0, 1.0, 1.0, 0.0));
}

TEST(ProductionFittingTest, RecoverySearchPreservesInputAndReportsQualification)
{
    auto fixture{ Fixture() };
    const auto before{ fixture.state };
    const std::vector<d::ClusterKey> keys{{0},{1}};
    const auto domain{ d::BuildObjectiveDomain(fixture.context,
        d::BuildSecondStageModelSnapshot(fixture.context, fixture.state), keys) };
    const auto objective{ d::EvaluateAuditObjective(domain, fixture.context,
        d::BuildSecondStageModelSnapshot(fixture.context, fixture.state)) };
    ASSERT_TRUE(objective);
    d::BestAuditState best{d::AuditedState{*objective, fixture.state}};
    d::TrustRegionStateSet radii; radii.Reconcile(keys);
    const d::SuspiciousBlockActivity activity{{0,0},{0,0},{0,0}};
    const auto recovery{ d::RunFixedPointRecovery(fixture.context, keys, fixture.state,
        second_stage_test::MakeSecondStageOptions(), {1.0,1.0}, domain, best, radii, activity) };
    EXPECT_TRUE(recovery.diagnostics.attempted);
    EXPECT_LE(recovery.diagnostics.trials.size(), 8U);
    EXPECT_GE(recovery.diagnostics.operator_evaluations, 1U);
    ASSERT_TRUE(recovery.state) << recovery.diagnostics.reason;
    EXPECT_TRUE(recovery.assessment->certificate.solver_qualified);
    for (std::size_t i = 0; i < before.size(); ++i)
        second_stage_test::ExpectGaussianModelsNear(before[i].mdpde.GetModel(), fixture.state[i].mdpde.GetModel(), 0.0);
}

TEST(ProductionFittingTest, RecoveryExhaustsEightTrialsWithoutChangingHistory)
{
    auto fixture{ Fixture() };
    const std::vector<d::ClusterKey> keys{{0},{1}};
    const auto domain{ d::BuildObjectiveDomain(fixture.context,
        d::BuildSecondStageModelSnapshot(fixture.context, fixture.state), keys) };
    auto better{ fixture.state };
    better[0].mdpde = second_stage_test::MakeGaussianResult({6.0, 0.5, 0.1}).mdpde;
    better[1].mdpde = second_stage_test::MakeGaussianResult({7.0, 0.5, -0.1}).mdpde;
    const auto objective{ d::EvaluateAuditObjective(domain, fixture.context,
        d::BuildSecondStageModelSnapshot(fixture.context, better)) };
    ASSERT_TRUE(objective);
    d::BestAuditState best{d::AuditedState{*objective, better}};
    d::TrustRegionStateSet radii; radii.Reconcile(keys); radii.ResetToMinimum(keys);
    const d::SuspiciousBlockActivity activity{{0,0},{0,0},{0,0}};
    const auto recovery{ d::RunFixedPointRecovery(fixture.context, keys, fixture.state,
        second_stage_test::MakeSecondStageOptions(), {1.0,1.0}, domain, best, radii, activity) };
    EXPECT_FALSE(recovery.state);
    ASSERT_EQ(recovery.diagnostics.trials.size(), 8U) << recovery.diagnostics.reason;
    EXPECT_EQ(recovery.diagnostics.reason, "search-exhausted");
    for (std::size_t i = 0; i < 8; ++i)
        EXPECT_DOUBLE_EQ(recovery.diagnostics.trials[i].factor, std::ldexp(1.0,-static_cast<int>(i)));
    for (std::size_t i = 0; i < better.size(); ++i)
        second_stage_test::ExpectGaussianModelsNear(better[i].mdpde.GetModel(), best->state[i].mdpde.GetModel(), 0.0);
}

TEST(ProductionFittingTest, ReplaysCapturedSolverFailureWithoutSimulationInputs)
{
    const auto * directory{ std::getenv("RHBM_TEST_REPLAY_DIR") };
    if (!directory) GTEST_SKIP() << "External solver fixtures are supplied by the validation runner.";
    std::size_t count{0};
    for (const auto & file : std::filesystem::directory_iterator(directory))
    {
        if (file.path().extension() != ".txt") continue;
        EXPECT_TRUE(second_stage_test::ReplaySolverFailure(file.path().string())) << file.path();
        ++count;
    }
    EXPECT_GT(count, 0U);
}
