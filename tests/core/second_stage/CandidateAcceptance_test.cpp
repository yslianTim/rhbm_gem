#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/CandidateState.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"
#include "core/detail/second_stage/CouplingGraph.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/JointFitting.hpp"
#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/Quarantine.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
namespace alg = rhbm_gem::algorithm;
using rhbm_gem::FittingStage;

using second_stage_test::BuildDefenseModel;
using second_stage_test::BuildJointPolishFixture;
using second_stage_test::BuildSeparatedRollbackDefenseModel;
using second_stage_test::CalculateSelectedAtomResponseMeanSquaredError;
using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::ExpectSelectedAtomEstimatesAreFinite;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeGaussianResult;
using second_stage_test::MakeSecondStageOptions;

std::unique_ptr<rg::ModelObject> BuildBoundaryComponentConflictDefenseModel(
    double intensity_scale = 1.0)
{
    std::vector<std::array<double, 3>> position_list;
    std::vector<Spot> spot_list;
    std::vector<Element> element_list;
    std::vector<rg::GaussianModel3D> truth_model_list;
    for (std::size_t i = 0; i < 103; i++)
    {
        // Break repeated edge-weight ties while retaining a connected 101-atom chain.
        const auto atom_position{ static_cast<double>(i) };
        const auto x_position{
            i < 101 ? 0.60 * atom_position + 0.002 * atom_position * atom_position :
                110.0 + 0.45 * static_cast<double>(i - 101)
        };
        position_list.push_back({ x_position, 0.0, 0.0 });
        spot_list.push_back(i % 2 == 0 ? Spot::C : Spot::O);
        element_list.push_back(
            i % 2 == 0 ? Element::CARBON : Element::OXYGEN);
        const auto truth_model{
            i >= 101 ? rg::GaussianModel3D{ 7.5, 0.70, 0.10 } :
            i % 2 == 0 ?
                rg::GaussianModel3D{ 10.0, 0.85, 0.25 } :
                rg::GaussianModel3D{ 2.0, 0.40, -0.15 }
        };
        truth_model_list.emplace_back(rg::GaussianModel3D{
            truth_model.GetAmplitude() * intensity_scale,
            truth_model.GetWidth(),
            truth_model.GetOffset() * intensity_scale
        });
    }
    return BuildDefenseModel(
        position_list,
        spot_list,
        element_list,
        truth_model_list,
        rg::GaussianModel3D{ 5.5 * intensity_scale, 0.55, 0.0 });
}


} // namespace


TEST(EstimatorSecondStageDefenseTest, AuditObjectiveProgressGuardChecksPreviousAndBest)
{
    const detail::ObjectiveTolerance tolerance{
        1.0e-8,
        1.0e-3
    };
    const detail::ObjectiveBreakdown best_one{ 1.0, 0.0, 0.0 };
    const detail::ObjectiveBreakdown best_below{ 0.99, 0.0, 0.0 };
    EXPECT_TRUE(detail::IsAuditObjectiveAcceptableForProgress(
        1.0005, 1.0, &best_one, tolerance));
    EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(
        1.002, 1.0, &best_one, tolerance));
    EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(
        1.0, 1.0, &best_below, tolerance));
}


TEST(EstimatorSecondStageDefenseTest, GlobalObjectiveWeightsClustersByAtomCount)
{
    EXPECT_DOUBLE_EQ(
        detail::CalculateClusterAtomWeight(1, 4),
        0.25);
    EXPECT_DOUBLE_EQ(
        detail::CalculateClusterAtomWeight(3, 4),
        0.75);
    EXPECT_DOUBLE_EQ(
        0.25 * 2.0 + 0.75 * 6.0,
        5.0);
}


TEST(EstimatorSecondStageDefenseTest, LocalSearchChecksTrustAndGuardBeforeStartingObjectiveBacktracking)
{
    const auto previous_model{ rg::GaussianModel3D::FromTransformedCoordinates({ 1.0, std::log(0.8), 0.0 }) };
    const auto proposal_model{ rg::GaussianModel3D::FromTransformedCoordinates({ 1.0, std::log(0.8) + 0.6, 0.0 }) };
    ASSERT_TRUE(previous_model);
    ASSERT_TRUE(proposal_model);
    for (const double target_factor : { 0.125, 0.04 })
    {
        SCOPED_TRACE(target_factor);
        const auto target_model{ rg::GaussianModel3D::FromTransformedCoordinates(
            { 1.0, std::log(0.8) + 0.6 * target_factor, 0.0 }) };
        ASSERT_TRUE(target_model);
        auto fixture{ BuildJointPolishFixture({ *previous_model }, { *target_model }) };
        const detail::ClusterKey key{ 0 };
        detail::CouplingGraphPartition partition;
        partition.sample_id_list_by_key[key] = fixture.sample_ref_list;
        const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
        const auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { key }) };
        const auto previous{ detail::BuildObjectiveByKey(partition, domain, baseline) };
        const detail::FitState proposal{ MakeGaussianResult(*proposal_model) };
        const detail::PolishProvenance provenance{ 0 };
        const detail::SuspiciousBlockActivity activity{ { 0 }, { 1 }, { 0 } };
        const std::vector<double> ridge{ 1.0 };
        const detail::ClusterHealthMap health{ { key, detail::ClusterHealth{ detail::JointOffsetSolveStatus::SystemBuildFailed } } };
        const detail::BestAuditState audit;
        detail::TrustRegionStateSet trust;
        trust.Reconcile({ key });
        detail::ClusterSolverWorkspaceMap workspaces;
        workspaces.try_emplace(key);
        detail::BoundaryJointCorrectionWorkspaceMap corrections;
        detail::PerformanceCounters counters{ true };
        const auto options{ MakeSecondStageOptions() };
        const detail::CandidateSelectionInputs inputs{
            fixture.context, options, baseline, partition, health, fixture.state, provenance,
            proposal, activity, ridge, domain, previous, audit, trust, workspaces, corrections, counters };
        detail::CandidateTransactionBuilder builder;
        builder.Select(inputs);
        const auto & selection{ builder.View() };
        ASSERT_EQ(selection.accepted_cluster_evidence_list.size(), 1U);
        const auto & evidence{ selection.accepted_cluster_evidence_list.front().evidence };
        ASSERT_TRUE(evidence.accepted_factor);
        EXPECT_EQ(evidence.invalid_trial_count, 0U);
        EXPECT_EQ(evidence.guard_rejected_trial_count, 2U);
        EXPECT_TRUE(evidence.terminal_evidence_list.empty());
        const bool objective_backtracked{ target_factor < 0.125 };
        EXPECT_DOUBLE_EQ(*evidence.accepted_factor, objective_backtracked ? 0.0625 : 0.125);
        EXPECT_EQ(evidence.objective_rejected_trial_count, objective_backtracked ? 1U : 0U);
        EXPECT_EQ(selection.shrink_trust_region_key_list,
            objective_backtracked ? std::vector<detail::ClusterKey>{ key } : std::vector<detail::ClusterKey>{});
    }
}


TEST(EstimatorSecondStageDefenseTest, CandidateCommitPublishesAcceptedRejectedAndMixedTransactions)
{
    using Commit = decltype(&detail::CandidateTransaction::Commit);
    static_assert(!std::is_copy_constructible_v<detail::CandidateTransaction>);
    static_assert(std::is_invocable_v<Commit, detail::CandidateTransaction &&,
        detail::FitState &, detail::FitState &, detail::PolishProvenance &,
        detail::QuarantineState &, detail::TrustRegionStateSet &>);
    static_assert(!std::is_invocable_v<Commit, detail::CandidateTransaction &,
        detail::FitState &, detail::FitState &, detail::PolishProvenance &,
        detail::QuarantineState &, detail::TrustRegionStateSet &>);

    for (const std::size_t accepted_count : { 0U, 1U, 2U })
    {
        SCOPED_TRACE(testing::Message() << "accepted=" << accepted_count);
        const std::vector<rg::GaussianModel3D> models{ { 6.0, 0.5, 0.0 }, { 7.0, 0.6, 0.1 } };
        auto fixture{ BuildJointPolishFixture(models, models) };
        const auto original{ fixture.state };
        const std::vector<detail::ClusterKey> keys{ { 0 }, { 1 } };
        const std::vector<detail::ClusterKey> accepted_keys(keys.begin(), keys.begin() + static_cast<std::ptrdiff_t>(accepted_count));
        const std::vector<detail::ClusterKey> rejected_keys(keys.begin() + static_cast<std::ptrdiff_t>(accepted_count), keys.end());
        detail::CouplingGraphPartition partition;
        for (const auto & key : keys) partition.sample_id_list_by_key[key] = fixture.sample_ref_list;
        const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
        const auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, keys) };
        const auto previous_objectives{ detail::BuildObjectiveByKey(partition, domain, baseline) };
        auto candidate{ fixture.state };
        for (std::size_t atom = 0; atom < accepted_count; ++atom)
        {
            candidate.at(atom).alpha_r = 0.25;
            candidate.at(atom).ols = rg::GaussianModel3DWithUncertainty{
                rg::GaussianModel3D{ 8.0, 0.7, 0.2 }, rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 } };
            candidate.at(atom).mdpde = rg::GaussianModel3DWithUncertainty{
                models.at(atom), rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 } };
        }
        const detail::PolishProvenance original_provenance{ 0, 1 };
        auto provenance{ original_provenance };
        auto candidate_provenance{ original_provenance };
        for (std::size_t atom = 0; atom < accepted_count; ++atom) candidate_provenance.at(atom) = 1 - provenance.at(atom);
        const detail::SuspiciousBlockActivity activity{ { 1, 0 }, { 0, 0 }, { 0, 0 } };
        std::vector<detail::SuspiciousGaussianAssessment> assessments(2);
        assessments.at(0).reason = detail::SuspiciousGaussianReason::WidthGrowth;
        const detail::QuarantineTarget target{ detail::QuarantineTargetKind::ShapeAtom, { 0 } };
        const detail::StabilizationTerminalFailure failure{
            detail::StabilizationTerminalReason::GuardInfeasible, detail::SuspiciousGaussianReason::WidthGrowth };
        detail::QuarantineState quarantine(2);
        quarantine.state_by_target[target] = { failure, detail::kPersistentQuarantineFailureIterationLimit - 1, 0,
            detail::QuarantineLifecycle::Active };
        detail::TrustRegionStateSet trust;
        trust.Reconcile(keys);
        trust.ResetToMinimum({ keys.at(1) });
        const std::vector<double> ridge(2, 1.0);
        const detail::ClusterHealthMap health;
        const detail::BestAuditState audit;
        detail::ClusterSolverWorkspaceMap workspaces;
        detail::BoundaryJointCorrectionWorkspaceMap corrections;
        detail::PerformanceCounters counters{ true };
        const auto options{ MakeSecondStageOptions() };
        const auto saved_level{ Logger::GetLogLevel() };
        Logger::SetLogLevel(LogLevel::Debug);
        const detail::CandidateSelectionInputs inputs{
            fixture.context, options, baseline, partition, health, fixture.state, provenance,
            candidate, activity, ridge, domain, previous_objectives, audit, trust, workspaces, corrections, counters };
        detail::CandidateSelection selection;
        selection.block_activity = activity;
        selection.assembled_state = candidate;
        selection.assembled_polish_provenance = candidate_provenance;
        selection.accepted_key_list = accepted_keys;
        selection.rejected_key_list = rejected_keys;
        selection.shrink_trust_region_key_list = accepted_keys;
        if (accepted_count == 0) selection.exhausted_key_list = { keys.at(1) };
        selection.final_audit_objective = detail::ObjectiveBreakdown{ 1.0, 2.0, 3.0 };
        selection.polish_progress = { 2, accepted_count, 2 - accepted_count, 0 };
        for (std::size_t atom = keys.size(); atom-- > 0;)
        {
            detail::ClusterCandidateDecision decision{ keys.at(atom), {} };
            if (atom == 0) decision.evidence.terminal_evidence_list.push_back({
                detail::StabilizationTerminalReason::GuardInfeasible, 0,
                detail::SuspiciousUpdateMode::PostRefit, detail::SuspiciousGaussianReason::WidthGrowth });
            (atom < accepted_count ? selection.accepted_cluster_evidence_list : selection.rejected_cluster_evidence_list).push_back(decision);
        }
        // Use the existing builder path to materialize keys and preserve rejection event order.
        detail::CandidateTransactionBuilder builder(std::move(selection));
        builder.ReconcileSelectedBoundaries(inputs);
        ASSERT_EQ(builder.View().accepted_key_list, accepted_keys);
        ASSERT_EQ(builder.View().rejected_key_list, rejected_keys);
        detail::FitState published{ MakeGaussianResult({ 99.0, 0.9, 0.9 }) };
        testing::internal::CaptureStdout();
        auto transaction{ std::move(builder).Finish(inputs, quarantine, assessments, health, {}, 0) };
        const auto finish_output{ testing::internal::GetCapturedStdout() };
        EXPECT_EQ(finish_output.find("Cluster best publication:"), std::string::npos);
        EXPECT_EQ(quarantine.TargetCount(), 0U);
        EXPECT_EQ(quarantine.state_by_target.at(target).stable_iteration_count, detail::kPersistentQuarantineFailureIterationLimit - 1);
        EXPECT_DOUBLE_EQ(trust.GetRadius(keys.at(0)), 1.0);
        EXPECT_EQ(provenance, original_provenance);
        EXPECT_EQ(published.size(), 1U);
        EXPECT_DOUBLE_EQ(published.front().mdpde.GetModel().GetAmplitude(), 99.0);
        EXPECT_EQ(fixture.state.size(), original.size());
        const auto committed{ std::move(transaction).Commit(fixture.state, published, provenance, quarantine, trust) };
        Logger::SetLogLevel(saved_level);
        EXPECT_EQ(committed.accepted_key_list, accepted_keys);
        EXPECT_EQ(committed.rejected_key_list, rejected_keys);
        EXPECT_EQ(committed.trust_region_update.changed_key_list, (std::vector<detail::ClusterKey>{ keys.at(0) }));
        EXPECT_EQ(committed.trust_region_update.saturated_key_list,
            accepted_count == 0 ? std::vector<detail::ClusterKey>{} : std::vector<detail::ClusterKey>{ keys.at(1) });
        EXPECT_DOUBLE_EQ(trust.GetRadius(keys.at(0)), 0.5);
        EXPECT_DOUBLE_EQ(trust.GetRadius(keys.at(1)), 0.0625);
        EXPECT_EQ(committed.accepted, accepted_count != 0);
        EXPECT_EQ(committed.rejected_cluster, accepted_count != keys.size());
        EXPECT_EQ(committed.suspicious_atom_count, 1U);
        EXPECT_TRUE(committed.quarantine_transition);
        EXPECT_EQ(quarantine.TargetCount(), 1U);
        EXPECT_EQ(quarantine.entered_target_count, 1U);
        EXPECT_EQ(quarantine.state_by_target.at(target).lifecycle, detail::QuarantineLifecycle::Frozen);
        EXPECT_EQ(committed.block_activity.shape_fixed_atom_mask, activity.shape_fixed_atom_mask);
        EXPECT_EQ(committed.block_activity.offset_fixed_atom_mask, activity.offset_fixed_atom_mask);
        EXPECT_EQ(committed.block_activity.hard_failure_atom_mask, activity.hard_failure_atom_mask);
        ASSERT_TRUE(committed.final_audit_objective);
        EXPECT_DOUBLE_EQ(committed.final_audit_objective->GetTotalObjective(), 4.5);
        EXPECT_EQ(committed.polish_progress.eligible_count, 2U);
        EXPECT_EQ(committed.polish_progress.accepted_count, accepted_count);
        EXPECT_EQ(committed.polish_progress.rejected_count, 2 - accepted_count);
        EXPECT_EQ(committed.polish_progress.skipped_count, 0U);
        EXPECT_EQ(provenance, accepted_count ? candidate_provenance : original_provenance);
        const auto & expected{ accepted_count ? candidate : original };
        ASSERT_EQ(published.size(), expected.size());
        for (std::size_t atom = 0; atom < expected.size(); ++atom)
        {
            EXPECT_DOUBLE_EQ(published.at(atom).alpha_r, expected.at(atom).alpha_r);
            EXPECT_EQ(published.at(atom).fit_result.has_value(), expected.at(atom).fit_result.has_value());
            for (const auto member : { &rg::LocalGaussianResult::ols, &rg::LocalGaussianResult::mdpde })
            {
                const auto & actual_model{ published.at(atom).*member };
                const auto & expected_model{ expected.at(atom).*member };
                ExpectGaussianModelsNear(actual_model.GetModel(), expected_model.GetModel(), 0.0);
                EXPECT_DOUBLE_EQ(actual_model.GetStandardDeviationModel().GetAmplitude(), expected_model.GetStandardDeviationModel().GetAmplitude());
                EXPECT_DOUBLE_EQ(actual_model.GetStandardDeviationModel().GetWidth(), expected_model.GetStandardDeviationModel().GetWidth());
                EXPECT_DOUBLE_EQ(actual_model.GetStandardDeviationModel().GetOffset(), expected_model.GetStandardDeviationModel().GetOffset());
            }
        }
    }
}


TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationBacktracksAndPreservesRemoteCluster)
{
    auto model{ BuildBoundaryComponentConflictDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 101, 103)
    };
    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 101, 103),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
