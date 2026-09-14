#include "support/SecondStageNumericalProbe.hpp"
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
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
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

std::unique_ptr<rg::ModelObject> BuildBoundaryJointCorrectionDefenseModel(
    double intensity_scale = 1.0)
{
    return BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 3.4, 0.0, 0.0 }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{
                10.0 * intensity_scale,
                0.65,
                0.2 * intensity_scale
            },
            rg::GaussianModel3D{
                2.0 * intensity_scale,
                0.65,
                -0.1 * intensity_scale
            }
        },
        rg::GaussianModel3D{
            5.5 * intensity_scale,
            0.6,
            0.0
        });
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveKeepsEarlierBestOnTie)
{
    const detail::ObjectiveTolerance tolerance{
        1.0e-10,
        1.0e-8
    };
    EXPECT_TRUE(detail::IsBetterAuditObjective(
        0.8, 1.0, tolerance));
    EXPECT_FALSE(detail::IsBetterAuditObjective(
        1.2, 1.0, tolerance));
    EXPECT_FALSE(detail::IsBetterAuditObjective(
        1.0 - 0.5e-8, 1.0, tolerance));
}

TEST(EstimatorSecondStageDefenseTest, BestAuditStateUpdateUsesPrecomputedObjective)
{
    detail::BestAuditState audit_state;
    const detail::FitState initial_state;
    const auto initial_objective{
        detail::BuildObjectiveBreakdown(2.0, 0.0, 1.0)
    };
    const auto tied_objective{
        detail::BuildObjectiveBreakdown(3.0, 0.0, 0.0)
    };
    const auto worse_objective{
        detail::BuildObjectiveBreakdown(3.5, 0.0, 0.0)
    };
    const auto improved_objective{
        detail::BuildObjectiveBreakdown(1.0, 0.0, 0.0)
    };
    ASSERT_TRUE(initial_objective.has_value());
    ASSERT_TRUE(tied_objective.has_value());
    ASSERT_TRUE(worse_objective.has_value());
    ASSERT_TRUE(improved_objective.has_value());

    EXPECT_TRUE(detail::TryUpdateBestAuditState(
        initial_state,
        false,
        0,
        *initial_objective,
        audit_state));
    ASSERT_TRUE(audit_state.has_value());
    EXPECT_DOUBLE_EQ(
        audit_state->objective.GetTotalObjective(),
        initial_objective->GetTotalObjective());
    EXPECT_FALSE(audit_state->uses_polish);
    EXPECT_EQ(audit_state->source_iteration, 0U);

    EXPECT_FALSE(detail::TryUpdateBestAuditState(
        initial_state,
        false,
        0,
        *tied_objective,
        audit_state));
    EXPECT_FALSE(detail::TryUpdateBestAuditState(
        initial_state,
        false,
        0,
        *worse_objective,
        audit_state));

    EXPECT_TRUE(detail::TryUpdateBestAuditState(
        initial_state,
        true,
        7,
        *improved_objective,
        audit_state));
    ASSERT_TRUE(audit_state.has_value());
    EXPECT_DOUBLE_EQ(
        audit_state->objective.GetTotalObjective(),
        improved_objective->GetTotalObjective());
    EXPECT_TRUE(audit_state->uses_polish);
    EXPECT_EQ(audit_state->source_iteration, 7U);
    // A background refresh must rescore both retained best and previous under
    // the same fixed domain; an old zero score must not block a new exact fit.
    detail::SecondStageContext context;
    context.atom_list.resize(1);
    context.atom_list.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0, 0 };
    const detail::FitState seed{ MakeGaussianResult({ 4.0, 0.5, 0.0 }) };
    const detail::FitState earlier_best{ MakeGaussianResult({ 6.0, 0.5, 0.0 }) };
    const detail::FitState previous{ MakeGaussianResult({ 5.0, 0.5, 0.0 }) };
    for (const double distance : { 0.15, 0.35, 0.60 })
    {
        context.atom_list.at(0).raw_sampling_entries.emplace_back(LocalPotentialSample{
            2.0 * previous.front().mdpde.GetModel().ResponseAtDistance(distance),
            SamplingPoint{ distance } });
        context.atom_list.at(0).unselected_distance_list_by_sample.push_back({ distance });
    }
    context.frozen_background = detail::BuildFrozenBackground(context, seed);
    ASSERT_TRUE(context.frozen_background);
    const auto old_snapshot{ detail::BuildSecondStageModelSnapshot(context, earlier_best) };
    const auto domain{ detail::BuildObjectiveDomain(context, old_snapshot, { { 0 } }) };
    const auto old_score{ detail::EvaluateAuditObjective(
        domain, context, old_snapshot) };
    ASSERT_TRUE(old_score.has_value());
    audit_state.reset();
    ASSERT_TRUE(detail::TryUpdateBestAuditState(earlier_best, true, 3, *old_score, audit_state));
    context.frozen_background = detail::BuildFrozenBackground(context, previous);
    ASSERT_TRUE(context.frozen_background);
    detail::ReevaluateBestAuditState(context, domain, audit_state);
    ASSERT_TRUE(audit_state.has_value());
    EXPECT_GT(audit_state->objective.GetTotalObjective(), old_score->GetTotalObjective());
    EXPECT_EQ(audit_state->source_iteration, 3U);
    EXPECT_TRUE(audit_state->uses_polish);
    ExpectGaussianModelsNear(audit_state->state.front().mdpde.GetModel(), earlier_best.front().mdpde.GetModel(), 0.0);
    const auto refreshed_snapshot{ detail::BuildSecondStageModelSnapshot(context, previous) };
    const auto refreshed_score{ detail::EvaluateAuditObjective(
        domain, context, refreshed_snapshot) };
    ASSERT_TRUE(refreshed_score.has_value());
    EXPECT_TRUE(detail::TryUpdateBestAuditState(previous, false, 4, *refreshed_score, audit_state));
    EXPECT_EQ(audit_state->source_iteration, 4U);
    EXPECT_NEAR(audit_state->objective.GetTotalObjective(), 0.0, 1.0e-12);

}

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

TEST(EstimatorSecondStageDefenseTest, AuditToleranceUsesAbsolutePlusRelativeReference)
{
    const detail::ObjectiveTolerance tolerance{
        1.0e-8,
        1.0e-3
    };
    EXPECT_TRUE(detail::IsAuditObjectiveAcceptableForProgress(
        1.0e-8,
        0.0,
        nullptr,
        tolerance));
    EXPECT_TRUE(detail::IsAuditObjectiveAcceptableForProgress(
        -2.0 + 1.0e-8 + 2.0e-3,
        -2.0,
        nullptr,
        tolerance));
    const auto boundary{ 2.0 + 1.0e-8 + 2.0e-3 };
    EXPECT_TRUE(detail::IsAuditObjectiveAcceptableForProgress(
        boundary,
        2.0,
        nullptr,
        tolerance));
    EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(
        boundary + 1.0e-9,
        2.0,
        nullptr,
        tolerance));
}

TEST(EstimatorSecondStageDefenseTest, ProgressGateValidatesToleranceBeforeTouchingEvidence)
{
    const auto infinity{ std::numeric_limits<double>::infinity() };
    const auto nan{ std::numeric_limits<double>::quiet_NaN() };
    for (const auto invalid : { -1.0, infinity, nan })
    for (const auto tolerance : { detail::ObjectiveTolerance{invalid, 0.0}, detail::ObjectiveTolerance{0.0, invalid} })
    {
        detail::ObjectiveProgressGateEvidence evidence{true, true, "untouched"};
        EXPECT_THROW(detail::IsAuditObjectiveAcceptableForProgress(nan, nan, nullptr, tolerance, &evidence),
            std::invalid_argument);
        EXPECT_TRUE(evidence.previous_checked); EXPECT_TRUE(evidence.best_checked);
        EXPECT_EQ(evidence.reason, "untouched");
    }
}

TEST(EstimatorSecondStageDefenseTest, ProgressGateEvidencePreservesShortCircuitAndResetsBetweenCalls)
{
    const auto tolerance{ detail::kObjectiveProgressTolerance };
    const detail::ObjectiveBreakdown best{0.5, 0.0, 0.0};
    detail::ObjectiveProgressGateEvidence evidence;
    EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(2.0, 1.0, &best, tolerance, &evidence));
    EXPECT_TRUE(evidence.previous_checked); EXPECT_FALSE(evidence.best_checked);
    EXPECT_EQ(evidence.reason, "previous-gate");
    EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(1.0, 1.0, &best, tolerance, &evidence));
    EXPECT_TRUE(evidence.previous_checked); EXPECT_TRUE(evidence.best_checked);
    EXPECT_EQ(evidence.reason, "best-gate");
    EXPECT_TRUE(detail::IsAuditObjectiveAcceptableForProgress(1.0, 1.0, nullptr, tolerance, &evidence));
    EXPECT_TRUE(evidence.previous_checked); EXPECT_FALSE(evidence.best_checked);
    EXPECT_TRUE(evidence.reason.empty());
    for (const auto nonfinite : { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN() })
    {
        EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(nonfinite, 1.0, &best, tolerance, &evidence));
        EXPECT_FALSE(evidence.previous_checked); EXPECT_FALSE(evidence.best_checked);
        EXPECT_EQ(evidence.reason, "objective-nonfinite");
        EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(1.0, nonfinite, &best, tolerance, &evidence));
        EXPECT_FALSE(evidence.previous_checked); EXPECT_FALSE(evidence.best_checked);
        const detail::ObjectiveBreakdown invalid_best{nonfinite, 0.0, 0.0};
        EXPECT_FALSE(detail::IsAuditObjectiveAcceptableForProgress(1.0, 1.0, &invalid_best, tolerance, &evidence));
        EXPECT_TRUE(evidence.previous_checked); EXPECT_TRUE(evidence.best_checked);
        EXPECT_EQ(evidence.reason, "best-gate");
    }
}

TEST(EstimatorSecondStageDefenseTest, GlobalCandidateRequiresPreviousAndAvailableDelta)
{
    auto fixture{ BuildJointPolishFixture({ { 6.0, 0.5, 0.0 } }, { { 6.4, 0.5, 0.0 } }) };
    const detail::ClusterKey key{ 0 };
    const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
    auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { key }) };
    const auto previous{ detail::EvaluateAuditObjective(domain, baseline) };
    ASSERT_TRUE(previous);
    const auto patch{ detail::FitStatePatch::FromState(fixture.state, key) };
    const detail::CandidateEvaluationOverlay candidate{ fixture.context, baseline, fixture.state, patch };
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    const auto saved_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    second_stage_test::BeginNumericalCapture();
    testing::internal::CaptureStdout();
    {
        detail::PerformanceCounters counters{ false, fixture.context, workspaces, corrections };
        EXPECT_FALSE(detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
            fixture.sample_ref_list, domain, nullptr, nullptr, counters }));
    }
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(saved_level);
    EXPECT_EQ(second_stage_test::EndNumericalCapture().work[2], 0U);

    domain.cluster_by_key.at(key).scale.reset();
    detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
    EXPECT_FALSE(detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
        fixture.sample_ref_list, domain, nullptr, &*previous, counters }));
}

TEST(EstimatorSecondStageDefenseTest, GlobalCandidatePreservesGatesAndOverlappingFitTailDelta)
{
    for (const bool overlap : { false, true })
    {
        SCOPED_TRACE(overlap);
        auto fixture{ BuildJointPolishFixture({ { 6.0, 0.5, 0.0 } }, { { 6.4, 0.5, 0.0 } }) };
        const detail::ClusterKey key{ 0 };
        const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
        auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { key }) };
        if (overlap)
        {
            auto & cluster{ domain.cluster_by_key.at(key) };
            cluster.tail_sample_ref_list = { { 0, 0 } };
            cluster.scale->tail = 2.0 * cluster.scale->fit;
            domain.tail_sample_mask_by_atom[0][0] = 1;
            domain.tail_sample_count = 1;
        }
        const auto previous{ detail::EvaluateAuditObjective(domain, baseline) };
        ASSERT_TRUE(previous);
        detail::ClusterSolverWorkspaceMap workspaces;
        detail::BoundaryJointCorrectionWorkspaceMap corrections;
        detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
        const detail::FitState improved{ MakeGaussianResult({ 6.2, 0.5, 0.0 }) };
        const auto patch{ detail::FitStatePatch::FromState(improved, key) };
        const detail::CandidateEvaluationOverlay candidate{ fixture.context, baseline, fixture.state, patch };
        const auto full{ detail::EvaluateAuditObjective(domain, detail::BuildResidualBaseline(fixture.context, improved)) };
        ASSERT_TRUE(full);
        ASSERT_LT(full->GetTotalObjective(), previous->GetTotalObjective());
        const auto accepted{ detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
            fixture.sample_ref_list, domain, nullptr, &*previous, counters }) };
        ASSERT_TRUE(accepted);
        EXPECT_NEAR(accepted->fit_range_residual_objective, full->fit_range_residual_objective, 1.0e-12);
        EXPECT_NEAR(accepted->tail_validation_loss, full->tail_validation_loss, 1.0e-12);
        EXPECT_NEAR(accepted->offset_plausibility_penalty, full->offset_plausibility_penalty, 1.0e-12);
        EXPECT_TRUE(detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
            fixture.sample_ref_list, domain, &*previous, &*previous, counters }));
        const detail::ObjectiveBreakdown best{ 0.0, 0.0, 0.0 };
        EXPECT_FALSE(detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
            fixture.sample_ref_list, domain, &best, &*previous, counters }));

        const detail::FitState worse{ MakeGaussianResult({ 20.0, 0.5, 0.0 }) };
        const auto worse_patch{ detail::FitStatePatch::FromState(worse, key) };
        const detail::CandidateEvaluationOverlay worse_candidate{ fixture.context, baseline, fixture.state, worse_patch };
        EXPECT_FALSE(detail::EvaluateGlobalCandidate(worse_candidate, detail::GlobalCandidateReference{
            fixture.sample_ref_list, domain, nullptr, &*previous, counters }));
    }
}

TEST(EstimatorSecondStageDefenseTest, GlobalCandidateKeepsInclusiveProgressTolerance)
{
    auto fixture{ BuildJointPolishFixture({ { 6.0, 0.5, 0.0 } }, { { 6.4, 0.5, 0.0 } }) };
    const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
    const auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { { 0 } }) };
    const detail::FitStatePatch patch;
    const detail::CandidateEvaluationOverlay candidate{ fixture.context, baseline, fixture.state, patch };
    const std::vector<detail::SampleRef> samples;
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
    const detail::ObjectiveBreakdown best{ 2.0, 0.0, 0.0 };
    // An empty delta preserves the supplied baseline exactly, isolating the global gate.
    const detail::ObjectiveBreakdown boundary{ 2.0 + 1.0e-8 + 2.0e-3, 0.0, 0.0 };
    const auto accepted{ detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
        samples, domain, &best, &boundary, counters }) };
    ASSERT_TRUE(accepted);
    EXPECT_DOUBLE_EQ(accepted->GetTotalObjective(), boundary.GetTotalObjective());
    const detail::ObjectiveBreakdown beyond{ boundary.GetTotalObjective() + 1.0e-9, 0.0, 0.0 };
    EXPECT_FALSE(detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
        samples, domain, &best, &beyond, counters }));
    const detail::ObjectiveBreakdown nonfinite{ std::numeric_limits<double>::infinity(), 0.0, 0.0 };
    EXPECT_FALSE(detail::EvaluateGlobalCandidate(candidate, detail::GlobalCandidateReference{
        samples, domain, nullptr, &nonfinite, counters }));
}

TEST(EstimatorSecondStageDefenseTest, BoundaryCandidateKeepsPolicyReferencesAndMissingEvidence)
{
    auto fixture{ BuildJointPolishFixture({ { 6.0, 0.5, 0.0 } }, { { 6.4, 0.5, 0.0 } }) };
    const detail::ClusterKey key{ 0 };
    detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key[key] = fixture.sample_ref_list;
    const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
    auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { key }) };
    auto previous_by_key{ detail::BuildObjectiveByKey(partition, domain, baseline) };
    const auto previous{ detail::EvaluateAuditObjective(domain, baseline) };
    ASSERT_TRUE(previous);
    const detail::BoundaryReconciliationComponent component{ .key_list = { key },
        .affected_sample_ref_list = fixture.sample_ref_list, .halo_atom_index_list = key };
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
    const detail::FitState improved{ MakeGaussianResult({ 6.2, 0.5, 0.0 }) };
    const auto patch{ detail::FitStatePatch::FromState(improved, key) };
    const detail::CandidateEvaluationOverlay candidate{ fixture.context, baseline, fixture.state, patch };
    const detail::FitStatePatch unchanged_patch;
    const detail::CandidateEvaluationOverlay unchanged{ fixture.context, baseline, fixture.state, unchanged_patch };
    const detail::ObjectiveBreakdown best{ 0.0, 0.0, 0.0 };
    const detail::ObjectiveBreakdown nonfinite{ std::numeric_limits<double>::infinity(), 0.0, 0.0 };
    for (const auto policy : { detail::BoundaryAcceptancePolicy::Ordinary, detail::BoundaryAcceptancePolicy::CooperativeRescue })
    {
        SCOPED_TRACE(static_cast<int>(policy));
        detail::BoundaryCandidateReference reference{
            .policy = policy, .samples_by_key = partition.sample_id_list_by_key, .domain = domain,
            .previous_objective_by_key = previous_by_key, .best_audit = &best, .counters = counters,
            .component = component };
        EXPECT_EQ(detail::EvaluateBoundaryCandidate(candidate, reference, &*previous).has_value(), policy == detail::BoundaryAcceptancePolicy::Ordinary);
        reference.best_audit = &*previous;
        const auto accepted{ detail::EvaluateBoundaryCandidate(candidate, reference, &*previous) };
        ASSERT_TRUE(accepted);
        const auto full{ detail::EvaluateAuditObjective(domain, detail::BuildResidualBaseline(fixture.context, improved)) };
        ASSERT_TRUE(full);
        EXPECT_NEAR(accepted->GetTotalObjective(), full->GetTotalObjective(), 1.0e-12);
        EXPECT_EQ(detail::EvaluateBoundaryCandidate(unchanged, reference, &*previous).has_value(), policy == detail::BoundaryAcceptancePolicy::Ordinary);
        EXPECT_FALSE(detail::EvaluateBoundaryCandidate(candidate, reference, nullptr));
        EXPECT_FALSE(detail::EvaluateBoundaryCandidate(candidate, reference, &nonfinite));
        const auto saved_member{ previous_by_key.at(key) };
        previous_by_key.at(key).reset();
        EXPECT_FALSE(detail::EvaluateBoundaryCandidate(candidate, reference, &*previous));
        previous_by_key.at(key) = nonfinite;
        EXPECT_FALSE(detail::EvaluateBoundaryCandidate(candidate, reference, &*previous));
        previous_by_key.at(key) = detail::ObjectiveBreakdown{ -1.0, 0.0, 0.0 };
        EXPECT_FALSE(detail::EvaluateBoundaryCandidate(candidate, reference, &*previous));
        previous_by_key.at(key) = saved_member;
        const auto saved_scale{ domain.cluster_by_key.at(key).scale };
        domain.cluster_by_key.at(key).scale.reset();
        EXPECT_FALSE(detail::EvaluateBoundaryCandidate(candidate, reference, &*previous));
        domain.cluster_by_key.at(key).scale = saved_scale;
    }
}

TEST(EstimatorSecondStageDefenseTest, BoundaryCorrectionKeepsStrictReferenceAndSingleDeltaEvaluation)
{
    auto fixture{ BuildJointPolishFixture({ { 6.0, 0.5, 0.0 } }, { { 6.4, 0.5, 0.0 } }) };
    const detail::ClusterKey key{ 0 };
    detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key[key] = fixture.sample_ref_list;
    const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
    const auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { key }) };
    const auto previous_by_key{ detail::BuildObjectiveByKey(partition, domain, baseline) };
    const auto previous{ detail::EvaluateAuditObjective(domain, baseline) };
    ASSERT_TRUE(previous);
    const detail::BoundaryReconciliationComponent component{ .key_list = { key },
        .affected_sample_ref_list = fixture.sample_ref_list, .halo_atom_index_list = key };
    const detail::FitState improved{ MakeGaussianResult({ 6.2, 0.5, 0.0 }) };
    const auto patch{ detail::FitStatePatch::FromState(improved, key) };
    const detail::CandidateEvaluationOverlay candidate{ fixture.context, baseline, fixture.state, patch };
    const detail::FitStatePatch empty_patch;
    const detail::FitStateView endpoint{ fixture.state, empty_patch };
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    const auto saved_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    for (const bool unavailable : { false, true })
    {
        SCOPED_TRACE(unavailable);
        second_stage_test::BeginNumericalCapture();
        testing::internal::CaptureStdout();
        {
            detail::PerformanceCounters counters{ false, fixture.context, workspaces, corrections };
            const detail::ObjectiveBreakdown audit{ unavailable ? std::numeric_limits<double>::infinity() :
                previous->fit_range_residual_objective, previous->tail_validation_loss, previous->offset_plausibility_penalty };
            const detail::BoundaryCandidateReference reference{
                .policy = detail::BoundaryAcceptancePolicy::Ordinary,
                .samples_by_key = partition.sample_id_list_by_key, .domain = domain,
                .previous_objective_by_key = previous_by_key, .best_audit = nullptr, .counters = counters,
                .component = component };
            const auto result{ detail::EvaluateBoundaryCorrection(candidate, reference, endpoint, audit, *previous) };
            EXPECT_EQ(result, !unavailable);
        }
        const auto output{ testing::internal::GetCapturedStdout() };
        // One delta, its two residual contributions, and the member contribution.
        EXPECT_EQ(second_stage_test::EndNumericalCapture().work[2], 4U);
    }
    Logger::SetLogLevel(saved_level);
    detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
    const auto objective{ detail::EvaluateObjectiveDelta(candidate, fixture.sample_ref_list, domain, *previous, counters) };
    ASSERT_TRUE(objective);
    const detail::BoundaryCandidateReference tied_reference{
        .policy = detail::BoundaryAcceptancePolicy::Ordinary,
        .samples_by_key = partition.sample_id_list_by_key, .domain = domain,
        .previous_objective_by_key = previous_by_key, .best_audit = nullptr, .counters = counters,
        .component = component };
    const auto tied{ detail::EvaluateBoundaryCorrection(candidate, tied_reference, endpoint, *previous, *objective) };
    EXPECT_FALSE(tied);
    for (const auto [margin, expected] : { std::pair{1.0e-9, false}, std::pair{1.0e-6, true} })
    {
        const detail::ObjectiveBreakdown improvement{ objective->GetTotalObjective() + margin, 0.0, 0.0 };
        const detail::BoundaryCandidateReference reference{
            .policy = detail::BoundaryAcceptancePolicy::Ordinary,
            .samples_by_key = partition.sample_id_list_by_key, .domain = domain,
            .previous_objective_by_key = previous_by_key, .best_audit = nullptr, .counters = counters,
            .component = component };
        const auto result{ detail::EvaluateBoundaryCorrection(candidate, reference, endpoint, *previous, improvement) };
        EXPECT_EQ(result, expected);
    }
}

TEST(EstimatorSecondStageDefenseTest, ScientificObjectiveUsesFitTailAndOffsetOnly)
{
    const auto objective{
        detail::BuildObjectiveBreakdown(
            0.4,
            2.0,
            0.18)
    };
    const auto previous{
        detail::BuildObjectiveBreakdown(
            1.4,
            0.0,
            0.0)
    };

    ASSERT_TRUE(objective.has_value());
    ASSERT_TRUE(previous.has_value());
    const auto empty_tail{
        detail::BuildObjectiveBreakdown(
            0.4,
            0.0,
            0.0)
    };
    ASSERT_TRUE(empty_tail.has_value());
    EXPECT_DOUBLE_EQ(empty_tail->tail_validation_loss, 0.0);
    EXPECT_DOUBLE_EQ(empty_tail->GetTailValidationPenalty(), 0.0);
    EXPECT_DOUBLE_EQ(objective->tail_validation_loss, 2.0);
    EXPECT_DOUBLE_EQ(objective->GetTailValidationPenalty(), 0.5);
    EXPECT_DOUBLE_EQ(objective->offset_plausibility_penalty, 0.18);
    EXPECT_DOUBLE_EQ(objective->GetTotalObjective(), 1.08);
    EXPECT_DOUBLE_EQ(
        objective->GetTotalObjective(),
        objective->fit_range_residual_objective +
            objective->GetTailValidationPenalty() +
            objective->offset_plausibility_penalty);
    EXPECT_TRUE(detail::IsBetterAuditObjective(
        objective->GetTotalObjective(),
        previous->GetTotalObjective(),
        detail::ObjectiveTolerance{ 1.0e-10, 1.0e-8 }));
    EXPECT_FALSE(
        detail::BuildObjectiveBreakdown(
            std::numeric_limits<double>::infinity(),
            2.0,
            0.18).has_value());
    EXPECT_FALSE(
        detail::BuildObjectiveBreakdown(
            std::numeric_limits<double>::max(),
            0.0,
            std::numeric_limits<double>::max()).has_value());
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

TEST(EstimatorSecondStageDefenseTest, AcceptedTrustRegionRadiusShrinksOnlyAfterObjectiveBacktracking)
{
    EXPECT_FALSE(detail::ShouldShrinkAcceptedTrustRegionRadius(0.5, 0.5));
    EXPECT_TRUE(detail::ShouldShrinkAcceptedTrustRegionRadius(0.5, 0.25));
    EXPECT_FALSE(detail::ShouldShrinkAcceptedTrustRegionRadius(0.5, 1.0));
    EXPECT_FALSE(detail::ShouldShrinkAcceptedTrustRegionRadius(std::nullopt, 0.5));
    EXPECT_FALSE(detail::ShouldShrinkAcceptedTrustRegionRadius(0.5, std::nullopt));
    EXPECT_FALSE(detail::ShouldShrinkAcceptedTrustRegionRadius(std::nullopt, std::nullopt));
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionStateReconcilesKeepsShrinksAndSaturates)
{
    detail::TrustRegionStateSet state;
    const detail::ClusterKey key{ 0 };
    state.Reconcile({ key });
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 1.0);

    for (const auto expected : { 0.5, 0.25, 0.125, 0.0625 })
    {
        const auto update{
            state.ApplyRadiusUpdates({ key }, {}, {})
        };
        EXPECT_EQ(
            update.changed_key_list,
            std::vector<detail::ClusterKey>{ key });
        EXPECT_TRUE(update.saturated_key_list.empty());
        EXPECT_DOUBLE_EQ(state.GetRadius(key), expected);
    }
    const auto saturated{
        state.ApplyRadiusUpdates({ key }, {}, {})
    };
    EXPECT_TRUE(saturated.changed_key_list.empty());
    EXPECT_EQ(
        saturated.saturated_key_list,
        std::vector<detail::ClusterKey>{ key });

    state.ApplyRadiusUpdates({}, {}, {});
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 0.0625);

    const detail::ClusterKey replacement_key{ 1 };
    state.Reconcile({ replacement_key });
    EXPECT_DOUBLE_EQ(state.GetRadius(replacement_key), 1.0);
}

TEST(EstimatorSecondStageDefenseTest, ExhaustedRejectionsAreExcludedFromRadiusShrink)
{
    using Key = detail::ClusterKey;
    const Key keep_key{ 0 };
    const Key shrink_key{ 1 };

    detail::TrustRegionStateSet state;
    state.Reconcile({ keep_key, shrink_key });
    const auto update{ state.ApplyRadiusUpdates(
        {},
        { shrink_key },
        { shrink_key }) };

    EXPECT_DOUBLE_EQ(state.GetRadius(keep_key), 1.0);
    EXPECT_DOUBLE_EQ(state.GetRadius(shrink_key), 1.0);
    EXPECT_TRUE(update.changed_key_list.empty());
    EXPECT_TRUE(update.saturated_key_list.empty());
}

TEST(EstimatorSecondStageDefenseTest, TerminalRejectionShrinksRadiusOncePerControllerUpdate)
{
    using Key = detail::ClusterKey;
    const Key first_key{ 0 };
    const Key second_key{ 1 };
    detail::TrustRegionStateSet state;
    state.Reconcile({ first_key, second_key });

    const auto update{ state.ApplyRadiusUpdates(
        {},
        { first_key, second_key },
        {}) };

    EXPECT_EQ(
        update.changed_key_list,
        (std::vector<Key>{ first_key, second_key }));
    EXPECT_DOUBLE_EQ(state.GetRadius(first_key), 0.5);
    EXPECT_DOUBLE_EQ(state.GetRadius(second_key), 0.5);
}

TEST(EstimatorSecondStageDefenseTest,
    BacktrackingWorkspaceGeneratesCandidatesAndMergesProvenance)
{
    const std::vector<rg::GaussianModel3D> previous_model_list{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 },
        rg::GaussianModel3D{ 10.0, 0.60, 0.20 }
    };
    const std::vector<rg::GaussianModel3D> endpoint_model_list{
        rg::GaussianModel3D{ 12.0, 0.75, 0.40 },
        rg::GaussianModel3D{ 14.0, 0.90, 0.80 }
    };
    const std::vector<rg::GaussianModel3DUncertainty> endpoint_uncertainty_list{
        rg::GaussianModel3DUncertainty{ 0.10, 0.02, 0.03 },
        rg::GaussianModel3DUncertainty{ 0.20, 0.04, 0.05 }
    };

    detail::FitState previous_state;
    detail::FitState endpoint_state;
    previous_state.resize(previous_model_list.size());
    endpoint_state.resize(endpoint_model_list.size());
    for (std::size_t atom_index = 0;
        atom_index < previous_model_list.size();
        atom_index++)
    {
        previous_state.at(atom_index).mdpde =
            rg::GaussianModel3DWithUncertainty{
                previous_model_list.at(atom_index),
                rg::GaussianModel3DUncertainty{ 0.01, 0.02, 0.03 }
            };
        endpoint_state.at(atom_index).mdpde =
            rg::GaussianModel3DWithUncertainty{
                endpoint_model_list.at(atom_index),
                endpoint_uncertainty_list.at(atom_index)
            };
    }

    const auto endpoint_patch{
        detail::FitStatePatch::FromState(
            endpoint_state,
            std::vector<std::size_t>{ 1, 0 })
    };
    detail::BacktrackingWorkspace workspace{
        previous_state,
        endpoint_patch,
        1.0e-4
    };
    const auto step{ workspace.BuildNextCandidate() };
    ASSERT_EQ(
        step.status,
        detail::BacktrackingStepStatus::CandidateReady);
    EXPECT_DOUBLE_EQ(step.factor, 0.5);
    EXPECT_EQ(step.trial_number, 2U);
    const auto & candidate_patch{ workspace.GetCandidatePatch() };
    EXPECT_EQ(
        candidate_patch.atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_DOUBLE_EQ(
        candidate_patch.mdpde_list.at(0)
            .GetStandardDeviationModel().GetAmplitude(),
        endpoint_uncertainty_list.at(0).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        candidate_patch.mdpde_list.at(1)
            .GetStandardDeviationModel().GetWidth(),
        endpoint_uncertainty_list.at(1).GetWidth());

    auto candidate_state{ previous_state };
    candidate_patch.ApplyTo(candidate_state);
    for (const auto atom_index : candidate_patch.atom_index_list)
    {
        ExpectGaussianModelsNear(
            candidate_state.at(atom_index).mdpde.GetModel(),
            candidate_patch.mdpde_list.at(atom_index).GetModel(),
            1.0e-12);
    }
    EXPECT_DOUBLE_EQ(candidate_state.at(0).mdpde.GetModel().GetOffset(),
        std::lerp(-0.10, 0.40, 0.5));
    EXPECT_DOUBLE_EQ(candidate_state.at(1).mdpde.GetModel().GetOffset(),
        std::lerp(0.20, 0.80, 0.5));

    const auto merged_provenance{
        workspace.BuildCandidatePolishProvenance(
            std::vector<char>{ 0, 1 },
            std::vector<char>{ 1, 0 })
    };
    EXPECT_EQ(merged_provenance, (std::vector<char>{ 1, 0 }));

    detail::SecondStageContext median_context;
    median_context.atom_list.resize(3);
    detail::FitState median_previous;
    detail::FitState median_endpoint;
    const std::array previous_offsets{ 0.0, 10.0, 20.0 };
    const std::array endpoint_offsets{ 30.0, 5.0, 0.0 };
    for (std::size_t node = 0; node < 3; node++)
    {
        median_context.atom_list.at(node).raw_sampling_entries.resize(1);
        median_context.atom_list.at(node).unselected_distance_list_by_sample = { { 0.3 } };
        median_previous.emplace_back(MakeGaussianResult({ 8.0 + static_cast<double>(node), 0.5, previous_offsets.at(node) }));
        median_endpoint.emplace_back(MakeGaussianResult({ 17.0 - static_cast<double>(node), 0.7, endpoint_offsets.at(node) }));
    }
    median_context.frozen_background = detail::BuildFrozenBackground(median_context, median_previous);
    ASSERT_TRUE(median_context.frozen_background);
    const auto frozen{ median_context.frozen_background };
    detail::BacktrackingWorkspace median_workspace{ median_previous,
        detail::FitStatePatch::FromState(median_endpoint, { 0, 1, 2 }), 1.0e-4 };
    ASSERT_EQ(median_workspace.BuildNextCandidate().status,
        detail::BacktrackingStepStatus::CandidateReady);
    EXPECT_EQ(median_workspace.GetCandidatePatch().atom_index_list.size(), 3U);
    EXPECT_EQ(median_context.frozen_background, frozen);
    ExpectGaussianModelsNear(frozen->model_by_atom.front(), { 9.0, 0.5, 10.0 }, 1.0e-12);
    const auto next{ detail::BuildFrozenBackground(median_context, median_endpoint) };
    ASSERT_TRUE(next);
    ExpectGaussianModelsNear(next->model_by_atom.front(), { 16.0, 0.7, 5.0 }, 1.0e-12);
    EXPECT_NE(next->response_by_atom, frozen->response_by_atom);

}

TEST(EstimatorSecondStageDefenseTest,
    BacktrackingWorkspaceStopsWhenChangeBecomesNonmaterial)
{
    detail::SecondStageContext context;
    context.atom_list.resize(1);
    detail::FitState previous_state(1);
    previous_state.at(0).mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
    };

    detail::FitState endpoint_state(1);
    endpoint_state.at(0).mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 12.0, 0.75, 0.40 },
        rg::GaussianModel3DUncertainty{ 0.2, 0.04, 0.05 }
    };
    const auto endpoint_patch{
        detail::FitStatePatch::FromState(
            endpoint_state,
            std::vector<std::size_t>{ 0 })
    };
    detail::BacktrackingWorkspace change_exhausted_workspace{
        previous_state,
        endpoint_patch,
        1.0e6
    };
    const auto change_exhausted_step{
        change_exhausted_workspace.BuildNextCandidate()
    };
    EXPECT_EQ(
        change_exhausted_step.status,
        detail::BacktrackingStepStatus::Exhausted);
    EXPECT_EQ(change_exhausted_step.trial_number, 1U);
}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveSourcesAgreeAcrossTailPartitions)
{
    const rg::GaussianModel3D model{ 8.0, 0.50, -0.10 };
    const detail::ClusterKey key{ 0 };
    for (const bool include_tail : { false, true })
    {
        detail::SecondStageContext context;
        context.atom_list.resize(1);
        auto & atom{ context.atom_list.front() };
        const std::vector<double> distances{ include_tail ?
            std::vector<double>{ 0.0, 1.0, 1.1, 1.2, 2.0, 2.1 } :
            std::vector<double>{ 0.0, 1.0, 1.1, 2.1 }
        };
        atom.neighbor_atom_sample_offset_list.assign(distances.size() + 1, 0);
        std::vector<detail::SampleRef> all_samples;
        for (const auto distance : distances)
        {
            all_samples.push_back({ 0, atom.raw_sampling_entries.size() });
            atom.raw_sampling_entries.push_back({ model.ResponseAtDistance(distance) + 0.1,
                SamplingPoint{ distance } });
        }
        const detail::FitState state{ MakeGaussianResult(model) };
        auto baseline{ detail::BuildResidualBaseline(context, state) };
        auto domain{ detail::BuildObjectiveDomain(context, baseline.model_snapshot, { key }) };
        EXPECT_EQ(domain.fit_sample_count, 2U);
        EXPECT_EQ(domain.tail_sample_count, include_tail ? 2U : 0U);
        EXPECT_EQ(domain.unique_sample_count, include_tail ? 4U : 2U);
        EXPECT_EQ(domain.fit_sample_mask_by_atom.front(),
            (include_tail ? std::vector<char>{ 1, 1, 0, 0, 0, 0 } : std::vector<char>{ 1, 1, 0, 0 }));
        EXPECT_EQ(domain.tail_sample_mask_by_atom.front(),
            (include_tail ? std::vector<char>{ 0, 0, 0, 1, 1, 0 } : std::vector<char>{ 0, 0, 0, 0 }));

        // Excluded samples must be skipped even when their cached residual is unavailable.
        baseline.sample_list.front().at(2).reset();
        baseline.sample_list.front().back().reset();
        detail::ClusterSolverWorkspaceMap workspaces;
        detail::BoundaryJointCorrectionWorkspaceMap corrections;
        detail::PerformanceCounters counters{ true, context, workspaces, corrections };
        auto candidate_state{ state };
        candidate_state.front() = MakeGaussianResult({ 8.1, 0.51, -0.10 });
        const auto patch{ detail::FitStatePatch::FromState(candidate_state, key) };
        const detail::CandidateEvaluationOverlay overlay{ context, baseline, state, patch };
        const auto candidate_snapshot{ detail::BuildSecondStageModelSnapshot(context, candidate_state) };
        for (const bool overlap : { false, true })
        {
            auto & cluster{ domain.cluster_by_key.at(key) };
            if (overlap)
            {
                // Keep the physical sample union unchanged; give r=1 both roles.
                cluster.tail_sample_ref_list.push_back({ 0, 1 });
                domain.tail_sample_mask_by_atom.front().at(1) = 1;
                domain.tail_sample_count++;
                cluster.scale->tail = 2.0 * cluster.scale->fit;
            }
            const auto snapshot_objective{ detail::EvaluateAuditObjective(domain, context, baseline.model_snapshot) };
            const auto baseline_objective{ detail::EvaluateAuditObjective(domain, baseline) };
            const auto contribution{ detail::EvaluateObjectiveContribution(baseline, key, all_samples, domain) };
            ASSERT_TRUE(snapshot_objective.has_value());
            ASSERT_TRUE(baseline_objective.has_value());
            ASSERT_TRUE(contribution.has_value());
            EXPECT_DOUBLE_EQ(snapshot_objective->GetTotalObjective(), baseline_objective->GetTotalObjective());
            EXPECT_DOUBLE_EQ(contribution->GetTotalObjective(), baseline_objective->GetTotalObjective());
            if (!include_tail && !overlap) EXPECT_DOUBLE_EQ(baseline_objective->tail_validation_loss, 0.0);
            double expected_tail{ 0.0 };
            for (const auto & ref : cluster.tail_sample_ref_list)
            {
                expected_tail += alg::CalculateCauchyLoss(
                    baseline(ref)->residual / cluster.scale->tail,
                    detail::kObjectiveRobustLossCutoffMultiplier) /
                    static_cast<double>(cluster.tail_sample_ref_list.size());
            }
            EXPECT_NEAR(baseline_objective->tail_validation_loss, expected_tail, 1.0e-12);
            const auto delta{ detail::EvaluateObjectiveDelta(
                overlay, all_samples, domain, *baseline_objective, counters) };
            const auto full{ detail::EvaluateAuditObjective(domain, context, candidate_snapshot) };
            ASSERT_TRUE(delta.has_value());
            ASSERT_TRUE(full.has_value());
            EXPECT_NEAR(delta->fit_range_residual_objective, full->fit_range_residual_objective, 1.0e-12);
            EXPECT_NEAR(delta->tail_validation_loss, full->tail_validation_loss, 1.0e-12);
            EXPECT_EQ(domain.unique_sample_count, include_tail ? 4U : 2U);
        }
    }

    detail::SecondStageContext tail_only_context;
    tail_only_context.atom_list.resize(1);
    auto & tail_only_atom{ tail_only_context.atom_list.front() };
    tail_only_atom.neighbor_atom_sample_offset_list = { 0, 0 };
    tail_only_atom.raw_sampling_entries.push_back({ model.ResponseAtDistance(1.2), SamplingPoint{ 1.2 } });
    const detail::FitState tail_only_state{ MakeGaussianResult(model) };
    const auto tail_only_baseline{ detail::BuildResidualBaseline(tail_only_context, tail_only_state) };
    const auto tail_only_domain{
        detail::BuildObjectiveDomain(tail_only_context, tail_only_baseline.model_snapshot, { key })
    };
    EXPECT_FALSE(tail_only_domain.cluster_by_key.at(key).scale.has_value());
    EXPECT_FALSE(detail::EvaluateAuditObjective(tail_only_domain, tail_only_baseline).has_value());

}

TEST(EstimatorSecondStageDefenseTest, LocalCandidateUsesOnlyPreviousAndRejectsUnavailableObjective)
{
    auto fixture{ BuildJointPolishFixture({ { 6.0, 0.5, 0.0 } }, { { 6.4, 0.5, 0.0 } }) };
    const detail::ClusterKey key{ 0 };
    const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
    auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, { key }) };
    const auto previous{ detail::EvaluateObjectiveContribution(baseline, key, fixture.sample_ref_list, domain) };
    ASSERT_TRUE(previous);
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
    for (double amplitude : { 6.4, 20.0 })
    {
        const detail::FitState candidate{ MakeGaussianResult({ amplitude, 0.5, 0.0 }) };
        const auto patch{ detail::FitStatePatch::FromState(candidate, key) };
        const detail::CandidateEvaluationOverlay overlay{ fixture.context, baseline, fixture.state, patch };
        const auto current{ detail::EvaluateObjectiveContribution(overlay, key, fixture.sample_ref_list, domain) };
        ASSERT_TRUE(current);
        detail::CandidateDecisionEvidence evidence;
        const auto result{ detail::EvaluateLocalCandidate(overlay,
            detail::LocalCandidateReference{detail::LocalObjectivePolicy::PreviousNonRegression, key, fixture.sample_ref_list,
                amplitude == 6.4 ? &*previous : &*current, domain, evidence, counters}) };
        EXPECT_TRUE(result.accepted);
        EXPECT_FALSE(result.evidence.rejected_by_previous);
        const auto saved_scale{ domain.cluster_by_key.at(key).scale };
        domain.cluster_by_key.at(key).scale.reset();
        EXPECT_FALSE(detail::EvaluateLocalCandidate(overlay,
            detail::LocalCandidateReference{detail::LocalObjectivePolicy::PreviousNonRegression, key, fixture.sample_ref_list,
                &*current, domain, evidence, counters}).accepted);
        domain.cluster_by_key.at(key).scale = saved_scale;
    }
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
        detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
        const auto options{ MakeSecondStageOptions() };
        const detail::CandidateSelectionInputs inputs{
            fixture.context, options, baseline, partition, health, fixture.state, provenance,
            proposal, activity, ridge, domain, previous, audit, trust, workspaces, corrections, counters };
        detail::CandidateTransactionBuilder builder;
        second_stage_test::BeginNumericalCapture();
        builder.Select(inputs);
        const auto capture{ second_stage_test::EndNumericalCapture() };
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
        // The factor-1 trial fails trust without building a guard snapshot. The
        // factor-1/2 and factor-1/4 trials fail guard without evaluating the
        // local objective. Each guard enters both snapshot overloads.
        EXPECT_EQ(capture.work[static_cast<std::size_t>(second_stage_test::Work::Snapshot)],
            objective_backtracked ? 8U : 6U);
        EXPECT_EQ(capture.work[static_cast<std::size_t>(second_stage_test::Work::Objective)],
            objective_backtracked ? 3U : 2U);
    }
}

TEST(EstimatorSecondStageDefenseTest, BoundaryRejectionRestoresPreviousModels)
{
    const std::vector<rg::GaussianModel3D> models{ { 6.0, 0.5, 0.0 }, { 7.0, 0.5, 0.0 } };
    auto fixture{ BuildJointPolishFixture(models, models) };
    const std::vector<detail::ClusterKey> keys{ { 0 }, { 1 } };
    detail::CouplingGraphPartition partition;
    for (const auto & key : keys) partition.sample_id_list_by_key[key] = fixture.sample_ref_list;
    partition.boundary_sample_dependency_list = { { { 0, 0 }, keys, { 0, 1 } } };
    const auto baseline{ detail::BuildResidualBaseline(fixture.context, fixture.state) };
    const auto domain{ detail::BuildObjectiveDomain(fixture.context, baseline.model_snapshot, keys) };
    const auto previous{ detail::BuildObjectiveByKey(partition, domain, baseline) };
    auto candidate{ fixture.state };
    for (auto & result : candidate) result = MakeGaussianResult({ 20.0, 0.5, 0.0 });
    const detail::PolishProvenance provenance(2, 0);
    const detail::SuspiciousBlockActivity fixed{ { 1, 1 }, { 1, 1 }, { 1, 1 } };
    const std::vector<double> ridge(2, 1.0);
    const detail::ClusterHealthMap health;
    const detail::BestAuditState audit;
    detail::TrustRegionStateSet trust;
    trust.Reconcile(keys);
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
    const auto options{ MakeSecondStageOptions() };
    detail::SecondStageObservationSession observation;
    const detail::CandidateSelectionInputs inputs{
        fixture.context, options, baseline, partition, health, fixture.state, provenance,
        candidate, fixed, ridge, domain, previous, audit, trust, workspaces, corrections, counters, &observation };
    detail::CandidateSelection selection;
    selection.block_activity = fixed;
    selection.assembled_state = candidate;
    selection.assembled_polish_provenance = provenance;
    selection.accepted_key_list = keys;
    detail::CandidateTransactionBuilder builder(std::move(selection));
    builder.ReconcileSelectedBoundaries(inputs);
    selection = builder.View();
    EXPECT_TRUE(selection.accepted_key_list.empty());
    EXPECT_EQ(selection.rejected_key_list.size(), 2U);
    for (const auto & key : keys)
    {
        ExpectGaussianModelsNear(selection.assembled_state.at(key.front()).mdpde.GetModel(), models.at(key.front()), 0.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, CandidateCommitPublishesAcceptedRejectedAndMixedTransactions)
{
    using Commit = decltype(&detail::CandidateTransaction::Commit);
    static_assert(!std::is_copy_constructible_v<detail::CandidateTransaction>);
    static_assert(std::is_invocable_v<Commit, detail::CandidateTransaction &&,
        detail::FitState &, detail::FitState &, detail::PolishProvenance &,
        detail::QuarantineState &, detail::TrustRegionStateSet &,
        detail::SecondStageObservationSession *>);
    static_assert(!std::is_invocable_v<Commit, detail::CandidateTransaction &,
        detail::FitState &, detail::FitState &, detail::PolishProvenance &,
        detail::QuarantineState &, detail::TrustRegionStateSet &,
        detail::SecondStageObservationSession *>);

    for (const std::size_t accepted_count : { 0U, 1U, 2U })
    for (const bool observe : { false, true })
    {
        SCOPED_TRACE(testing::Message() << "accepted=" << accepted_count << ", observation=" << observe);
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
        detail::PerformanceCounters counters{ true, fixture.context, workspaces, corrections };
        const auto options{ MakeSecondStageOptions() };
        const auto saved_level{ Logger::GetLogLevel() };
        Logger::SetLogLevel(LogLevel::Debug);
        detail::SecondStageObservationSession observation(!observe);
        const detail::CandidateSelectionInputs inputs{
            fixture.context, options, baseline, partition, health, fixture.state, provenance,
            candidate, activity, ridge, domain, previous_objectives, audit, trust, workspaces, corrections, counters,
            observe ? &observation : nullptr };
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
        testing::internal::CaptureStdout();
        const auto committed{ std::move(transaction).Commit(fixture.state, published, provenance, quarantine, trust,
            observe ? &observation : nullptr) };
        const auto commit_output{ testing::internal::GetCapturedStdout() };
        Logger::SetLogLevel(saved_level);
        if (observation.Enabled())
            EXPECT_EQ(observation.Audit()->batch.stages[static_cast<std::size_t>(detail::AuditStage::Commit)].total,
                keys.size() + committed.trust_region_update.changed_key_list.size() + (accepted_count == 0 ? 1U : 0U));
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

TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsMatchesSerialAndParallelSelection)
{
    auto serial_model{ BuildSeparatedRollbackDefenseModel() };
    auto parallel_model{ BuildSeparatedRollbackDefenseModel() };
    auto serial_options{ MakeSecondStageOptions() };
    auto parallel_options{ MakeSecondStageOptions() };
    serial_options.thread_size = 1;
    parallel_options.thread_size = 2;

    serial_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*serial_model, serial_options);
    parallel_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*parallel_model, parallel_options);

    const auto & serial_atoms{ serial_model->GetSelectedAtoms() };
    const auto & parallel_atoms{ parallel_model->GetSelectedAtoms() };
    ASSERT_EQ(serial_atoms.size(), parallel_atoms.size());
    for (std::size_t i = 0; i < serial_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*serial_atoms.at(i)),
            GetEstimateModel(*parallel_atoms.at(i)),
            1.0e-12);
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationMatchesSerialAndParallelSelection)
{
    auto serial_model{ BuildBoundaryComponentConflictDefenseModel() };
    auto parallel_model{ BuildBoundaryComponentConflictDefenseModel() };
    auto serial_options{ MakeSecondStageOptions() };
    auto parallel_options{ MakeSecondStageOptions() };
    serial_options.thread_size = 1;
    parallel_options.thread_size = 2;

    const auto previous_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    serial_options.quiet_mode = false;
    testing::internal::CaptureStdout();
    serial_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*serial_model, serial_options);
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_level);
    parallel_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*parallel_model, parallel_options);
    const auto cutoff_position{ output.find("Local-fitting atom cutoff: atoms=103, limit=100, clusters=") };
    ASSERT_NE(cutoff_position, std::string::npos);
    const auto maximum_position{ output.find(", max-atoms=", cutoff_position) };
    const auto cuts_position{ output.find(", cutoff-edges=", cutoff_position) };
    ASSERT_NE(maximum_position, std::string::npos);
    ASSERT_NE(cuts_position, std::string::npos);
    EXPECT_LE(std::stoull(output.substr(maximum_position + 12)), 100U);
    EXPECT_GT(std::stoull(output.substr(cuts_position + 15)), 0U);
    const auto & serial_atoms{ serial_model->GetSelectedAtoms() };
    const auto & parallel_atoms{ parallel_model->GetSelectedAtoms() };
    ASSERT_EQ(serial_atoms.size(), parallel_atoms.size());
    for (std::size_t i = 0; i < serial_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*serial_atoms.at(i)),
            GetEstimateModel(*parallel_atoms.at(i)),
            1.0e-12);
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationIsIntensityScaleInvariant)
{
    constexpr double intensity_scale{ 100.0 };
    auto base_model{ BuildBoundaryComponentConflictDefenseModel() };
    auto scaled_model{
        BuildBoundaryComponentConflictDefenseModel(intensity_scale)
    };

    base_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*base_model, MakeSecondStageOptions());
    scaled_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*scaled_model, MakeSecondStageOptions());
    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    for (std::size_t i = 0; i < base_atoms.size(); i++)
    {
        const auto base{ GetEstimateModel(*base_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        EXPECT_NEAR(
            base.GetAmplitude() * intensity_scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 5.0e-5));
        EXPECT_NEAR(base.GetWidth(), scaled.GetWidth(), 1.0e-5);
        EXPECT_NEAR(
            base.GetOffset() * intensity_scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 1.0e-4));
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryJointCorrectionMatchesSerialParallelAndIntensityScaling)
{
    constexpr double intensity_scale{ 100.0 };
    auto serial_model{ BuildBoundaryJointCorrectionDefenseModel() };
    auto parallel_model{ BuildBoundaryJointCorrectionDefenseModel() };
    auto scaled_model{
        BuildBoundaryJointCorrectionDefenseModel(intensity_scale)
    };
    auto serial_options{ MakeSecondStageOptions() };
    auto parallel_options{ MakeSecondStageOptions() };
    serial_options.thread_size = 1;
    parallel_options.thread_size = 2;

    serial_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*serial_model, serial_options);
    parallel_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*parallel_model, parallel_options);
    scaled_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*scaled_model, MakeSecondStageOptions());

    const auto & serial_atoms{ serial_model->GetSelectedAtoms() };
    const auto & parallel_atoms{ parallel_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(serial_atoms.size(), parallel_atoms.size());
    ASSERT_EQ(serial_atoms.size(), scaled_atoms.size());
    for (std::size_t i = 0; i < serial_atoms.size(); i++)
    {
        const auto serial{ GetEstimateModel(*serial_atoms.at(i)) };
        const auto parallel{ GetEstimateModel(*parallel_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        ExpectGaussianModelsNear(serial, parallel, 1.0e-12);
        EXPECT_NEAR(
            serial.GetAmplitude() * intensity_scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 5.0e-5));
        EXPECT_NEAR(serial.GetWidth(), scaled.GetWidth(), 5.0e-6);
        EXPECT_NEAR(
            serial.GetOffset() * intensity_scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 1.0e-4));
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

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageIterationsDampsOffsetStepIntoInitialTrustRadius)
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto truth_coordinates{ initial_model.ToTransformedCoordinates() };
    ASSERT_TRUE(truth_coordinates.has_value());
    (*truth_coordinates)(static_cast<Eigen::Index>(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) = 1.25;
    const auto truth_model{
        rg::GaussianModel3D::FromTransformedCoordinates(
            *truth_coordinates)
    };
    ASSERT_TRUE(truth_model.has_value());

    auto model{
        BuildDefenseModel(
            { std::array<double, 3>{ 0.0, 0.0, 0.0 } },
            { Spot::O },
            { Element::OXYGEN },
            { *truth_model },
            initial_model)
    };
    const auto previous_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fitted_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    EXPECT_NE(fitted_model.GetOffset(), previous_model.GetOffset());
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
