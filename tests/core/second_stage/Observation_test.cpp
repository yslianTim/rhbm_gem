#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <future>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/JointFitting.hpp"
#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include "core/detail/second_stage/observation/PhaseAudit.hpp"
#include "core/detail/second_stage/observation/SecondStageLogging.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/observation/TrustModelAudit.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
namespace alg = rhbm_gem::algorithm;
using rhbm_gem::FittingStage;

using second_stage_test::BuildJointPolishDefenseModel;
using second_stage_test::BuildJointPolishFixture;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeGaussianResult;
using second_stage_test::MakeSecondStageOptions;

} // namespace

TEST(EstimatorSecondStageDefenseTest, PerformanceCountersRespectLogLevelQuietAndDestruction)
{
    const auto saved_level{ Logger::GetLogLevel() };
    for (const auto level : { LogLevel::Info, LogLevel::Debug })
    {
        for (const bool quiet : { false, true })
        {
            for (const bool unwind : { false, true })
            {
                SCOPED_TRACE(::testing::Message() << static_cast<int>(level) << "/" << quiet << "/" << unwind);
                detail::SecondStageContext context;
                detail::ClusterSolverWorkspaceMap workspaces;
                detail::BoundaryJointCorrectionWorkspaceMap boundary_workspaces;
                Logger::SetLogLevel(level);
                testing::internal::CaptureStdout();
                try
                {
                    detail::PerformanceCounters counters{ quiet, context, workspaces, boundary_workspaces };
                    counters.RecordFullStateMaterialization();
                    const auto before_destruction{ testing::internal::GetCapturedStdout() };
                    testing::internal::CaptureStdout();
                    EXPECT_TRUE(before_destruction.empty());
                    if (unwind) throw std::runtime_error("counter lifetime test");
                }
                catch (const std::runtime_error &)
                {
                }
                const auto output{ testing::internal::GetCapturedStdout() };
                if (quiet)
                {
                    EXPECT_TRUE(output.empty());
                    continue;
                }
                const std::string heading{ " Second-Stage Local Fitting Performance :\n" };
                EXPECT_EQ(output.find(heading), 0);
                EXPECT_EQ(output.find(heading, heading.size()), std::string::npos);
                EXPECT_NE(output.find(" - boundary_reconciliation_ms = 0.000\n"
                    " - boundary_joint_correction_ms = 0.000\n"
                    " - dependency_polish_ms = 0.000\n"
                    " - iteration/candidate/topology/total_ms = 0.000/0.000/0.000/"), std::string::npos);
                EXPECT_EQ(output.find("full_state_materializations = 1") != std::string::npos,
                    level == LogLevel::Debug);
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                EXPECT_EQ(output.find("Trust-model performance: schema=1") != std::string::npos,
                    level == LogLevel::Debug);
#else
                EXPECT_EQ(output.find("Trust-model performance:"), std::string::npos);
#endif
            }
        }
    }
    Logger::SetLogLevel(saved_level);
}

TEST(EstimatorSecondStageDefenseTest, PerformanceCountersAccumulateParallelUpdatesAndRetiredWorkspaces)
{
    detail::SecondStageContext context;
    context.atom_list.resize(2);
    context.atom_list[0].raw_sampling_entries.resize(2);
    context.atom_list[1].raw_sampling_entries.resize(3);
    detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap boundary_workspaces;
    alg::WeightedRidgeSystem system;
    system.design_matrix.resize(1, 1);
    system.design_matrix.insert(0, 0) = 1.0;
    system.response = Eigen::VectorXd::Ones(1);
    system.previous_parameter = Eigen::VectorXd::Zero(1);
    system.ridge_diagonal = Eigen::VectorXd::Ones(1);
    const auto saved_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    {
        detail::PerformanceCounters counters{ false, context, workspaces, boundary_workspaces };
        std::vector<std::future<void>> workers;
        for (std::size_t worker = 0; worker < 4; worker++)
        {
            workers.push_back(std::async(std::launch::async, [&counters] {
                for (std::size_t i = 0; i < 100; i++)
                {
                    counters.RecordFullStateMaterialization();
                    counters.RecordGaussianCacheHits();
                    counters.RecordGaussianCacheMisses();
                    counters.RecordObjectiveSampleEvaluation(2, 5);
                    counters.RecordObjectiveSampleEvaluation(7, 5);
                }
            }));
        }
        for (auto & worker : workers) worker.get();
        for (std::size_t generation = 0; generation < 2; generation++)
        {
            auto & workspace{ workspaces[{ 0 }] };
            EXPECT_TRUE(workspace.joint_offset.AnalyzePattern(system));
            EXPECT_TRUE(workspace.joint_polish.AnalyzePattern(system));
            EXPECT_TRUE(boundary_workspaces[{}].AnalyzePattern(system));
            if (generation == 0)
            {
                counters.RecordSolverWorkspaceReset();
                workspaces.clear();
                boundary_workspaces.clear();
            }
        }
        counters.RecordTopologyRebuild(1.25, true);
        counters.RecordTopologyRebuild(2.5, false);
        counters.RecordBoundaryReconciliation(3, 2, 1, 4.25);
        counters.RecordBoundaryJointCorrection(true, 1.0);
        counters.RecordBoundaryJointCorrection(false, 2.0);
        counters.RecordBoundaryRescue(true, false);
        counters.RecordBoundaryRescue(false, true);
        counters.RecordBoundaryRescueExclusions(1, 2, 3);
        counters.RecordDependencyPolish(3, 2, 1, 1, 5, 6, 7, 8.5);
        counters.FinishIterationPhase(std::chrono::steady_clock::now());
        counters.FinishCandidatePhase(std::chrono::steady_clock::now());
    }
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(saved_level);
    EXPECT_NE(output.find(" - boundary_reconciliation_ms = 4.250\n"
        " - boundary_joint_correction_ms = 3.000\n"
        " - dependency_polish_ms = 8.500\n"), std::string::npos);
    EXPECT_NE(output.find("[Debug]  - full_state_materializations = 400\n"
        " - gaussian_cache_hit/miss = 2000/2000\n"
        " - objective_recomputed/reused_samples = 3600/1200\n"
        " - solver_symbolic_analyses = 6\n"
        " - topology_rebuilds/partition_changes = 2/1\n"
        " - boundary_reconciliations/backtracked/rejected = 3/2/1\n"
        " - boundary_joint_correction_attempts/accepted/fallback = 2/1/1\n"
        " - boundary_rescues/accepted/fallback/rejected = 2/1/1/1\n"
        " - boundary_rescue_exclusions_hard/invalid/no-objective = 1/2/3\n"
        " - dependency_polish_components/attempted/accepted/fallback = 3/2/1/1\n"
        " - dependency_polish_atoms/parameters/rounds = 5/6/7\n"), std::string::npos);
}

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
TEST(EstimatorSecondStageDefenseTest, TrustModelShadowUsesFrozenIrlsDirectionalPrediction)
{
    using Action = detail::TrustRegionRadiusAction;
    const rg::GaussianModel3D previous_model{ 5.4, 0.52, 0.05 };
    const rg::GaussianModel3D target_model{ 6.0, 0.55, 0.10 };
    auto fixture{
        BuildJointPolishFixture( { previous_model }, { target_model })
    };
    constexpr double unselected_distance{ 0.25 };
    auto & atom_context{ fixture.context.atom_list.at(0) };
    atom_context.unselected_distance_list_by_sample.assign(atom_context.raw_sampling_entries.size(),
        { unselected_distance });
    for (auto & sample : atom_context.raw_sampling_entries)
        sample.response += previous_model.ResponseAtDistance(unselected_distance);
    const detail::ClusterKey key{ 0 };
    fixture.context.frozen_background = detail::BuildFrozenBackground(fixture.context, fixture.state);
    ASSERT_TRUE(fixture.context.frozen_background);
    const auto previous_snapshot{
        detail::BuildSecondStageModelSnapshot(
            fixture.context,
            fixture.state)
    };
    const auto objective_domain{
        detail::BuildObjectiveDomain(
            fixture.context,
            previous_snapshot,
            { key })
    };
    const auto residual_baseline{
        detail::BuildResidualBaseline(fixture.context, fixture.state)
    };
    const auto previous_objective{
        detail::EvaluateAuditObjective(
            objective_domain,
            residual_baseline)
    };
    ASSERT_TRUE(previous_objective.has_value());

    const auto previous_coordinates{
        previous_model.ToTransformedCoordinates()
    };
    const auto target_coordinates{
        target_model.ToTransformedCoordinates()
    };
    ASSERT_TRUE(previous_coordinates.has_value());
    ASSERT_TRUE(target_coordinates.has_value());
    const auto candidate_model{
        rg::GaussianModel3D::FromTransformedCoordinates(
            *previous_coordinates + 0.05 *
                (*target_coordinates - *previous_coordinates))
    };
    ASSERT_TRUE(candidate_model.has_value());
    auto candidate_state{ fixture.state };
    candidate_state.at(0) = MakeGaussianResult(*candidate_model);
    const auto candidate_patch{
        detail::FitStatePatch::FromState(candidate_state, key)
    };
    const auto candidate_snapshot{
        detail::BuildSecondStageModelSnapshot(
            fixture.context,
            candidate_state)
    };
    const auto candidate_objective{
        detail::EvaluateAuditObjective(
            objective_domain,
            fixture.context, candidate_snapshot)
    };
    ASSERT_TRUE(candidate_objective.has_value());

    const auto diagnostic{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            candidate_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        diagnostic.status,
        detail::TrustModelPredictionStatus::Available);
    ASSERT_TRUE(diagnostic.actual_reduction.has_value());
    ASSERT_TRUE(diagnostic.predicted_reduction.has_value());
    ASSERT_TRUE(diagnostic.rho.has_value());
    EXPECT_GT(*diagnostic.actual_reduction, 0.0);
    EXPECT_GT(*diagnostic.predicted_reduction, 0.0);
    EXPECT_NEAR(*diagnostic.rho, 1.0, 0.10);
    EXPECT_EQ(candidate_state.size(), 1U);

    constexpr double intensity_scale{ 1.0e4 };
    auto scaled_fixture{ fixture };
    auto scaled_background{ std::make_shared<detail::FrozenBackground>(*fixture.context.frozen_background) };
    for (auto & responses : scaled_background->response_by_atom)
        for (auto & response : responses) response *= intensity_scale;
    for (auto & model : scaled_background->model_by_atom)
        model = { model.GetAmplitude() * intensity_scale, model.GetWidth(), model.GetOffset() * intensity_scale };
    scaled_fixture.context.frozen_background = std::move(scaled_background);

    for (auto & sample : scaled_fixture.context.atom_list.at(0).raw_sampling_entries)
    {
        sample.response *= intensity_scale;
    }
    const auto scale_model = [](const rg::GaussianModel3D & model)
    {
        return rg::GaussianModel3D{
            model.GetAmplitude() * intensity_scale,
            model.GetWidth(),
            model.GetOffset() * intensity_scale
        };
    };
    scaled_fixture.state.at(0) = MakeGaussianResult(scale_model(previous_model));
    auto scaled_candidate_state{ scaled_fixture.state };
    scaled_candidate_state.at(0) = MakeGaussianResult(scale_model(*candidate_model));
    const auto scaled_previous_snapshot{
        detail::BuildSecondStageModelSnapshot(
            scaled_fixture.context,
            scaled_fixture.state)
    };
    const auto scaled_domain{
        detail::BuildObjectiveDomain(
            scaled_fixture.context,
            scaled_previous_snapshot,
            { key })
    };
    const auto scaled_baseline{
        detail::BuildResidualBaseline(
            scaled_fixture.context,
            scaled_fixture.state)
    };
    const auto scaled_previous_objective{
        detail::EvaluateAuditObjective(scaled_domain, scaled_baseline)
    };
    const auto scaled_candidate_snapshot{
        detail::BuildSecondStageModelSnapshot(
            scaled_fixture.context,
            scaled_candidate_state)
    };
    const auto scaled_candidate_objective{
        detail::EvaluateAuditObjective(
            scaled_domain,
            scaled_fixture.context, scaled_candidate_snapshot)
    };
    ASSERT_TRUE(scaled_previous_objective.has_value());
    ASSERT_TRUE(scaled_candidate_objective.has_value());
    const auto scaled_diagnostic{
        detail::EvaluateTrustModelShadow(
            scaled_fixture.context,
            scaled_baseline,
            scaled_fixture.state,
            detail::FitStatePatch::FromState(
                scaled_candidate_state,
                key),
            key,
            scaled_fixture.sample_ref_list,
            scaled_domain,
            scaled_previous_objective,
            scaled_candidate_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    ASSERT_TRUE(scaled_diagnostic.predicted_reduction.has_value());
    ASSERT_TRUE(scaled_diagnostic.rho.has_value());
    EXPECT_NEAR(
        *scaled_diagnostic.predicted_reduction,
        *diagnostic.predicted_reduction,
        1.0e-6);
    EXPECT_NEAR(*scaled_diagnostic.rho, *diagnostic.rho, 1.0e-6);

    auto previous_with_penalty{ *previous_objective };
    auto candidate_with_penalty{ *candidate_objective };
    previous_with_penalty.offset_plausibility_penalty = 0.20;
    candidate_with_penalty.offset_plausibility_penalty = 0.05;
    const auto penalty_diagnostic{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_with_penalty,
            candidate_with_penalty,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Polish,
            false)
    };
    ASSERT_TRUE(penalty_diagnostic.predicted_reduction.has_value());
    EXPECT_NEAR(
        *penalty_diagnostic.predicted_reduction -
            *diagnostic.predicted_reduction,
        0.15,
        1.0e-12);

    // Overlap is an additional objective role, not another residual evaluation.
    auto tail_domain{ objective_domain };
    auto & tail_cluster{ tail_domain.cluster_by_key.at(key) };
    tail_cluster.tail_sample_ref_list = tail_cluster.fit_sample_ref_list;
    tail_cluster.scale->tail = tail_cluster.scale->fit;
    tail_domain.tail_sample_mask_by_atom = tail_domain.fit_sample_mask_by_atom;
    tail_domain.tail_sample_count = tail_domain.fit_sample_count;
    const auto tail_previous_objective{
        detail::EvaluateAuditObjective(
            tail_domain,
            residual_baseline)
    };
    const auto tail_candidate_objective{
        detail::EvaluateAuditObjective(
            tail_domain,
            fixture.context, candidate_snapshot)
    };
    ASSERT_TRUE(tail_previous_objective.has_value());
    ASSERT_TRUE(tail_candidate_objective.has_value());
    const auto tail_diagnostic{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            tail_domain,
            tail_previous_objective,
            tail_candidate_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        tail_diagnostic.status,
        detail::TrustModelPredictionStatus::Available);
    ASSERT_TRUE(tail_diagnostic.rho.has_value());
    EXPECT_NEAR(*tail_diagnostic.rho, 1.0, 0.10);
    ASSERT_TRUE(tail_diagnostic.predicted_residual_reduction.has_value());
    ASSERT_TRUE(diagnostic.predicted_residual_reduction.has_value());
    EXPECT_NEAR(*tail_diagnostic.predicted_residual_reduction,
        (1.0 + detail::kTailValidationWeight) * *diagnostic.predicted_residual_reduction,
        1.0e-12);
    EXPECT_EQ(tail_domain.unique_sample_count, objective_domain.unique_sample_count);


    auto nonmaterial_prediction_domain{ objective_domain };
    nonmaterial_prediction_domain.cluster_by_key.at(key).scale->fit *= 1.0e6;
    const auto nonmaterial_prediction{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            nonmaterial_prediction_domain,
            previous_objective,
            candidate_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonmaterial_prediction.status,
        detail::TrustModelPredictionStatus::NonmaterialPrediction);
    EXPECT_FALSE(nonmaterial_prediction.rho.has_value());

    auto nonfinite_baseline{ residual_baseline };
    nonfinite_baseline.sample_list.at(0).at(0)->residual =
        std::numeric_limits<double>::infinity();
    const auto nonfinite_prediction{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            nonfinite_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            candidate_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonfinite_prediction.status,
        detail::TrustModelPredictionStatus::Nonfinite);
    EXPECT_FALSE(nonfinite_prediction.rho.has_value());

    const auto opposite_model{
        rg::GaussianModel3D::FromTransformedCoordinates(
            *previous_coordinates - 0.05 *
                (*target_coordinates - *previous_coordinates))
    };
    ASSERT_TRUE(opposite_model.has_value());
    auto opposite_state{ fixture.state };
    opposite_state.at(0) = MakeGaussianResult(*opposite_model);
    const auto opposite_patch{
        detail::FitStatePatch::FromState(opposite_state, key)
    };
    const auto opposite_snapshot{
        detail::BuildSecondStageModelSnapshot(
            fixture.context,
            opposite_state)
    };
    const auto opposite_objective{
        detail::EvaluateAuditObjective(
            objective_domain,
            fixture.context, opposite_snapshot)
    };
    ASSERT_TRUE(opposite_objective.has_value());
    const auto nonpositive_prediction{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            opposite_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            opposite_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonpositive_prediction.status,
        detail::TrustModelPredictionStatus::NonpositivePrediction);
    EXPECT_FALSE(nonpositive_prediction.rho.has_value());

    const auto nonmaterial{
        detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            detail::FitStatePatch::FromState(fixture.state, key),
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            previous_objective,
            1.0,
            Action::Keep,
            detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonmaterial.status,
        detail::TrustModelPredictionStatus::NonmaterialStep);
    EXPECT_FALSE(nonmaterial.rho.has_value());
}
#endif

TEST(EstimatorSecondStageDefenseTest, BestReferenceDiagnosticsReportTheEvaluatedGate)
{
    detail::JointCandidateObjectiveDiagnostic record;
    record.source = "endpoint";
    record.stored_best = detail::ObjectiveBreakdown{ 1.0, 0.0, 0.0 };
    detail::RecordJointMemberRejection(&record, { 0 },
        detail::ObjectiveBreakdown{ 4.0, 0.0, 0.0 },
        detail::ObjectiveBreakdown{ 2.0, 0.0, 0.0 },
        detail::ObjectiveBreakdown{ 3.0, 0.0, 0.0 }, true);
    EXPECT_EQ(record.outcome, "best");
    detail::IterationObservation observation;
    auto & boundary{ observation.boundary_reconciliation_diagnostic_list.emplace_back() };
    boundary.key_list = { { 0 }, { 1 } };
    boundary.objective_diagnostic_list.push_back(record);
    const auto saved_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    detail::LogAcceptedCandidateSearchDiagnostics(false, observation);
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(saved_level);
    EXPECT_NE(output.find("Joint candidate objective rejection: schema=2"), std::string::npos);
    EXPECT_NE(output.find("reference-environment=candidate"), std::string::npos);
    const auto reference_position{ output.find("best-reference=") };
    const auto stored_position{ output.find("stored-best=") };
    ASSERT_NE(reference_position, std::string::npos);
    ASSERT_NE(stored_position, std::string::npos);
    EXPECT_DOUBLE_EQ(std::stod(output.substr(reference_position + std::string("best-reference=").size())), 2.0);
    EXPECT_DOUBLE_EQ(std::stod(output.substr(stored_position + std::string("stored-best=").size())), 1.0);
    detail::RecordJointMemberRejection(&record, { 0 },
        detail::ObjectiveBreakdown{ 4.0, 0.0, 0.0 }, std::nullopt,
        detail::ObjectiveBreakdown{ 3.0, 0.0, 0.0 }, true);
    EXPECT_EQ(record.outcome, "best-reference-unavailable");
}

TEST(EstimatorSecondStageDefenseTest, PhaseAuditSnapshotsUseWholeFrozenDomainAndDoNotPublishCandidates)
{
    detail::SecondStageContext context;
    context.atom_list.resize(3);
    const detail::FitState baseline{
        MakeGaussianResult({ 6.0, 0.5, 0.05 }), MakeGaussianResult({ 7.0, 0.6, -0.02 }),
        MakeGaussianResult({ 8.0, 0.5, 0.03 }) };
    for (std::size_t atom = 0; atom < 3; atom++)
    {
        auto & target{ context.atom_list.at(atom) };
        for (const auto distance : { 0.0, 0.25, 0.5, 0.75, 1.5 })
        {
            target.neighbor_atom_sample_offset_list.push_back(target.neighbor_atom_sample_list.size());
            double response{ baseline.at(atom).mdpde.GetModel().ResponseAtDistance(distance) + 0.1 };
            for (std::size_t neighbor = 0; neighbor < 3; neighbor++)
                if (neighbor != atom)
                {
                    target.neighbor_atom_sample_list.push_back({ neighbor, 1.0 });
                    response += baseline.at(neighbor).mdpde.GetModel().ResponseAtDistance(1.0);
                }
            target.raw_sampling_entries.push_back({ response, SamplingPoint{ distance } });
        }
        target.neighbor_atom_sample_offset_list.push_back(target.neighbor_atom_sample_list.size());
        target.refit_design = detail::PreparedLocalGaussianDesign(target.raw_sampling_entries, 0.0, 1.0);
    }
    auto background{ std::make_shared<detail::FrozenBackground>() };
    background->response_by_atom.assign(3, std::vector<double>(5, 0.025));
    context.frozen_background = background;
    const std::vector<detail::ClusterKey> keys{ { 0, 1 }, { 2 } };
    const auto domain{ detail::BuildObjectiveDomain(context,
        detail::BuildSecondStageModelSnapshot(context, baseline), keys) };
    const detail::SuspiciousBlockActivity activity{ std::vector<char>(3, 0), std::vector<char>(3, 0), std::vector<char>(3, 0) };
    detail::ClusterSolverWorkspaceMap workspaces;
    for (const auto & key : keys) workspaces.try_emplace(key);
    auto options{ MakeSecondStageOptions() };
    const auto production{ detail::BuildIterationProposal(context, keys, baseline, options,
        std::vector<double>(3, 1.0), activity, workspaces) };
    const auto analyses{ workspaces.at(keys.front()).joint_offset.GetSymbolicAnalysisCount() };
    auto candidate{ baseline };
    candidate[0] = MakeGaussianResult({ 6.1, 0.51, 0.06 });
    candidate[2] = MakeGaussianResult({ 7.9, 0.5, 0.03 });
    const auto patch{ detail::FitStatePatch::FromState(candidate, { 0, 2 }) };
    const detail::FitStateView view{ baseline, patch };
    const auto saved_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    detail::PhaseAudit collector(context, domain, baseline, keys, 6, 12);
    collector.Capture("local-polish", { 0, 2 }, view, nullptr, 1.0, "rejected", "best", true);
    collector.CaptureState("assembly-after-polish", candidate, true);
    collector.Missing("local-search", { 2 }, "search-exhausted");
    auto incomplete{ production.fixed_point_operator };
    incomplete.shape_available_atom_mask[0] = 0;
    collector.CaptureOperator(incomplete);
    testing::internal::CaptureStdout();
    collector.Finish(options, std::vector<double>(3, 1.0), activity, production, baseline);
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(saved_level);
    EXPECT_EQ(workspaces.at(keys.front()).joint_offset.GetSymbolicAnalysisCount(), analyses);
    EXPECT_DOUBLE_EQ(baseline[0].mdpde.GetModel().GetAmplitude(), 6.0);
    const auto independent{ detail::EvaluateAuditObjective(domain,
        detail::BuildResidualBaseline(context, candidate)) };
    ASSERT_TRUE(independent);
    // Compare every serialized component with an independently materialized residual baseline.
    for (const auto & [field, value] : std::vector<std::pair<std::string, double>>{
        { "fit", independent->fit_range_residual_objective },
        { "tail_weighted", independent->GetTailValidationPenalty() },
        { "offset", independent->offset_plausibility_penalty },
        { "total", independent->GetTotalObjective() } })
    {
        const auto line_begin{ output.find("\"stage\":\"local-polish\"") };
        ASSERT_NE(line_begin, std::string::npos);
        const auto value_begin{ output.find("\"" + field + "\":", line_begin) };
        ASSERT_NE(value_begin, std::string::npos);
        EXPECT_NEAR(std::stod(output.substr(value_begin + field.size() + 3)), value, 1e-12);
    }
    EXPECT_NE(output.find("\"disposition\":\"rejected\",\"reason\":\"best\",\"final_retained\":false"), std::string::npos);
    EXPECT_NE(output.find("\"alpha\":0.001"), std::string::npos);
    // Independently rerun T(candidate), then compare its residual (not the candidate movement).
    detail::ClusterSolverWorkspaceMap candidate_workspaces;
    for (const auto & key : keys) candidate_workspaces.try_emplace(key);
    const auto candidate_operator{ detail::BuildIterationProposal(context, keys, candidate, options,
        std::vector<double>(3, 1.0), activity, candidate_workspaces) };
    std::vector<detail::TransformedChange> residuals;
    for (std::size_t i = 0; i < candidate.size(); i++)
        residuals.push_back(detail::CalculateTransformedChange(
            candidate_operator.fixed_point_operator.state[i], candidate[i].mdpde.GetModel()));
    const auto expected_residual{ detail::SummarizeActiveDofChanges(residuals, { { 0, 1, 2 }, { 0, 1, 2 } }) };
    const auto polish_position{ output.find("\"stage\":\"local-polish\"") };
    auto residual_position{ output.find("\"p99\":[", polish_position) };
    ASSERT_NE(residual_position, std::string::npos);
    residual_position += 7;
    for (std::size_t coordinate = 0; coordinate < 3; coordinate++)
    {
        std::size_t consumed{ 0 };
        EXPECT_NEAR(std::stod(output.substr(residual_position), &consumed), expected_residual.percentile_list[coordinate], 1e-12);
        residual_position += consumed + 1;
    }
    auto half{ baseline };
    for (std::size_t i = 0; i < half.size(); i++)
    {
        const auto a{ baseline[i].mdpde.GetModel() }, b{ candidate[i].mdpde.GetModel() };
        half[i] = MakeGaussianResult({ std::sqrt(a.GetAmplitude() * b.GetAmplitude()),
            std::sqrt(a.GetWidth() * b.GetWidth()), 0.5 * (a.GetOffset() + b.GetOffset()) });
    }
    const auto half_objective{ detail::EvaluateAuditObjective(domain, detail::BuildResidualBaseline(context, half)) };
    ASSERT_TRUE(half_objective);
    const auto half_position{ output.find("\"alpha\":0.5", polish_position) };
    ASSERT_NE(half_position, std::string::npos);
    const auto half_total{ output.find("\"total\":", half_position) };
    ASSERT_NE(half_total, std::string::npos);
    EXPECT_NEAR(std::stod(output.substr(half_total + 8)), half_objective->GetTotalObjective(), 1e-12);

    EXPECT_NE(output.find("search-exhausted"), std::string::npos);
    EXPECT_NE(output.find("incomplete-production-operator"), std::string::npos);
    EXPECT_NE(output.find("\"population\":[3,3,3]"), std::string::npos);
    EXPECT_NE(output.find("\"production_operator\":"), std::string::npos);
    EXPECT_NE(output.find("\"domain_id\":12"), std::string::npos);
    const auto collect_workers = [&](bool parallel)
    {
        detail::PhaseAudit workers(context, domain, baseline, keys, 7, 12);
        const auto capture = [&](const detail::ClusterKey & key)
        {
            const auto worker_patch{ detail::FitStatePatch::FromState(candidate, key) };
            const detail::FitStateView worker_view{ baseline, worker_patch };
            workers.Capture("local-search", key, worker_view, nullptr, 1.0, "accepted");
            workers.Capture("local-polish", key, worker_view, &worker_view, 1.0, "rejected", "strict-improvement", true);
        };
        if (parallel)
        {
            auto first{ std::async(std::launch::async, capture, keys[0]) };
            auto second{ std::async(std::launch::async, capture, keys[1]) };
            first.get(); second.get();
        }
        else for (const auto & key : keys) capture(key);
        workers.CaptureSearchAssembly();
        Logger::SetLogLevel(LogLevel::Debug);
        testing::internal::CaptureStdout();
        workers.Finish(options, std::vector<double>(3, 1.0), activity, production, baseline);
        auto log{ testing::internal::GetCapturedStdout() };
        Logger::SetLogLevel(saved_level);
        log.resize(log.find("[Debug] Second-stage phase audit counters:"));
        return log;
    };
    EXPECT_EQ(collect_workers(false), collect_workers(true));
    auto incomplete_context{ context };
    incomplete_context.atom_list[0].refit_design = {};
    detail::PhaseAudit failed_operator(incomplete_context, domain, baseline, keys, 8, 12);
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    failed_operator.Finish(options, std::vector<double>(3, 1.0), activity, production, baseline);
    const auto failed_output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(saved_level);
    EXPECT_NE(failed_output.find("\"reason\":\"incomplete-operator\""), std::string::npos);
    EXPECT_NE(failed_output.find("\"p99\":null"), std::string::npos);
    EXPECT_NE(failed_output.find("\"operator_reproduced\":false"), std::string::npos);
    EXPECT_EQ(workspaces.at(keys.front()).joint_offset.GetSymbolicAnalysisCount(), analyses);

    EXPECT_FALSE(detail::BeginPhaseAudit(context, true, domain, baseline, keys, 1, 1));
    Logger::SetLogLevel(LogLevel::Info);
    EXPECT_FALSE(detail::BeginPhaseAudit(context, false, domain, baseline, keys, 1, 1));
    Logger::SetLogLevel(saved_level);
}

TEST(EstimatorSecondStageDefenseTest, PhaseAuditQuietAndEnabledRunsPreserveFinalParameters)
{
    auto quiet_model{ BuildJointPolishDefenseModel() };
    auto traced_model{ BuildJointPolishDefenseModel() };
    quiet_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    traced_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    const auto saved_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    auto options{ MakeSecondStageOptions() };
    testing::internal::CaptureStdout();
    detail::RunSecondStageIterations(*quiet_model, options);
    const auto quiet_output{ testing::internal::GetCapturedStdout() };
    options.quiet_mode = false;
    testing::internal::CaptureStdout();
    detail::RunSecondStageIterations(*traced_model, options);
    const auto traced_output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(saved_level);
    EXPECT_EQ(quiet_output.find("Second-stage phase audit:"), std::string::npos);
    EXPECT_EQ(traced_output.find("Second-stage phase audit error:"), std::string::npos);
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE
    EXPECT_NE(traced_output.find("Second-stage phase audit:"), std::string::npos);
#else
    EXPECT_EQ(traced_output.find("Second-stage phase audit:"), std::string::npos);
#endif
    EXPECT_EQ(quiet_output.find("Trust-model funnel:"), std::string::npos);
    EXPECT_EQ(quiet_output.find("Trust-model shadow:"), std::string::npos);
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    EXPECT_NE(traced_output.find("Trust-model funnel:"), std::string::npos);
    EXPECT_NE(traced_output.find("Trust-model shadow:"), std::string::npos);
#else
    EXPECT_EQ(traced_output.find("Trust-model funnel:"), std::string::npos);
    EXPECT_EQ(traced_output.find("Trust-model shadow:"), std::string::npos);
#endif
    const auto & quiet_atoms{ quiet_model->GetSelectedAtoms() };
    const auto & traced_atoms{ traced_model->GetSelectedAtoms() };
    ASSERT_EQ(quiet_atoms.size(), traced_atoms.size());
    for (std::size_t i = 0; i < quiet_atoms.size(); i++)
    {
        const auto a{ GetEstimateModel(*quiet_atoms[i]) }, b{ GetEstimateModel(*traced_atoms[i]) };
        EXPECT_DOUBLE_EQ(a.GetAmplitude(), b.GetAmplitude());
        EXPECT_DOUBLE_EQ(a.GetWidth(), b.GetWidth());
        EXPECT_DOUBLE_EQ(a.GetOffset(), b.GetOffset());
    }
}
