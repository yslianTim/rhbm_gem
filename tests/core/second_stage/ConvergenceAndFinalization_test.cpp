#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include "core/detail/second_stage/CouplingGraph.hpp"
#include "core/detail/second_stage/DependencyPolish.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/JointFitting.hpp"
#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
namespace rt = rhbm_gem::core;
using rhbm_gem::FittingStage;

using second_stage_test::BuildIndependentOffsetDefenseModel;
using second_stage_test::BuildJointPolishDefenseModel;
using second_stage_test::BuildJointPolishFixture;
using second_stage_test::BuildNearCollinearDefenseModel;
using second_stage_test::CalculateSelectedAtomResponseMeanSquaredError;
using second_stage_test::Distance;
using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::ExpectSelectedAtomEstimatesAreFinite;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeSecondStageOptions;

void ExpectPeelingSamplingEntriesMatchFinalModels(
    const rg::ModelObject & model)
{
    constexpr double neighbor_atom_search_range{ 5.0 };
    constexpr double neighbor_contribution_distance_max{ 2.5 };
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    for (const auto * target_atom : selected_atoms)
    {
        const auto target_view{
            rg::AtomLocalPotentialView::For(*target_atom)
        };
        const auto raw_sampling_entries{
            target_view.GetRawSamplingEntries(false)
        };
        const auto peeling_sampling_entries{
            target_view.GetPeelingSamplingEntries(false)
        };
        ASSERT_EQ(peeling_sampling_entries.size(), raw_sampling_entries.size());
        for (std::size_t sample_index = 0;
            sample_index < raw_sampling_entries.size();
            sample_index++)
        {
            const auto & raw_sample{ raw_sampling_entries.at(sample_index) };
            const auto & peeling_sample{
                peeling_sampling_entries.at(sample_index)
            };
            auto expected_response{ raw_sample.response };
            for (const auto * neighbor_atom : selected_atoms)
            {
                if (neighbor_atom == target_atom) continue;
                if (Distance(
                        target_atom->GetPosition(),
                        neighbor_atom->GetPosition()) >
                    neighbor_atom_search_range)
                {
                    continue;
                }
                const auto sample_distance{
                    Distance(
                        raw_sample.point.position,
                        neighbor_atom->GetPosition())
                };
                if (sample_distance > neighbor_contribution_distance_max)
                {
                    continue;
                }
                expected_response -= GetEstimateModel(*neighbor_atom)
                    .ResponseAtDistance(sample_distance);
            }

            EXPECT_DOUBLE_EQ(
                peeling_sample.response,
                expected_response);
            EXPECT_DOUBLE_EQ(
                peeling_sample.point.distance,
                raw_sample.point.distance);
            EXPECT_EQ(
                peeling_sample.point.position,
                raw_sample.point.position);
            EXPECT_EQ(
                peeling_sample.point.is_selected,
                raw_sample.point.is_selected);
        }
    }
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, DependencyPolishDefaultsAndIterationValidation)
{
    const rt::FitOptions defaults;
    EXPECT_EQ(defaults.second_stage_boundary_halo_depth, 1U);
    EXPECT_TRUE(defaults.enable_second_stage_dependency_polish);
    EXPECT_EQ(defaults.second_stage_dependency_polish_max_iterations, 10U);

    auto model{ BuildJointPolishDefenseModel() };
    std::vector<rg::GaussianModel3D> original_model_list;
    for (const auto * atom : model->GetSelectedAtoms())
    {
        original_model_list.emplace_back(GetEstimateModel(*atom));
    }
    auto options{ MakeSecondStageOptions() };
    options.second_stage_dependency_polish_max_iterations = 0;
    EXPECT_THROW(
        detail::RunSecondStageIterations(*model, options),
        std::invalid_argument);
    for (std::size_t atom_index = 0;
        atom_index < model->GetSelectedAtoms().size();
        atom_index++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*model->GetSelectedAtoms().at(atom_index)),
            original_model_list.at(atom_index),
            0.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, FinalDependencyPolishImprovesUncutComponent)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        { 6.0, 0.55, 0.0 },
        { 4.5, 0.70, 0.0 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        { 6.5, 0.60, 0.15 },
        { 4.0, 0.65, -0.15 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    detail::GraphTopology topology;
    topology.adjacency_list.resize(2);
    for (const auto & sample_ref : fixture.sample_ref_list)
    {
        topology.sample_dependency_list.emplace_back(
            detail::GraphSampleDependency{
                sample_ref,
                sample_ref == fixture.sample_ref_list.front() ?
                    std::vector<std::size_t>{ 0, 1 } :
                    std::vector<std::size_t>{ sample_ref.atom_index }
            });
    }
    const auto partition{
        detail::BuildGraphPartition(topology, { 0, 1 })
    };
    const auto base_snapshot{
        detail::BuildSecondStageModelSnapshot(
            fixture.context,
            fixture.state)
    };
    const auto objective_domain{
        detail::BuildObjectiveDomain(
            fixture.context,
            base_snapshot,
            detail::BuildGraphClusterKeyList(partition))
    };
    detail::ClusterSolverWorkspaceMap solver_workspace_by_key;
    detail::BoundaryJointCorrectionWorkspaceMap correction_workspace_by_key;
    detail::PerformanceCounters performance_counters{
        true,
        fixture.context,
        solver_workspace_by_key,
        correction_workspace_by_key
    };
    auto options{ MakeSecondStageOptions() };
    const detail::SuspiciousBlockActivity all_active{
        std::vector<char>(fixture.context.atom_list.size(), 0),
        std::vector<char>(fixture.context.atom_list.size(), 0),
        std::vector<char>(fixture.context.atom_list.size(), 0)
    };
    detail::SecondStageObservationSession observation;
    const auto polish_result{
        detail::RunFinalDependencyPolish(
            fixture.context,
            options,
            topology,
            partition,
            objective_domain,
            all_active,
            fixture.state,
            correction_workspace_by_key,
            performance_counters, &observation)
    };
    ASSERT_TRUE(polish_result.accepted);
    ASSERT_TRUE(polish_result.objective.has_value());
    ASSERT_TRUE(observation.final_polish.objective_before.has_value());
    ASSERT_TRUE(observation.final_polish.objective_after.has_value());
    EXPECT_LT(
        *observation.final_polish.objective_after,
        *observation.final_polish.objective_before);
    ASSERT_EQ(observation.final_polish.component_list.size(), 1U);
    EXPECT_GE(observation.final_polish.component_list.front().round_count, 1U);
    EXPECT_LE(
        observation.final_polish.component_list.front().round_count,
        options.second_stage_dependency_polish_max_iterations);
    EXPECT_EQ(observation.final_polish.component_list.front().parameter_count, 6U);
    EXPECT_NE(
        polish_result.state.at(0).mdpde.GetModel().GetOffset(),
        polish_result.state.at(1).mdpde.GetModel().GetOffset());
    EXPECT_NE(
        polish_result.state.at(0).mdpde.GetModel().GetAmplitude(),
        fixture.state.at(0).mdpde.GetModel().GetAmplitude());
}

TEST(EstimatorSecondStageDefenseTest, TransformedConvergenceIgnoresHiddenMaximumTail)
{
    std::vector<detail::TransformedChange> change_list(
        1000,
        detail::TransformedChange{});
    change_list.back().at(
        rg::GaussianModel3D::LogPeakHeightCoordinateIndex()) =
        2.0e-3;
    const auto summary{ detail::SummarizeTransformedChanges(change_list) };
    EXPECT_LT(
        summary.percentile_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1.0e-4);
    EXPECT_GT(
        summary.maximum_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1.0e-3);
    EXPECT_TRUE(detail::IsTransformedPercentileConverged(summary));
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateKeepsAcceptedResidualIndependent)
{
    const auto make_summary = [](double percentile, double maximum)
    {
        detail::TransformedChangeSummary summary;
        summary.percentile_list.fill(percentile);
        summary.maximum_list.fill(maximum);
        summary.population_size_list.fill(100);
        return summary;
    };
    const auto small{ make_summary(5.0e-5, 2.0e-3) };
    const auto large{ make_summary(2.0e-4, 2.0e-3) };
    detail::ConvergenceCertificate certificate;
    certificate.accepted_active_p99 = small.percentile_list;
    certificate.operator_nominal_p99 = large.percentile_list;
    certificate.solver_qualified = true;
    EXPECT_FALSE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest,
    ActiveCoordinatePopulationKeepsMovingSelectedNodesWithoutFixedDilution)
{
    constexpr std::size_t atom_count{ 1000 };
    constexpr std::size_t fixed_atom_count{ 990 };
    std::vector<detail::TransformedChange> change_list(
        atom_count,
        detail::TransformedChange{});
    for (std::size_t i = fixed_atom_count; i < atom_count; i++)
    {
        change_list.at(i).fill(5.0e-4);
    }
    std::vector<std::size_t> atom_index_list(atom_count);
    for (std::size_t i = 0; i < atom_count; i++) atom_index_list.at(i) = i;

    detail::SuspiciousBlockActivity block_activity{
        detail::SuspiciousUpdateMask(atom_count, 0),
        detail::SuspiciousUpdateMask(atom_count, 0),
        detail::SuspiciousUpdateMask(atom_count, 0)
    };
    for (std::size_t i = 0; i < fixed_atom_count; i++)
    {
        block_activity.shape_fixed_atom_mask.at(i) = 1;
        block_activity.offset_fixed_atom_mask.at(i) = 1;
    }
    const auto all_selected{
        detail::SummarizeTransformedChanges(change_list)
    };

    EXPECT_TRUE(detail::IsTransformedPercentileConverged(all_selected));

    const auto population{
        detail::BuildActiveCoordinatePopulation(
            atom_index_list,
            block_activity)
    };
    const auto dual_shadow{
        detail::SummarizeActiveDofChanges(change_list, population)
    };
    EXPECT_FALSE(detail::IsTransformedPercentileConverged(dual_shadow));
    EXPECT_EQ(
        dual_shadow.population_size_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        10U);
    EXPECT_EQ(
        dual_shadow.population_size_list.at(
            rg::GaussianModel3D::LogWidthCoordinateIndex()),
        10U);
    EXPECT_EQ(
        dual_shadow.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        10U);

    std::vector<detail::SuspiciousGaussianAssessment> assessment_by_atom(
        atom_count);
    const auto benign_fixed_mask{
        detail::BuildSuspiciousFailureAtomMask(
            block_activity,
            assessment_by_atom)
    };
    EXPECT_EQ(
        std::ranges::count(benign_fixed_mask, 1),
        0);
    assessment_by_atom.front().reason =
        detail::SuspiciousGaussianReason::WidthGrowth;
    const auto suspicious_fixed_mask{
        detail::BuildSuspiciousFailureAtomMask(
            block_activity,
            assessment_by_atom)
    };
    EXPECT_EQ(
        std::ranges::count(suspicious_fixed_mask, 1),
        1);
}

TEST(EstimatorSecondStageDefenseTest,
    ActiveCoordinatePopulationCountsEverySelectedOffset)
{
    constexpr std::size_t stable_atom_count{ 100 };
    constexpr std::size_t atom_count{ stable_atom_count + 1 };
    std::vector<detail::TransformedChange> change_list(
        atom_count,
        detail::TransformedChange{});
    change_list.back().at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        5.0e-4;
    std::vector<std::size_t> atom_index_list(atom_count);
    std::iota(atom_index_list.begin(), atom_index_list.end(), 0);
    detail::SuspiciousBlockActivity activity{
        detail::SuspiciousUpdateMask(atom_count, 1),
        detail::SuspiciousUpdateMask(atom_count, 0),
        detail::SuspiciousUpdateMask(atom_count, 0)
    };
    const auto population{
        detail::BuildActiveCoordinatePopulation(
            atom_index_list,
            activity)
    };
    const auto audit{
        detail::SummarizeActiveDofChanges(change_list, population)
    };

    EXPECT_EQ(
        audit.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        atom_count);
    EXPECT_TRUE(detail::IsTransformedPercentileConverged(audit));
    EXPECT_TRUE(detail::IsTransformedChangeMaterial(change_list.back(), 1.0e-4));
}

TEST(EstimatorSecondStageDefenseTest, ActiveCoordinatePopulationPreservesExtremeAndNonFiniteMembers)
{
    const std::vector<std::size_t> atom_index_list{ 0, 1, 2 };
    detail::SuspiciousBlockActivity activity{
        detail::SuspiciousUpdateMask(3, 1),
        detail::SuspiciousUpdateMask(3, 0),
        detail::SuspiciousUpdateMask(3, 0)
    };
    const auto population{
        detail::BuildActiveCoordinatePopulation(
            atom_index_list,
            activity)
    };
    std::vector<detail::TransformedChange> change_list(
        3,
        detail::TransformedChange{});
    change_list.at(1).at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        2.0e-3;
    const auto extreme{
        detail::SummarizeActiveDofChanges(change_list, population)
    };
    EXPECT_DOUBLE_EQ(
        extreme.maximum_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        2.0e-3);
    EXPECT_FALSE(detail::IsTransformedPercentileConverged(extreme));

    change_list.at(1).at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite{
        detail::SummarizeActiveDofChanges(change_list, population)
    };
    EXPECT_TRUE(std::isinf(nonfinite.maximum_list.at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())));
    EXPECT_FALSE(detail::IsTransformedPercentileConverged(nonfinite));

    activity.shape_fixed_atom_mask = { 1, 0, 1 };
    activity.offset_fixed_atom_mask.assign(3, 1);
    const auto shape_population{ detail::BuildActiveCoordinatePopulation(
        atom_index_list, activity) };
    EXPECT_EQ(shape_population.active_shape_atom_index_list, (std::vector<std::size_t>{ 1 }));
    for (const auto parameter_index : std::array<std::size_t, 2>{
        rg::GaussianModel3D::LogPeakHeightCoordinateIndex(),
        rg::GaussianModel3D::LogWidthCoordinateIndex() })
    {
        change_list.assign(3, {});
        change_list.at(1).at(parameter_index) = std::numeric_limits<double>::quiet_NaN();
        const auto shape_nonfinite{ detail::SummarizeActiveDofChanges(change_list, shape_population) };
        EXPECT_EQ(shape_nonfinite.population_size_list.at(parameter_index), 1U);
        EXPECT_FALSE(std::isfinite(shape_nonfinite.percentile_list.at(parameter_index)));
        EXPECT_FALSE(detail::IsTransformedPercentileConverged(shape_nonfinite));
        EXPECT_EQ(shape_nonfinite.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()), 0U);
    }
}

TEST(EstimatorSecondStageDefenseTest, IndependentOffsetActivityRequiresItsOwnClusterQualification)
{
    const std::vector<std::size_t> atoms{ 0, 1 };
    detail::SuspiciousBlockActivity activity{
        { 1, 1 }, { 0, 1 }, { 0, 0 } };
    const auto population{ detail::BuildActiveCoordinatePopulation(atoms, activity) };
    EXPECT_EQ(population.active_offset_atom_index_list, (std::vector<std::size_t>{ 0 }));
    const std::vector<detail::TransformedChange> changes(2);
    const auto summary{ detail::SummarizeActiveDofChanges(changes, population) };
    EXPECT_TRUE(detail::IsTransformedPercentileConverged(summary));
    const std::vector<std::optional<rg::RHBMEstimationStatus>> status(2);
    EXPECT_FALSE(detail::AreActiveCoordinatesSolverQualified(
        atoms, { atoms }, activity, status, {}));
    detail::ClusterHealthMap health;
    health.emplace(atoms, detail::ClusterHealth{
        detail::JointOffsetSolveStatus::Converged });
    EXPECT_TRUE(detail::AreActiveCoordinatesSolverQualified(
        atoms, { atoms }, activity, status, health));
    EXPECT_FALSE(detail::AreActiveCoordinatesSolverQualified(
        atoms, {}, activity, status, health));
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateSeparatesAcceptedAndNominalPopulations)
{
    const std::vector<std::size_t> atom_index_list{ 0, 1, 2 };
    detail::SuspiciousBlockActivity activity{
        detail::SuspiciousUpdateMask{ 1, 1, 0 },
        detail::SuspiciousUpdateMask{ 1, 1, 0 },
        detail::SuspiciousUpdateMask(3, 0)
    };
    detail::SuspiciousBlockActivity nominal_activity{
        detail::SuspiciousUpdateMask(3, 0),
        detail::SuspiciousUpdateMask(3, 0),
        detail::SuspiciousUpdateMask(3, 0)
    };
    const auto accepted_population{ detail::BuildActiveCoordinatePopulation(
        atom_index_list,
        activity) };
    const auto nominal_population{ detail::BuildActiveCoordinatePopulation(
        atom_index_list,
        nominal_activity) };
    std::vector<detail::TransformedChange> changes(
        3, detail::TransformedChange{});
    changes.at(0).fill(2.0e-3);
    changes.at(1).fill(2.0e-3);
    changes.at(2).fill(5.0e-5);

    detail::ConvergenceCertificate certificate;
    detail::ConvergenceDiagnostics diagnostics;
    diagnostics.accepted_active_movement =
        detail::SummarizeActiveDofChanges(changes, accepted_population);
    certificate.accepted_active_p99 = diagnostics.accepted_active_movement.percentile_list;
    diagnostics.operator_nominal_residual =
        detail::SummarizeActiveDofChanges(changes, nominal_population);
    certificate.operator_nominal_p99 = diagnostics.operator_nominal_residual.percentile_list;
    certificate.solver_qualified = true;

    EXPECT_EQ(
        diagnostics.accepted_active_movement.population_size_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1U);
    EXPECT_EQ(
        diagnostics.accepted_active_movement.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        1U);
    EXPECT_EQ(
        diagnostics.operator_nominal_residual.population_size_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        3U);
    EXPECT_EQ(
        diagnostics.operator_nominal_residual.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        3U);
    EXPECT_TRUE(detail::IsTransformedPercentileConverged(
        diagnostics.accepted_active_movement));
    EXPECT_FALSE(detail::IsTransformedPercentileConverged(
        diagnostics.operator_nominal_residual));
    EXPECT_FALSE(certificate.ProductionConverged());

    changes.at(0).fill(5.0e-5);
    changes.at(1).fill(5.0e-5);
    diagnostics.operator_nominal_residual =
        detail::SummarizeActiveDofChanges(changes, nominal_population);
    certificate.operator_nominal_p99 = diagnostics.operator_nominal_residual.percentile_list;
    EXPECT_TRUE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateQualifiesIndependentOffsetActivity)
{
    const std::vector<std::size_t> atom_index_list{ 0, 1 };
    detail::SuspiciousBlockActivity activity{
        detail::SuspiciousUpdateMask(2, 0),
        detail::SuspiciousUpdateMask{ 0, 1 },
        detail::SuspiciousUpdateMask(2, 0)
    };
    const std::vector<detail::ClusterKey> cluster_key_list{
        atom_index_list
    };
    const auto population{ detail::BuildActiveCoordinatePopulation(
        atom_index_list,
        activity) };
    std::vector<detail::TransformedChange> changes(
        2, detail::TransformedChange{});

    detail::ConvergenceCertificate certificate;
    certificate.accepted_active_p99 =
        detail::SummarizeActiveDofChanges(changes, population).percentile_list;
    const detail::SuspiciousBlockActivity nominal{ { 0, 0 }, { 0, 0 }, { 0, 0 } };
    const auto nominal_population{ detail::BuildActiveCoordinatePopulation(atom_index_list, nominal) };
    certificate.operator_nominal_p99 =
        detail::SummarizeActiveDofChanges(changes, nominal_population).percentile_list;
    const std::vector<std::optional<rg::RHBMEstimationStatus>>
        local_refit_status_by_atom(2, rg::RHBMEstimationStatus::SUCCESS);
    detail::ClusterHealthMap health_by_key;
    health_by_key.emplace(
        atom_index_list,
        detail::ClusterHealth{
            detail::JointOffsetSolveStatus::Converged
        });
    certificate.solver_qualified =
        detail::AreActiveCoordinatesSolverQualified(
            atom_index_list,
            cluster_key_list,
            activity,
            local_refit_status_by_atom,
            health_by_key);

    EXPECT_TRUE(certificate.solver_qualified);
    EXPECT_TRUE(certificate.StrictOperatorPassed());
    EXPECT_TRUE(certificate.ProductionConverged());
    changes.at(1).at(rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) = 2.0e-3;
    certificate.operator_nominal_p99 =
        detail::SummarizeActiveDofChanges(changes, nominal_population).percentile_list;
    EXPECT_FALSE(certificate.StrictOperatorPassed());
    EXPECT_FALSE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateAllFixedStillRequiresCompleteOperator)
{
    const auto make_summary = [](double value, std::size_t population)
    {
        detail::TransformedChangeSummary summary;
        summary.percentile_list.fill(value);
        summary.maximum_list.fill(value);
        summary.population_size_list.fill(population);
        return summary;
    };

    detail::ConvergenceCertificate certificate;
    certificate.accepted_active_p99 =
        make_summary(0.0, 0).percentile_list;
    certificate.operator_nominal_p99 =
        make_summary(5.0e-5, 4).percentile_list;
    certificate.solver_qualified = true;

    EXPECT_TRUE(detail::IsTransformedPercentileConverged(
        certificate.accepted_active_p99));
    EXPECT_TRUE(certificate.StrictOperatorPassed());
    EXPECT_TRUE(certificate.ProductionConverged());

    certificate.operator_complete = false;
    EXPECT_FALSE(certificate.StrictOperatorPassed());
    EXPECT_FALSE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest, NonFiniteChangeFailsPercentilePredicate)
{
    detail::TransformedChangeSummary summary;
    summary.percentile_list.fill(0.0);
    summary.maximum_list.fill(0.0);
    summary.population_size_list.fill(1);
    summary.percentile_list.at(
        rg::GaussianModel3D::LogPeakHeightCoordinateIndex()) =
        std::numeric_limits<double>::infinity();

    EXPECT_FALSE(detail::IsTransformedPercentileConverged(summary));
}

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageIterationsPersistsFinalModelAndPeelingWithoutGroupFitting)
{
    auto model{ BuildIndependentOffsetDefenseModel() };
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto options{ MakeSecondStageOptions() };
    const auto initial_analysis_view{ model->GetAnalysisView() };
    const auto group_key_list{ initial_analysis_view.CollectAtomGroupKeys() };
    std::vector<rg::GaussianModel3D> initial_group_mean_list;
    std::vector<rg::GaussianModel3D> initial_group_mdpde_list;
    std::vector<rg::GaussianModel3D> initial_group_prior_list;
    initial_group_mean_list.reserve(group_key_list.size());
    initial_group_mdpde_list.reserve(group_key_list.size());
    initial_group_prior_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        initial_group_mean_list.emplace_back(
            initial_analysis_view.GetAtomGroupMean(group_key));
        initial_group_mdpde_list.emplace_back(
            initial_analysis_view.GetAtomGroupMDPDE(group_key));
        initial_group_prior_list.emplace_back(
            initial_analysis_view.GetAtomGroupPrior(group_key));
    }

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, options);

    ExpectPeelingSamplingEntriesMatchFinalModels(*model);
    const auto analysis_view{ model->GetAnalysisView() };
    ASSERT_EQ(
        analysis_view.CollectAtomGroupKeys(),
        group_key_list);
    for (std::size_t group_index = 0;
        group_index < group_key_list.size();
        group_index++)
    {
        const auto group_key{ group_key_list.at(group_index) };
        const auto & atom_list{
            analysis_view.GetAtomObjectList(group_key)
        };
        EXPECT_DOUBLE_EQ(
            analysis_view.GetAtomAlphaG(group_key),
            0.0);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupMean(group_key),
            initial_group_mean_list.at(group_index),
            0.0);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupMDPDE(group_key),
            initial_group_mdpde_list.at(group_index),
            0.0);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupPrior(group_key),
            initial_group_prior_list.at(group_index),
            0.0);
        for (const auto * atom : atom_list)
        {
            EXPECT_FALSE(
                rg::AtomLocalPotentialView::For(*atom)
                    .GetGroupMemberResult().has_value());
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, SameChemicalKeyAtomsKeepIndependentOffsetsAndGroupKeyInvariantResults)
{
    auto original{ BuildIndependentOffsetDefenseModel() };
    auto relabeled{ BuildIndependentOffsetDefenseModel(1.0, true) };
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*original) };
    const auto run_logged = [](rg::ModelObject & model)
    {
        auto options{ MakeSecondStageOptions() };
        options.quiet_mode = false;
        const auto previous_level{ Logger::GetLogLevel() };
        Logger::SetLogLevel(LogLevel::Debug);
        testing::internal::CaptureStdout();
        model.EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
        detail::RunSecondStageIterations(model, options);
        const auto output{ testing::internal::GetCapturedStdout() };
        Logger::SetLogLevel(previous_level);
        std::vector<std::string> evidence;
        std::istringstream lines{ output };
        std::string line;
        while (std::getline(lines, line))
        {
            for (const std::string marker : {
                "Convergence safeguard audit:", "Second-stage audit terminal:",
                "Second-stage audit terminal atom:" })
            {
                const auto position{ line.find(marker) };
                if (position != std::string::npos) evidence.emplace_back(line.substr(position));
            }
        }
        EXPECT_FALSE(evidence.empty());
        return evidence;
    };
    const auto original_trace{ run_logged(*original) };
    EXPECT_EQ(original_trace, run_logged(*relabeled));
    const auto & atoms{ original->GetSelectedAtoms() };
    const auto & changed{ relabeled->GetSelectedAtoms() };
    ASSERT_EQ(atoms.size(), 2U);
    EXPECT_EQ(rg::data_internal::GetGroupKey(atoms.at(0)), rg::data_internal::GetGroupKey(atoms.at(1)));
    EXPECT_NE(rg::data_internal::GetGroupKey(changed.at(0)), rg::data_internal::GetGroupKey(changed.at(1)));
    EXPECT_GT(std::abs(GetEstimateModel(*atoms.at(0)).GetOffset() -
        GetEstimateModel(*atoms.at(1)).GetOffset()), 1.0e-5);
    for (std::size_t atom = 0; atom < atoms.size(); atom++)
    {
        ExpectGaussianModelsNear(GetEstimateModel(*atoms.at(atom)), GetEstimateModel(*changed.at(atom)), 1.0e-12);
        const auto peeled{ rg::AtomLocalPotentialView::For(*atoms.at(atom)).GetPeelingSamplingEntries(false) };
        const auto other{ rg::AtomLocalPotentialView::For(*changed.at(atom)).GetPeelingSamplingEntries(false) };
        ASSERT_EQ(peeled.size(), other.size());
        for (std::size_t row = 0; row < peeled.size(); row++)
            EXPECT_DOUBLE_EQ(peeled.at(row).response, other.at(row).response);
    }
    ExpectPeelingSamplingEntriesMatchFinalModels(*original);
    EXPECT_LT(CalculateSelectedAtomResponseMeanSquaredError(*original), initial_error);
    ExpectSelectedAtomEstimatesAreFinite(*original);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    auto base_model{ BuildNearCollinearDefenseModel() };
    auto scaled_model{ BuildNearCollinearDefenseModel(scale) };

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
            base.GetAmplitude() * scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 1.0e-5));
        EXPECT_NEAR(base.GetWidth(), scaled.GetWidth(), 1.0e-6);
        EXPECT_NEAR(
            base.GetOffset() * scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 5.0e-5));
    }
}
