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
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/JointFitting.hpp"
#include "core/detail/second_stage/ObjectiveEvaluation.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
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
using second_stage_test::BuildJointPolishFixture;
using second_stage_test::CalculateSelectedAtomResponseMeanSquaredError;
using second_stage_test::Distance;
using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::ExpectSelectedAtomEstimatesAreFinite;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeGaussianResult;
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
    detail::BoundaryJointCorrectionWorkspaceMap correction_workspace_by_key;
    detail::PerformanceCounters performance_counters{ true };
    auto options{ MakeSecondStageOptions() };
    const detail::SuspiciousBlockActivity all_active{
        std::vector<char>(fixture.context.atom_list.size(), 0),
        std::vector<char>(fixture.context.atom_list.size(), 0),
        std::vector<char>(fixture.context.atom_list.size(), 0)
    };
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
            performance_counters)
    };
    ASSERT_TRUE(polish_result.accepted);
    ASSERT_TRUE(polish_result.objective.has_value());
    ASSERT_TRUE(polish_result.diagnostic.objective_before.has_value());
    ASSERT_TRUE(polish_result.diagnostic.objective_after.has_value());
    EXPECT_LT(
        *polish_result.diagnostic.objective_after,
        *polish_result.diagnostic.objective_before);
    ASSERT_EQ(polish_result.diagnostic.component_count, 1U);
    EXPECT_GE(polish_result.diagnostic.round_count, 1U);
    EXPECT_LE(
        polish_result.diagnostic.round_count,
        options.second_stage_dependency_polish_max_iterations);
    EXPECT_EQ(polish_result.diagnostic.parameter_count, 6U);
    EXPECT_NE(
        polish_result.state.at(0).mdpde.GetModel().GetOffset(),
        polish_result.state.at(1).mdpde.GetModel().GetOffset());
    EXPECT_NE(
        polish_result.state.at(0).mdpde.GetModel().GetAmplitude(),
        fixture.state.at(0).mdpde.GetModel().GetAmplitude());
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
