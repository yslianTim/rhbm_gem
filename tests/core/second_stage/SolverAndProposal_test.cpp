#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/JointFitting.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
namespace rt = rhbm_gem::core;
namespace alg = rhbm_gem::algorithm;
using rhbm_gem::FittingStage;

using second_stage_test::BuildIndependentOffsetDefenseModel;
using second_stage_test::BuildJointPolishDefenseModel;
using second_stage_test::BuildJointPolishFixture;
using second_stage_test::CalculateSelectedAtomResponseMeanSquaredError;
using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::ExpectSelectedAtomEstimatesAreFinite;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeGaussianResult;
using second_stage_test::MakeSecondStageOptions;

static_assert(!std::is_default_constructible_v<detail::ClusterHealth>);

std::pair<detail::SecondStageContext,
    detail::SecondStageModelSnapshot>
BuildJointOffsetEstimationFixture(
    const std::vector<rg::GaussianModel3D> & model_list,
    const std::vector<double> & target_offset_list)
{
    if (model_list.size() != target_offset_list.size())
    {
        throw std::invalid_argument(
            "Joint offset fixture input sizes are inconsistent.");
    }

    detail::SecondStageContext context;
    context.atom_list.resize(model_list.size());
    for (std::size_t atom_index = 0;
        atom_index < model_list.size();
        atom_index++)
    {
        auto & atom_context{ context.atom_list.at(atom_index) };
        SamplingPoint point;
        point.distance = 0.35;
        atom_context.raw_sampling_entries.emplace_back(LocalPotentialSample{
            model_list.at(atom_index).SignalAtDistance(point.distance) +
                target_offset_list.at(atom_index) *
                    model_list.at(atom_index).OffsetBasisAtDistance(point.distance),
            point
        });
        atom_context.neighbor_atom_sample_offset_list = { 0, 0 };
    }
    return {
        std::move(context),
        detail::SecondStageModelSnapshot{
            model_list
        }
    };
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, JointFittingConditioningDetectsJointDependence)
{
    Eigen::SparseMatrix<double> design_matrix{ 3, 3 };
    const std::vector<Eigen::Triplet<double>> entries{
        { 0, 0, 1.0 }, { 0, 2, 1.0 },
        { 1, 1, 1.0 }, { 1, 2, 1.0 },
        { 2, 0, 1.0 }, { 2, 1, 1.0 }, { 2, 2, 2.0 }
    };
    design_matrix.setFromTriplets(entries.begin(), entries.end());

    const Eigen::MatrixXd dense{ design_matrix };
    for (Eigen::Index left = 0; left < dense.cols(); left++)
    {
        for (Eigen::Index right = left + 1; right < dense.cols(); right++)
        {
            const auto overlap{
                std::abs(dense.col(left).dot(dense.col(right))) /
                (dense.col(left).norm() * dense.col(right).norm())
            };
            EXPECT_LT(overlap, 0.98);
        }
    }

    const auto diagnostics{
        detail::EvaluateJointFittingConditioning(
            design_matrix,
            1.0e-8)
    };
    EXPECT_TRUE(diagnostics.guard_required);
    EXPECT_LE(diagnostics.pivot_ratio, 1.0e-8);
}

TEST(EstimatorSecondStageDefenseTest, JointFittingConditioningKeepsIndependentColumns)
{
    Eigen::SparseMatrix<double> design_matrix{ 3, 3 };
    design_matrix.setIdentity();

    const auto diagnostics{
        detail::EvaluateJointFittingConditioning(
            design_matrix,
            1.0e-8)
    };
    EXPECT_FALSE(diagnostics.guard_required);
    EXPECT_NEAR(diagnostics.pivot_ratio, 1.0, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, JointFittingConditioningFailsClosedAndIncludesThreshold)
{
    Eigen::SparseMatrix<double> zero_column{ 2, 2 };
    zero_column.insert(0, 0) = 1.0;
    auto result{ detail::EvaluateJointFittingConditioning(zero_column, 1.0e-8) };
    EXPECT_TRUE(result.guard_required);
    EXPECT_DOUBLE_EQ(result.pivot_ratio, 0.0);
    zero_column.insert(1, 1) = std::numeric_limits<double>::infinity();
    result = detail::EvaluateJointFittingConditioning(zero_column, 1.0e-8);
    EXPECT_TRUE(result.guard_required);
    EXPECT_DOUBLE_EQ(result.pivot_ratio, 0.0);

    Eigen::SparseMatrix<double> identity{ 2, 2 };
    identity.setIdentity();
    result = detail::EvaluateJointFittingConditioning(identity, 1.0);
    EXPECT_TRUE(result.guard_required);
    EXPECT_DOUBLE_EQ(result.pivot_ratio, 1.0);
    result = detail::EvaluateJointFittingConditioning(
        identity, std::nextafter(1.0, 0.0));
    EXPECT_FALSE(result.guard_required);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorPreservesIndividualRidgeAnchors)
{
    const std::vector<rg::GaussianModel3D> models{
        { 6.0, 0.55, 1.0 }, { 6.0, 0.55, 6.0 } };
    auto fixture{ BuildJointOffsetEstimationFixture(models, { 1.0, 6.0 }) };
    // Identical columns constrain only the sum; ridge must retain each atom's own seed.
    for (std::size_t atom = 0; atom < 2; atom++)
    {
        fixture.first.atom_list.at(atom).raw_sampling_entries.front().response +=
            models.at(1 - atom).ResponseAtDistance(0.35);
        fixture.first.atom_list.at(atom).neighbor_atom_sample_list = { { 1 - atom, 0.35 } };
        fixture.first.atom_list.at(atom).neighbor_atom_sample_offset_list = { 0, 1 };
    }
    alg::WeightedRidgeSolver solver;
    const auto result{ detail::EstimateJointOffsets(
        fixture.first, { 0, 1 }, fixture.second, { 1.0, 1.0 }, solver) };
    ASSERT_EQ(result.status, detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(result.offset.size(), 2);
    EXPECT_NEAR(result.offset(0), 1.0, 1.0e-10);
    EXPECT_NEAR(result.offset(1), 6.0, 1.0e-10);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorMapsPermutedAtomColumns)
{
    const std::vector<double> target_offsets{ 1.15, 3.9, 3.2, -0.2 };
    auto fixture{ BuildJointOffsetEstimationFixture(
        { { 6.0, 0.55, 1.0 }, { 7.0, 0.60, 4.0 },
            { 8.0, 0.65, 3.0 }, { 4.0, 0.50, -0.2 } },
        target_offsets) };
    // Each active target sees two active neighbors and one fixed, outside-cluster atom.
    for (std::size_t atom_index = 0; atom_index < 3; atom_index++)
    {
        auto & atom{ fixture.first.atom_list.at(atom_index) };
        for (std::size_t neighbor_index = 0; neighbor_index < 4; neighbor_index++)
        {
            if (neighbor_index == atom_index) continue;
            const auto distance{ 0.7 + 0.2 * static_cast<double>(neighbor_index) +
                0.1 * static_cast<double>(atom_index) };
            atom.neighbor_atom_sample_list.push_back({ neighbor_index, distance });
            atom.raw_sampling_entries.front().response += fixture.second.node.at(neighbor_index)
                .WithOffset(target_offsets.at(neighbor_index)).ResponseAtDistance(distance);
        }
        atom.neighbor_atom_sample_offset_list = { 0, atom.neighbor_atom_sample_list.size() };
    }
    alg::WeightedRidgeSolver solver;
    const auto original{ detail::EstimateJointOffsets(
        fixture.first, { 0, 1, 2 }, fixture.second, { 1.0, 2.0, 4.0, 9.0 }, solver) };
    const auto reordered{ detail::EstimateJointOffsets(
        fixture.first, { 1, 2, 0 }, fixture.second, { 1.0, 2.0, 4.0, 9.0 }, solver) };
    ASSERT_EQ(original.status, detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(reordered.status, detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(original.offset.size(), 3);
    ASSERT_EQ(reordered.offset.size(), 3);
    EXPECT_NEAR(reordered.offset(0), original.offset(1), 1.0e-12);
    EXPECT_NEAR(reordered.offset(1), original.offset(2), 1.0e-12);
    EXPECT_NEAR(reordered.offset(2), original.offset(0), 1.0e-12);
    for (auto & atom : fixture.first.atom_list)
    {
        std::ranges::reverse(atom.neighbor_atom_sample_list);
    }
    const auto reversed_neighbors{ detail::EstimateJointOffsets(
        fixture.first, { 0, 1, 2 }, fixture.second, { 1.0, 2.0, 4.0, 9.0 }, solver) };
    ASSERT_EQ(reversed_neighbors.status, original.status);
    ASSERT_EQ(reversed_neighbors.offset.size(), original.offset.size());
    for (Eigen::Index column = 0; column < original.offset.size(); column++)
    {
        EXPECT_NEAR(reversed_neighbors.offset(column), original.offset(column), 1.0e-12);
        const auto atom_index{ static_cast<std::size_t>(column) };
        EXPECT_LT(std::abs(original.offset(column) - target_offsets.at(atom_index)),
            std::abs(fixture.second.node.at(atom_index).GetOffset() - target_offsets.at(atom_index)));
    }
}

TEST(EstimatorSecondStageDefenseTest, LocalRefitHealthTracksSolverQualification)
{
    EXPECT_TRUE(rg::GetSolveQualification(rg::RHBMEstimationStatus::SUCCESS, {}) != rg::RHBMSolveQualification::Unqualified);
    EXPECT_FALSE(rg::GetSolveQualification(rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED, {}) != rg::RHBMSolveQualification::Unqualified);
    for (const auto status : {
        rg::RHBMEstimationStatus::NUMERICAL_FALLBACK,
        rg::RHBMEstimationStatus::INSUFFICIENT_DATA,
        rg::RHBMEstimationStatus::SINGLE_MEMBER })
    {
        EXPECT_FALSE(rg::GetSolveQualification(status, {}) != rg::RHBMSolveQualification::Unqualified);
    }
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetHealthSeparatesHardFailureFromSolverQualification)
{
    using Status = detail::JointOffsetSolveStatus;

    EXPECT_TRUE(detail::ClusterHealth{ Status::Converged }.IsSolverQualified());
    EXPECT_FALSE(detail::IsJointOffsetSolveHardFailure(Status::Converged));

    for (const auto status : {
        Status::IrlsObjectiveDeteriorated,
        Status::IrlsMaximumIterationsReached })
    {
        EXPECT_FALSE(detail::ClusterHealth{ status }.IsSolverQualified());
        EXPECT_FALSE(detail::IsJointOffsetSolveHardFailure(status));
    }

    for (const auto status : {
        Status::SystemBuildFailed,
        Status::EmptySystem,
        Status::InitialSolveFailed,
        Status::IrlsSolveFailed })
    {
        EXPECT_FALSE(detail::ClusterHealth{ status }.IsSolverQualified());
        EXPECT_TRUE(detail::IsJointOffsetSolveHardFailure(status));
    }
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorFitsIndependentAtomOffsets)
{
    auto fixture{
        BuildJointOffsetEstimationFixture(
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.0 },
                rg::GaussianModel3D{ 7.0, 0.60, 4.0 }
            },
            { 1.0, 3.0 })
    };
    alg::WeightedRidgeSolver solver;
    const auto result{
        detail::EstimateJointOffsets(
            fixture.first,
            { 0, 1 },
            fixture.second,
            { 1.0, 1.0 },
            solver)
    };

    EXPECT_EQ(
        result.status,
        detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(result.offset.size(), 2);
    EXPECT_NEAR(result.offset(0), 1.0, 0.01);
    EXPECT_NEAR(result.offset(1), 3.0, 0.01);
    EXPECT_GT(result.offset(1) - result.offset(0), 1.9);
    for (const auto multiplier : { 0.0, 0.25, std::numeric_limits<double>::quiet_NaN() })
    {
        const auto clamped{ detail::EstimateJointOffsets(
            fixture.first, { 0, 1 }, fixture.second,
            { multiplier, multiplier }, solver) };
        ASSERT_EQ(clamped.status, result.status);
        ASSERT_EQ(clamped.offset.size(), result.offset.size());
        EXPECT_DOUBLE_EQ(clamped.offset(0), result.offset(0));
        EXPECT_DOUBLE_EQ(clamped.offset(1), result.offset(1));
    }
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorKeepsFrozenBackgroundInRhs)
{
    auto fixture{
        BuildJointOffsetEstimationFixture(
            {
                rg::GaussianModel3D{ 6.0, 0.55, 1.0 },
                rg::GaussianModel3D{ 7.0, 0.60, 3.0 }
            },
            { 1.0, 3.0 })
    };
    alg::WeightedRidgeSolver solver;
    const auto result{
        detail::EstimateJointOffsets(
            fixture.first,
            { 0, 1 },
            fixture.second,
            { 1.0, 1.0 },
            solver)
    };

    EXPECT_EQ(
        result.status,
        detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(result.offset.size(), 2);
    EXPECT_NEAR(result.offset(0), 1.0, 1.0e-5);
    EXPECT_NEAR(result.offset(1), 3.0, 1.0e-5);
    EXPECT_GT(std::abs(result.offset(0) - result.offset(1)), 1.0);


    detail::FitState state;
    for (const auto & gaussian : fixture.second.node) state.emplace_back(MakeGaussianResult(gaussian));
    for (auto & atom : fixture.first.atom_list)
        atom.unselected_distance_list_by_sample.assign(atom.raw_sampling_entries.size(), { 0.4 });
    fixture.first.frozen_background = detail::BuildFrozenBackground(fixture.first, state);
    ASSERT_TRUE(fixture.first.frozen_background);
    for (std::size_t node = 0; node < state.size(); node++)
        for (std::size_t row = 0; row < fixture.first.atom_list.at(node).raw_sampling_entries.size(); row++)
            fixture.first.atom_list.at(node).raw_sampling_entries.at(row).response +=
                fixture.first.frozen_background->response_by_atom.at(node).at(row);
    fixture.second = detail::BuildSecondStageModelSnapshot(fixture.first, state);
    rg::algorithm::WeightedRidgeSolver background_solver;
    const auto with_background{ detail::EstimateJointOffsets(fixture.first, { 0, 1 },
        fixture.second, { 1.0, 1.0 }, background_solver) };
    ASSERT_EQ(with_background.status, detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(with_background.offset.size(), 2);
    EXPECT_NEAR(with_background.offset(0), result.offset(0), 1.0e-12);
    EXPECT_NEAR(with_background.offset(1), result.offset(1), 1.0e-12);

}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorReportsBuildAndEmptyFailures)
{
    auto empty_fixture{
        BuildJointOffsetEstimationFixture(
            { rg::GaussianModel3D{ 6.0, 0.55, 2.0 } },
            { 2.0 })
    };
    empty_fixture.first.atom_list.at(0).raw_sampling_entries.clear();
    empty_fixture.first.atom_list.at(0).neighbor_atom_sample_offset_list.clear();
    alg::WeightedRidgeSolver empty_solver;
    const auto empty_result{
        detail::EstimateJointOffsets(
            empty_fixture.first,
            { 0 },
            empty_fixture.second,
            { 1.0 },
            empty_solver)
    };
    EXPECT_EQ(
        empty_result.status,
        detail::JointOffsetSolveStatus::EmptySystem);
    ASSERT_EQ(empty_result.offset.size(), 1);
    EXPECT_DOUBLE_EQ(empty_result.offset(0), 2.0);

    auto invalid_fixture{
        BuildJointOffsetEstimationFixture(
            { rg::GaussianModel3D{ 6.0, 0.55, 2.0 } },
            { 2.0 })
    };
    invalid_fixture.first.atom_list.at(0).raw_sampling_entries.at(0).response =
        std::numeric_limits<double>::infinity();
    alg::WeightedRidgeSolver invalid_solver;
    const auto invalid_result{
        detail::EstimateJointOffsets(
            invalid_fixture.first,
            { 0 },
            invalid_fixture.second,
            { 1.0 },
            invalid_solver)
    };
    EXPECT_EQ(
        invalid_result.status,
        detail::JointOffsetSolveStatus::SystemBuildFailed);
    ASSERT_EQ(invalid_result.offset.size(), 1);
    EXPECT_DOUBLE_EQ(invalid_result.offset(0), 2.0);

    auto negligible_basis_fixture{ BuildJointOffsetEstimationFixture(
        { { 6.0, 0.55, 2.0 } }, { 2.0 }) };
    negligible_basis_fixture.first.atom_list.at(0).raw_sampling_entries.front().point.distance = 1.0e20;
    const auto negligible_basis_result{ detail::EstimateJointOffsets(
        negligible_basis_fixture.first, { 0 }, negligible_basis_fixture.second,
        { 1.0 }, empty_solver) };
    EXPECT_EQ(negligible_basis_result.status, detail::JointOffsetSolveStatus::EmptySystem);
    ASSERT_EQ(negligible_basis_result.offset.size(), 1);
    EXPECT_DOUBLE_EQ(negligible_basis_result.offset(0), 2.0);

    for (const bool active_neighbor : { false, true })
    {
        auto non_finite_neighbor_fixture{ BuildJointOffsetEstimationFixture(
            { { 6.0, 0.55, 2.0 }, { 7.0, 0.60, 3.0 } }, { 2.0, 3.0 }) };
        auto & target{ non_finite_neighbor_fixture.first.atom_list.at(0) };
        target.neighbor_atom_sample_list = { { 1, std::numeric_limits<double>::quiet_NaN() } };
        target.neighbor_atom_sample_offset_list = { 0, 1 };
        const auto non_finite_neighbor_result{ detail::EstimateJointOffsets(
            non_finite_neighbor_fixture.first,
            active_neighbor ? detail::ClusterKey{ 0, 1 } : detail::ClusterKey{ 0 },
            non_finite_neighbor_fixture.second, { 1.0, 1.0 }, invalid_solver) };
        EXPECT_EQ(non_finite_neighbor_result.status, detail::JointOffsetSolveStatus::SystemBuildFailed);
        ASSERT_EQ(non_finite_neighbor_result.offset.size(), active_neighbor ? 2 : 1);
        EXPECT_DOUBLE_EQ(non_finite_neighbor_result.offset(0), 2.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, JointPolishJacobianMatchesFiniteDifference)
{
    constexpr double step{ 1.0e-6 };
    const std::array<rg::GaussianModel3D, 2> model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.70, -0.15 }
    };
    const std::array<double, 4> distance_list{
        0.0,
        5.0e-6,
        1.0e-5,
        0.35
    };
    for (const auto & model : model_list)
    {
        const auto transformed{ model.ToTransformedCoordinates() };
        ASSERT_TRUE(transformed.has_value());
        const auto invariants{
            detail::BuildTransformedModelInvariants(model)
        };
        ASSERT_TRUE(invariants.has_value());
        for (const auto distance : distance_list)
        {
            const auto evaluation{
                detail::EvaluateTransformedJacobian(
                    *invariants,
                    distance)
            };
            ASSERT_TRUE(evaluation.has_value());
            EXPECT_TRUE(evaluation->allFinite());

            for (Eigen::Index parameter_index = 0;
                parameter_index < transformed->size();
                parameter_index++)
            {
                auto lower{ *transformed };
                auto upper{ *transformed };
                lower(parameter_index) -= step;
                upper(parameter_index) += step;
                const auto lower_model{
                    rg::GaussianModel3D::FromTransformedCoordinates(lower)
                };
                const auto upper_model{
                    rg::GaussianModel3D::FromTransformedCoordinates(upper)
                };
                ASSERT_TRUE(lower_model.has_value());
                ASSERT_TRUE(upper_model.has_value());
                const auto finite_difference{
                    (upper_model->ResponseAtDistance(distance) -
                        lower_model->ResponseAtDistance(distance)) /
                    (2.0 * step)
                };
                const auto analytic{
                    (*evaluation)(parameter_index)
                };
                const auto tolerance{
                    1.0e-6 * std::max(std::abs(finite_difference), 1.0)
                };
                EXPECT_NEAR(analytic, finite_difference, tolerance);
            }
        }
    }
}

TEST(EstimatorSecondStageDefenseTest,
    IndividualRefitUsesOwnModelForTargetOffsetResponse)
{
    const rg::GaussianModel3D offset_model{ 3.0, 0.70, 0.2 };
    const rg::GaussianModel3D truth_shape{ 6.0, 0.55, 0.0 };
    LocalPotentialSampleList sample_list;
    for (const auto distance : std::array<double, 5>{
        0.0, 0.15, 0.30, 0.45, 0.60 })
    {
        sample_list.emplace_back(LocalPotentialSample{
            truth_shape.SignalAtDistance(distance) +
                offset_model.GetOffset() *
                    offset_model.OffsetBasisAtDistance(distance),
            SamplingPoint{ distance }
        });
    }

    const auto result{
        rt::EstimateLocalGaussian(
            sample_list,
            0.0,
            MakeSecondStageOptions(),
            offset_model)
    };
    EXPECT_NEAR(
        result.mdpde.GetModel().GetAmplitude(),
        truth_shape.GetAmplitude(),
        1.0e-4);
    EXPECT_NEAR(
        result.mdpde.GetModel().GetWidth(),
        truth_shape.GetWidth(),
        1.0e-6);
    EXPECT_DOUBLE_EQ(
        result.mdpde.GetModel().GetOffset(),
        offset_model.GetOffset());
}

TEST(EstimatorSecondStageDefenseTest,
    JointPolishParameterizationKeepsGlobalBackgroundFrozen)
{
    const std::vector<rg::GaussianModel3D> models{
        { 6.0, 0.55, 1.0 }, { 7.0, 0.60, 4.0 }, { 8.0, 0.65, 3.0 } };
    detail::SecondStageContext context;
    context.atom_list.resize(models.size());
    detail::FitState state;
    for (std::size_t i = 0; i < models.size(); i++)
    {
        state.emplace_back(MakeGaussianResult(models.at(i)));
        context.atom_list.at(i).raw_sampling_entries.resize(1);
        context.atom_list.at(i).unselected_distance_list_by_sample = { { 0.35 } };
    }
    context.frozen_background = detail::BuildFrozenBackground(context, state);
    ASSERT_TRUE(context.frozen_background);
    ExpectGaussianModelsNear(context.frozen_background->model_by_atom.front(), { 7.0, 0.60, 3.0 }, 1.0e-12);
    const auto parameterization{ detail::JointPolishParameterization::Build( models) };
    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->AtomCount(), 3U);
    EXPECT_EQ(parameterization->ParameterCount(), 9);
    EXPECT_NE(parameterization->OffsetColumn(0), parameterization->OffsetColumn(2));
    const auto frozen{ context.frozen_background };
    auto direction{ Eigen::VectorXd::Zero(9).eval() };
    direction(parameterization->OffsetColumn(0)) = 6.0;
    const auto candidate{ parameterization->DecodeModels(direction, 1.0) };
    ASSERT_TRUE(candidate.has_value());
    EXPECT_EQ(context.frozen_background, frozen);
    ExpectGaussianModelsNear(frozen->model_by_atom.front(), { 7.0, 0.60, 3.0 }, 1.0e-12);

    // A fixed selected member still contributes equally to the even median.
    auto even_context{ context };
    even_context.atom_list.resize(2);
    const detail::FitState even_state{ state.at(0), state.at(1) };
    const auto even{ detail::BuildFrozenBackground(even_context, even_state) };
    ASSERT_TRUE(even);
    ExpectGaussianModelsNear(even->model_by_atom.front(), { 6.5, 0.575, 2.5 }, 1.0e-12);
    const auto fixed{ detail::JointPolishParameterization::BuildActiveSet( models, { 1, 0, 0 }, { 1, 0, 0 }) };
    ASSERT_TRUE(fixed.has_value());
    EXPECT_EQ(fixed->ParameterCount(), 3);
    constexpr double step{ 1.0e-6 };
    auto perturbation{ Eigen::VectorXd::Zero(3).eval() };
    perturbation(fixed->OffsetColumn(0)) = step;
    const auto upper{ fixed->DecodeModels(perturbation, 1.0) };
    const auto lower{ fixed->DecodeModels(-perturbation, 1.0) };
    ASSERT_TRUE(upper.has_value());
    ASSERT_TRUE(lower.has_value());
    const auto response{ detail::EvaluatePhysicalOffsetResponse(models.front(), 0.35) };
    ASSERT_TRUE(response.has_value());
    const double background{ frozen->response_by_atom.front().front() };
    EXPECT_NEAR(((upper->front().ResponseAtDistance(0.35) + background) -
        (lower->front().ResponseAtDistance(0.35) + background)) / (2.0 * step),
        response->offset_jacobian, 1.0e-8);

}

TEST(EstimatorSecondStageDefenseTest, JointPolishSeedsAndDecodesIndividualOffsets)
{
    const std::vector<rg::GaussianModel3D> models{
        { 6.0, 0.55, 1.0 }, { 7.0, 0.60, 4.0 }, { 8.0, 0.65, 3.0 },
        { 9.0, 0.70, 2.0 }, { 10.0, 0.75, 6.0 } };
    const auto parameterization{ detail::JointPolishParameterization::Build(models) };
    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->ParameterCount(), 15);
    const auto seed{ parameterization->DecodeSeedModels() };
    ASSERT_TRUE(seed.has_value());
    auto direction{ Eigen::VectorXd::Zero(parameterization->ParameterCount()).eval() };
    direction(parameterization->OffsetColumn(0)) = -1.0;
    direction(parameterization->OffsetColumn(1)) = 2.0;
    const auto candidate{ parameterization->DecodeModels(direction, 0.5) };
    const auto zero{ parameterization->DecodeModels(direction, 0.0) };
    ASSERT_TRUE(candidate.has_value());
    ASSERT_TRUE(zero.has_value());
    for (std::size_t atom = 0; atom < models.size(); atom++)
    {
        ExpectGaussianModelsNear(seed->at(atom), models.at(atom), 1.0e-12);
        ExpectGaussianModelsNear(zero->at(atom), models.at(atom), 1.0e-12);
        const auto delta{ atom == 0 ? -0.5 : atom == 1 ? 1.0 : 0.0 };
        EXPECT_DOUBLE_EQ(candidate->at(atom).GetOffset(), models.at(atom).GetOffset() + delta);
        for (std::size_t other = atom + 1; other < models.size(); other++)
            EXPECT_NE(parameterization->OffsetColumn(atom), parameterization->OffsetColumn(other));
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    ActiveSetJointPolishKeepsInactiveCoordinatesAndIndependentOffsets)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.20 }
    };
    const auto parameterization{
        detail::JointPolishParameterization::BuildActiveSet(
            base_model_list,
            { 1, 0, 1 },
            { 1, 1, 1 })
    };
    ASSERT_TRUE(parameterization.has_value());
    EXPECT_TRUE(parameterization->HasShapeColumn(0));
    EXPECT_FALSE(parameterization->HasShapeColumn(1));
    EXPECT_TRUE(parameterization->HasShapeColumn(2));
    EXPECT_EQ(parameterization->ParameterCount(), 7);

    Eigen::VectorXd direction{
        Eigen::VectorXd::Zero(parameterization->ParameterCount())
    };
    direction(parameterization->ShapeColumn(0, 0)) = 0.2;
    direction(parameterization->ShapeColumn(2, 1)) = -0.1;
    direction(parameterization->OffsetColumn(0)) = 0.4;
    direction(parameterization->OffsetColumn(2)) = -0.2;
    const auto candidate_model_list{
        parameterization->DecodeModels(direction, 1.0)
    };
    ASSERT_TRUE(candidate_model_list.has_value());
    EXPECT_NE(
        candidate_model_list->at(0).GetAmplitude(),
        base_model_list.at(0).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        candidate_model_list->at(1).GetAmplitude(),
        base_model_list.at(1).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        candidate_model_list->at(1).GetWidth(),
        base_model_list.at(1).GetWidth());
    EXPECT_NE(
        candidate_model_list->at(2).GetWidth(),
        base_model_list.at(2).GetWidth());
    EXPECT_NE(
        candidate_model_list->at(0).GetOffset(),
        candidate_model_list->at(1).GetOffset());

    const auto fixed_offset_parameterization{
        detail::JointPolishParameterization::BuildActiveSet(
            base_model_list,
            { 1, 0, 1 },
            { 1, 1, 0 })
    };
    ASSERT_TRUE(fixed_offset_parameterization.has_value());
    EXPECT_TRUE(fixed_offset_parameterization->HasOffsetColumn(0));
    EXPECT_FALSE(fixed_offset_parameterization->HasOffsetColumn(2));
    EXPECT_EQ(fixed_offset_parameterization->ParameterCount(), 6);
    const auto fixed_offset_models{
        fixed_offset_parameterization->DecodeModels(
            Eigen::VectorXd::Zero(fixed_offset_parameterization->ParameterCount()),
            1.0)
    };
    ASSERT_TRUE(fixed_offset_models.has_value());
    EXPECT_DOUBLE_EQ(
        fixed_offset_models->at(2).GetOffset(),
        base_model_list.at(2).GetOffset());

    auto fixture{ BuildJointPolishFixture(
        { { 6.0, 0.5, 0.1 }, { 7.0, 0.6, 0.3 } },
        { { 6.4, 0.55, 0.2 }, { 7.2, 0.62, 0.4 } }) };
    for (auto & atom : fixture.context.atom_list)
        atom.unselected_distance_list_by_sample.assign(atom.raw_sampling_entries.size(), { 0.35 });
    fixture.context.frozen_background = detail::BuildFrozenBackground(fixture.context, fixture.state);
    ASSERT_TRUE(fixture.context.frozen_background);
    const auto frozen{ fixture.context.frozen_background };
    for (std::size_t node = 0; node < fixture.state.size(); node++)
        for (std::size_t row = 0; row < fixture.context.atom_list.at(node).raw_sampling_entries.size(); row++)
            fixture.context.atom_list.at(node).raw_sampling_entries.at(row).response += frozen->response_by_atom.at(node).at(row);
    const detail::FitStatePatch empty_patch;
    const detail::FitStateView endpoint{ fixture.state, empty_patch };
    for (const bool freeze_second_shape : { false, true })
    {
        rg::algorithm::WeightedRidgeSolver solver;
        const auto correction{ detail::BuildBoundaryJointCorrection(fixture.context, endpoint,
            freeze_second_shape ? detail::ClusterKey{ 0 } : detail::ClusterKey{ 0, 1 },
            { 0, 1 }, fixture.sample_ref_list, { 1.0, 1.0 }, { { { 0, 1 }, 100.0 } }, solver) };
        ASSERT_EQ(correction.status, detail::BoundaryJointCorrectionStatus::CandidateReady);
        ASSERT_TRUE(correction.patch.has_value());
        EXPECT_EQ(correction.parameter_count, freeze_second_shape ? 4U : 6U);
        EXPECT_EQ(correction.patch->atom_index_list, (detail::ClusterKey{ 0, 1 }));
        EXPECT_EQ(fixture.context.frozen_background, frozen);
        ExpectGaussianModelsNear(frozen->model_by_atom.at(0), { 6.5, 0.55, 0.2 }, 1.0e-12);
        ExpectGaussianModelsNear(frozen->model_by_atom.at(1), { 6.5, 0.55, 0.2 }, 1.0e-12);
    }

}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryJointCorrectionUsesIndependentActiveCoordinatesAndPerClusterTrust)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.20 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.8, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.20 },
        rg::GaussianModel3D{ 6.4, 0.56, -0.15 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    const detail::FitStatePatch endpoint_patch;
    const detail::FitStateView endpoint_state{
        fixture.state,
        endpoint_patch
    };
    alg::WeightedRidgeSolver solver;
    const auto result{
        detail::BuildBoundaryJointCorrection(
            fixture.context,
            endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 4.0 },
                { { 2 }, 0.5 }
            },
            solver)
    };
    ASSERT_EQ(
        result.status,
        detail::BoundaryJointCorrectionStatus::CandidateReady);
    ASSERT_TRUE(result.patch.has_value());
    EXPECT_EQ(result.parameter_count, 7U);
    EXPECT_LE(result.maximum_normalized_trust_step, 1.0 + 1.0e-12);
    ASSERT_EQ(result.patch->mdpde_list.size(), 3U);
    const auto & offset_only_candidate{
        result.patch->mdpde_list.at(1).GetModel()
    };
    EXPECT_DOUBLE_EQ(
        offset_only_candidate.GetAmplitude(),
        base_model_list.at(1).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        offset_only_candidate.GetWidth(),
        base_model_list.at(1).GetWidth());
    EXPECT_NE(
        result.patch->mdpde_list.at(0).GetModel().GetOffset(),
        offset_only_candidate.GetOffset());
    EXPECT_DOUBLE_EQ(
        result.patch->mdpde_list.at(1)
            .GetStandardDeviationModel().GetAmplitude(),
        fixture.state.at(1).mdpde
            .GetStandardDeviationModel().GetAmplitude());

    // Atom 1 has only an offset column: the ridge floor must not affect shape columns.
    for (const auto multiplier : { 0.25, std::numeric_limits<double>::quiet_NaN() })
    {
        const auto clamped{ detail::BuildBoundaryJointCorrection(
            fixture.context, endpoint_state, { 0, 2 }, { 0, 1, 2 },
            fixture.sample_ref_list, { 1.0, multiplier, 1.0 },
            { { { 0, 1 }, 4.0 }, { { 2 }, 0.5 } }, solver) };
        ASSERT_EQ(clamped.status, result.status);
        ASSERT_TRUE(clamped.patch.has_value());
        ASSERT_EQ(clamped.patch->mdpde_list.size(), result.patch->mdpde_list.size());
        EXPECT_EQ(clamped.parameter_count, result.parameter_count);
        for (std::size_t atom_index = 0; atom_index < result.patch->mdpde_list.size(); atom_index++)
        {
            ExpectGaussianModelsNear(clamped.patch->mdpde_list.at(atom_index).GetModel(),
                result.patch->mdpde_list.at(atom_index).GetModel(), 0.0);
        }
    }

    alg::WeightedRidgeSolver invalid_solver;
    EXPECT_EQ(
        detail::BuildBoundaryJointCorrection(
            fixture.context,
            endpoint_state,
            {},
            {},
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            { { { 0, 1, 2 }, 4.0 } },
            invalid_solver).status,
        detail::BoundaryJointCorrectionStatus::InvalidInput);

    alg::WeightedRidgeSolver small_step_solver;
    EXPECT_EQ(
        detail::BuildBoundaryJointCorrection(
            fixture.context,
            endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 1.0e-8 },
                { { 2 }, 1.0e-8 }
            },
            small_step_solver).status,
        detail::BoundaryJointCorrectionStatus::NoMaterialChange);

    auto outside_patch{ detail::FitStatePatch::FromState(fixture.state, { 0 }) };
    outside_patch.mdpde_list.front() = MakeGaussianResult({ 10.0, 0.55, 0.10 }).mdpde;
    const detail::FitStateView outside_state{ fixture.state, outside_patch };
    alg::WeightedRidgeSolver unavailable_solver;
    EXPECT_EQ(
        detail::BuildBoundaryJointCorrection(
            fixture.context,
            outside_state,
            { 0, 2 },
            { 0, 1, 2 },
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 1.0e-8 },
                { { 2 }, 1.0e-8 }
            },
            unavailable_solver).status,
        detail::BoundaryJointCorrectionStatus::TrustRegionUnavailable);

    auto non_finite_fixture{ fixture };
    non_finite_fixture.context.atom_list.at(0).raw_sampling_entries.at(0).response =
        std::numeric_limits<double>::infinity();
    alg::WeightedRidgeSolver non_finite_solver;
    EXPECT_EQ(
        detail::BuildBoundaryJointCorrection(
            non_finite_fixture.context,
            endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            non_finite_fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 4.0 },
                { { 2 }, 0.5 }
            },
            non_finite_solver).status,
        detail::BoundaryJointCorrectionStatus::SystemBuildFailed);

    const std::vector<rg::GaussianModel3D> stationary_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.0 },
        rg::GaussianModel3D{ 4.5, 0.70, 0.0 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.0 }
    };
    auto stationary_fixture{
        BuildJointPolishFixture(
            stationary_model_list,
            stationary_model_list)
    };
    const detail::FitStateView stationary_endpoint_state{
        stationary_fixture.state,
        endpoint_patch
    };
    alg::WeightedRidgeSolver stationary_solver;
    EXPECT_EQ(
        detail::BuildBoundaryJointCorrection(
            stationary_fixture.context,
            stationary_endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            stationary_fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 4.0 },
                { { 2 }, 0.5 }
            },
            stationary_solver).status,
        detail::BoundaryJointCorrectionStatus::NoMaterialChange);
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalKeepsIndependentOffsets)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, -0.20 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    const detail::FitStatePatch base_patch;
    const detail::FitStateView base_state_view{ fixture.state, base_patch };
    alg::WeightedRidgeSolver proposal_solver;
    const auto proposal{
        detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            detail::ClusterKey{ 0, 1 },
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            proposal_solver,
            4.0)
    };
    ASSERT_TRUE(proposal.has_value());
    EXPECT_EQ(proposal->patch.atom_index_list,
        (detail::ClusterKey{ 0, 1 }));
    ASSERT_EQ(proposal->patch.mdpde_list.size(), 2U);
    EXPECT_NE(
        proposal->patch.mdpde_list.at(0).GetModel().GetOffset(),
        proposal->patch.mdpde_list.at(1).GetModel().GetOffset());
    bool shape_changed{ false };
    bool offset_changed{ false };
    for (std::size_t atom_index = 0; atom_index < base_model_list.size(); atom_index++)
    {
        const auto & candidate{
            proposal->patch.mdpde_list.at(atom_index).GetModel()
        };
        const auto & base{ base_model_list.at(atom_index) };
        shape_changed = shape_changed ||
            std::abs(candidate.GetAmplitude() - base.GetAmplitude()) > 1.0e-8 ||
            std::abs(candidate.GetWidth() - base.GetWidth()) > 1.0e-8;
        offset_changed = offset_changed ||
            std::abs(candidate.GetOffset() - base.GetOffset()) > 1.0e-8;
    }
    EXPECT_TRUE(shape_changed);
    EXPECT_TRUE(offset_changed);
    EXPECT_LE(proposal->step_norm, 4.0 + 1.0e-12);
    EXPECT_DOUBLE_EQ(
        proposal->patch.mdpde_list.at(0)
            .GetStandardDeviationModel().GetAmplitude(),
        0.0);
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalRejectsUnchangedAndRespectsSmallTrustRegion)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, -0.20 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    const detail::FitStatePatch base_patch;
    const detail::FitStateView base_state_view{ fixture.state, base_patch };
    const auto key{ detail::ClusterKey{ 0, 1 } };

    const std::vector<rg::GaussianModel3D> unchanged_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, 0.10 }
    };
    auto unchanged_fixture{
        BuildJointPolishFixture(
            unchanged_model_list,
            unchanged_model_list)
    };
    const detail::FitStatePatch unchanged_patch;
    const detail::FitStateView unchanged_state_view{
        unchanged_fixture.state,
        unchanged_patch
    };
    alg::WeightedRidgeSolver unchanged_solver;
    EXPECT_FALSE(
        detail::BuildJointPolishProposal(
            unchanged_fixture.context,
            unchanged_state_view,
            key,
            unchanged_fixture.sample_ref_list,
            { 1.0, 1.0 },
            unchanged_solver,
            4.0).has_value());

    alg::WeightedRidgeSolver trust_region_solver;
    const auto small_proposal{
        detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            key,
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            trust_region_solver,
            0.01)
    };
    ASSERT_TRUE(small_proposal.has_value());
    EXPECT_LE(small_proposal->step_norm, 0.01 + 1.0e-12);
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalUsesFitStateViewBaseForTrustRegionOrigin)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, -0.20 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    auto patch{
        detail::FitStatePatch::FromState(
            fixture.state,
            detail::ClusterKey{ 0, 1 })
    };
    patch.mdpde_list.at(0) = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 10.0, 0.55, 0.10 },
        rg::GaussianModel3DUncertainty{}
    };
    const detail::FitStateView base_state_view{ fixture.state, patch };
    EXPECT_DOUBLE_EQ(base_state_view.GetBaseModel(0).GetAmplitude(), 6.0);
    EXPECT_DOUBLE_EQ(base_state_view.GetModel(0).GetAmplitude(), 10.0);
    const auto patched_parameterization{
        detail::JointPolishParameterization::Build(
            std::vector<rg::GaussianModel3D>{
                base_state_view.GetModel(0),
                base_state_view.GetModel(1) })
    };
    ASSERT_TRUE(patched_parameterization.has_value());
    const auto patched_seed{ patched_parameterization->DecodeSeedModels() };
    ASSERT_TRUE(patched_seed.has_value());
    EXPECT_DOUBLE_EQ(patched_seed->at(0).GetAmplitude(), 10.0);

    alg::WeightedRidgeSolver solver;
    EXPECT_FALSE(
        detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            detail::ClusterKey{ 0, 1 },
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            solver,
            0.2).has_value());
}

TEST(EstimatorSecondStageDefenseTest, PhysicalOffsetJacobianMatchesFiniteDifference)
{
    constexpr double step{ 1.0e-6 };
    const std::array<rg::GaussianModel3D, 2> model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.70, -0.15 }
    };
    const std::array<double, 4> distance_list{
        0.0,
        5.0e-6,
        1.0e-5,
        0.35
    };
    for (const auto & model : model_list)
    {
        const auto transformed{ model.ToTransformedCoordinates() };
        ASSERT_TRUE(transformed.has_value());
        for (const auto distance : distance_list)
        {
            const auto evaluation{
                detail::EvaluatePhysicalOffsetResponse(
                    model,
                    distance)
            };
            ASSERT_TRUE(evaluation.has_value());
            EXPECT_TRUE(std::isfinite(evaluation->response));
            EXPECT_TRUE(evaluation->shape_jacobian.allFinite());
            EXPECT_TRUE(std::isfinite(evaluation->offset_jacobian));

            for (Eigen::Index parameter_index = 0;
                parameter_index < evaluation->shape_jacobian.size();
                parameter_index++)
            {
                auto lower{ *transformed };
                auto upper{ *transformed };
                lower(parameter_index) -= step;
                upper(parameter_index) += step;
                lower(static_cast<Eigen::Index>(
                    rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) = 0.0;
                upper(static_cast<Eigen::Index>(
                    rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) = 0.0;
                const auto lower_shape{
                    rg::GaussianModel3D::FromTransformedCoordinates(lower)
                };
                const auto upper_shape{
                    rg::GaussianModel3D::FromTransformedCoordinates(upper)
                };
                ASSERT_TRUE(lower_shape.has_value());
                ASSERT_TRUE(upper_shape.has_value());
                const auto finite_difference{
                    (upper_shape->WithOffset(model.GetOffset())
                            .ResponseAtDistance(distance) -
                        lower_shape->WithOffset(model.GetOffset())
                            .ResponseAtDistance(distance)) /
                    (2.0 * step)
                };
                const auto tolerance{
                    1.0e-6 * std::max(std::abs(finite_difference), 1.0)
                };
                EXPECT_NEAR(
                    evaluation->shape_jacobian(parameter_index),
                    finite_difference,
                    tolerance);
            }

            const auto offset_finite_difference{
                (model.WithOffset(model.GetOffset() + step)
                        .ResponseAtDistance(distance) -
                    model.WithOffset(model.GetOffset() - step)
                        .ResponseAtDistance(distance)) /
                (2.0 * step)
            };
            EXPECT_NEAR(
                evaluation->offset_jacobian,
                offset_finite_difference,
                1.0e-8 * std::max(std::abs(offset_finite_difference), 1.0));
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsJointlyPolishesClusterParameters)
{
    auto model{ BuildJointPolishDefenseModel() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    std::vector<rg::GaussianModel3D> initial_model_list;
    for (const auto * atom : atom_list)
    {
        initial_model_list.emplace_back(GetEstimateModel(*atom));
    }
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fitted_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };
    EXPECT_LT(fitted_error, initial_error);
    bool amplitude_changed{ false };
    bool width_changed{ false };
    bool offset_changed{ false };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto fitted{ GetEstimateModel(*atom_list.at(i)) };
        const auto & initial{ initial_model_list.at(i) };
        amplitude_changed = amplitude_changed ||
            std::abs(fitted.GetAmplitude() - initial.GetAmplitude()) > 1.0e-6;
        width_changed = width_changed ||
            std::abs(fitted.GetWidth() - initial.GetWidth()) > 1.0e-6;
        offset_changed = offset_changed ||
            std::abs(fitted.GetOffset() - initial.GetOffset()) > 1.0e-6;
    }
    EXPECT_TRUE(amplitude_changed);
    EXPECT_TRUE(width_changed);
    EXPECT_TRUE(offset_changed);
    EXPECT_GT(
        std::abs(
            GetEstimateModel(*atom_list.at(0)).GetOffset() -
            GetEstimateModel(*atom_list.at(1)).GetOffset()),
        1.0e-6);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, IndependentOffsetJointPolishIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    auto base_model{ BuildIndependentOffsetDefenseModel() };
    auto scaled_model{ BuildIndependentOffsetDefenseModel(scale) };

    base_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*base_model, MakeSecondStageOptions());
    scaled_model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*scaled_model, MakeSecondStageOptions());

    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    ASSERT_EQ(base_atoms.size(), 2U);
    EXPECT_NE(
        GetEstimateModel(*base_atoms.at(0)).GetOffset(),
        GetEstimateModel(*base_atoms.at(1)).GetOffset());
    EXPECT_NE(
        GetEstimateModel(*scaled_atoms.at(0)).GetOffset(),
        GetEstimateModel(*scaled_atoms.at(1)).GetOffset());
    for (std::size_t i = 0; i < base_atoms.size(); i++)
    {
        const auto base{ GetEstimateModel(*base_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        EXPECT_NEAR(
            base.GetAmplitude() * scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 5.0e-5));
        EXPECT_NEAR(
            base.GetWidth(),
            scaled.GetWidth(),
            std::max(1.0e-8, std::abs(scaled.GetWidth()) * 5.0e-5));
        EXPECT_NEAR(
            base.GetOffset() * scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 5.0e-5));
    }
}
