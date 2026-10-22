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
