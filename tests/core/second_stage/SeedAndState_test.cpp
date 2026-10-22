#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
using rhbm_gem::FittingStage;

using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeSecondStageOptions;

} // namespace

TEST(EstimatorSecondStageDefenseTest, SeedSelectionUsesLocalMdpdeThenGlobalMedian)
{
    const auto make_candidate = [](double amplitude)
    {
        return rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ amplitude, 0.5, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
        };
    };
    const auto invalid_candidate{
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
            rg::GaussianModel3DUncertainty{}
        }
    };
    auto local_mdpde{ make_candidate(1.0) };
    std::optional<rg::GaussianModel3D> global_median{
        make_candidate(4.0).GetModel()
    };

    const auto expect_source = [&](detail::SecondStageSeedSource source)
    {
        const auto selection{
            detail::SelectSecondStageSeed(local_mdpde, global_median)
        };
        ASSERT_TRUE(selection.has_value());
        EXPECT_EQ(selection->source, source);
    };
    expect_source(detail::SecondStageSeedSource::LocalMdpde);
    local_mdpde = invalid_candidate;
    expect_source(detail::SecondStageSeedSource::GlobalMedian);
    global_median = invalid_candidate.GetModel();
    EXPECT_FALSE(detail::SelectSecondStageSeed(
        local_mdpde,
        global_median).has_value());
}

TEST(EstimatorSecondStageDefenseTest, SeedSelectionReturnsCompleteSourceModelAndUncertainty)
{
    auto local_mdpde{
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 6.0, 0.55, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
        }
    };
    const std::optional<rg::GaussianModel3D> global_median{
        rg::GaussianModel3D{ 4.0, 0.65, 0.4 }
    };

    const auto local_selection{
        detail::SelectSecondStageSeed(local_mdpde, global_median)
    };
    ASSERT_TRUE(local_selection.has_value());
    EXPECT_EQ(local_selection->source, detail::SecondStageSeedSource::LocalMdpde);
    EXPECT_DOUBLE_EQ(local_selection->model.GetModel().GetAmplitude(), 6.0);
    EXPECT_DOUBLE_EQ(local_selection->model.GetModel().GetWidth(), 0.55);
    EXPECT_DOUBLE_EQ(local_selection->model.GetModel().GetOffset(), -0.2);
    EXPECT_DOUBLE_EQ(
        local_selection->model.GetStandardDeviationModel().GetAmplitude(),
        0.1);
    EXPECT_DOUBLE_EQ(
        local_selection->model.GetStandardDeviationModel().GetWidth(),
        0.02);
    EXPECT_DOUBLE_EQ(
        local_selection->model.GetStandardDeviationModel().GetOffset(),
        0.03);

    local_mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
        rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 }
    };
    const auto median_selection{
        detail::SelectSecondStageSeed(local_mdpde, global_median)
    };
    ASSERT_TRUE(median_selection.has_value());
    EXPECT_EQ(median_selection->source, detail::SecondStageSeedSource::GlobalMedian);
    EXPECT_DOUBLE_EQ(median_selection->model.GetModel().GetAmplitude(), 4.0);
    EXPECT_DOUBLE_EQ(median_selection->model.GetModel().GetWidth(), 0.65);
    EXPECT_DOUBLE_EQ(median_selection->model.GetModel().GetOffset(), 0.4);
    EXPECT_DOUBLE_EQ(
        median_selection->model.GetStandardDeviationModel().GetAmplitude(),
        0.0);
    EXPECT_DOUBLE_EQ(
        median_selection->model.GetStandardDeviationModel().GetWidth(),
        0.0);
    EXPECT_DOUBLE_EQ(
        median_selection->model.GetStandardDeviationModel().GetOffset(),
        0.0);
}


TEST(EstimatorSecondStageDefenseTest,
    DampedModelsInterpolateIndividualPhysicalOffsets)
{
    const std::vector<rg::GaussianModel3D> previous_model_list{
        rg::GaussianModel3D{ 4.0, 0.40, 0.1 },
        rg::GaussianModel3D{ 6.0, 0.60, 0.9 },
        rg::GaussianModel3D{ 8.0, 0.80, -1.0 }
    };
    const std::vector<rg::GaussianModel3D> raw_model_list{
        rg::GaussianModel3D{ 5.0, 0.50, 0.5 },
        rg::GaussianModel3D{ 9.0, 0.90, 0.7 },
        rg::GaussianModel3D{ 7.0, 0.70, 2.0 }
    };

    for (const auto damping : std::array<double, 3>{ 0.0, 0.25, 1.0 })
    {
        const auto candidate_model_list{
            detail::BuildDampedModelList(
                previous_model_list,
                raw_model_list,
                damping)
        };
        ASSERT_TRUE(candidate_model_list.has_value());
        for (std::size_t atom_position = 0;
            atom_position < candidate_model_list->size();
            atom_position++)
        {
            EXPECT_NEAR(candidate_model_list->at(atom_position).GetOffset(),
                std::lerp(previous_model_list.at(atom_position).GetOffset(),
                    raw_model_list.at(atom_position).GetOffset(), damping), 1.0e-12);
            const auto previous_coordinates{
                previous_model_list.at(atom_position).ToTransformedCoordinates()
            };
            const auto raw_coordinates{
                raw_model_list.at(atom_position).ToTransformedCoordinates()
            };
            const auto candidate_coordinates{
                candidate_model_list->at(atom_position).ToTransformedCoordinates()
            };
            ASSERT_TRUE(previous_coordinates.has_value());
            ASSERT_TRUE(raw_coordinates.has_value());
            ASSERT_TRUE(candidate_coordinates.has_value());
            for (const auto parameter_index : std::array<int, 2>{
                rg::GaussianModel3D::LogPeakHeightCoordinateIndex(),
                rg::GaussianModel3D::LogWidthCoordinateIndex() })
            {
                const auto eigen_index{
                    static_cast<Eigen::Index>(parameter_index)
                };
                EXPECT_NEAR(
                    (*candidate_coordinates)(eigen_index),
                    (*previous_coordinates)(eigen_index) +
                        damping * ((*raw_coordinates)(eigen_index) -
                            (*previous_coordinates)(eigen_index)),
                    1.0e-12);
            }
        }
    }
}
