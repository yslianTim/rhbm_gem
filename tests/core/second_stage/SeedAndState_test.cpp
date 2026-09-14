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

using second_stage_test::BuildNearCollinearDefenseModel;
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

TEST(EstimatorSecondStageDefenseTest, TransformedChangeIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 8.8, 0.55, -0.12 };
    const auto base_change{
        detail::CalculateTransformedChange(current, previous)
    };

    for (const auto scale : { 1.0e-2, 1.0e2 })
    {
        const auto scaled_change{
            detail::CalculateTransformedChange(
                rg::GaussianModel3D{
                    current.GetAmplitude() * scale,
                    current.GetWidth(),
                    current.GetOffset() * scale
                },
                rg::GaussianModel3D{
                    previous.GetAmplitude() * scale,
                    previous.GetWidth(),
                    previous.GetOffset() * scale
                })
        };
        ASSERT_EQ(base_change.size(), scaled_change.size());
        for (std::size_t i = 0; i < base_change.size(); i++)
        {
            EXPECT_NEAR(
                base_change.at(i),
                scaled_change.at(i),
                1.0e-12);
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, GaussianParameterMedianUsesValidComponentMedians)
{
    const std::vector<rg::GaussianModel3D> models{
        { 1.0, 0.6, 1.0 }, { 9.0, 0.2, 3.0 }, { 5.0, 0.4, 2.0 }, {} };
    const auto odd{ detail::BuildGaussianParameterMedian(models) };
    ASSERT_TRUE(odd.has_value());
    ExpectGaussianModelsNear(*odd, { 5.0, 0.4, 2.0 }, 1.0e-12);
    const auto even{ detail::BuildGaussianParameterMedian(
        { { 4.0, 0.5, -2.0 }, { 8.0, 0.9, 2.0 } }) };
    ASSERT_TRUE(even.has_value());
    ExpectGaussianModelsNear(*even, { 6.0, 0.7, 0.0 }, 1.0e-12);
    EXPECT_FALSE(detail::BuildGaussianParameterMedian({}).has_value());
    EXPECT_FALSE(detail::BuildGaussianParameterMedian({ { } }).has_value());
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

TEST(EstimatorSecondStageDefenseTest, GaussianParameterMedianIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    const std::vector<rg::GaussianModel3D> models{
        { 3.0, 0.4, -0.2 }, { 5.0, 0.6, 0.1 }, { 7.0, 0.8, 0.4 } };
    std::vector<rg::GaussianModel3D> scaled_models;
    for (const auto & model : models)
    {
        scaled_models.emplace_back(scale * model.GetAmplitude(), model.GetWidth(),
            scale * model.GetOffset());
    }
    const auto median{ detail::BuildGaussianParameterMedian(models) };
    const auto scaled{ detail::BuildGaussianParameterMedian(scaled_models) };
    ASSERT_TRUE(median.has_value());
    ASSERT_TRUE(scaled.has_value());
    EXPECT_DOUBLE_EQ(scale * median->GetAmplitude(), scaled->GetAmplitude());
    EXPECT_DOUBLE_EQ(median->GetWidth(), scaled->GetWidth());
    EXPECT_DOUBLE_EQ(scale * median->GetOffset(), scaled->GetOffset());
}

TEST(EstimatorSecondStageDefenseTest, TransformedDampingIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 9.0, 0.60, -0.15 };
    constexpr double damping{ 0.25 };

    const auto damp = [&](const rg::GaussianModel3D & lhs,
                          const rg::GaussianModel3D & rhs)
    {
        const auto lhs_coordinates{ lhs.ToTransformedCoordinates() };
        const auto rhs_coordinates{ rhs.ToTransformedCoordinates() };
        EXPECT_TRUE(lhs_coordinates.has_value());
        EXPECT_TRUE(rhs_coordinates.has_value());
        return rg::GaussianModel3D::FromTransformedCoordinates(
            *lhs_coordinates + damping * (*rhs_coordinates - *lhs_coordinates));
    };

    const auto base{ damp(previous, current) };
    ASSERT_TRUE(base.has_value());
    for (const auto scale : { 1.0e-2, 1.0e2 })
    {
        const auto scaled{
            damp(
                rg::GaussianModel3D{
                    previous.GetAmplitude() * scale,
                    previous.GetWidth(),
                    previous.GetOffset() * scale
                },
                rg::GaussianModel3D{
                    current.GetAmplitude() * scale,
                    current.GetWidth(),
                    current.GetOffset() * scale
                })
        };
        ASSERT_TRUE(scaled.has_value());
        EXPECT_NEAR(base->GetAmplitude() * scale, scaled->GetAmplitude(), 1.0e-10);
        EXPECT_NEAR(base->GetWidth(), scaled->GetWidth(), 1.0e-12);
        EXPECT_NEAR(base->GetOffset() * scale, scaled->GetOffset(), 1.0e-12);
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedBacktrackingIncludesOffset)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D endpoint{ 12.0, 0.75, 0.40 };
    const auto previous_coordinates{ previous.ToTransformedCoordinates() };
    const auto endpoint_coordinates{ endpoint.ToTransformedCoordinates() };
    ASSERT_TRUE(previous_coordinates.has_value());
    ASSERT_TRUE(endpoint_coordinates.has_value());

    double previous_offset_distance{
        std::abs(endpoint.GetOffset() - previous.GetOffset())
    };
    for (const auto factor : { 0.5, 0.25, 0.125 })
    {
        const auto candidate{
            rg::GaussianModel3D::FromTransformedCoordinates(
                *previous_coordinates +
                factor * (*endpoint_coordinates - *previous_coordinates))
        };
        ASSERT_TRUE(candidate.has_value());
        const auto offset_distance{
            std::abs(candidate->GetOffset() - previous.GetOffset())
        };
        EXPECT_LT(offset_distance, previous_offset_distance);
        previous_offset_distance = offset_distance;
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedExtrapolationKeepsPositiveShape)
{
    const auto left{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 }.ToTransformedCoordinates()
    };
    const auto right{
        rg::GaussianModel3D{ 9.0, 0.60, -0.15 }.ToTransformedCoordinates()
    };
    ASSERT_TRUE(left.has_value());
    ASSERT_TRUE(right.has_value());

    const auto extrapolated{
        rg::GaussianModel3D::FromTransformedCoordinates(
            2.0 * *right - *left)
    };
    ASSERT_TRUE(extrapolated.has_value());
    EXPECT_GT(extrapolated->GetAmplitude(), 0.0);
    EXPECT_GT(extrapolated->GetWidth(), 0.0);
    EXPECT_TRUE(std::isfinite(extrapolated->GetOffset()));
}

TEST(EstimatorSecondStageDefenseTest, TransformedChangeSeparatesPeakHeightAndWidth)
{
    const auto change{
        detail::CalculateTransformedChange(
            rg::GaussianModel3D{ 8.0, 1.0, 0.0 },
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 })
    };

    EXPECT_NEAR(
        0.0,
        change.at(rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1.0e-12);
    EXPECT_NEAR(
        std::log(2.0),
        change.at(rg::GaussianModel3D::LogWidthCoordinateIndex()),
        1.0e-12);
    EXPECT_DOUBLE_EQ(
        0.0,
        change.at(rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()));
}

TEST(
    EstimatorSecondStageDefenseTest,
    MissingGroupSeedsUsesLocalMdpdeAndGlobalFallbackWithoutGroupFitting)
{
    auto model{ BuildNearCollinearDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    for (auto * atom : model->GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::For(*atom).GetGaussianResult(
                FittingStage::Second)
        };
        ASSERT_TRUE(detail::IsValidSecondStageGaussianModel(
            result.mdpde.GetModel()));
        ASSERT_TRUE(detail::IsValidSecondStageGaussianModel(
            result.ols.GetModel()));
        if (atom == model->GetSelectedAtoms().front())
        {
            result.mdpde = rg::GaussianModel3DWithUncertainty{
                rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
                rg::GaussianModel3DUncertainty{}
            };
        }
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            std::move(result));
        LocalPotentialSampleList sentinel_peeling_sampling_entries{
            LocalPotentialSample{
                100.0 + static_cast<double>(atom->GetSerialID()),
                SamplingPoint{ 0.5, atom->GetPosition(), true }
            }
        };
        analysis.SetAtomLocalPeelingSamplingEntries(
            *atom, sentinel_peeling_sampling_entries);
        analysis.SetAtomLocalNeighborCountForPeeling(*atom, 99);
    }
    const auto previous_analysis_view{ model->GetAnalysisView() };
    std::vector<rg::GaussianModel3D> previous_group_prior_list;
    for (const auto group_key : previous_analysis_view.CollectAtomGroupKeys())
    {
        previous_group_prior_list.emplace_back(
            previous_analysis_view.GetAtomGroupPrior(group_key));
    }

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, options);

    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        const auto * atom{ model->GetSelectedAtoms().at(i) };
        EXPECT_TRUE(detail::IsValidSecondStageGaussianModel(
            GetEstimateModel(*atom)));
        const auto peeling_sampling_entries{
            rg::AtomLocalPotentialView::For(
                *atom)
                .GetPeelingSamplingEntries(false)
        };
        EXPECT_GT(peeling_sampling_entries.size(), 1U);
        EXPECT_NE(
            rg::AtomLocalPotentialView::For(
                *atom)
                .GetNeighborCountForPeeling(),
            99);
        EXPECT_FALSE(
            rg::AtomLocalPotentialView::For(*atom)
                .GetGroupMemberResult().has_value());
    }
    const auto final_analysis_view{ model->GetAnalysisView() };
    const auto group_key_list{ final_analysis_view.CollectAtomGroupKeys() };
    ASSERT_EQ(group_key_list.size(), previous_group_prior_list.size());
    for (std::size_t i = 0; i < group_key_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            final_analysis_view.GetAtomGroupPrior(group_key_list.at(i)),
            previous_group_prior_list.at(i),
            0.0);
    }
}
