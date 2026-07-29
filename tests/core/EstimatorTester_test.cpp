#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include "core/detail/GaussianEstimatorStages.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {
namespace rt = rhbm_gem::core;
namespace tdf = rhbm_gem::core;
namespace rg = rhbm_gem;

tdf::GaussianParameterDistribution MakeDistribution(
    const rg::GaussianModel3D & mean,
    const rg::GaussianModel3DUncertainty & sigma = rg::GaussianModel3DUncertainty{ 0.05, 0.025, 0.01 })
{
    return tdf::GaussianParameterDistribution{ mean, sigma };
}

rt::LocalTestOptions MakeLocalOptions(
    double alpha_r,
    bool alpha_training)
{
    rt::LocalTestOptions options;
    options.requested_alpha_r = alpha_r;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

rt::GroupTestOptions MakeGroupOptions(
    double alpha_g,
    bool alpha_training)
{
    rt::GroupTestOptions options;
    options.requested_alpha_g = alpha_g;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

rt::FitOptions MakeSecondStageOptions()
{
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    return options;
}

void ExpectBiasStatisticSize(const rt::BiasStatistics & bias)
{
    EXPECT_EQ(bias.mean.size(), rg::GaussianModel3D::ParameterSize());
    EXPECT_EQ(bias.sigma.size(), rg::GaussianModel3D::ParameterSize());
}

double Distance(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    const auto dx{ static_cast<double>(lhs.at(0) - rhs.at(0)) };
    const auto dy{ static_cast<double>(lhs.at(1) - rhs.at(1)) };
    const auto dz{ static_cast<double>(lhs.at(2) - rhs.at(2)) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double CalculateSelectedAtomResponseMeanSquaredError(const rg::ModelObject & model_object)
{
    double squared_error_sum{ 0.0 };
    std::size_t sample_count{ 0 };
    const auto & selected_atoms{ model_object.GetSelectedAtoms() };
    for (const auto * atom : selected_atoms)
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        for (const auto & sample : local_view.GetSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{ rg::AtomLocalPotentialView::RequireFor(*fitted_atom) };
                fitted_response += fitted_view.GetEstimateMDPDE().ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ static_cast<double>(sample.response) - fitted_response };
            squared_error_sum += residual * residual;
            sample_count++;
        }
    }
    return squared_error_sum / static_cast<double>(sample_count);
}

void SetSelectedAtomPosteriorFromMdpde(rg::ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult()
        };
        result.posterior = result.mdpde;
        analysis.EnsureAtomLocalPotential(*atom).SetGaussianResult(
            std::move(result));
    }
}

std::unique_ptr<rg::ModelObject> BuildSecondStageScaleDiagnosticModel()
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };
    const auto options{ MakeSecondStageOptions() };
    rt::RunLocalAlphaTraining(*model, options, rt::LocalFittingPass::FirstStage);
    rt::RunFixedOffsetLocalFitting(*model, options, rt::LocalFittingPass::FirstStage);
    SetSelectedAtomPosteriorFromMdpde(*model);
    return model;
}

void SetSelectedAtomEstimateModel(
    rg::ModelObject & model_object,
    const rg::GaussianModel3D & model)
{
    auto analysis{ model_object.EditAnalysis() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult()
        };
        result.ols = rg::GaussianModel3DWithUncertainty{
            model,
            result.ols.GetStandardDeviationModel()
        };
        result.mdpde = rg::GaussianModel3DWithUncertainty{
            model,
            result.mdpde.GetStandardDeviationModel()
        };
        analysis.EnsureAtomLocalPotential(*atom).SetGaussianResult(std::move(result));
    }
}

void RewriteSamplingResponsesFromSelectedAtomEstimates(rg::ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & selected_atoms{ model_object.GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto sampling_entries{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
        };
        for (auto & sample : sampling_entries)
        {
            double response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::RequireFor(*fitted_atom)
                };
                response += fitted_view.GetEstimateMDPDE().ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            sample.response = static_cast<float>(response);
        }
        analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(
            std::move(sampling_entries));
    }
}

std::unique_ptr<rg::ModelObject> BuildSecondStageSuspiciousOffsetDiagnosticModel()
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    rt::RunLocalAlphaTraining(*model, options, rt::LocalFittingPass::FirstStage);
    rt::RunFixedOffsetLocalFitting(*model, options, rt::LocalFittingPass::FirstStage);
    SetSelectedAtomPosteriorFromMdpde(*model);

    const auto & atom_list{ model->GetSelectedAtoms() };
    auto * target_atom{ atom_list.at(0) };
    auto target_position{ target_atom->GetPosition() };

    auto analysis{ model->EditAnalysis() };
    auto target_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetSamplingEntries(false)
    };
    target_sampling_entries.resize(256);
    target_sampling_entries.front().response = 0.0F;
    target_sampling_entries.front().point.distance = 0.0F;
    target_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < target_sampling_entries.size(); i++)
    {
        auto & sample{ target_sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(target_sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }
    analysis.EnsureAtomLocalPotential(*target_atom).SetSamplingEntries(
        std::move(target_sampling_entries));
    return model;
}

void ExpectSelectedAtomEstimatesAreFinite(const rg::ModelObject & model_object)
{
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        const auto model{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE()
        };
        EXPECT_TRUE(std::isfinite(model.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(model.GetWidth()));
        EXPECT_TRUE(std::isfinite(model.GetOffset()));
    }
}

} // namespace

TEST(EstimatorTesterTest, RunLocalEstimationTestPopulatesBiasOutputs)
{
    constexpr double alpha_r{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(alpha_r, true))
    };

    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_median.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_median.value(), 0.0);
}

TEST(EstimatorTesterTest, RunGroupEstimationTestPopulatesBiasOutputs)
{
    constexpr double alpha_g{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunGroupEstimationTest(test_input, MakeGroupOptions(alpha_g, true))
    };

    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_median.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_median.value(), 0.0);
}

TEST(EstimatorTesterTest, RunLocalEstimationTestSkipsTrainedAlphaWhenDisabled)
{
    constexpr double alpha_r{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(alpha_r, false))
    };

    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_median.has_value());
}

TEST(EstimatorTesterTest, RunGroupEstimationTestSkipsTrainedAlphaWhenDisabled)
{
    constexpr double alpha_g{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunGroupEstimationTest(test_input, MakeGroupOptions(alpha_g, false))
    };

    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_median.has_value());
}

TEST(
    EstimatorTesterTest,
    RunPotentialFittingWorkflowBootstrapsGroupPosteriorFromAdjustedSamples)
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };
    ASSERT_EQ(model->GetSelectedAtoms().size(), 1u);
    const auto initial_view{
        rg::AtomLocalPotentialView::RequireFor(
            *model->GetSelectedAtoms().front())
    };
    ASSERT_FALSE(initial_view.GetSamplingEntries(false).empty());
    ASSERT_TRUE(initial_view.GetSamplingEntries(false, true).empty());

    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    testing::internal::CaptureStdout();
    rt::RunPotentialFittingWorkflow(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };

    const auto fitted_view{
        rg::AtomLocalPotentialView::RequireFor(
            *model->GetSelectedAtoms().front())
    };
    EXPECT_FALSE(fitted_view.GetSamplingEntries(false, true).empty());
    EXPECT_NE(
        out.find(
            "Selected second-stage initial seeds = 1, sources = "
            "group-posterior:1, group-prior:0, group-median:0, "
            "global-median:0."),
        std::string::npos);
    EXPECT_EQ(out.find("stop_reason=no-valid-seed"), std::string::npos);
}

TEST(EstimatorTesterTest, RunLocalEstimationTestRejectsNonFiniteTruth)
{
    auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = rg::GaussianModel3D{
        std::numeric_limits<double>::quiet_NaN(),
        0.5,
        0.0
    };

    EXPECT_THROW(
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(0.5, true)),
        std::invalid_argument
    );
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingImprovesBadFiniteEntryScale)
{
    auto model{ BuildSecondStageScaleDiagnosticModel() };
    SetSelectedAtomEstimateModel(
        *model,
        rg::GaussianModel3D{ 1.0e4, 0.25, 1.0e3 });
    const auto entry_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_LT(fitted_error, entry_error);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingHandlesNearPerfectEntryScale)
{
    auto model{ BuildSecondStageScaleDiagnosticModel() };
    RewriteSamplingResponsesFromSelectedAtomEstimates(*model);
    const auto entry_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_LE(fitted_error, entry_error + 1.0e-8);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingRollsBackSuspiciousJointOffset)
{
    auto model{ BuildSecondStageSuspiciousOffsetDiagnosticModel() };
    auto * target_atom{ model->GetSelectedAtoms().front() };
    const auto previous_offset{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetEstimateMDPDE().GetOffset()
    };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    rt::RunSecondStageLocalFitting(*model, options);

    const auto fitted_offset{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetEstimateMDPDE().GetOffset()
    };
    EXPECT_NEAR(fitted_offset, previous_offset, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
