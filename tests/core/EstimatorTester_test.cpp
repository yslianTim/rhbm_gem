#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

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

std::unique_ptr<rg::ModelObject> BuildSecondStageFallbackDiagnosticModel()
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
    rt::RunLocalAlphaTraining(*model, options);
    rt::RunFirstStageLocalFitting(*model, options);

    auto analysis{ model->EditAnalysis() };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.front().response = std::numeric_limits<float>::quiet_NaN();
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(std::move(sampling_entries));
    return model;
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

TEST(EstimatorTesterTest, RunLocalPotentialFittingDoesNotWorsenCoupledResponseResidual)
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    const auto input{
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
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    rg::ModelObject first_stage_model{ *input.replica_model_objects.front() };
    rt::RunLocalAlphaTraining(first_stage_model, options);
    rt::RunFirstStageLocalFitting(first_stage_model, options);
    const auto first_stage_error{
        CalculateSelectedAtomResponseMeanSquaredError(first_stage_model)
    };

    rg::ModelObject full_model{ *input.replica_model_objects.front() };
    rt::RunLocalAlphaTraining(full_model, options);
    rt::RunLocalPotentialFitting(full_model, options);
    const auto full_error{
        CalculateSelectedAtomResponseMeanSquaredError(full_model)
    };

    const auto tolerance{ 1.0e-3 * std::max(first_stage_error, 1.0) };
    EXPECT_LE(full_error, first_stage_error + tolerance);
}

TEST(EstimatorTesterTest, RunLocalPotentialFittingHandlesCoupledNonzeroOffsetNeighbors)
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    const auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::C,
            Element::CARBON,
            0.15,
            rg::GaussianModel3D{ 6.0, 0.55, -0.2 },
            potential_model,
            0.0,
            1,
            77
        })
    };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    rg::ModelObject first_stage_model{ *input.replica_model_objects.front() };
    rt::RunLocalAlphaTraining(first_stage_model, options);
    rt::RunFirstStageLocalFitting(first_stage_model, options);
    const auto first_stage_error{
        CalculateSelectedAtomResponseMeanSquaredError(first_stage_model)
    };

    rg::ModelObject full_model{ *input.replica_model_objects.front() };
    rt::RunLocalAlphaTraining(full_model, options);
    rt::RunLocalPotentialFitting(full_model, options);
    const auto full_error{
        CalculateSelectedAtomResponseMeanSquaredError(full_model)
    };

    const auto tolerance{ 1.0e-3 * std::max(first_stage_error, 1.0) };
    EXPECT_LE(full_error, first_stage_error + tolerance);
}

TEST(EstimatorTesterTest, RunLocalPotentialFittingCoupledNonzeroOffsetDoesNotWarnFallback)
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    const auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::C,
            Element::CARBON,
            0.15,
            rg::GaussianModel3D{ 6.0, 0.55, -0.2 },
            potential_model,
            0.0,
            1,
            77
        })
    };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = false;

    rg::ModelObject model{ *input.replica_model_objects.front() };
    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunLocalAlphaTraining(model, options);
    rt::RunLocalPotentialFitting(model, options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingLogsFallbackSummary)
{
    auto model{ BuildSecondStageFallbackDiagnosticModel() };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = false;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, model->GetSelectedAtoms(), options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("joint offset fallback iterations = 1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("refit fallback atom-events = 1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("refit fallback distinct atoms = 1"),
        std::string::npos);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingQuietModeSuppressesFallbackSummary)
{
    auto model{ BuildSecondStageFallbackDiagnosticModel() };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, model->GetSelectedAtoms(), options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
}
