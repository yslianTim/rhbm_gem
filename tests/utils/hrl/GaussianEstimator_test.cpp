#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimator.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rg = rhbm_gem;
namespace ls = rhbm_gem::linearization_service;

namespace
{

Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (double value : values)
    {
        result(index++) = value;
    }
    return result;
}

rg::RHBMMemberDataset MakeLinearDataset(double slope)
{
    rg::RHBMMemberDataset dataset;
    dataset.X = rg::RHBMDesignMatrix::Zero(6, 2);
    dataset.y = rg::RHBMResponseVector::Zero(6);
    for (Eigen::Index i = 0; i < dataset.X.rows(); i++)
    {
        dataset.X(i, 0) = 1.0;
        dataset.X(i, 1) = static_cast<double>(i);
        dataset.y(i) = 1.0 + static_cast<double>(i) * slope;
    }
    return dataset;
}

LocalPotentialSampleList MakeSampleEntries()
{
    LocalPotentialSampleList sample_entries;
    sample_entries.reserve(6);
    for (int i = 0; i < 6; i++)
    {
        const auto distance{ static_cast<float>(0.1 * static_cast<double>(i)) };
        const auto response{ static_cast<float>(std::exp(1.0 - 0.5 * distance * distance)) };
        sample_entries.emplace_back(LocalPotentialSample{ distance, response, std::nullopt });
    }
    return sample_entries;
}

rg::gaussian_estimator::CrossValidationOptions MakeOptions()
{
    rg::gaussian_estimator::CrossValidationOptions options;
    options.alpha_min = 0.0;
    options.alpha_max = 0.5;
    options.alpha_step = 0.5;
    options.thread_size = 1;
    return options;
}

} // namespace

TEST(GaussianEstimatorTest, SampleListAlphaRReturnsFiniteAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{ MakeSampleEntries() };
    const auto alpha_r{
        rg::gaussian_estimator::CrossValidationAlphaR(sample_entries_list, options)
    };

    EXPECT_TRUE(std::isfinite(alpha_r));
    EXPECT_GE(alpha_r, options.alpha_min);
    EXPECT_LE(alpha_r, options.alpha_max);
}

TEST(GaussianEstimatorTest, AlphaRMatchesAlphaTrainerBestAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<rg::RHBMMemberDataset> dataset_list{
        MakeLinearDataset(2.0),
        MakeLinearDataset(-0.5)
    };
    const rg::rhbm_trainer::AlphaTrainer trainer{
        options.alpha_min,
        options.alpha_max,
        options.alpha_step
    };
    rg::rhbm_trainer::AlphaTrainer::AlphaTrainingOptions trainer_options;
    trainer_options.subset_size = 5;
    trainer_options.execution_options.quiet_mode = true;
    trainer_options.execution_options.thread_size = options.thread_size;

    const auto expected{ trainer.TrainAlphaR(dataset_list, trainer_options).best_alpha };
    const auto actual{
        rg::gaussian_estimator::CrossValidationAlphaR(dataset_list, options)
    };

    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, AlphaGMatchesAlphaTrainerBestAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<std::vector<rg::RHBMParameterVector>> beta_group_list{
        std::vector<rg::RHBMParameterVector>(10, MakeVector({ 1.5, -0.5 }))
    };
    const rg::rhbm_trainer::AlphaTrainer trainer{
        options.alpha_min,
        options.alpha_max,
        options.alpha_step
    };
    rg::rhbm_trainer::AlphaTrainer::AlphaTrainingOptions trainer_options;
    trainer_options.subset_size = 10;
    trainer_options.execution_options.quiet_mode = true;
    trainer_options.execution_options.thread_size = options.thread_size;

    const auto expected{ trainer.TrainAlphaG(beta_group_list, trainer_options).best_alpha };
    const auto actual{
        rg::gaussian_estimator::CrossValidationAlphaG(beta_group_list, options)
    };

    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, QuietAlphaGOptionsDoNotChangeBestAlpha)
{
    auto quiet_options{ MakeOptions() };
    quiet_options.output_summary_log = false;
    quiet_options.output_progress = false;
    auto verbose_options{ quiet_options };
    verbose_options.output_summary_log = true;
    verbose_options.output_progress = true;
    const std::vector<std::vector<rg::RHBMParameterVector>> beta_group_list{
        std::vector<rg::RHBMParameterVector>(10, MakeVector({ 1.5, -0.5 }))
    };

    const auto quiet_alpha{
        rg::gaussian_estimator::CrossValidationAlphaG(beta_group_list, quiet_options)
    };
    const auto verbose_alpha{
        rg::gaussian_estimator::CrossValidationAlphaG(beta_group_list, verbose_options)
    };

    EXPECT_DOUBLE_EQ(quiet_alpha, verbose_alpha);
}

TEST(GaussianEstimatorTest, RejectsEmptyAlphaRTrainingInputs)
{
    const auto options{ MakeOptions() };
    const std::vector<rg::RHBMMemberDataset> empty_dataset_list;

    EXPECT_THROW(
        rg::gaussian_estimator::CrossValidationAlphaR(empty_dataset_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EmptyAlphaGTrainingInputReturnsFallbackAlpha)
{
    auto options{ MakeOptions() };
    options.output_summary_log = false;
    options.output_progress = false;
    const std::vector<std::vector<rg::RHBMParameterVector>> empty_beta_group_list;

    const auto alpha_g{
        rg::gaussian_estimator::CrossValidationAlphaG(empty_beta_group_list, options)
    };

    EXPECT_DOUBLE_EQ(alpha_g, options.alpha_min);
}

TEST(GaussianEstimatorTest, PlotRequestWithEmptyDirectoryDoesNotRequireRoot)
{
    const auto options{ MakeOptions() };
    const std::vector<std::vector<rg::RHBMParameterVector>> beta_group_list{
        std::vector<rg::RHBMParameterVector>(10, MakeVector({ 1.5, -0.5 }))
    };

    EXPECT_NO_THROW({
        const auto alpha_g{
            rg::gaussian_estimator::CrossValidationAlphaG(beta_group_list, options, true)
        };
        EXPECT_TRUE(std::isfinite(alpha_g));
    });
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianMatchesHelperPath)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };
    constexpr double alpha_r{ 0.2 };
    const auto dataset{
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries, options.fit_range_min, options.fit_range_max)
    };
    const auto expected_fit{
        rg::rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset)
    };

    const auto actual{
        rg::gaussian_estimator::EstimateLocalGaussian(sample_entries, alpha_r, options)
    };

    const auto expected_ols{ ls::DecodeParameterVector(expected_fit.beta_ols) };
    const auto expected_mdpde{ ls::DecodeParameterVector(expected_fit.beta_mdpde) };
    EXPECT_NEAR(expected_ols.GetAmplitude(), actual.ols.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_ols.GetWidth(), actual.ols.GetModel().GetWidth(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetAmplitude(), actual.mdpde.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetWidth(), actual.mdpde.GetModel().GetWidth(), 1e-12);
    EXPECT_DOUBLE_EQ(alpha_r, actual.alpha_r);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianMatchesHelperPath)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<double> alpha_r_list{ 0.0, 0.1, 0.2 };
    constexpr double alpha_g{ 0.3 };

    std::vector<rg::RHBMMemberDataset> dataset_list;
    std::vector<rg::RHBMBetaEstimateResult> fit_result_list;
    dataset_list.reserve(sample_entries_list.size());
    fit_result_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        auto dataset{
            rg::rhbm_helper::BuildMemberDataset(
                sample_entries_list.at(i), options.fit_range_min, options.fit_range_max)
        };
        fit_result_list.emplace_back(
            rg::rhbm_helper::EstimateBetaMDPDE(alpha_r_list.at(i), dataset));
        dataset_list.emplace_back(std::move(dataset));
    }
    const auto expected_raw{
        rg::rhbm_helper::EstimateGroup(
            alpha_g,
            rg::rhbm_helper::BuildGroupInput(dataset_list, fit_result_list))
    };

    const auto actual{
        rg::gaussian_estimator::EstimateGroupGaussian(
            sample_entries_list, alpha_r_list, alpha_g, options)
    };

    const auto expected_mean{ ls::DecodeParameterVector(expected_raw.mu_mean) };
    const auto expected_mdpde{ ls::DecodeParameterVector(expected_raw.mu_mdpde) };
    const auto expected_prior{
        ls::DecodeParameterVector(expected_raw.mu_prior, expected_raw.capital_lambda)
    };
    EXPECT_NEAR(expected_mean.GetAmplitude(), actual.mean.GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mean.GetWidth(), actual.mean.GetWidth(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetAmplitude(), actual.mdpde.GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetWidth(), actual.mdpde.GetWidth(), 1e-12);
    EXPECT_NEAR(expected_prior.GetModel().GetAmplitude(), actual.prior.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_prior.GetModel().GetWidth(), actual.prior.GetModel().GetWidth(), 1e-12);
    ASSERT_EQ(sample_entries_list.size(), actual.member_results.size());
}
