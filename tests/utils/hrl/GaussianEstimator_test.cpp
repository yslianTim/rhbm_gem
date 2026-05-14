#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <filesystem>
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

LocalPotentialSampleList MakeSampleEntries(double log_response_shift = 0.0)
{
    LocalPotentialSampleList sample_entries;
    sample_entries.reserve(6);
    for (int i = 0; i < 6; i++)
    {
        const auto distance{ static_cast<float>(0.1 * static_cast<double>(i)) };
        const auto response{
            static_cast<float>(std::exp(1.0 + log_response_shift - 0.5 * distance * distance))
        };
        sample_entries.emplace_back(LocalPotentialSample{ distance, response, std::nullopt });
    }
    return sample_entries;
}

rg::gaussian_estimator::TrainingOptions MakeOptions()
{
    rg::gaussian_estimator::TrainingOptions options;
    options.alpha_min = 0.0;
    options.alpha_max = 0.5;
    options.alpha_step = 0.5;
    options.thread_size = 1;
    return options;
}

std::vector<LocalPotentialSampleList> MakeSampleGroup(std::size_t member_size)
{
    std::vector<LocalPotentialSampleList> sample_group;
    sample_group.reserve(member_size);
    for (std::size_t i = 0; i < member_size; i++)
    {
        sample_group.emplace_back(MakeSampleEntries(0.01 * static_cast<double>(i)));
    }
    return sample_group;
}

std::vector<LocalPotentialSampleList> MakeIdenticalSampleGroup(std::size_t member_size)
{
    std::vector<LocalPotentialSampleList> sample_group;
    sample_group.reserve(member_size);
    for (std::size_t i = 0; i < member_size; i++)
    {
        sample_group.emplace_back(MakeSampleEntries());
    }
    return sample_group;
}

std::vector<rg::LocalGaussianResult> EstimateMemberResults(
    const std::vector<LocalPotentialSampleList> & sample_group,
    const rg::gaussian_estimator::TrainingOptions & options)
{
    std::vector<rg::LocalGaussianResult> member_results;
    member_results.reserve(sample_group.size());
    for (const auto & sample_entries : sample_group)
    {
        member_results.emplace_back(
            rg::gaussian_estimator::EstimateLocalGaussian(sample_entries, 0.0, options));
    }
    return member_results;
}

std::vector<std::vector<rg::RHBMParameterVector>> BuildBetaGroupList(
    const std::vector<std::vector<rg::LocalGaussianResult>> & member_result_list)
{
    std::vector<std::vector<rg::RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(member_result_list.size());
    for (const auto & group_results : member_result_list)
    {
        std::vector<rg::RHBMParameterVector> beta_list;
        beta_list.reserve(group_results.size());
        for (const auto & member_result : group_results)
        {
            beta_list.emplace_back(
                ls::EncodeGaussianToParameterVector(member_result.mdpde.GetModel()));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }
    return beta_group_list;
}

} // namespace

TEST(GaussianEstimatorTest, SampleListAlphaRReturnsFiniteAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{ MakeSampleEntries() };
    const auto alpha_r{
        rg::gaussian_estimator::TrainAlphaR(sample_entries_list, options)
    };

    EXPECT_TRUE(std::isfinite(alpha_r));
    EXPECT_GE(alpha_r, options.alpha_min);
    EXPECT_LE(alpha_r, options.alpha_max);
}

TEST(GaussianEstimatorTest, AlphaRMatchesTrainingFunctionBestAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries(0.2)
    };
    const std::vector<rg::RHBMMemberDataset> dataset_list{
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries_list.at(0), options.fit_range_min, options.fit_range_max),
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries_list.at(1), options.fit_range_min, options.fit_range_max)
    };
    rg::rhbm_trainer::RHBMTrainingOptions trainer_options;
    trainer_options.subset_size = 5;
    trainer_options.execution_options.quiet_mode = true;
    trainer_options.execution_options.thread_size = options.thread_size;

    const auto expected{
        rg::rhbm_trainer::CrossValidationAlphaR(
            dataset_list,
            options.alpha_min,
            options.alpha_max,
            options.alpha_step,
            trainer_options).best_alpha
    };
    const auto actual{
        rg::gaussian_estimator::TrainAlphaR(sample_entries_list, options)
    };

    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, AlphaGMatchesTrainingFunctionBestAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<std::vector<LocalPotentialSampleList>> sample_group_list{
        MakeIdenticalSampleGroup(10)
    };
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        EstimateMemberResults(sample_group_list.front(), options)
    };
    const auto beta_group_list{ BuildBetaGroupList(member_result_list) };
    rg::rhbm_trainer::RHBMTrainingOptions trainer_options;
    trainer_options.subset_size = 10;
    trainer_options.execution_options.quiet_mode = true;
    trainer_options.execution_options.thread_size = options.thread_size;

    const auto expected{
        rg::rhbm_trainer::CrossValidationAlphaG(
            beta_group_list,
            options.alpha_min,
            options.alpha_max,
            options.alpha_step,
            trainer_options).best_alpha
    };
    const auto actual{
        rg::gaussian_estimator::TrainAlphaG(
            sample_group_list, member_result_list, options)
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
    const std::vector<std::vector<LocalPotentialSampleList>> sample_group_list{
        MakeIdenticalSampleGroup(10)
    };
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        EstimateMemberResults(sample_group_list.front(), quiet_options)
    };

    const auto quiet_alpha{
        rg::gaussian_estimator::TrainAlphaG(
            sample_group_list, member_result_list, quiet_options)
    };
    const auto verbose_alpha{
        rg::gaussian_estimator::TrainAlphaG(
            sample_group_list, member_result_list, verbose_options)
    };

    EXPECT_DOUBLE_EQ(quiet_alpha, verbose_alpha);
}

TEST(GaussianEstimatorTest, AlphaRStudyPlotPathDoesNotChangeBestAlpha)
{
    auto options{ MakeOptions() };
    options.study_plot_dir = std::filesystem::path{ testing::TempDir() } / "alpha_r_study";
    std::filesystem::create_directories(options.study_plot_dir);
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries(0.2)
    };

    const auto expected{
        rg::gaussian_estimator::TrainAlphaR(sample_entries_list, options)
    };
    const auto actual{
        rg::gaussian_estimator::TrainAlphaR(sample_entries_list, options, true)
    };

    EXPECT_TRUE(std::isfinite(actual));
    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, AlphaGStudyPlotPathDoesNotChangeBestAlpha)
{
    auto options{ MakeOptions() };
    options.study_plot_dir = std::filesystem::path{ testing::TempDir() } / "alpha_g_study";
    std::filesystem::create_directories(options.study_plot_dir);
    const std::vector<std::vector<LocalPotentialSampleList>> sample_group_list{
        MakeIdenticalSampleGroup(10)
    };
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        EstimateMemberResults(sample_group_list.front(), options)
    };

    const auto expected{
        rg::gaussian_estimator::TrainAlphaG(
            sample_group_list, member_result_list, options)
    };
    const auto actual{
        rg::gaussian_estimator::TrainAlphaG(
            sample_group_list, member_result_list, options, true)
    };

    EXPECT_TRUE(std::isfinite(actual));
    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, RejectsEmptyAlphaRTrainingInputs)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> empty_sample_entries_list;

    EXPECT_THROW(
        rg::gaussian_estimator::TrainAlphaR(empty_sample_entries_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EmptyAlphaGTrainingInputReturnsFallbackAlpha)
{
    auto options{ MakeOptions() };
    options.output_summary_log = false;
    options.output_progress = false;
    const std::vector<std::vector<LocalPotentialSampleList>> empty_sample_group_list;
    const std::vector<std::vector<rg::LocalGaussianResult>> empty_member_result_list;

    const auto alpha_g{
        rg::gaussian_estimator::TrainAlphaG(
            empty_sample_group_list, empty_member_result_list, options)
    };

    EXPECT_DOUBLE_EQ(alpha_g, options.alpha_min);
}

TEST(GaussianEstimatorTest, PlotRequestWithEmptyDirectoryDoesNotRequireRoot)
{
    const auto options{ MakeOptions() };
    const std::vector<std::vector<LocalPotentialSampleList>> sample_group_list{
        MakeSampleGroup(10)
    };
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        EstimateMemberResults(sample_group_list.front(), options)
    };

    EXPECT_NO_THROW({
        const auto alpha_g{
            rg::gaussian_estimator::TrainAlphaG(
                sample_group_list, member_result_list, options, true)
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
    ASSERT_TRUE(actual.fit_result.has_value());
    EXPECT_TRUE(actual.fit_result->beta_mdpde.isApprox(expected_fit.beta_mdpde, 1e-12));
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
    std::vector<rg::LocalGaussianResult> member_result_list;
    constexpr double alpha_g{ 0.3 };

    std::vector<rg::RHBMMemberDataset> dataset_list;
    std::vector<rg::RHBMBetaEstimateResult> fit_result_list;
    dataset_list.reserve(sample_entries_list.size());
    fit_result_list.reserve(sample_entries_list.size());
    member_result_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        auto dataset{
            rg::rhbm_helper::BuildMemberDataset(
                sample_entries_list.at(i), options.fit_range_min, options.fit_range_max)
        };
        auto member_result{
            rg::gaussian_estimator::EstimateLocalGaussian(
                sample_entries_list.at(i), alpha_r_list.at(i), options)
        };
        ASSERT_TRUE(member_result.fit_result.has_value());
        fit_result_list.emplace_back(*member_result.fit_result);
        member_result_list.emplace_back(std::move(member_result));
        dataset_list.emplace_back(std::move(dataset));
    }
    const auto expected_raw{
        rg::rhbm_helper::EstimateGroup(
            alpha_g,
            rg::rhbm_helper::BuildGroupInput(dataset_list, fit_result_list))
    };

    const auto actual{
        rg::gaussian_estimator::EstimateGroupGaussian(
            sample_entries_list, member_result_list, alpha_g, options)
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

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsInconsistentMemberCount)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<rg::LocalGaussianResult> member_result_list{
        rg::gaussian_estimator::EstimateLocalGaussian(sample_entries_list.front(), 0.0, options)
    };

    EXPECT_THROW(
        rg::gaussian_estimator::EstimateGroupGaussian(
            sample_entries_list, member_result_list, 0.0, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsMissingTransientFitState)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<rg::LocalGaussianResult> member_result_list(sample_entries_list.size());

    EXPECT_THROW(
        rg::gaussian_estimator::EstimateGroupGaussian(
            sample_entries_list, member_result_list, 0.0, options),
        std::invalid_argument);
}
