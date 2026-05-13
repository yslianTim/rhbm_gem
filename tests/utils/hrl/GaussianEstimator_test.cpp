#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimator.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rg = rhbm_gem;

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

TEST(GaussianEstimatorTest, SimpleAlphaRReturnsFiniteDefaultAlpha)
{
    const auto alpha_r{ rg::gaussian_estimator::CrossValidationAlphaR(MakeSampleEntries()) };

    EXPECT_TRUE(std::isfinite(alpha_r));
    EXPECT_GE(alpha_r, 0.0);
    EXPECT_LE(alpha_r, 2.0);
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
    const auto options{ MakeOptions() };
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
