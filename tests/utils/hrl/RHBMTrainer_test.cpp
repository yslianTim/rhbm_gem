#include <gtest/gtest.h>

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rg = rhbm_gem;

namespace
{
constexpr double kAlphaMin{ 0.0 };
constexpr double kAlphaMax{ 1.0 };
constexpr double kAlphaStep{ 0.5 };

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

double BestAlphaForErrors(
    const std::vector<double> & alpha_list,
    const Eigen::VectorXd & error_sum_list)
{
    int error_min_id{ 0 };
    error_sum_list.minCoeff(&error_min_id);
    return alpha_list.at(static_cast<std::size_t>(error_min_id));
}
} // namespace

TEST(RHBMTrainerTest, TrainAlphaRBuildsAlphaGridWithExactMax)
{
    const auto dataset{ MakeLinearDataset(2.0) };
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;

    const auto result{
        rg::rhbm_trainer::TrainAlphaR({ dataset }, 0.2, 0.8, 0.2, options)
    };

    ASSERT_EQ(result.alpha_grid.size(), 4U);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(0), 0.2);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(1), 0.4);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(2), 0.6);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(3), 0.8);
}

TEST(RHBMTrainerTest, TrainAlphaRBuildsAlphaGridWithoutExceedingMax)
{
    const auto dataset{ MakeLinearDataset(2.0) };
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;

    const auto result{
        rg::rhbm_trainer::TrainAlphaR({ dataset }, 0.0, 1.0, 0.3, options)
    };

    ASSERT_EQ(result.alpha_grid.size(), 4U);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(0), 0.0);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(1), 0.3);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(2), 0.6);
    EXPECT_DOUBLE_EQ(result.alpha_grid.at(3), 0.9);
}

TEST(RHBMTrainerTest, TrainAlphaRejectsInvalidAlphaGrid)
{
    const auto dataset{ MakeLinearDataset(1.0) };
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;

    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR({ dataset }, -0.1, 1.0, 0.1, options),
        std::invalid_argument);
    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR({ dataset }, 1.0, 0.0, 0.1, options),
        std::invalid_argument);
    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR({ dataset }, 0.0, 1.0, 0.0, options),
        std::invalid_argument);
}

TEST(RHBMTrainerTest, TrainAlphaRSingleDatasetWithExactLinearDataReturnsZeroError)
{
    const auto dataset{ MakeLinearDataset(2.0) };
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;

    const auto result{
        rg::rhbm_trainer::TrainAlphaR({ dataset }, 0.0, 0.5, 0.5, options)
    };

    ASSERT_EQ(result.error_sum_list.size(), 2);
    EXPECT_NEAR(0.0, result.error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, result.error_sum_list(1), 1e-12);
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(result.alpha_grid, result.error_sum_list));
}

TEST(RHBMTrainerTest, TrainAlphaGSingleGroupWithIdenticalBetasReturnsZeroError)
{
    const std::vector<Eigen::VectorXd> beta_list(6, MakeVector({ 1.5, -0.5 }));
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 4;

    const auto result{
        rg::rhbm_trainer::TrainAlphaG({ beta_list }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };

    ASSERT_EQ(result.error_sum_list.size(), 3);
    EXPECT_NEAR(0.0, result.error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, result.error_sum_list(1), 1e-12);
    EXPECT_NEAR(0.0, result.error_sum_list(2), 1e-12);
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(result.alpha_grid, result.error_sum_list));
}

TEST(RHBMTrainerTest, TrainAlphaRAggregatesSingleBatchResultsAndSelectsMinimum)
{
    const std::vector<rg::RHBMMemberDataset> dataset_list{
        MakeLinearDataset(2.0),
        MakeLinearDataset(-0.5)
    };

    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;
    const auto result{
        rg::rhbm_trainer::TrainAlphaR(dataset_list, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto first_single{
        rg::rhbm_trainer::TrainAlphaR({ dataset_list.at(0) }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto second_single{
        rg::rhbm_trainer::TrainAlphaR({ dataset_list.at(1) }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto expected_error_sum{ (first_single.error_sum_list + second_single.error_sum_list).eval() };

    EXPECT_TRUE(result.error_sum_list.isApprox(expected_error_sum, 1e-12));
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(result.alpha_grid, result.error_sum_list));
}

TEST(RHBMTrainerTest, TrainAlphaGAggregatesSingleBatchResultsAndSelectsMinimum)
{
    const std::vector<std::vector<Eigen::VectorXd>> beta_group_list{
        {
            MakeVector({ 1.0, 2.0 }),
            MakeVector({ 1.5, 2.5 }),
            MakeVector({ 2.0, 3.0 }),
            MakeVector({ 2.5, 3.5 }),
            MakeVector({ 3.0, 4.0 }),
            MakeVector({ 3.5, 4.5 })
        },
        {
            MakeVector({ -1.0, 0.5 }),
            MakeVector({ -1.5, 1.0 }),
            MakeVector({ -2.0, 1.5 }),
            MakeVector({ -2.5, 2.0 }),
            MakeVector({ -3.0, 2.5 }),
            MakeVector({ -3.5, 3.0 })
        }
    };
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;
    options.execution_options.random_seed = 11U;

    const auto result{
        rg::rhbm_trainer::TrainAlphaG(beta_group_list, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto first_single{
        rg::rhbm_trainer::TrainAlphaG({ beta_group_list.at(0) }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto second_single{
        rg::rhbm_trainer::TrainAlphaG({ beta_group_list.at(1) }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto expected_error_sum{ (first_single.error_sum_list + second_single.error_sum_list).eval() };

    EXPECT_TRUE(result.error_sum_list.isApprox(expected_error_sum, 1e-12));
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(result.alpha_grid, result.error_sum_list));
}

TEST(RHBMTrainerTest, TrainAlphaGWithSeedIsDeterministic)
{
    const std::vector<std::vector<Eigen::VectorXd>> beta_group_list{
        {
            MakeVector({ 1.0, 2.0 }),
            MakeVector({ 1.5, 2.5 }),
            MakeVector({ 2.0, 3.0 }),
            MakeVector({ 2.5, 3.5 }),
            MakeVector({ 3.0, 4.0 }),
            MakeVector({ 3.5, 4.5 })
        },
        {
            MakeVector({ -1.0, 0.5 }),
            MakeVector({ -1.5, 1.0 }),
            MakeVector({ -2.0, 1.5 }),
            MakeVector({ -2.5, 2.0 }),
            MakeVector({ -3.0, 2.5 }),
            MakeVector({ -3.5, 3.0 })
        }
    };
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;
    options.execution_options.random_seed = 7U;

    const auto first{
        rg::rhbm_trainer::TrainAlphaG(beta_group_list, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto second{
        rg::rhbm_trainer::TrainAlphaG(beta_group_list, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto single_first{
        rg::rhbm_trainer::TrainAlphaG({ beta_group_list.at(0) }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };
    const auto single_second{
        rg::rhbm_trainer::TrainAlphaG({ beta_group_list.at(0) }, kAlphaMin, kAlphaMax, kAlphaStep, options)
    };

    EXPECT_TRUE(first.error_sum_list.isApprox(second.error_sum_list, 1e-12));
    EXPECT_TRUE(single_first.error_sum_list.isApprox(single_second.error_sum_list, 1e-12));
    EXPECT_DOUBLE_EQ(first.best_alpha, second.best_alpha);
    EXPECT_DOUBLE_EQ(single_first.best_alpha, single_second.best_alpha);
    EXPECT_DOUBLE_EQ(first.best_alpha, BestAlphaForErrors(first.alpha_grid, first.error_sum_list));
}

TEST(RHBMTrainerTest, TrainAlphaRequiresExplicitSubsetSizeAndRejectsInvalidInputs)
{
    const auto dataset{ MakeLinearDataset(1.0) };
    const std::vector<rg::RHBMMemberDataset> empty_dataset_list;
    const std::vector<rg::RHBMMemberDataset> dataset_list{ dataset };
    rg::rhbm_trainer::AlphaTrainingOptions default_options;

    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR(dataset_list, kAlphaMin, kAlphaMax, kAlphaStep, default_options),
        std::invalid_argument);

    rg::rhbm_trainer::AlphaTrainingOptions alpha_r_empty_options;
    alpha_r_empty_options.subset_size = 3;

    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR(
            empty_dataset_list, kAlphaMin, kAlphaMax, kAlphaStep, alpha_r_empty_options),
        std::invalid_argument);

    rg::rhbm_trainer::AlphaTrainingOptions alpha_r_options;
    alpha_r_options.subset_size = 7;
    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR(dataset_list, kAlphaMin, kAlphaMax, kAlphaStep, alpha_r_options),
        std::invalid_argument);

    rg::rhbm_trainer::AlphaTrainingOptions alpha_g_options;
    alpha_g_options.subset_size = 3;
    const std::vector<std::vector<Eigen::VectorXd>> beta_group_list{
        {
            MakeVector({ 1.0, 2.0 }),
            MakeVector({ 1.5, 2.5 })
        }
    };
    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaG(
            beta_group_list, kAlphaMin, kAlphaMax, kAlphaStep, alpha_g_options),
        std::invalid_argument);
}

TEST(RHBMTrainerTest, TrainAlphaRRejectsInconsistentDatasetShape)
{
    auto dataset{ MakeLinearDataset(1.0) };
    dataset.y = rg::RHBMResponseVector::Zero(dataset.X.rows() - 1);
    rg::rhbm_trainer::AlphaTrainingOptions options;
    options.subset_size = 3;

    EXPECT_THROW(
        rg::rhbm_trainer::TrainAlphaR({ dataset }, kAlphaMin, kAlphaMax, kAlphaStep, options),
        std::invalid_argument);
}
