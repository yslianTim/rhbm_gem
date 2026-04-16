#include <gtest/gtest.h>

#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLAlphaTrainer.hpp>

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

HRLMemberDataset MakeLinearDataset(double slope)
{
    const SeriesPointList data_list{
        SeriesPoint({ 1.0, 0.0 }, 1.0),
        SeriesPoint({ 1.0, 1.0 }, 1.0 + slope),
        SeriesPoint({ 1.0, 2.0 }, 1.0 + 2.0 * slope),
        SeriesPoint({ 1.0, 3.0 }, 1.0 + 3.0 * slope),
        SeriesPoint({ 1.0, 4.0 }, 1.0 + 4.0 * slope),
        SeriesPoint({ 1.0, 5.0 }, 1.0 + 5.0 * slope)
    };
    return HRLDataTransform::BuildMemberDataset(data_list);
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

TEST(HRLAlphaTrainerTest, ConstructorBuildsAlphaGridWithExactMax)
{
    const HRLAlphaTrainer trainer{ 0.2, 0.8, 0.2 };
    const auto & alpha_grid{ trainer.AlphaGrid() };
    const auto summary{ trainer.GetAlphaGridSummary() };

    ASSERT_EQ(alpha_grid.size(), 4U);
    EXPECT_DOUBLE_EQ(alpha_grid.at(0), 0.2);
    EXPECT_DOUBLE_EQ(alpha_grid.at(1), 0.4);
    EXPECT_DOUBLE_EQ(alpha_grid.at(2), 0.6);
    EXPECT_DOUBLE_EQ(alpha_grid.at(3), 0.8);
    EXPECT_EQ(
        summary.str(),
        "Alpha training search grid: min = 0.200000, max = 0.800000, "
        "step = 0.200000, count = 4");
}

TEST(HRLAlphaTrainerTest, ConstructorBuildsAlphaGridWithoutExceedingMax)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.3 };
    const auto & alpha_grid{ trainer.AlphaGrid() };
    const auto summary{ trainer.GetAlphaGridSummary() };

    ASSERT_EQ(alpha_grid.size(), 4U);
    EXPECT_DOUBLE_EQ(alpha_grid.at(0), 0.0);
    EXPECT_DOUBLE_EQ(alpha_grid.at(1), 0.3);
    EXPECT_DOUBLE_EQ(alpha_grid.at(2), 0.6);
    EXPECT_DOUBLE_EQ(alpha_grid.at(3), 0.9);
    EXPECT_EQ(
        summary.str(),
        "Alpha training search grid: min = 0.000000, max = 1.000000, "
        "step = 0.300000, count = 4");
}

TEST(HRLAlphaTrainerTest, ConstructorRejectsInvalidAlphaGrid)
{
    EXPECT_THROW(HRLAlphaTrainer(-0.1, 1.0, 0.1), std::invalid_argument);
    EXPECT_THROW(HRLAlphaTrainer(1.0, 0.0, 0.1), std::invalid_argument);
    EXPECT_THROW(HRLAlphaTrainer(0.0, 1.0, 0.0), std::invalid_argument);
}

TEST(HRLAlphaTrainerTest, TrainAlphaRSingleDatasetWithExactLinearDataReturnsZeroError)
{
    const SeriesPointList data_list{
        SeriesPoint({ 1.0, 0.0 }, 1.0),
        SeriesPoint({ 1.0, 1.0 }, 3.0),
        SeriesPoint({ 1.0, 2.0 }, 5.0),
        SeriesPoint({ 1.0, 3.0 }, 7.0),
        SeriesPoint({ 1.0, 4.0 }, 9.0),
        SeriesPoint({ 1.0, 5.0 }, 11.0)
    };
    const auto dataset{ HRLDataTransform::BuildMemberDataset(data_list) };
    const HRLAlphaTrainer trainer{ 0.0, 0.5, 0.5 };
    HRLAlphaTrainer::AlphaRTrainingOptions options;
    options.subset_size = 3;

    const auto result{ trainer.TrainAlphaR({ dataset }, options) };

    ASSERT_EQ(result.error_sum_list.size(), 2);
    EXPECT_NEAR(0.0, result.error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, result.error_sum_list(1), 1e-12);
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(trainer.AlphaGrid(), result.error_sum_list));
}

TEST(HRLAlphaTrainerTest, TrainAlphaGSingleGroupWithIdenticalBetasReturnsZeroError)
{
    const std::vector<Eigen::VectorXd> beta_list(6, MakeVector({ 1.5, -0.5 }));
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    HRLAlphaTrainer::AlphaGTrainingOptions options;
    options.subset_size = 4;

    const auto result{ trainer.TrainAlphaG({ beta_list }, options) };

    ASSERT_EQ(result.error_sum_list.size(), 3);
    EXPECT_NEAR(0.0, result.error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, result.error_sum_list(1), 1e-12);
    EXPECT_NEAR(0.0, result.error_sum_list(2), 1e-12);
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(trainer.AlphaGrid(), result.error_sum_list));
}

TEST(HRLAlphaTrainerTest, TrainAlphaRAggregatesSingleBatchResultsAndSelectsMinimum)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    const auto & alpha_list{ trainer.AlphaGrid() };
    const std::vector<HRLMemberDataset> dataset_list{
        MakeLinearDataset(2.0),
        MakeLinearDataset(-0.5)
    };

    HRLAlphaTrainer::AlphaRTrainingOptions options;
    options.subset_size = 3;
    const auto result{ trainer.TrainAlphaR(dataset_list, options) };
    const auto first_single{ trainer.TrainAlphaR({ dataset_list.at(0) }, options) };
    const auto second_single{ trainer.TrainAlphaR({ dataset_list.at(1) }, options) };
    const auto expected_error_sum{ (first_single.error_sum_list + second_single.error_sum_list).eval() };

    EXPECT_TRUE(result.error_sum_list.isApprox(expected_error_sum, 1e-12));
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(alpha_list, result.error_sum_list));
}

TEST(HRLAlphaTrainerTest, TrainAlphaGAggregatesSingleBatchResultsAndSelectsMinimum)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    const auto & alpha_list{ trainer.AlphaGrid() };
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
    HRLAlphaTrainer::AlphaGTrainingOptions options;
    options.subset_size = 3;
    options.execution_options.random_seed = 11U;

    const auto result{ trainer.TrainAlphaG(beta_group_list, options) };
    const auto first_single{ trainer.TrainAlphaG({ beta_group_list.at(0) }, options) };
    const auto second_single{ trainer.TrainAlphaG({ beta_group_list.at(1) }, options) };
    const auto expected_error_sum{ (first_single.error_sum_list + second_single.error_sum_list).eval() };

    EXPECT_TRUE(result.error_sum_list.isApprox(expected_error_sum, 1e-12));
    EXPECT_DOUBLE_EQ(result.best_alpha, BestAlphaForErrors(alpha_list, result.error_sum_list));
}

TEST(HRLAlphaTrainerTest, TrainAlphaGWithSeedIsDeterministic)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
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
    const auto & alpha_list{ trainer.AlphaGrid() };
    HRLAlphaTrainer::AlphaGTrainingOptions options;
    options.subset_size = 3;
    options.execution_options.random_seed = 7U;

    const auto first{ trainer.TrainAlphaG(beta_group_list, options) };
    const auto second{ trainer.TrainAlphaG(beta_group_list, options) };
    const auto single_first{ trainer.TrainAlphaG({ beta_group_list.at(0) }, options) };
    const auto single_second{ trainer.TrainAlphaG({ beta_group_list.at(0) }, options) };

    EXPECT_TRUE(first.error_sum_list.isApprox(second.error_sum_list, 1e-12));
    EXPECT_TRUE(single_first.error_sum_list.isApprox(single_second.error_sum_list, 1e-12));
    EXPECT_DOUBLE_EQ(first.best_alpha, second.best_alpha);
    EXPECT_DOUBLE_EQ(single_first.best_alpha, single_second.best_alpha);
    EXPECT_DOUBLE_EQ(first.best_alpha, BestAlphaForErrors(alpha_list, first.error_sum_list));
}

TEST(HRLAlphaTrainerTest, TrainAlphaRejectsInvalidInputs)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    const auto dataset{ MakeLinearDataset(1.0) };
    const std::vector<HRLMemberDataset> empty_dataset_list;
    const std::vector<HRLMemberDataset> dataset_list{ dataset };

    EXPECT_THROW(
        trainer.TrainAlphaR(empty_dataset_list),
        std::invalid_argument);

    HRLAlphaTrainer::AlphaRTrainingOptions alpha_r_options;
    alpha_r_options.subset_size = 7;
    EXPECT_THROW(
        trainer.TrainAlphaR(dataset_list, alpha_r_options),
        std::invalid_argument);

    HRLAlphaTrainer::AlphaGTrainingOptions alpha_g_options;
    alpha_g_options.subset_size = 3;
    const std::vector<std::vector<Eigen::VectorXd>> beta_group_list{
        {
            MakeVector({ 1.0, 2.0 }),
            MakeVector({ 1.5, 2.5 })
        }
    };
    EXPECT_THROW(
        trainer.TrainAlphaG(beta_group_list, alpha_g_options),
        std::invalid_argument);
}

TEST(HRLAlphaTrainerTest, StudyAlphaRBiasReturnsFiniteMatrixAndReportsProgress)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    const std::vector<HRLMemberDataset> dataset_list{
        MakeLinearDataset(2.0),
        MakeLinearDataset(-0.5)
    };
    std::size_t progress_count{ 0 };
    HRLAlphaTrainer::AlphaBiasStudyOptions options;
    options.execution_options.thread_size = 1;
    options.progress_callback =
        [&progress_count](std::size_t, std::size_t)
        {
            progress_count++;
        };

    const auto bias_matrix{ trainer.StudyAlphaRBias(dataset_list, options) };

    EXPECT_EQ(bias_matrix.rows(), 3);
    EXPECT_EQ(bias_matrix.cols(), static_cast<Eigen::Index>(trainer.AlphaGrid().size()));
    EXPECT_TRUE(bias_matrix.array().isFinite().all());
    EXPECT_EQ(progress_count, dataset_list.size());
}

TEST(HRLAlphaTrainerTest, StudyAlphaGBiasReturnsFiniteMatrixAndReportsProgress)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    const std::vector<std::vector<Eigen::VectorXd>> beta_group_list{
        {
            MakeVector({ 1.0, 2.0 }),
            MakeVector({ 1.5, 2.5 }),
            MakeVector({ 2.0, 3.0 })
        },
        {
            MakeVector({ -1.0, 0.5 }),
            MakeVector({ -1.5, 1.0 }),
            MakeVector({ -2.0, 1.5 })
        }
    };
    std::size_t progress_count{ 0 };
    HRLAlphaTrainer::AlphaBiasStudyOptions options;
    options.execution_options.thread_size = 1;
    options.progress_callback =
        [&progress_count](std::size_t, std::size_t)
        {
            progress_count++;
        };

    const auto bias_matrix{ trainer.StudyAlphaGBias(beta_group_list, options) };

    EXPECT_EQ(bias_matrix.rows(), 3);
    EXPECT_EQ(bias_matrix.cols(), static_cast<Eigen::Index>(trainer.AlphaGrid().size()));
    EXPECT_TRUE(bias_matrix.array().isFinite().all());
    EXPECT_EQ(progress_count, beta_group_list.size());
}

TEST(HRLAlphaTrainerTest, StudyAlphaBiasRejectsEmptyInputs)
{
    const HRLAlphaTrainer trainer{ 0.0, 1.0, 0.5 };
    const std::vector<HRLMemberDataset> empty_dataset_list;
    const std::vector<std::vector<Eigen::VectorXd>> empty_beta_group_list;
    const std::vector<std::vector<Eigen::VectorXd>> beta_group_with_empty_member_list{
        {}
    };

    EXPECT_THROW(trainer.StudyAlphaRBias(empty_dataset_list), std::invalid_argument);
    EXPECT_THROW(trainer.StudyAlphaGBias(empty_beta_group_list), std::invalid_argument);
    EXPECT_THROW(trainer.StudyAlphaGBias(beta_group_with_empty_member_list), std::invalid_argument);
}
