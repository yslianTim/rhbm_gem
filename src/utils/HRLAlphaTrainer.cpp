#include "HRLAlphaTrainer.hpp"

#include "HRLDataTransform.hpp"
#include "HRLModelAlgorithms.hpp"

#include <algorithm>
#include <random>
#include <stdexcept>
#include <vector>

namespace
{
void ValidateTrainingInputs(
    std::size_t data_size,
    std::size_t subset_size,
    const std::vector<double> & alpha_list)
{
    if (data_size == 0)
    {
        throw std::invalid_argument("training data must not be empty.");
    }
    if (subset_size == 0)
    {
        throw std::invalid_argument("subset_size must be greater than zero.");
    }
    if (subset_size > data_size)
    {
        throw std::invalid_argument("subset_size must not exceed data size.");
    }
    if (alpha_list.empty())
    {
        throw std::invalid_argument("alpha_list must not be empty.");
    }
}

std::mt19937 BuildGenerator(
    const std::optional<std::uint32_t> & seed,
    std::size_t offset)
{
    if (seed.has_value())
    {
        return std::mt19937(static_cast<std::uint32_t>(*seed + offset));
    }
    return std::mt19937(std::random_device{}());
}
} // namespace

Eigen::VectorXd HRLAlphaTrainer::EvaluateAlphaR(
    const std::vector<Eigen::VectorXd> & data_list,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const HRLExecutionOptions & options)
{
    ValidateTrainingInputs(data_list.size(), subset_size, alpha_list);

    std::vector<std::vector<Eigen::VectorXd>> data_subset_list(subset_size);
    const auto total_entry_size{ data_list.size() };
    const auto entries_in_subset_size{ total_entry_size / subset_size + 1 };
    for (std::size_t i = 0; i < subset_size; i++)
    {
        data_subset_list[i].reserve(entries_in_subset_size);
    }

    std::size_t count{ 0 };
    for (const auto & entry : data_list)
    {
        data_subset_list[count % subset_size].emplace_back(entry);
        count++;
    }

    std::vector<std::vector<Eigen::VectorXd>> data_set_test(subset_size);
    std::vector<std::vector<Eigen::VectorXd>> data_set_training(subset_size);
    for (std::size_t i = 0; i < subset_size; i++)
    {
        const auto test_set_size{ data_subset_list[i].size() };
        const auto training_set_size{ data_list.size() - test_set_size };
        data_set_test[i].reserve(test_set_size);
        data_set_training[i].reserve(training_set_size);
        data_set_test[i].insert(
            data_set_test[i].end(),
            data_subset_list[i].begin(),
            data_subset_list[i].end()
        );
        for (std::size_t j = 0; j < subset_size; j++)
        {
            if (i == j)
            {
                continue;
            }
            data_set_training[i].insert(
                data_set_training[i].end(),
                data_subset_list[j].begin(),
                data_subset_list[j].end()
            );
        }
    }

    Eigen::VectorXd error_sum_list{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(alpha_list.size()))
    };
    auto algorithm_options{ options };
    algorithm_options.quiet_mode = true;
    for (int p = 0; p < error_sum_list.size(); p++)
    {
        const auto alpha{ alpha_list.at(static_cast<std::size_t>(p)) };
        auto beta_error_sum{ 0.0 };
        for (std::size_t i = 0; i < subset_size; i++)
        {
            const auto dataset_test{
                HRLDataTransform::BuildMemberDataset(data_set_test.at(i), true)
            };
            const auto beta_result_test{
                HRLModelAlgorithms::EstimateBetaMDPDE(
                    alpha,
                    dataset_test.X,
                    dataset_test.y,
                    algorithm_options
                )
            };

            const auto dataset_training{
                HRLDataTransform::BuildMemberDataset(data_set_training.at(i), true)
            };
            const auto beta_result_training{
                HRLModelAlgorithms::EstimateBetaMDPDE(
                    alpha,
                    dataset_training.X,
                    dataset_training.y,
                    algorithm_options
                )
            };

            beta_error_sum +=
                (beta_result_test.beta_mdpde - beta_result_training.beta_mdpde).norm();
        }
        error_sum_list(p) = beta_error_sum;
    }
    return error_sum_list;
}

Eigen::VectorXd HRLAlphaTrainer::EvaluateAlphaG(
    const std::vector<Eigen::VectorXd> & beta_list,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const HRLExecutionOptions & options)
{
    ValidateTrainingInputs(beta_list.size(), subset_size, alpha_list);

    const auto data_size_in_half{ beta_list.size() / 2 };
    std::vector<std::vector<Eigen::VectorXd>> data_set_test(subset_size);
    std::vector<std::vector<Eigen::VectorXd>> data_set_training(subset_size);
    for (std::size_t i = 0; i < subset_size; i++)
    {
        data_set_test[i].reserve(data_size_in_half);
        data_set_training[i].reserve(beta_list.size() - data_size_in_half);
        std::vector<Eigen::VectorXd> shuffled_data{ beta_list };
        auto generator{ BuildGenerator(options.random_seed, i) };
        std::shuffle(shuffled_data.begin(), shuffled_data.end(), generator);
        const auto diff{
            static_cast<std::vector<Eigen::VectorXd>::difference_type>(data_size_in_half)
        };
        data_set_test[i].assign(shuffled_data.begin(), shuffled_data.begin() + diff);
        data_set_training[i].assign(shuffled_data.begin() + diff, shuffled_data.end());
    }

    Eigen::VectorXd error_sum_list{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(alpha_list.size()))
    };
    auto algorithm_options{ options };
    algorithm_options.quiet_mode = true;
    for (int p = 0; p < error_sum_list.size(); p++)
    {
        const auto alpha{ alpha_list.at(static_cast<std::size_t>(p)) };
        auto mu_error_sum{ 0.0 };
        for (std::size_t i = 0; i < subset_size; i++)
        {
            const auto beta_matrix_test{
                HRLDataTransform::BuildBetaMatrix(data_set_test.at(i), true)
            };
            const auto mu_result_test{
                HRLModelAlgorithms::EstimateMuMDPDE(alpha, beta_matrix_test, algorithm_options)
            };

            const auto beta_matrix_training{
                HRLDataTransform::BuildBetaMatrix(data_set_training.at(i), true)
            };
            const auto mu_result_training{
                HRLModelAlgorithms::EstimateMuMDPDE(alpha, beta_matrix_training, algorithm_options)
            };

            mu_error_sum +=
                (mu_result_test.mu_mdpde - mu_result_training.mu_mdpde).norm();
        }
        error_sum_list(p) = mu_error_sum;
    }

    return error_sum_list;
}
