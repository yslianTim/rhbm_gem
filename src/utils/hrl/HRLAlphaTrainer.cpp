#include <rhbm_gem/utils/hrl/HRLAlphaTrainer.hpp>

#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>

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

void ValidateMemberDataset(const HRLMemberDataset & dataset)
{
    if (dataset.X.rows() != dataset.y.size())
    {
        throw std::invalid_argument("dataset shape is inconsistent.");
    }
    if (dataset.X.rows() == 0 || dataset.X.cols() == 0)
    {
        throw std::invalid_argument("dataset must not be empty.");
    }
}

HRLMemberDataset BuildDatasetSlice(
    const HRLMemberDataset & dataset,
    const std::vector<int> & row_indices)
{
    if (row_indices.empty())
    {
        throw std::invalid_argument("dataset slice must not be empty.");
    }

    const auto row_count{ static_cast<Eigen::Index>(row_indices.size()) };
    HRLMemberDataset slice;
    slice.X = Eigen::MatrixXd::Zero(row_count, dataset.X.cols());
    slice.y = Eigen::VectorXd::Zero(row_count);
    for (Eigen::Index i = 0; i < row_count; i++)
    {
        const auto row_index{ static_cast<Eigen::Index>(row_indices.at(static_cast<std::size_t>(i))) };
        if (row_index < 0 || row_index >= dataset.X.rows())
        {
            throw std::out_of_range("dataset slice index is out of range.");
        }
        slice.X.row(i) = dataset.X.row(row_index);
        slice.y(i) = dataset.y(row_index);
    }
    return slice;
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
    const HRLMemberDataset & dataset,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const HRLExecutionOptions & options)
{
    ValidateMemberDataset(dataset);
    ValidateTrainingInputs(static_cast<std::size_t>(dataset.y.size()), subset_size, alpha_list);

    std::vector<std::vector<int>> data_subset_rows(subset_size);
    const auto total_entry_size{ static_cast<std::size_t>(dataset.y.size()) };
    const auto entries_in_subset_size{ total_entry_size / subset_size + 1 };
    for (std::size_t i = 0; i < subset_size; i++)
    {
        data_subset_rows[i].reserve(entries_in_subset_size);
    }

    for (std::size_t row = 0; row < total_entry_size; row++)
    {
        data_subset_rows[row % subset_size].emplace_back(static_cast<int>(row));
    }

    std::vector<HRLMemberDataset> data_set_test;
    std::vector<HRLMemberDataset> data_set_training;
    data_set_test.reserve(subset_size);
    data_set_training.reserve(subset_size);
    for (std::size_t i = 0; i < subset_size; i++)
    {
        const auto test_set_size{ data_subset_rows[i].size() };
        std::vector<int> training_rows;
        training_rows.reserve(total_entry_size - test_set_size);
        for (std::size_t j = 0; j < subset_size; j++)
        {
            if (i == j)
            {
                continue;
            }
            training_rows.insert(
                training_rows.end(),
                data_subset_rows[j].begin(),
                data_subset_rows[j].end()
            );
        }
        data_set_test.emplace_back(BuildDatasetSlice(dataset, data_subset_rows[i]));
        data_set_training.emplace_back(BuildDatasetSlice(dataset, training_rows));
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
            const auto beta_result_test{
                HRLModelAlgorithms::EstimateBetaMDPDE(
                    alpha,
                    data_set_test.at(i),
                    algorithm_options
                )
            };

            const auto beta_result_training{
                HRLModelAlgorithms::EstimateBetaMDPDE(
                    alpha,
                    data_set_training.at(i),
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
