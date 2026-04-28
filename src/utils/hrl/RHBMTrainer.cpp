#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>

#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <random>
#include <string>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::rhbm_trainer
{

namespace
{
constexpr double kAlphaGridTolerance{ 1.0e-12 };

const linearization_service::LinearizationSpec & MetricDecodeSpec()
{
    static const auto spec{ linearization_service::LinearizationSpec::DefaultMetricModel() };
    return spec;
}

Eigen::VectorXd CalculateAbsoluteGaussianDifference(
    const RHBMParameterVector & linear_a,
    const RHBMParameterVector & linear_b)
{
    const auto gaussian_a{ linearization_service::DecodeParameterVector(
        MetricDecodeSpec(),
        linear_a).ToVector() };
    const auto gaussian_b{ linearization_service::DecodeParameterVector(
        MetricDecodeSpec(),
        linear_b).ToVector() };
    eigen_validation::RequireVectorSize(gaussian_a, gaussian_b.rows(), "gaussian");
    return (gaussian_a - gaussian_b).array().abs().matrix();
}

void ValidateTrainingInputs(
    std::size_t data_size,
    std::size_t subset_size,
    const std::vector<double> & alpha_list)
{
    if (data_size == 0)
    {
        throw std::invalid_argument("training data must not be empty.");
    }
    numeric_validation::RequirePositive(subset_size, "subset_size");
    if (subset_size > data_size)
    {
        throw std::invalid_argument("subset_size must not exceed data size.");
    }
    if (alpha_list.empty())
    {
        throw std::invalid_argument("alpha_list must not be empty.");
    }
}

void ValidateTrainingBatch(
    std::size_t batch_size,
    std::size_t subset_size,
    const std::vector<double> & alpha_list)
{
    if (batch_size == 0)
    {
        throw std::invalid_argument("training data must not be empty.");
    }
    numeric_validation::RequirePositive(subset_size, "subset_size");
    if (alpha_list.empty())
    {
        throw std::invalid_argument("alpha_list must not be empty.");
    }
}

void ValidateMemberDataset(const RHBMMemberDataset & dataset)
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

RHBMMemberDataset BuildDatasetSlice(
    const RHBMMemberDataset & dataset,
    const std::vector<int> & row_indices)
{
    if (row_indices.empty())
    {
        throw std::invalid_argument("dataset slice must not be empty.");
    }

    const auto row_count{ static_cast<Eigen::Index>(row_indices.size()) };
    RHBMMemberDataset slice;
    slice.X = RHBMDesignMatrix::Zero(row_count, dataset.X.cols());
    slice.y = RHBMResponseVector::Zero(row_count);
    slice.score = RHBMScoreVector::Zero(row_count);
    for (Eigen::Index i = 0; i < row_count; i++)
    {
        const auto row_index{ static_cast<Eigen::Index>(row_indices.at(static_cast<std::size_t>(i))) };
        if (row_index < 0 || row_index >= dataset.X.rows())
        {
            throw std::out_of_range("dataset slice index is out of range.");
        }
        slice.X.row(i) = dataset.X.row(row_index);
        slice.y(i) = dataset.y(row_index);
        slice.score(i) = dataset.score(row_index);
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

AlphaTrainer::AlphaTrainingResult BuildTrainingResult(
    const std::vector<double> & alpha_list,
    const Eigen::VectorXd & error_sum_list)
{
    int error_min_id{ 0 };
    error_sum_list.minCoeff(&error_min_id);

    AlphaTrainer::AlphaTrainingResult result;
    result.best_alpha = alpha_list.at(static_cast<std::size_t>(error_min_id));
    result.error_sum_list = error_sum_list;
    return result;
}

Eigen::VectorXd EvaluateAlphaRForDataset(
    const RHBMMemberDataset & dataset,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const RHBMExecutionOptions & options)
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

    std::vector<RHBMMemberDataset> data_set_test;
    std::vector<RHBMMemberDataset> data_set_training;
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
                rhbm_helper::EstimateBetaMDPDE(
                    alpha,
                    data_set_test.at(i),
                    algorithm_options
                )
            };

            const auto beta_result_training{
                rhbm_helper::EstimateBetaMDPDE(
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

Eigen::VectorXd EvaluateAlphaGForGroup(
    const std::vector<RHBMParameterVector> & beta_list,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const RHBMExecutionOptions & options)
{
    ValidateTrainingInputs(beta_list.size(), subset_size, alpha_list);

    const auto data_size_in_half{ beta_list.size() / 2 };
    std::vector<std::vector<RHBMParameterVector>> data_set_test(subset_size);
    std::vector<std::vector<RHBMParameterVector>> data_set_training(subset_size);
    for (std::size_t i = 0; i < subset_size; i++)
    {
        data_set_test[i].reserve(data_size_in_half);
        data_set_training[i].reserve(beta_list.size() - data_size_in_half);
        std::vector<RHBMParameterVector> shuffled_data{ beta_list };
        auto generator{ BuildGenerator(options.random_seed, i) };
        std::shuffle(shuffled_data.begin(), shuffled_data.end(), generator);
        const auto diff{
            static_cast<std::vector<RHBMParameterVector>::difference_type>(data_size_in_half)
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
                rhbm_helper::BuildBetaMatrix(data_set_test.at(i))
            };
            const auto mu_result_test{
                rhbm_helper::EstimateMuMDPDE(alpha, beta_matrix_test, algorithm_options)
            };

            const auto beta_matrix_training{
                rhbm_helper::BuildBetaMatrix(data_set_training.at(i))
            };
            const auto mu_result_training{
                rhbm_helper::EstimateMuMDPDE(alpha, beta_matrix_training, algorithm_options)
            };

            mu_error_sum +=
                (mu_result_test.mu_mdpde - mu_result_training.mu_mdpde).norm();
        }
        error_sum_list(p) = mu_error_sum;
    }

    return error_sum_list;
}
} // namespace

AlphaTrainer::AlphaTrainer(
    double alpha_min,
    double alpha_max,
    double alpha_step) :
    m_alpha_min{ alpha_min },
    m_alpha_max{ alpha_max },
    m_alpha_step{ alpha_step },
    m_alpha_grid{ BuildAlphaGrid(alpha_min, alpha_max, alpha_step) }
{
}

std::ostringstream AlphaTrainer::GetAlphaGridSummary() const
{
    std::ostringstream summary;
    summary
        << "Alpha training search grid: min = " << std::to_string(m_alpha_min)
        << ", max = " << std::to_string(m_alpha_max)
        << ", step = " << std::to_string(m_alpha_step)
        << ", count = " << std::to_string(m_alpha_grid.size());
    return summary;
}

std::vector<double> AlphaTrainer::BuildAlphaGrid(
    double alpha_min,
    double alpha_max,
    double alpha_step)
{
    numeric_validation::RequireFiniteNonNegativeRange(alpha_min, alpha_max, "alpha training range");
    numeric_validation::RequireFinitePositive(alpha_step, "alpha training step");

    std::vector<double> alpha_list;
    for (double alpha{ alpha_min };
         alpha <= alpha_max + kAlphaGridTolerance;
         alpha += alpha_step)
    {
        auto alpha_value{ alpha };
        if (std::abs(alpha_value - alpha_max) <= kAlphaGridTolerance)
        {
            alpha_value = alpha_max;
        }
        if (alpha_value > alpha_max)
        {
            break;
        }

        alpha_list.emplace_back(alpha_value);
        if (alpha_value == alpha_max)
        {
            break;
        }
    }
    return alpha_list;
}

AlphaTrainer::AlphaTrainingResult AlphaTrainer::TrainAlphaR(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const AlphaTrainingOptions & options) const
{
    ValidateTrainingBatch(dataset_list.size(), options.subset_size, m_alpha_grid);
    for (const auto & dataset : dataset_list)
    {
        ValidateMemberDataset(dataset);
        ValidateTrainingInputs(
            static_cast<std::size_t>(dataset.y.size()),
            options.subset_size,
            m_alpha_grid);
    }

    const auto dataset_size{ dataset_list.size() };
    std::atomic<std::size_t> completed_count{ 0 };
    Eigen::ArrayXd error_sum_array{
        Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(m_alpha_grid.size()))
    };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.execution_options.thread_size)
#endif
    for (std::size_t i = 0; i < dataset_size; i++)
    {
        const auto error_array{
            EvaluateAlphaRForDataset(
                dataset_list.at(i),
                options.subset_size,
                m_alpha_grid,
                options.execution_options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            error_sum_array += error_array.array();
            const auto completed{ ++completed_count };
            if (options.progress_callback)
            {
                options.progress_callback(completed, dataset_size);
            }
        }
    }

    return BuildTrainingResult(m_alpha_grid, error_sum_array.matrix());
}

AlphaTrainer::AlphaTrainingResult AlphaTrainer::TrainAlphaG(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    const AlphaTrainingOptions & options) const
{
    ValidateTrainingBatch(beta_group_list.size(), options.subset_size, m_alpha_grid);
    for (const auto & beta_list : beta_group_list)
    {
        ValidateTrainingInputs(beta_list.size(), options.subset_size, m_alpha_grid);
    }

    const auto group_size{ beta_group_list.size() };
    std::atomic<std::size_t> completed_count{ 0 };
    Eigen::ArrayXd error_sum_array{
        Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(m_alpha_grid.size()))
    };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.execution_options.thread_size)
#endif
    for (std::size_t i = 0; i < group_size; i++)
    {
        const auto error_array{
            EvaluateAlphaGForGroup(
                beta_group_list.at(i),
                options.subset_size,
                m_alpha_grid,
                options.execution_options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            error_sum_array += error_array.array();
            const auto completed{ ++completed_count };
            if (options.progress_callback)
            {
                options.progress_callback(completed, group_size);
            }
        }
    }

    return BuildTrainingResult(m_alpha_grid, error_sum_array.matrix());
}

Eigen::MatrixXd AlphaTrainer::StudyAlphaRBias(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const AlphaRunOptions & options) const
{
    ValidateTrainingBatch(dataset_list.size(), 1, m_alpha_grid);
    for (const auto & dataset : dataset_list)
    {
        ValidateMemberDataset(dataset);
    }

    const auto dataset_size{ dataset_list.size() };
    const auto alpha_size{ static_cast<int>(m_alpha_grid.size()) };
    std::atomic<std::size_t> completed_count{ 0 };
    Eigen::MatrixXd gaus_bias_matrix{ Eigen::MatrixXd::Zero(3, alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.execution_options.thread_size)
#endif
    for (std::size_t i = 0; i < dataset_size; i++)
    {
        auto algorithm_options{ options.execution_options };
        algorithm_options.quiet_mode = true;

        Eigen::MatrixXd local_bias_array{ Eigen::MatrixXd::Zero(3, alpha_size) };
        for (int j = 0; j < alpha_size; j++)
        {
            const auto alpha_r{ m_alpha_grid.at(static_cast<std::size_t>(j)) };
            const auto result{
                rhbm_helper::EstimateBetaMDPDE(
                    alpha_r,
                    dataset_list.at(i),
                    algorithm_options)
            };
            local_bias_array.col(j) = CalculateAbsoluteGaussianDifference(
                result.beta_mdpde,
                result.beta_ols
            );
        }

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            gaus_bias_matrix += local_bias_array;
            const auto completed{ ++completed_count };
            if (options.progress_callback)
            {
                options.progress_callback(completed, dataset_size);
            }
        }
    }

    gaus_bias_matrix /= static_cast<double>(dataset_size);
    return gaus_bias_matrix;
}

Eigen::MatrixXd AlphaTrainer::StudyAlphaGBias(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    const AlphaRunOptions & options) const
{
    ValidateTrainingBatch(beta_group_list.size(), 1, m_alpha_grid);
    for (const auto & beta_list : beta_group_list)
    {
        if (beta_list.empty())
        {
            throw std::invalid_argument("training data must not be empty.");
        }
    }

    const auto group_size{ beta_group_list.size() };
    const auto alpha_size{ static_cast<int>(m_alpha_grid.size()) };
    std::atomic<std::size_t> completed_count{ 0 };
    Eigen::MatrixXd gaus_bias_matrix{ Eigen::MatrixXd::Zero(3, alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.execution_options.thread_size)
#endif
    for (std::size_t i = 0; i < group_size; i++)
    {
        auto algorithm_options{ options.execution_options };
        algorithm_options.quiet_mode = true;
        const auto beta_matrix{ rhbm_helper::BuildBetaMatrix(beta_group_list.at(i)) };

        Eigen::MatrixXd local_bias_array{ Eigen::MatrixXd::Zero(3, alpha_size) };
        for (int j = 0; j < alpha_size; j++)
        {
            const auto alpha_g{ m_alpha_grid.at(static_cast<std::size_t>(j)) };
            const auto result{
                rhbm_helper::EstimateMuMDPDE(alpha_g, beta_matrix, algorithm_options)
            };
            local_bias_array.col(j) = CalculateAbsoluteGaussianDifference(
                result.mu_mdpde,
                result.mu_mean
            );
        }

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            gaus_bias_matrix += local_bias_array;
            const auto completed{ ++completed_count };
            if (options.progress_callback)
            {
                options.progress_callback(completed, group_size);
            }
        }
    }

    gaus_bias_matrix /= static_cast<double>(group_size);
    return gaus_bias_matrix;
}

} // namespace rhbm_gem::rhbm_trainer
