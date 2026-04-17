#include <rhbm_gem/utils/hrl/HRLModelTester.hpp>
#include <rhbm_gem/utils/hrl/HRLAlphaTrainer.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>
#include <rhbm_gem/utils/hrl/HRLModelMetrics.hpp>

#include <cmath>
#include <stdexcept>

namespace
{
constexpr std::size_t kAlphaRSubsetSize{ 5 };
constexpr double kAlphaRMin{ 0.0 };
constexpr double kAlphaRMax{ 2.0 };
constexpr double kAlphaRStep{ 0.1 };

constexpr std::size_t kAlphaGSubsetSize{ 10 };
constexpr double kAlphaGMin{ 0.0 };
constexpr double kAlphaGMax{ 1.0 };
constexpr double kAlphaGStep{ 0.1 };

HRLExecutionOptions MakeTesterExecutionOptions()
{
    HRLExecutionOptions options;
    options.quiet_mode = true;
    options.thread_size = 1;
    return options;
}

Eigen::VectorXd CalculateReplicaSigma(
    const Eigen::MatrixXd & residual_matrix,
    const Eigen::VectorXd & residual_mean)
{
    if (residual_matrix.cols() <= 1)
    {
        return Eigen::VectorXd::Zero(residual_matrix.rows());
    }

    return (residual_matrix.colwise() - residual_mean).rowwise().norm()
        / std::sqrt(static_cast<double>(residual_matrix.cols() - 1));
}
} // namespace

HRLModelTester::HRLModelTester(int gaus_par_size) :
    m_gaus_par_size{ gaus_par_size }
{
    if (m_gaus_par_size <= 0)
    {
        throw std::invalid_argument("gaus_par_size must be positive value");
    }
}

bool HRLModelTester::RunBetaMDPDETest(
    const std::vector<double> & alpha_r_list,
    std::vector<Eigen::VectorXd> & residual_mean_ols_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    const HRLBetaTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    CheckGausParametersDimension(test_input.gaus_true);
    const auto replica_size{ static_cast<int>(test_input.replica_datasets.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.replica_datasets must not be empty");
    }

    const auto local_alpha_r_list{ alpha_r_list };
    const auto alpha_size{ local_alpha_r_list.size() + 1 }; // add one for training alpha_r
    std::vector<Eigen::MatrixXd> residual_matrix_ols_list(alpha_size);
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_ols_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, replica_size)
    );
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, replica_size)
    );
    residual_mean_ols_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_mean_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_ols_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));

    const HRLAlphaTrainer alpha_r_trainer{ kAlphaRMin, kAlphaRMax, kAlphaRStep };
    HRLAlphaTrainer::AlphaTrainingOptions alpha_r_training_options;
    alpha_r_training_options.subset_size = kAlphaRSubsetSize;
    alpha_r_training_options.execution_options = MakeTesterExecutionOptions();

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & dataset{ test_input.replica_datasets.at(static_cast<size_t>(i)) };

        const auto alpha_r_training_result{
            alpha_r_trainer.TrainAlphaR(
                std::vector<HRLMemberDataset>{ dataset }, alpha_r_training_options)
        };
        const auto trained_alpha_r{ alpha_r_training_result.best_alpha };
        const auto options{ MakeTesterExecutionOptions() };

        for (size_t j = 0; j < alpha_size; j++)
        {
            const auto alpha{ (j < local_alpha_r_list.size()) ?
                local_alpha_r_list.at(j) : trained_alpha_r
            };
            const auto beta_result{ HRLModelAlgorithms::EstimateBetaMDPDE(alpha, dataset, options) };

            residual_matrix_ols_list.at(j).col(i) =
                HRLModelMetrics::CalculateNormalizedResidual(
                    beta_result.beta_ols,
                    test_input.gaus_true
                );
            residual_matrix_mdpde_list.at(j).col(i) =
                HRLModelMetrics::CalculateNormalizedResidual(
                    beta_result.beta_mdpde,
                    test_input.gaus_true
                );
        }
    }

    for (size_t j = 0; j < alpha_size; j++)
    {
        residual_mean_ols_list.at(j) = residual_matrix_ols_list.at(j).rowwise().mean();
        residual_mean_mdpde_list.at(j) = residual_matrix_mdpde_list.at(j).rowwise().mean();
        residual_sigma_ols_list.at(j) = CalculateReplicaSigma(
            residual_matrix_ols_list.at(j),
            residual_mean_ols_list.at(j)
        );
        residual_sigma_mdpde_list.at(j) = CalculateReplicaSigma(
            residual_matrix_mdpde_list.at(j),
            residual_mean_mdpde_list.at(j)
        );
    }

    return true;
}

bool HRLModelTester::RunMuMDPDETest(
    const std::vector<double> & alpha_g_list,
    std::vector<Eigen::VectorXd> & residual_mean_median_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_median_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    const HRLMuTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    CheckGausParametersDimension(test_input.gaus_true);
    const auto replica_size{ static_cast<int>(test_input.replica_beta_matrices.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.replica_beta_matrices must not be empty");
    }

    const auto local_alpha_g_list{ alpha_g_list };
    const auto alpha_size{ local_alpha_g_list.size() + 1 }; // add one for training alpha_g
    std::vector<Eigen::MatrixXd> residual_matrix_median_list(alpha_size);
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_median_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, replica_size)
    );
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, replica_size)
    );
    residual_mean_median_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_mean_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_median_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));

    const HRLAlphaTrainer alpha_g_trainer{ kAlphaGMin, kAlphaGMax, kAlphaGStep };
    HRLAlphaTrainer::AlphaTrainingOptions alpha_g_training_options;
    alpha_g_training_options.subset_size = kAlphaGSubsetSize;
    alpha_g_training_options.execution_options = MakeTesterExecutionOptions();

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & beta_matrix{ test_input.replica_beta_matrices.at(static_cast<size_t>(i)) };
        std::vector<Eigen::VectorXd> train_data_entry_list;
        train_data_entry_list.reserve(static_cast<size_t>(beta_matrix.cols()));
        for (int m = 0; m < beta_matrix.cols(); m++)
        {
            train_data_entry_list.emplace_back(beta_matrix.col(m));
        }

        const auto alpha_g_training_result{
            alpha_g_trainer.TrainAlphaG(
                std::vector<std::vector<Eigen::VectorXd>>{ train_data_entry_list },
                alpha_g_training_options
            )
        };
        const auto trained_alpha_g{ alpha_g_training_result.best_alpha };
        const auto options{ MakeTesterExecutionOptions() };

        for (size_t j = 0; j < alpha_size; j++)
        {
            const auto alpha{
                j < local_alpha_g_list.size()
                    ? local_alpha_g_list.at(j)
                    : trained_alpha_g
            };
            const auto mdpde_result = HRLModelAlgorithms::EstimateMuMDPDE(
                alpha,
                beta_matrix,
                options
            );
            const auto ols_result = HRLModelAlgorithms::EstimateMuMDPDE(
                0.0,
                beta_matrix,
                options
            );

            residual_matrix_median_list.at(j).col(i) =
                HRLModelMetrics::CalculateNormalizedResidual(
                    ols_result.mu_mdpde,
                    test_input.gaus_true
                );
            residual_matrix_mdpde_list.at(j).col(i) =
                HRLModelMetrics::CalculateNormalizedResidual(
                    mdpde_result.mu_mdpde,
                    test_input.gaus_true
                );
        }
    }

    for (size_t j = 0; j < alpha_size; j++)
    {
        residual_mean_median_list.at(j) = residual_matrix_median_list.at(j).rowwise().mean();
        residual_mean_mdpde_list.at(j) = residual_matrix_mdpde_list.at(j).rowwise().mean();
        residual_sigma_median_list.at(j) = CalculateReplicaSigma(
            residual_matrix_median_list.at(j),
            residual_mean_median_list.at(j)
        );
        residual_sigma_mdpde_list.at(j) = CalculateReplicaSigma(
            residual_matrix_mdpde_list.at(j),
            residual_mean_mdpde_list.at(j)
        );
    }

    return true;
}

bool HRLModelTester::RunBetaMDPDEWithNeighborhoodTest(
    std::vector<Eigen::VectorXd> & residual_mean_list,
    std::vector<Eigen::VectorXd> & residual_sigma_list,
    const HRLNeighborhoodTestInput & test_input,
    double & training_alpha_r_average,
    int thread_size,
    double angle)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif
    (void)angle;

    CheckGausParametersDimension(test_input.gaus_true);
    const auto replica_size{ static_cast<int>(test_input.no_cut_datasets.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.no_cut_datasets must not be empty");
    }
    if (test_input.cut_datasets.size() != test_input.no_cut_datasets.size())
    {
        throw std::invalid_argument("test_input cut/no_cut dataset sizes must match");
    }

    const size_t method_size{ 3 }; // OLS (no cut), MDPDE (no cut and cut)
    std::vector<Eigen::MatrixXd> replica_residual_list(method_size);
    replica_residual_list.assign(method_size, Eigen::MatrixXd::Zero(m_gaus_par_size, replica_size));
    residual_mean_list.assign(method_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_list.assign(method_size, Eigen::VectorXd::Zero(m_gaus_par_size));

    const HRLAlphaTrainer alpha_r_trainer{ kAlphaRMin, kAlphaRMax, kAlphaRStep };
    HRLAlphaTrainer::AlphaTrainingOptions alpha_r_training_options;
    alpha_r_training_options.subset_size = kAlphaRSubsetSize;
    alpha_r_training_options.execution_options = MakeTesterExecutionOptions();

    training_alpha_r_average = 0.0;

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & no_cut_dataset{ test_input.no_cut_datasets.at(static_cast<size_t>(i)) };
        const auto & cut_dataset{ test_input.cut_datasets.at(static_cast<size_t>(i)) };

        const auto no_cut_training_result{
            alpha_r_trainer.TrainAlphaR(
                std::vector<HRLMemberDataset>{ no_cut_dataset }, alpha_r_training_options)
        };
        const auto cut_training_result{
            alpha_r_trainer.TrainAlphaR(
                std::vector<HRLMemberDataset>{ cut_dataset }, alpha_r_training_options)
        };
        const auto no_cut_alpha_r_train{ no_cut_training_result.best_alpha };
        const auto cut_alpha_r_train{ cut_training_result.best_alpha };
        const auto options{ MakeTesterExecutionOptions() };

        auto no_cut_result{
            HRLModelAlgorithms::EstimateBetaMDPDE(no_cut_alpha_r_train, no_cut_dataset, options)
        };
        auto cut_result{
            HRLModelAlgorithms::EstimateBetaMDPDE(cut_alpha_r_train, cut_dataset, options)
        };

        replica_residual_list.at(0).col(i) =
            HRLModelMetrics::CalculateNormalizedResidual(
                no_cut_result.beta_ols,
                test_input.gaus_true
            );
        replica_residual_list.at(1).col(i) =
            HRLModelMetrics::CalculateNormalizedResidual(
                no_cut_result.beta_mdpde,
                test_input.gaus_true
            );
        replica_residual_list.at(2).col(i) =
            HRLModelMetrics::CalculateNormalizedResidual(
                cut_result.beta_mdpde,
                test_input.gaus_true
            );

        #pragma omp critical
        {
            training_alpha_r_average += cut_alpha_r_train;
        }
    }
    training_alpha_r_average /= replica_size;

    for (size_t j = 0; j < method_size; j++)
    {
        residual_mean_list.at(j) = replica_residual_list.at(j).rowwise().mean();
        residual_sigma_list.at(j) = CalculateReplicaSigma(
            replica_residual_list.at(j),
            residual_mean_list.at(j)
        );
    }

    return true;
}

bool HRLModelTester::CheckGausParametersDimension(const Eigen::VectorXd & gaus_par)
{
    if (gaus_par.rows() != m_gaus_par_size)
    {
        throw std::invalid_argument("model parameters size invalid, must be : " +
            std::to_string(m_gaus_par_size));
    }
    return true;
}
