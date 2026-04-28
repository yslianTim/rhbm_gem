#include <rhbm_gem/utils/hrl/RHBMTester.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace rhbm_gem::rhbm_tester
{

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

struct BetaReplicaBias
{
    Eigen::VectorXd ols_bias;
    Eigen::VectorXd mdpde_bias;
};

struct MuReplicaBias
{
    Eigen::VectorXd median_bias;
    Eigen::VectorXd mdpde_bias;
};

RHBMExecutionOptions MakeTesterExecutionOptions()
{
    RHBMExecutionOptions options;
    options.quiet_mode = true;
    options.thread_size = 1;
    return options;
}

Eigen::VectorXd CalculateNormalizedBias(
    const RHBMBetaVector & linear_estimate,
    const Eigen::VectorXd & gaussian_truth)
{
    GaussianModel3D::RequireParameterVector(gaussian_truth, "gaussian_truth");
    const auto gaussian_estimate{
        linearization_service::DecodeGroupBeta(
            linearization_service::LinearizationSpec::DefaultMetricModel(), linear_estimate)
    };
    eigen_validation::RequireVectorSize(gaussian_estimate, gaussian_truth.rows(), "gaussian");
    return ((gaussian_estimate - gaussian_truth).array() / gaussian_truth.array()).matrix();
}

Eigen::VectorXd CalculateReplicaBiasSigma(
    const Eigen::MatrixXd & bias_matrix,
    const Eigen::VectorXd & bias_mean)
{
    if (bias_matrix.cols() <= 1)
    {
        return Eigen::VectorXd::Zero(bias_matrix.rows());
    }

    return (bias_matrix.colwise() - bias_mean).rowwise().norm()
        / std::sqrt(static_cast<double>(bias_matrix.cols() - 1));
}

BetaReplicaBias EstimateBetaReplicaBias(
    const RHBMMemberDataset & dataset,
    const Eigen::VectorXd & gaus_true,
    double alpha_r,
    const RHBMExecutionOptions & options)
{
    const auto beta_result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, options) };
    return BetaReplicaBias{
        CalculateNormalizedBias(beta_result.beta_ols, gaus_true),
        CalculateNormalizedBias(beta_result.beta_mdpde, gaus_true)
    };
}

MuReplicaBias EstimateMuReplicaBias(
    const RHBMBetaMatrix & beta_matrix,
    const Eigen::VectorXd & gaus_true,
    double alpha_g,
    const RHBMExecutionOptions & options)
{
    const auto mdpde_result{ rhbm_helper::EstimateMuMDPDE(alpha_g, beta_matrix, options) };
    const auto ols_result{ rhbm_helper::EstimateMuMDPDE(0.0, beta_matrix, options) };
    return MuReplicaBias{
        CalculateNormalizedBias(ols_result.mu_mdpde, gaus_true),
        CalculateNormalizedBias(mdpde_result.mu_mdpde, gaus_true)
    };
}

BiasStatistics FinalizeBiasStatistics(const Eigen::MatrixXd & bias_matrix)
{
    BiasStatistics result;
    result.mean = bias_matrix.rowwise().mean();
    result.sigma = CalculateReplicaBiasSigma(bias_matrix, result.mean);
    return result;
}

void FinalizeBiasStatisticsSeries(
    const std::vector<Eigen::MatrixXd> & bias_matrix_list,
    size_t requested_alpha_size,
    bool alpha_training,
    double trained_alpha_average,
    BiasStatisticsSeries & result)
{
    result.requested_alpha.assign(requested_alpha_size, BiasStatistics{});
    for (size_t i = 0; i < requested_alpha_size; i++)
    {
        result.requested_alpha.at(i) = FinalizeBiasStatistics(bias_matrix_list.at(i));
    }
    result.trained_alpha.reset();
    if (alpha_training)
    {
        result.trained_alpha = FinalizeBiasStatistics(bias_matrix_list.at(requested_alpha_size));
        result.trained_alpha_average = trained_alpha_average;
    }
    else
    {
        result.trained_alpha_average.reset();
    }
}
} // namespace

bool RunBetaMDPDETest(
    BetaMDPDETestBias & result,
    const test_data_factory::RHBMBetaTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    GaussianModel3D::RequireParameterVector(test_input.gaus_true, "test_input.gaus_true");
    const auto replica_size{ static_cast<int>(test_input.replica_datasets.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.replica_datasets must not be empty");
    }

    const auto local_alpha_r_list{ test_input.requested_alpha_r_list };
    const bool alpha_training{ test_input.alpha_training };
    const auto alpha_size{ local_alpha_r_list.size() + (alpha_training ? 1 : 0) };
    Eigen::MatrixXd bias_matrix_ols{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    std::vector<Eigen::MatrixXd> bias_matrix_mdpde_list(alpha_size);
    bias_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    );
    std::vector<double> trained_alpha_list(static_cast<size_t>(replica_size), 0.0);

    const rhbm_trainer::AlphaTrainer alpha_r_trainer{ kAlphaRMin, kAlphaRMax, kAlphaRStep };
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_r_training_options;
    alpha_r_training_options.subset_size = kAlphaRSubsetSize;
    alpha_r_training_options.execution_options = MakeTesterExecutionOptions();

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & dataset{ test_input.replica_datasets.at(static_cast<size_t>(i)) };

        const auto options{ MakeTesterExecutionOptions() };
        const auto baseline_result{
            EstimateBetaReplicaBias(dataset, test_input.gaus_true, 0.0, options)
        };
        bias_matrix_ols.col(i) = baseline_result.ols_bias;

        for (size_t j = 0; j < local_alpha_r_list.size(); j++)
        {
            const auto replica_result{
                EstimateBetaReplicaBias(
                    dataset,
                    test_input.gaus_true,
                    local_alpha_r_list.at(j),
                    options)
            };
            bias_matrix_mdpde_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (alpha_training)
        {
            const auto alpha_r_training_result{
                alpha_r_trainer.TrainAlphaR(
                    std::vector<RHBMMemberDataset>{ dataset }, alpha_r_training_options)
            };
            const auto trained_alpha_r{ alpha_r_training_result.best_alpha };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_r;
            const auto replica_result{
                EstimateBetaReplicaBias(dataset, test_input.gaus_true, trained_alpha_r, options)
            };
            bias_matrix_mdpde_list.at(local_alpha_r_list.size()).col(i) =
                replica_result.mdpde_bias;
        }
    }

    double trained_alpha_average{ 0.0 };
    if (alpha_training)
    {
        for (const auto trained_alpha : trained_alpha_list)
        {
            trained_alpha_average += trained_alpha;
        }
        trained_alpha_average /= static_cast<double>(replica_size);
    }

    result = BetaMDPDETestBias{};
    result.ols = FinalizeBiasStatistics(bias_matrix_ols);
    FinalizeBiasStatisticsSeries(
        bias_matrix_mdpde_list,
        local_alpha_r_list.size(),
        alpha_training,
        trained_alpha_average,
        result.mdpde);

    return true;
}

bool RunMuMDPDETest(
    MuMDPDETestBias & result,
    const test_data_factory::RHBMMuTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    GaussianModel3D::RequireParameterVector(test_input.gaus_true, "test_input.gaus_true");
    const auto replica_size{ static_cast<int>(test_input.replica_beta_matrices.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.replica_beta_matrices must not be empty");
    }

    const auto local_alpha_g_list{ test_input.requested_alpha_g_list };
    const bool alpha_training{ test_input.alpha_training };
    const auto alpha_size{ local_alpha_g_list.size() + (alpha_training ? 1 : 0) };
    Eigen::MatrixXd bias_matrix_median{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    std::vector<Eigen::MatrixXd> bias_matrix_mdpde_list(alpha_size);
    bias_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    );
    std::vector<double> trained_alpha_list(static_cast<size_t>(replica_size), 0.0);

    const rhbm_trainer::AlphaTrainer alpha_g_trainer{ kAlphaGMin, kAlphaGMax, kAlphaGStep };
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_g_training_options;
    alpha_g_training_options.subset_size = kAlphaGSubsetSize;
    alpha_g_training_options.execution_options = MakeTesterExecutionOptions();

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & beta_matrix{ test_input.replica_beta_matrices.at(static_cast<size_t>(i)) };
        const auto options{ MakeTesterExecutionOptions() };
        const auto baseline_result{
            EstimateMuReplicaBias(beta_matrix, test_input.gaus_true, 0.0, options)
        };
        bias_matrix_median.col(i) = baseline_result.median_bias;

        for (size_t j = 0; j < local_alpha_g_list.size(); j++)
        {
            const auto replica_result{
                EstimateMuReplicaBias(
                    beta_matrix,
                    test_input.gaus_true,
                    local_alpha_g_list.at(j),
                    options)
            };
            bias_matrix_mdpde_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (alpha_training)
        {
            std::vector<Eigen::VectorXd> train_data_entry_list;
            train_data_entry_list.reserve(static_cast<size_t>(beta_matrix.cols()));
            for (int m = 0; m < beta_matrix.cols(); m++)
            {
                train_data_entry_list.emplace_back(beta_matrix.col(m));
            }

            const auto alpha_g_training_result{
                alpha_g_trainer.TrainAlphaG(
                    std::vector<std::vector<Eigen::VectorXd>>{ train_data_entry_list },
                    alpha_g_training_options)
            };
            const auto trained_alpha_g{ alpha_g_training_result.best_alpha };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_g;
            const auto replica_result{
                EstimateMuReplicaBias(beta_matrix, test_input.gaus_true, trained_alpha_g, options)
            };
            bias_matrix_mdpde_list.at(local_alpha_g_list.size()).col(i) =
                replica_result.mdpde_bias;
        }
    }

    double trained_alpha_average{ 0.0 };
    if (alpha_training)
    {
        for (const auto trained_alpha : trained_alpha_list)
        {
            trained_alpha_average += trained_alpha;
        }
        trained_alpha_average /= static_cast<double>(replica_size);
    }

    result = MuMDPDETestBias{};
    result.median = FinalizeBiasStatistics(bias_matrix_median);
    FinalizeBiasStatisticsSeries(
        bias_matrix_mdpde_list,
        local_alpha_g_list.size(),
        alpha_training,
        trained_alpha_average,
        result.mdpde);

    return true;
}

} // namespace rhbm_gem::rhbm_tester
