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

struct BetaReplicaResidual
{
    Eigen::VectorXd ols_residual;
    Eigen::VectorXd mdpde_residual;
};

struct MuReplicaResidual
{
    Eigen::VectorXd median_residual;
    Eigen::VectorXd mdpde_residual;
};

RHBMExecutionOptions MakeTesterExecutionOptions()
{
    RHBMExecutionOptions options;
    options.quiet_mode = true;
    options.thread_size = 1;
    return options;
}

void ValidateGaussianTruthVector(
    const Eigen::VectorXd & gaussian_truth,
    const char * value_name)
{
    eigen_validation::RequireVectorSize(
        gaussian_truth,
        GaussianModel3D::kParameterSize,
        value_name);
}

GaussianParameterVector CalculateNormalizedResidual(
    const RHBMBetaVector & linear_estimate,
    const GaussianParameterVector & gaussian_truth)
{
    ValidateGaussianTruthVector(gaussian_truth, "gaussian_truth");
    const auto gaussian_estimate{
        linearization_service::DecodeGroupBeta(
            linearization_service::LinearizationSpec::DefaultMetricModel(), linear_estimate)
    };
    eigen_validation::RequireVectorSize(gaussian_estimate, gaussian_truth.rows(), "gaussian");
    return ((gaussian_estimate - gaussian_truth).array() / gaussian_truth.array()).matrix();
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

BetaReplicaResidual EstimateBetaReplicaResidual(
    const RHBMMemberDataset & dataset,
    const GaussianParameterVector & gaus_true,
    double alpha_r,
    const RHBMExecutionOptions & options)
{
    const auto beta_result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, options) };
    return BetaReplicaResidual{
        CalculateNormalizedResidual(beta_result.beta_ols, gaus_true),
        CalculateNormalizedResidual(beta_result.beta_mdpde, gaus_true)
    };
}

MuReplicaResidual EstimateMuReplicaResidual(
    const RHBMBetaMatrix & beta_matrix,
    const GaussianParameterVector & gaus_true,
    double alpha_g,
    const RHBMExecutionOptions & options)
{
    const auto mdpde_result{ rhbm_helper::EstimateMuMDPDE(alpha_g, beta_matrix, options) };
    const auto ols_result{ rhbm_helper::EstimateMuMDPDE(0.0, beta_matrix, options) };
    return MuReplicaResidual{
        CalculateNormalizedResidual(ols_result.mu_mdpde, gaus_true),
        CalculateNormalizedResidual(mdpde_result.mu_mdpde, gaus_true)
    };
}

ResidualStatistics FinalizeResidualStatistics(const Eigen::MatrixXd & residual_matrix)
{
    ResidualStatistics result;
    result.mean = residual_matrix.rowwise().mean();
    result.sigma = CalculateReplicaSigma(residual_matrix, result.mean);
    return result;
}

void FinalizeResidualStatisticsSeries(
    const std::vector<Eigen::MatrixXd> & residual_matrix_list,
    size_t requested_alpha_size,
    bool alpha_training,
    double trained_alpha_average,
    ResidualStatisticsSeries & result)
{
    result.requested_alpha.assign(requested_alpha_size, ResidualStatistics{});
    for (size_t i = 0; i < requested_alpha_size; i++)
    {
        result.requested_alpha.at(i) = FinalizeResidualStatistics(residual_matrix_list.at(i));
    }
    result.trained_alpha.reset();
    if (alpha_training)
    {
        result.trained_alpha = FinalizeResidualStatistics(residual_matrix_list.at(requested_alpha_size));
        result.trained_alpha_average = trained_alpha_average;
    }
    else
    {
        result.trained_alpha_average.reset();
    }
}
} // namespace

bool RunBetaMDPDETest(
    BetaMDPDETestResidual & result,
    const test_data_factory::RHBMBetaTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    ValidateGaussianTruthVector(test_input.gaus_true, "test_input.gaus_true");
    const auto replica_size{ static_cast<int>(test_input.replica_datasets.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.replica_datasets must not be empty");
    }

    const auto local_alpha_r_list{ test_input.requested_alpha_r_list };
    const bool alpha_training{ test_input.alpha_training };
    const auto alpha_size{ local_alpha_r_list.size() + (alpha_training ? 1 : 0) };
    Eigen::MatrixXd residual_matrix_ols{
        Eigen::MatrixXd::Zero(GaussianModel3D::kParameterSize, replica_size)
    };
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(GaussianModel3D::kParameterSize, replica_size)
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
            EstimateBetaReplicaResidual(dataset, test_input.gaus_true, 0.0, options)
        };
        residual_matrix_ols.col(i) = baseline_result.ols_residual;

        for (size_t j = 0; j < local_alpha_r_list.size(); j++)
        {
            const auto replica_result{
                EstimateBetaReplicaResidual(
                    dataset,
                    test_input.gaus_true,
                    local_alpha_r_list.at(j),
                    options)
            };
            residual_matrix_mdpde_list.at(j).col(i) = replica_result.mdpde_residual;
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
                EstimateBetaReplicaResidual(dataset, test_input.gaus_true, trained_alpha_r, options)
            };
            residual_matrix_mdpde_list.at(local_alpha_r_list.size()).col(i) =
                replica_result.mdpde_residual;
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

    result = BetaMDPDETestResidual{};
    result.ols = FinalizeResidualStatistics(residual_matrix_ols);
    FinalizeResidualStatisticsSeries(
        residual_matrix_mdpde_list,
        local_alpha_r_list.size(),
        alpha_training,
        trained_alpha_average,
        result.mdpde);

    return true;
}

bool RunMuMDPDETest(
    MuMDPDETestResidual & result,
    const test_data_factory::RHBMMuTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    ValidateGaussianTruthVector(test_input.gaus_true, "test_input.gaus_true");
    const auto replica_size{ static_cast<int>(test_input.replica_beta_matrices.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("test_input.replica_beta_matrices must not be empty");
    }

    const auto local_alpha_g_list{ test_input.requested_alpha_g_list };
    const bool alpha_training{ test_input.alpha_training };
    const auto alpha_size{ local_alpha_g_list.size() + (alpha_training ? 1 : 0) };
    Eigen::MatrixXd residual_matrix_median{
        Eigen::MatrixXd::Zero(GaussianModel3D::kParameterSize, replica_size)
    };
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(GaussianModel3D::kParameterSize, replica_size)
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
            EstimateMuReplicaResidual(beta_matrix, test_input.gaus_true, 0.0, options)
        };
        residual_matrix_median.col(i) = baseline_result.median_residual;

        for (size_t j = 0; j < local_alpha_g_list.size(); j++)
        {
            const auto replica_result{
                EstimateMuReplicaResidual(
                    beta_matrix,
                    test_input.gaus_true,
                    local_alpha_g_list.at(j),
                    options)
            };
            residual_matrix_mdpde_list.at(j).col(i) = replica_result.mdpde_residual;
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
                EstimateMuReplicaResidual(beta_matrix, test_input.gaus_true, trained_alpha_g, options)
            };
            residual_matrix_mdpde_list.at(local_alpha_g_list.size()).col(i) =
                replica_result.mdpde_residual;
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

    result = MuMDPDETestResidual{};
    result.median = FinalizeResidualStatistics(residual_matrix_median);
    FinalizeResidualStatisticsSeries(
        residual_matrix_mdpde_list,
        local_alpha_g_list.size(),
        alpha_training,
        trained_alpha_average,
        result.mdpde);

    return true;
}

} // namespace rhbm_gem::rhbm_tester
