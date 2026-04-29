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

Eigen::VectorXd CalculateNormalizedBias(
    const RHBMParameterVector & linear_estimate,
    const GaussianModel3D & gaussian_truth)
{
    GaussianModel3D::RequireFiniteModel(gaussian_truth, "gaussian_truth");
    const auto gaussian_truth_vector{ gaussian_truth.ToVector() };
    const auto gaussian_estimate{
        linearization_service::DecodeParameterVector(
            linearization_service::LinearizationSpec::AtomDecode(), linear_estimate).ToVector()
    };
    eigen_validation::RequireVectorSize(gaussian_estimate, gaussian_truth_vector.rows(), "gaussian");
    return ((gaussian_estimate - gaussian_truth_vector).array() / gaussian_truth_vector.array()).matrix();
}

BetaReplicaBias EstimateBetaReplicaBias(
    const RHBMMemberDataset & dataset,
    const GaussianModel3D & gaus_true,
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
    const GaussianModel3D & gaus_true,
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
    const auto replica_size{ bias_matrix.cols() };
    BiasStatistics result;
    result.mean = bias_matrix.rowwise().mean();
    if (replica_size <= 1)
    {
        result.sigma = Eigen::VectorXd::Zero(bias_matrix.rows());
        return result;
    }
    result.sigma = (bias_matrix.colwise() - result.mean).rowwise().norm()
        / std::sqrt(static_cast<double>(replica_size - 1));
    return result;
}

} // namespace

bool RunBetaMDPDETest(
    BetaMDPDETestBias & result,
    const test_data_factory::RHBMBetaTestInput & input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto replica_size{ static_cast<int>(input.replica_datasets.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_datasets must not be empty");
    }

    const auto local_alpha_r_list{ input.requested_alpha_r_list };
    Eigen::MatrixXd bias_matrix_ols{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    std::vector<Eigen::MatrixXd> bias_matrix_mdpde_requested_list(
        local_alpha_r_list.size(),
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    );
    Eigen::MatrixXd bias_matrix_mdpde_trained;
    std::vector<double> trained_alpha_list;
    if (input.alpha_training)
    {
        bias_matrix_mdpde_trained = Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size);
        trained_alpha_list.assign(static_cast<size_t>(replica_size), 0.0);
    }

    const rhbm_trainer::AlphaTrainer alpha_r_trainer{ kAlphaRMin, kAlphaRMax, kAlphaRStep };
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions trainer_options;
    trainer_options.subset_size = kAlphaRSubsetSize;
    trainer_options.execution_options = RHBMExecutionOptions{ true };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & dataset{ input.replica_datasets.at(static_cast<size_t>(i)) };
        const auto baseline_result{
            EstimateBetaReplicaBias(dataset, input.gaus_true, 0.0, trainer_options.execution_options)
        };
        bias_matrix_ols.col(i) = baseline_result.ols_bias;

        for (size_t j = 0; j < local_alpha_r_list.size(); j++)
        {
            const auto replica_result{
                EstimateBetaReplicaBias(
                    dataset,
                    input.gaus_true,
                    local_alpha_r_list.at(j),
                    trainer_options.execution_options)
            };
            bias_matrix_mdpde_requested_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (input.alpha_training)
        {
            const auto alpha_r_training_result{
                alpha_r_trainer.TrainAlphaR(
                    std::vector<RHBMMemberDataset>{ dataset }, trainer_options)
            };
            const auto trained_alpha_r{ alpha_r_training_result.best_alpha };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_r;
            const auto replica_result{
                EstimateBetaReplicaBias(dataset, input.gaus_true, trained_alpha_r, trainer_options.execution_options)
            };
            bias_matrix_mdpde_trained.col(i) = replica_result.mdpde_bias;
        }
    }

    double trained_alpha_average{ 0.0 };
    if (input.alpha_training)
    {
        for (const auto trained_alpha : trained_alpha_list)
        {
            trained_alpha_average += trained_alpha;
        }
        trained_alpha_average /= static_cast<double>(replica_size);
    }

    result = BetaMDPDETestBias{};
    result.ols = FinalizeBiasStatistics(bias_matrix_ols);
    result.mdpde.requested_alpha.assign(local_alpha_r_list.size(), BiasStatistics{});
    for (size_t i = 0; i < local_alpha_r_list.size(); i++)
    {
        result.mdpde.requested_alpha.at(i) =
            FinalizeBiasStatistics(bias_matrix_mdpde_requested_list.at(i));
    }
    result.mdpde.trained_alpha.reset();
    result.mdpde.trained_alpha_average.reset();
    if (input.alpha_training)
    {
        result.mdpde.trained_alpha = FinalizeBiasStatistics(bias_matrix_mdpde_trained);
        result.mdpde.trained_alpha_average = trained_alpha_average;
    }

    return true;
}

bool RunMuMDPDETest(
    MuMDPDETestBias & result,
    const test_data_factory::RHBMMuTestInput & input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto replica_size{ static_cast<int>(input.replica_beta_matrices.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_beta_matrices must not be empty");
    }

    const auto local_alpha_g_list{ input.requested_alpha_g_list };
    Eigen::MatrixXd bias_matrix_median{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    std::vector<Eigen::MatrixXd> bias_matrix_mdpde_requested_list(
        local_alpha_g_list.size(),
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    );
    Eigen::MatrixXd bias_matrix_mdpde_trained;
    std::vector<double> trained_alpha_list;
    if (input.alpha_training)
    {
        bias_matrix_mdpde_trained =
            Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size);
        trained_alpha_list.assign(static_cast<size_t>(replica_size), 0.0);
    }

    const rhbm_trainer::AlphaTrainer alpha_g_trainer{ kAlphaGMin, kAlphaGMax, kAlphaGStep };
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions trainer_options;
    trainer_options.subset_size = kAlphaGSubsetSize;
    trainer_options.execution_options = RHBMExecutionOptions{ true };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & beta_matrix{ input.replica_beta_matrices.at(static_cast<size_t>(i)) };
        const auto baseline_result{
            EstimateMuReplicaBias(beta_matrix, input.gaus_true, 0.0, trainer_options.execution_options)
        };
        bias_matrix_median.col(i) = baseline_result.median_bias;

        for (size_t j = 0; j < local_alpha_g_list.size(); j++)
        {
            const auto replica_result{
                EstimateMuReplicaBias(
                    beta_matrix,
                    input.gaus_true,
                    local_alpha_g_list.at(j),
                    trainer_options.execution_options)
            };
            bias_matrix_mdpde_requested_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (input.alpha_training)
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
                    trainer_options)
            };
            const auto trained_alpha_g{ alpha_g_training_result.best_alpha };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_g;
            const auto replica_result{
                EstimateMuReplicaBias(beta_matrix, input.gaus_true, trained_alpha_g, trainer_options.execution_options)
            };
            bias_matrix_mdpde_trained.col(i) = replica_result.mdpde_bias;
        }
    }

    double trained_alpha_average{ 0.0 };
    if (input.alpha_training)
    {
        for (const auto trained_alpha : trained_alpha_list)
        {
            trained_alpha_average += trained_alpha;
        }
        trained_alpha_average /= static_cast<double>(replica_size);
    }

    result = MuMDPDETestBias{};
    result.median = FinalizeBiasStatistics(bias_matrix_median);
    result.mdpde.requested_alpha.assign(local_alpha_g_list.size(), BiasStatistics{});
    for (size_t i = 0; i < local_alpha_g_list.size(); i++)
    {
        result.mdpde.requested_alpha.at(i) =
            FinalizeBiasStatistics(bias_matrix_mdpde_requested_list.at(i));
    }
    result.mdpde.trained_alpha.reset();
    result.mdpde.trained_alpha_average.reset();
    if (input.alpha_training)
    {
        result.mdpde.trained_alpha = FinalizeBiasStatistics(bias_matrix_mdpde_trained);
        result.mdpde.trained_alpha_average = trained_alpha_average;
    }

    return true;
}

} // namespace rhbm_gem::rhbm_tester
