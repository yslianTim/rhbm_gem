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

struct NeighborhoodReplicaResidual
{
    Eigen::VectorXd no_cut_ols_residual;
    Eigen::VectorXd no_cut_mdpde_residual;
    Eigen::VectorXd cut_mdpde_residual;
    double trained_alpha_r{ 0.0 };
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

NeighborhoodReplicaResidual EstimateNeighborhoodReplicaResidual(
    const RHBMMemberDataset & no_cut_dataset,
    const RHBMMemberDataset & cut_dataset,
    const GaussianParameterVector & gaus_true,
    const rhbm_trainer::AlphaTrainer & alpha_r_trainer,
    const rhbm_trainer::AlphaTrainer::AlphaTrainingOptions & alpha_r_training_options,
    const RHBMExecutionOptions & options)
{
    const auto no_cut_training_result{
        alpha_r_trainer.TrainAlphaR(
            std::vector<RHBMMemberDataset>{ no_cut_dataset },
            alpha_r_training_options)
    };
    const auto cut_training_result{
        alpha_r_trainer.TrainAlphaR(
            std::vector<RHBMMemberDataset>{ cut_dataset },
            alpha_r_training_options)
    };
    const auto no_cut_alpha_r_train{ no_cut_training_result.best_alpha };
    const auto cut_alpha_r_train{ cut_training_result.best_alpha };
    const auto no_cut_result{
        EstimateBetaReplicaResidual(
            no_cut_dataset,
            gaus_true,
            no_cut_alpha_r_train,
            options)
    };
    const auto cut_result{
        EstimateBetaReplicaResidual(
            cut_dataset,
            gaus_true,
            cut_alpha_r_train,
            options)
    };

    NeighborhoodReplicaResidual result;
    result.no_cut_ols_residual = no_cut_result.ols_residual;
    result.no_cut_mdpde_residual = no_cut_result.mdpde_residual;
    result.cut_mdpde_residual = cut_result.mdpde_residual;
    result.trained_alpha_r = cut_alpha_r_train;
    return result;
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
            const auto replica_result{
                EstimateBetaReplicaResidual(dataset, test_input.gaus_true, trained_alpha_r, options)
            };
            residual_matrix_mdpde_list.at(local_alpha_r_list.size()).col(i) =
                replica_result.mdpde_residual;
        }
    }

    result = BetaMDPDETestResidual{};
    result.ols = FinalizeResidualStatistics(residual_matrix_ols);
    FinalizeResidualStatisticsSeries(
        residual_matrix_mdpde_list,
        local_alpha_r_list.size(),
        alpha_training,
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
            const auto replica_result{
                EstimateMuReplicaResidual(beta_matrix, test_input.gaus_true, trained_alpha_g, options)
            };
            residual_matrix_mdpde_list.at(local_alpha_g_list.size()).col(i) =
                replica_result.mdpde_residual;
        }
    }

    result = MuMDPDETestResidual{};
    result.median = FinalizeResidualStatistics(residual_matrix_median);
    FinalizeResidualStatisticsSeries(
        residual_matrix_mdpde_list,
        local_alpha_g_list.size(),
        alpha_training,
        result.mdpde);

    return true;
}

bool RunBetaMDPDEWithNeighborhoodTest(
    NeighborhoodMDPDETestResidual & result,
    const test_data_factory::RHBMNeighborhoodTestInput & test_input,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    ValidateGaussianTruthVector(test_input.gaus_true, "test_input.gaus_true");
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
    replica_residual_list.assign(
        method_size, Eigen::MatrixXd::Zero(GaussianModel3D::kParameterSize, replica_size));

    const rhbm_trainer::AlphaTrainer alpha_r_trainer{ kAlphaRMin, kAlphaRMax, kAlphaRStep };
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_r_training_options;
    alpha_r_training_options.subset_size = kAlphaRSubsetSize;
    alpha_r_training_options.execution_options = MakeTesterExecutionOptions();

    result = NeighborhoodMDPDETestResidual{};

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & no_cut_dataset{ test_input.no_cut_datasets.at(static_cast<size_t>(i)) };
        const auto & cut_dataset{ test_input.cut_datasets.at(static_cast<size_t>(i)) };
        const auto options{ MakeTesterExecutionOptions() };
        const auto replica_result{
            EstimateNeighborhoodReplicaResidual(
                no_cut_dataset,
                cut_dataset,
                test_input.gaus_true,
                alpha_r_trainer,
                alpha_r_training_options,
                options)
        };

        replica_residual_list.at(0).col(i) = replica_result.no_cut_ols_residual;
        replica_residual_list.at(1).col(i) = replica_result.no_cut_mdpde_residual;
        replica_residual_list.at(2).col(i) = replica_result.cut_mdpde_residual;

        #pragma omp critical
        {
            result.trained_alpha_r_average += replica_result.trained_alpha_r;
        }
    }
    result.trained_alpha_r_average /= replica_size;

    result.no_cut_ols = FinalizeResidualStatistics(replica_residual_list.at(0));
    result.no_cut_mdpde = FinalizeResidualStatistics(replica_residual_list.at(1));
    result.cut_mdpde = FinalizeResidualStatistics(replica_residual_list.at(2));

    return true;
}

} // namespace rhbm_gem::rhbm_tester
