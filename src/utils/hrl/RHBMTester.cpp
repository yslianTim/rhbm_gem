#include <rhbm_gem/utils/hrl/RHBMTester.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimator.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace rhbm_gem::rhbm_tester
{

namespace
{

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

gaussian_estimator::CrossValidationOptions MakeAlphaOptions()
{
    gaussian_estimator::CrossValidationOptions options;
    options.thread_size = 1;
    return options;
}

Eigen::VectorXd CalculateNormalizedBias(
    const GaussianModel3D & gaussian_estimate,
    const GaussianModel3D & gaussian_truth)
{
    GaussianModel3D::RequireFiniteModel(gaussian_truth, "gaussian_truth");
    const auto gaussian_truth_vector{ gaussian_truth.ToVector() };
    const auto gaussian_estimate_vector{ gaussian_estimate.ToVector() };
    eigen_validation::RequireVectorSize(gaussian_estimate_vector, gaussian_truth_vector.rows(), "gaussian");
    return ((gaussian_estimate_vector - gaussian_truth_vector).array() / gaussian_truth_vector.array()).matrix();
}

BetaReplicaBias EstimateBetaReplicaBias(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & gaus_true,
    double alpha_r,
    const gaussian_estimator::CrossValidationOptions & options)
{
    const auto gaussian_result{
        gaussian_estimator::EstimateLocalGaussian(sample_entries, alpha_r, options)
    };
    return BetaReplicaBias{
        CalculateNormalizedBias(gaussian_result.ols.GetModel(), gaus_true),
        CalculateNormalizedBias(gaussian_result.mdpde.GetModel(), gaus_true)
    };
}

std::vector<double> BuildZeroAlphaRList(std::size_t member_size)
{
    return std::vector<double>(member_size, 0.0);
}

std::vector<LocalGaussianResult> EstimateLocalGaussianList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const gaussian_estimator::CrossValidationOptions & options)
{
    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        member_results.emplace_back(
            gaussian_estimator::EstimateLocalGaussian(sample_entries, 0.0, options));
    }
    return member_results;
}

MuReplicaBias EstimateMuReplicaBias(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const GaussianModel3D & gaus_true,
    double alpha_g,
    const gaussian_estimator::CrossValidationOptions & options)
{
    const auto alpha_r_list{ BuildZeroAlphaRList(sample_entries_list.size()) };
    const auto baseline_result{
        gaussian_estimator::EstimateGroupGaussian(sample_entries_list, alpha_r_list, 0.0, options)
    };
    const auto mdpde_result{
        gaussian_estimator::EstimateGroupGaussian(sample_entries_list, alpha_r_list, alpha_g, options)
    };
    return MuReplicaBias{
        CalculateNormalizedBias(baseline_result.mdpde, gaus_true),
        CalculateNormalizedBias(mdpde_result.mdpde, gaus_true)
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
    const auto replica_size{ static_cast<int>(input.replica_sampling_entries.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_sampling_entries must not be empty");
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

    const auto estimator_options{ MakeAlphaOptions() };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & sample_entries{ input.replica_sampling_entries.at(static_cast<size_t>(i)) };
        const auto baseline_result{
            EstimateBetaReplicaBias(sample_entries, input.gaus_true, 0.0, estimator_options)
        };
        bias_matrix_ols.col(i) = baseline_result.ols_bias;

        for (size_t j = 0; j < local_alpha_r_list.size(); j++)
        {
            const auto replica_result{
                EstimateBetaReplicaBias(
                    sample_entries,
                    input.gaus_true,
                    local_alpha_r_list.at(j),
                    estimator_options)
            };
            bias_matrix_mdpde_requested_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (input.alpha_training)
        {
            const auto trained_alpha_r{
                gaussian_estimator::CrossValidationAlphaR(
                    std::vector<LocalPotentialSampleList>{ sample_entries }, estimator_options)
            };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_r;
            const auto replica_result{
                EstimateBetaReplicaBias(
                    sample_entries, input.gaus_true, trained_alpha_r, estimator_options)
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
    const auto replica_size{ static_cast<int>(input.replica_member_sampling_entries.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_member_sampling_entries must not be empty");
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

    const auto estimator_options{ MakeAlphaOptions() };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & sample_entries_list{
            input.replica_member_sampling_entries.at(static_cast<size_t>(i))
        };
        const auto baseline_result{
            EstimateMuReplicaBias(sample_entries_list, input.gaus_true, 0.0, estimator_options)
        };
        bias_matrix_median.col(i) = baseline_result.median_bias;

        for (size_t j = 0; j < local_alpha_g_list.size(); j++)
        {
            const auto replica_result{
                EstimateMuReplicaBias(
                    sample_entries_list,
                    input.gaus_true,
                    local_alpha_g_list.at(j),
                    estimator_options)
            };
            bias_matrix_mdpde_requested_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (input.alpha_training)
        {
            const auto member_results{
                EstimateLocalGaussianList(sample_entries_list, estimator_options)
            };
            const auto trained_alpha_g{
                gaussian_estimator::CrossValidationAlphaG(
                    std::vector<std::vector<LocalPotentialSampleList>>{ sample_entries_list },
                    std::vector<std::vector<LocalGaussianResult>>{ member_results },
                    estimator_options)
            };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_g;
            const auto replica_result{
                EstimateMuReplicaBias(
                    sample_entries_list, input.gaus_true, trained_alpha_g, estimator_options)
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
