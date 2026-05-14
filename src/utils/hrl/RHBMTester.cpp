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

Eigen::VectorXd CalculateNormalizedBias(const GaussianModel3D & estimate, const GaussianModel3D & truth)
{
    GaussianModel3D::RequireFiniteModel(truth, "truth");
    const auto truth_vector{ truth.ToVector() };
    const auto estimate_vector{ estimate.ToVector() };
    eigen_validation::RequireVectorSize(estimate_vector, truth_vector.rows(), "gaussian");
    return ((estimate_vector - truth_vector).array() / truth_vector.array()).matrix();
}

BetaReplicaBias EstimateBetaReplicaBias(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & truth,
    double alpha_r,
    const gaussian_estimator::TrainingOptions & options)
{
    const auto estimate{
        gaussian_estimator::EstimateLocalGaussian(sample_entries, alpha_r, options)
    };
    return BetaReplicaBias{
        CalculateNormalizedBias(estimate.ols.GetModel(), truth),
        CalculateNormalizedBias(estimate.mdpde.GetModel(), truth)
    };
}

MuReplicaBias EstimateMuReplicaBias(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_results,
    const GaussianModel3D & truth,
    double alpha_g,
    const gaussian_estimator::TrainingOptions & options)
{
    const auto baseline_result{
        gaussian_estimator::EstimateGroupGaussian(sample_entries_list, member_results, 0.0, options)
    };
    const auto mdpde_result{
        gaussian_estimator::EstimateGroupGaussian(sample_entries_list, member_results, alpha_g, options)
    };
    return MuReplicaBias{
        CalculateNormalizedBias(baseline_result.mdpde, truth),
        CalculateNormalizedBias(mdpde_result.mdpde, truth)
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

BetaMDPDETestBias RunBetaMDPDETest(
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

    const gaussian_estimator::TrainingOptions estimator_options;

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
                gaussian_estimator::TrainAlphaR(
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

    BetaMDPDETestBias result;
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

    return result;
}

MuMDPDETestBias RunMuMDPDETest(
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

    const gaussian_estimator::TrainingOptions estimator_options;

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & sample_entries_list{
            input.replica_member_sampling_entries.at(static_cast<size_t>(i))
        };
        std::vector<LocalGaussianResult> member_results;
        member_results.reserve(sample_entries_list.size());
        for (const auto & sample_entries : sample_entries_list)
        {
            member_results.emplace_back(
                gaussian_estimator::EstimateLocalGaussian(sample_entries, 0.0, estimator_options));
        }
        const auto baseline_result{
            EstimateMuReplicaBias(
                sample_entries_list, member_results, input.gaus_true, 0.0, estimator_options)
        };
        bias_matrix_median.col(i) = baseline_result.median_bias;

        for (size_t j = 0; j < local_alpha_g_list.size(); j++)
        {
            const auto replica_result{
                EstimateMuReplicaBias(
                    sample_entries_list,
                    member_results,
                    input.gaus_true,
                    local_alpha_g_list.at(j),
                    estimator_options)
            };
            bias_matrix_mdpde_requested_list.at(j).col(i) = replica_result.mdpde_bias;
        }

        if (input.alpha_training)
        {
            const auto trained_alpha_g{
                gaussian_estimator::TrainAlphaG(
                    std::vector<std::vector<LocalPotentialSampleList>>{ sample_entries_list },
                    std::vector<std::vector<LocalGaussianResult>>{ member_results },
                    estimator_options)
            };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_g;
            const auto replica_result{
                EstimateMuReplicaBias(
                    sample_entries_list,
                    member_results,
                    input.gaus_true,
                    trained_alpha_g,
                    estimator_options)
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

    MuMDPDETestBias result;
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

    return result;
}

} // namespace rhbm_gem::rhbm_tester
