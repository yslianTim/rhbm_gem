#include <rhbm_gem/core/EstimatorTester.hpp>
#include "core/detail/GaussianEstimatorStages.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace rhbm_gem::core {

namespace {

struct LocalReplicaBias
{
    Eigen::VectorXd ols_bias;
    Eigen::VectorXd mdpde_bias;
};

struct GroupReplicaBias
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

    Eigen::VectorXd bias{ estimate_vector - truth_vector };
    for (Eigen::Index i = 0; i < truth_vector.rows(); i++)
    {
        if (truth_vector(i) != 0.0)
        {
            bias(i) /= truth_vector(i);
        }
    }
    return bias;
}

LocalReplicaBias EstimateLocalReplicaBias(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & truth,
    double alpha_r,
    const FitOptions & options)
{
    const auto estimate{ EstimateLocalGaussian(sample_entries, alpha_r, options) };
    return LocalReplicaBias{
        CalculateNormalizedBias(estimate.ols.GetModel(), truth),
        CalculateNormalizedBias(estimate.mdpde.GetModel(), truth)
    };
}

GroupReplicaBias EstimateGroupReplicaBias(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_results,
    const GaussianModel3D & truth,
    double alpha_g,
    const FitOptions & options)
{
    const auto estimate{
        EstimateGroupGaussian(sample_entries_list, member_results, alpha_g, options)
    };
    return GroupReplicaBias{
        CalculateNormalizedBias(estimate.mean, truth),
        CalculateNormalizedBias(estimate.mdpde, truth)
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

Eigen::MatrixXd EstimateAtomicModelFirstStageModels(
    const AtomicModelTestData & input,
    const FitOptions & options)
{
    const auto replica_size{ static_cast<int>(input.replica_model_objects.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_model_objects must not be empty");
    }

    Eigen::MatrixXd estimation_matrix{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };

    for (int i = 0; i < replica_size; i++)
    {
        ModelObject model_object{ *input.replica_model_objects.at(static_cast<size_t>(i)) };
        RunLocalAlphaTraining(model_object, options, LocalFittingPass::FirstStage);
        RunFixedOffsetLocalFitting(model_object, options, LocalFittingPass::FirstStage);

        const auto local_view{
            AtomLocalPotentialView::RequireFor(*model_object.GetSelectedAtoms().front())
        };
        const auto & gaussian_result{ local_view.GetGaussianResult() };
        estimation_matrix.col(i) = gaussian_result.mdpde.GetModel().ToVector();
    }

    return estimation_matrix;
}

} // namespace

LocalTestBias RunLocalEstimationTest(
    const LocalTestData & input,
    const LocalTestOptions & options)
{
    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto replica_size{ static_cast<int>(input.replica_sampling_entries.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_sampling_entries must not be empty");
    }

    Eigen::MatrixXd bias_matrix_ols{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    Eigen::MatrixXd bias_matrix_mdpde_requested{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    Eigen::MatrixXd bias_matrix_mdpde_trained{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    std::vector<double> trained_alpha_list;
    trained_alpha_list.assign(static_cast<size_t>(replica_size), 0.0);

    FitOptions estimator_options;
    estimator_options.quiet_mode = options.quiet_mode;

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
#endif
    for (int i = 0; i < replica_size; i++)
    {
        const auto & sample_entries{ input.replica_sampling_entries.at(static_cast<size_t>(i)) };
        const auto requested_result{
            EstimateLocalReplicaBias(
                sample_entries,
                input.gaus_true,
                options.requested_alpha_r,
                estimator_options)
        };
        bias_matrix_ols.col(i) = requested_result.ols_bias;
        bias_matrix_mdpde_requested.col(i) = requested_result.mdpde_bias;

        if (options.alpha_training)
        {
            const auto trained_alpha_r{
                TrainAlphaR(
                    std::vector<LocalPotentialSampleList>{ sample_entries }, estimator_options)
            };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_r;
            const auto replica_result{
                EstimateLocalReplicaBias(
                    sample_entries, input.gaus_true, trained_alpha_r, estimator_options)
            };
            bias_matrix_mdpde_trained.col(i) = replica_result.mdpde_bias;
        }
    }

    double trained_alpha_median{ 0.0 };
    if (options.alpha_training)
    {
        trained_alpha_median = array_helper::ComputeMedian(trained_alpha_list);
    }

    LocalTestBias result;
    result.ols = FinalizeBiasStatistics(bias_matrix_ols);
    result.mdpde.requested_alpha = FinalizeBiasStatistics(bias_matrix_mdpde_requested);
    result.mdpde.trained_alpha.reset();
    result.mdpde.trained_alpha_median.reset();
    if (options.alpha_training)
    {
        result.mdpde.trained_alpha = FinalizeBiasStatistics(bias_matrix_mdpde_trained);
        result.mdpde.trained_alpha_median = trained_alpha_median;
    }

    return result;
}

GroupTestBias RunGroupEstimationTest(
    const GroupTestData & input,
    const GroupTestOptions & options)
{
    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto replica_size{ static_cast<int>(input.replica_member_sampling_entries.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_member_sampling_entries must not be empty");
    }

    Eigen::MatrixXd bias_matrix_median{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    Eigen::MatrixXd bias_matrix_mdpde_requested{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    Eigen::MatrixXd bias_matrix_mdpde_trained{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    std::vector<double> trained_alpha_list;
    trained_alpha_list.assign(static_cast<size_t>(replica_size), 0.0);

    FitOptions estimator_options;
    estimator_options.quiet_mode = options.quiet_mode;

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
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
                EstimateLocalGaussian(sample_entries, 0.0, estimator_options));
        }

        const auto requested_result{
            EstimateGroupReplicaBias(
                sample_entries_list,
                member_results,
                input.gaus_true,
                options.requested_alpha_g,
                estimator_options)
        };
        bias_matrix_median.col(i) = requested_result.median_bias;
        bias_matrix_mdpde_requested.col(i) = requested_result.mdpde_bias;

        if (options.alpha_training)
        {
            const auto trained_alpha_g{
                TrainAlphaG(
                    std::vector<std::vector<LocalGaussianResult>>{ member_results },
                    estimator_options)
            };
            trained_alpha_list.at(static_cast<size_t>(i)) = trained_alpha_g;
            const auto replica_result{
                EstimateGroupReplicaBias(
                    sample_entries_list,
                    member_results,
                    input.gaus_true,
                    trained_alpha_g,
                    estimator_options)
            };
            bias_matrix_mdpde_trained.col(i) = replica_result.mdpde_bias;
        }
    }

    double trained_alpha_median{ 0.0 };
    if (options.alpha_training)
    {
        trained_alpha_median = array_helper::ComputeMedian(trained_alpha_list);
    }

    GroupTestBias result;
    result.median = FinalizeBiasStatistics(bias_matrix_median);
    result.mdpde.requested_alpha = FinalizeBiasStatistics(bias_matrix_mdpde_requested);
    result.mdpde.trained_alpha.reset();
    result.mdpde.trained_alpha_median.reset();
    if (options.alpha_training)
    {
        result.mdpde.trained_alpha = FinalizeBiasStatistics(bias_matrix_mdpde_trained);
        result.mdpde.trained_alpha_median = trained_alpha_median;
    }

    return result;
}

BiasStatistics RunAtomicModelFirstStageEstimationTest(
    const AtomicModelTestData & input,
    const FitOptions & options)
{
    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto estimation_matrix{ EstimateAtomicModelFirstStageModels(input, options) };
    const auto replica_size{ static_cast<int>(estimation_matrix.cols()) };
    Eigen::MatrixXd bias_matrix{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };

    for (int i = 0; i < replica_size; i++)
    {
        bias_matrix.col(i) =
            CalculateNormalizedBias(
                GaussianModel3D::FromVector(estimation_matrix.col(i)),
                input.gaus_true);
    }

    return FinalizeBiasStatistics(bias_matrix);
}

LocalGaussianEstimateBias RunAtomicModelLocalEstimationTest(
    const AtomicModelTestData & input,
    double alpha_r,
    const FitOptions & options)
{
    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto replica_size{ static_cast<int>(input.replica_model_objects.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_model_objects must not be empty");
    }

    Eigen::MatrixXd ols_bias{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };
    Eigen::MatrixXd mdpde_bias{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };

    for (int i = 0; i < replica_size; i++)
    {
        const auto & model_object{ *input.replica_model_objects.at(static_cast<size_t>(i)) };
        const auto local_view{
            AtomLocalPotentialView::RequireFor(*model_object.GetSelectedAtoms().front())
        };
        const auto estimate{
            EstimateLocalGaussian(
                local_view.GetSamplingEntries(false), alpha_r, options)
        };
        ols_bias.col(i) = CalculateNormalizedBias(estimate.ols.GetModel(), input.gaus_true);
        mdpde_bias.col(i) = CalculateNormalizedBias(estimate.mdpde.GetModel(), input.gaus_true);
    }

    return LocalGaussianEstimateBias{
        FinalizeBiasStatistics(ols_bias),
        FinalizeBiasStatistics(mdpde_bias)
    };
}

GaussianModel3D EstimateAtomicModelFirstStageMean(
    const AtomicModelTestData & input,
    const FitOptions & options)
{
    const auto estimation_matrix{ EstimateAtomicModelFirstStageModels(input, options) };
    return GaussianModel3D::FromVector(estimation_matrix.rowwise().mean());
}

BiasStatistics RunAtomicModelFullEstimationTest(
    const AtomicModelTestData & input,
    const FitOptions & options)
{
    GaussianModel3D::RequireFiniteModel(input.gaus_true, "input.gaus_true");
    const auto replica_size{ static_cast<int>(input.replica_model_objects.size()) };
    if (replica_size <= 0)
    {
        throw std::invalid_argument("input.replica_model_objects must not be empty");
    }

    Eigen::MatrixXd bias_matrix{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), replica_size)
    };

    for (int i = 0; i < replica_size; i++)
    {
        ModelObject model_object{ *input.replica_model_objects.at(static_cast<size_t>(i)) };
        RunPotentialFittingWorkflow(model_object, options);

        const auto local_view{
            AtomLocalPotentialView::RequireFor(*model_object.GetSelectedAtoms().front())
        };
        const auto & gaussian_result{ local_view.GetGaussianResult() };
        bias_matrix.col(i) =
            CalculateNormalizedBias(gaussian_result.mdpde.GetModel(), input.gaus_true);
    }

    return FinalizeBiasStatistics(bias_matrix);
}

} // namespace rhbm_gem::core
