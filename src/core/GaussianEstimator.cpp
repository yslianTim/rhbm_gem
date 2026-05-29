#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <atomic>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kAlphaRSubsetSize{ 5 };
constexpr std::size_t kAlphaGSubsetSize{ 10 };

struct AlphaStudyOptions
{
    RHBMExecutionOptions execution_options{};
    rhbm_trainer::ProgressCallback progress_callback{};
};

RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = true;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

double EstimateInitialIntercept(const LocalPotentialSampleList & sample_entries)
{
    float maximum_distance{ 0.0f };
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance > maximum_distance) maximum_distance = sample.point.distance;
    }

    std::vector<double> response_list;
    response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance != maximum_distance) continue;
        const auto response{ static_cast<double>(sample.response) };
        numeric_validation::RequireFinite(response, "maximum distance shell response");
        response_list.emplace_back(response);
    }
    if (response_list.empty()) return 0.0;

    const auto lowest_response_list{
        array_helper::ComputeSmallestProportionValues(response_list, 0.10)
    };
    return array_helper::ComputeMedian(lowest_response_list);
}

double EstimateResidualIntercept(const LocalPotentialSampleList & residual_sample_entries)
{
    std::vector<double> response_list;
    response_list.reserve(residual_sample_entries.size());
    for (const auto & sample : residual_sample_entries)
    {
        response_list.emplace_back(static_cast<double>(sample.response));
    }
    if (response_list.empty()) return 0.0;
    return array_helper::ComputeMean(response_list.data(), response_list.size());
}

LocalPotentialSampleList BuildShiftedSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    double intercept)
{
    LocalPotentialSampleList shifted_sample_entries;
    shifted_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        shifted_sample_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(static_cast<double>(sample.response) - intercept),
                sample.point
            }
        );
    }
    return shifted_sample_entries;
}

LocalPotentialSampleList BuildResidualSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    const RHBMParameterVector & beta,
    LocalGaussianFitModel fit_model,
    double range_min,
    double range_max)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(beta, fit_model) };
    LocalPotentialSampleList residual_sample_entries;
    residual_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ sample.point.distance };
        if (distance < static_cast<float>(range_min)) continue;
        if (distance > static_cast<float>(range_max)) continue;

        const auto residual{
            static_cast<double>(sample.response) -
                signal_model.ResponseAtDistance(static_cast<double>(distance))
        };
        residual_sample_entries.emplace_back(
            LocalPotentialSample{ static_cast<float>(residual), sample.point }
        );
    }
    return residual_sample_entries;
}

bool EmitTrainingReportIfRequested(
    const Eigen::MatrixXd & gaus_bias_matrix,
    const std::vector<double> & alpha_list,
    std::string_view x_label,
    std::string_view y_label,
    const std::filesystem::path & training_report_dir,
    std::string_view report_file_name)
{
    if (training_report_dir.empty())
    {
        Logger::Log(LogLevel::Debug,
            "Alpha training report was skipped because no study plot directory was provided.");
        return false;
    }

    const auto report_path{ training_report_dir / std::filesystem::path{ report_file_name } };
    if (gaus_bias_matrix.rows() < 2 ||
        gaus_bias_matrix.cols() != static_cast<Eigen::Index>(alpha_list.size()))
    {
        Logger::Log(LogLevel::Warning,
            "Skip training report '" + report_path.string() +
                "' because the alpha bias matrix shape does not match the alpha grid.");
        return false;
    }

    std::vector<double> amplitude_values;
    std::vector<double> width_values;
    std::vector<double> intercept_values;
    amplitude_values.reserve(alpha_list.size());
    width_values.reserve(alpha_list.size());
    intercept_values.reserve(alpha_list.size());
    for (Eigen::Index column = 0; column < gaus_bias_matrix.cols(); ++column)
    {
        amplitude_values.emplace_back(gaus_bias_matrix(0, column));
        width_values.emplace_back(gaus_bias_matrix(1, column));
        if (gaus_bias_matrix.rows() >= GaussianModel3D::ParameterSize())
        {
            intercept_values.emplace_back(gaus_bias_matrix(2, column));
        }
    }

    LinePlotRequest report_request;
    report_request.output_path = report_path;
    report_request.x_axis.title = std::string(x_label);
    report_request.shared_y_axis_title = std::string(y_label);
    report_request.panels = {
        LinePlotPanel{
            "Amplitude #font[2]{A}",
            AxisSpec{},
            { LineSeries{ "Amplitude", alpha_list, amplitude_values, std::nullopt } }
        },
        LinePlotPanel{
            "Width #tau",
            AxisSpec{},
            { LineSeries{ "Width", alpha_list, width_values, std::nullopt } }
        }
    };
    if (gaus_bias_matrix.rows() >= GaussianModel3D::ParameterSize())
    {
        report_request.panels.emplace_back(LinePlotPanel{
            "Intercept #font[2]{C}",
            AxisSpec{},
            { LineSeries{ "Intercept", alpha_list, intercept_values, std::nullopt } }
        });
    }

    const auto plot_result{ local_painter::SaveLinePlot(report_request) };
    if (!plot_result.Succeeded())
    {
        Logger::Log(LogLevel::Warning,
            "Failed to emit training report '" + report_path.string() + "': " + plot_result.message);
        return false;
    }

    if (!std::filesystem::exists(report_path))
    {
        Logger::Log(LogLevel::Warning,
            "Training report output was requested but no file was produced: " + report_path.string());
        return false;
    }

    return true;
}

bool ShouldOutputStudyPlot(
    bool output_study_plot,
    const FitOptions & options)
{
    if (!output_study_plot) return false;
    if (!options.study_plot_dir.empty()) return true;

    Logger::Log(LogLevel::Debug,
        "Alpha training report was skipped because no study plot directory was provided.");
    return false;
}

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(
    std::size_t subset_size,
    const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.subset_size = subset_size;
    training_options.execution_options = MakeExecutionOptions(options);
    if (options.output_progress)
    {
        training_options.progress_callback =
            [](std::size_t completed, std::size_t total)
            {
                Logger::ProgressPercent(completed, total);
            };
    }
    return training_options;
}

AlphaStudyOptions MakeStudyOptions(
    const FitOptions & options)
{
    AlphaStudyOptions study_options;
    study_options.execution_options = MakeExecutionOptions(options);
    if (options.output_progress)
    {
        study_options.progress_callback =
            [](std::size_t completed, std::size_t total)
            {
                Logger::ProgressPercent(completed, total);
            };
    }
    return study_options;
}

Eigen::VectorXd CalculateAbsoluteGaussianDifference(
    const RHBMParameterVector & linear_a,
    const RHBMParameterVector & linear_b,
    LocalGaussianFitModel fit_model)
{
    const auto gaussian_a{ linearization_service::DecodeParameterVector(linear_a, fit_model).ToVector() };
    const auto gaussian_b{ linearization_service::DecodeParameterVector(linear_b, fit_model).ToVector() };
    eigen_validation::RequireSameSize(gaussian_a, gaussian_b, "gaussian");
    return (gaussian_a - gaussian_b).array().abs().matrix();
}

void ValidateStudyBatch(
    std::size_t batch_size,
    const std::vector<double> & alpha_list)
{
    numeric_validation::RequirePositive(batch_size, "training data size");
    numeric_validation::RequirePositive(alpha_list.size(), "alpha_list");
    numeric_validation::RequireAllFinite(alpha_list, "alpha_list");
}

Eigen::MatrixXd StudyAlphaRBias(
    const std::vector<double> & alpha_grid,
    const std::vector<RHBMMemberDataset> & dataset_list,
    LocalGaussianFitModel fit_model,
    const AlphaStudyOptions & options)
{
    ValidateStudyBatch(dataset_list.size(), alpha_grid);
    for (const auto & dataset : dataset_list)
    {
        eigen_validation::RequireVectorSize(
            dataset.y, dataset.X.rows(), "dataset.y", "dataset shape is inconsistent.");
        eigen_validation::RequireNonEmpty(dataset.X, "dataset.X");
        numeric_validation::RequirePositive(dataset.X.cols(), "dataset.X column count");
        eigen_validation::RequireFinite(dataset.X, "dataset.X");
        eigen_validation::RequireFinite(dataset.y, "dataset.y");
    }

    const auto dataset_size{ dataset_list.size() };
    const auto alpha_size{ static_cast<int>(alpha_grid.size()) };
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
            const auto alpha_r{ alpha_grid.at(static_cast<std::size_t>(j)) };
            const auto result{
                rhbm_helper::EstimateBetaMDPDE(
                    alpha_r,
                    dataset_list.at(i),
                    algorithm_options)
            };
            local_bias_array.col(j) = CalculateAbsoluteGaussianDifference(
                result.beta_mdpde,
                result.beta_ols,
                fit_model
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

Eigen::MatrixXd StudyAlphaGBias(
    const std::vector<double> & alpha_grid,
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    LocalGaussianFitModel fit_model,
    const AlphaStudyOptions & options)
{
    ValidateStudyBatch(beta_group_list.size(), alpha_grid);
    for (const auto & beta_list : beta_group_list)
    {
        numeric_validation::RequirePositive(beta_list.size(), "training data size");
    }

    const auto group_size{ beta_group_list.size() };
    const auto alpha_size{ static_cast<int>(alpha_grid.size()) };
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
            const auto alpha_g{ alpha_grid.at(static_cast<std::size_t>(j)) };
            const auto result{
                rhbm_helper::EstimateMuMDPDE(alpha_g, beta_matrix, algorithm_options)
            };
            local_bias_array.col(j) = CalculateAbsoluteGaussianDifference(
                result.mu_mdpde, result.mu_mean, fit_model);
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

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        );
    }
    return dataset_list;
}

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    const FitOptions & options)
{
    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto intercept{ member_result_list.at(i).mdpde.GetModel().GetIntercept() };
        const auto shifted_sample_entries{
            BuildShiftedSampleEntries(sample_entries_list.at(i), intercept)
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        );
    }
    return dataset_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    LocalGaussianFitModel fit_model,
    double intercept = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols, fit_model)
            .WithIntercept(intercept)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde, fit_model)
            .WithIntercept(intercept)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        false,
        0.0,
        fit_model,
        fit_result
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    LocalGaussianFitModel fit_model)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean, fit_model),
        linearization_service::DecodeParameterVector(result.mu_mdpde, fit_model),
        linearization_service::DecodeParameterVector(
            result.mu_prior, result.capital_lambda, fit_model)
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    LocalGaussianFitModel fit_model)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    eigen_validation::RequireVectorSize(
        result.outlier_flag_array, result.beta_posterior_matrix.cols(),
        "outlier_flag_array", "Group Gaussian member result count is inconsistent.");
    eigen_validation::RequireVectorSize(
        result.statistical_distance_array, result.beta_posterior_matrix.cols(),
        "statistical_distance_array", "Group Gaussian member result count is inconsistent.");

    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(static_cast<std::size_t>(i)),
                fit_model)
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i),
            fit_model
        });
    }
    return member_results;
}

std::vector<RHBMBetaEstimateResult> BuildMemberFitResultList(
    const std::vector<LocalGaussianResult> & member_result_list,
    LocalGaussianFitModel fit_model)
{
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        if (!member_result.fit_result.has_value())
        {
            throw std::invalid_argument(
                "member_result_list contains a result without transient fit state.");
        }
        if (member_result.fit_model != fit_model)
        {
            throw std::invalid_argument(
                "group Gaussian estimation requires member results to match the requested local fit model.");
        }
        fit_result_list.emplace_back(*member_result.fit_result);
    }
    return fit_result_list;
}

} // namespace

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options,
    bool output_study_plot)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.alpha_min, options.alpha_max, "alpha training range");
    numeric_validation::RequireFinitePositive(options.alpha_step, "alpha training step");
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, options) };
    const auto training_result{
        rhbm_trainer::CrossValidationAlphaR(
            dataset_list,
            options.alpha_min,
            options.alpha_max,
            options.alpha_step,
            MakeTrainingOptions(kAlphaRSubsetSize, options))
    };

    if (ShouldOutputStudyPlot(output_study_plot, options))
    {
        const auto bias_matrix{
            StudyAlphaRBias(
                training_result.alpha_grid,
                dataset_list,
                options.local_fit_model,
                MakeStudyOptions(options))
        };
        (void)EmitTrainingReportIfRequested(
            bias_matrix, training_result.alpha_grid, "#alpha_{r}", "Deviation with OLS",
            options.study_plot_dir, "alpha_r_bias.pdf");
    }

    return training_result.best_alpha;
}

double TrainAlphaG(
    const std::vector<std::vector<LocalPotentialSampleList>> & sample_group_list,
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options,
    bool output_study_plot)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.alpha_min, options.alpha_max, "alpha training range");
    numeric_validation::RequireFinitePositive(options.alpha_step, "alpha training step");

    if (sample_group_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_group_list and member_result_list sizes are inconsistent.");
    }

    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(sample_group_list.size());
    for (std::size_t i = 0; i < sample_group_list.size(); i++)
    {
        if (sample_group_list.at(i).size() != member_result_list.at(i).size())
        {
            throw std::invalid_argument("sample group and member result group sizes are inconsistent.");
        }

        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(member_result_list.at(i).size());
        for (const auto & member_result : member_result_list.at(i))
        {
            if (member_result.fit_model != options.local_fit_model)
            {
                throw std::invalid_argument(
                    "Alpha_G training requires member results to match the requested local fit model.");
            }
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(
                    member_result.mdpde.GetModel(),
                    member_result.fit_model));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }
    
    if (options.output_summary_log)
    {
        Logger::Log(LogLevel::Info,
            "Run Alpha_G Training with " + std::to_string(beta_group_list.size()) + " groups.");
    }

    if (beta_group_list.empty())
    {
        if (options.output_summary_log)
        {
            Logger::Log(LogLevel::Warning,
                "Skip Alpha_G Training because no eligible groups were selected.");
            Logger::Log(LogLevel::Info,
                "Alpha_G Training Results Summary: minimum mu error sum alpha_g = "
                + std::to_string(options.alpha_min));
        }
        return options.alpha_min;
    }

    const auto training_result{
        rhbm_trainer::CrossValidationAlphaG(
            beta_group_list,
            options.alpha_min,
            options.alpha_max,
            options.alpha_step,
            MakeTrainingOptions(kAlphaGSubsetSize, options))
    };

    if (ShouldOutputStudyPlot(output_study_plot, options))
    {
        const auto bias_matrix{
            StudyAlphaGBias(
                training_result.alpha_grid,
                beta_group_list,
                options.local_fit_model,
                MakeStudyOptions(options))
        };
        (void)EmitTrainingReportIfRequested(
            bias_matrix, training_result.alpha_grid, "#alpha_{g}", "Deviation with Mean",
            options.study_plot_dir, "alpha_g_bias.pdf");
    }

    if (options.output_summary_log)
    {
        Logger::Log(LogLevel::Info,
            "Alpha_G Training Results Summary: minimum mu error sum alpha_g = "
            + std::to_string(training_result.best_alpha));
    }
    return training_result.best_alpha;
}

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");

    auto execution_options{ MakeExecutionOptions(options) };
    double intercept{ EstimateInitialIntercept(sample_entries) };
    RHBMBetaEstimateResult result;
    RHBMParameterVector beta_in_previous_iter;
    bool has_previous_beta{ false };

    for (int t = 0; t < execution_options.max_iterations; t++)
    {
        auto shifted_sample_entries{ BuildShiftedSampleEntries(sample_entries, intercept) };
        auto dataset{
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        };
        result = rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options);
        const auto residual_sample_entries{
            BuildResidualSampleEntries(
                sample_entries,
                result.beta_mdpde,
                options.local_fit_model,
                options.distance_min,
                options.distance_max)
        };
        const auto next_intercept{ EstimateResidualIntercept(residual_sample_entries) };
        if (has_previous_beta && std::fabs(next_intercept - intercept) < execution_options.tolerance)
        {
            break;
        }
        beta_in_previous_iter = result.beta_mdpde;
        has_previous_beta = true;
        if (t + 1 < execution_options.max_iterations)
        {
            intercept = next_intercept;
        }
    }
    return DecodeLocalGaussianResult(alpha_r, result, options.local_fit_model, intercept);
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    auto execution_options{ MakeExecutionOptions(options) };
    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, member_result_list, options) };
    const auto fit_result_list{ BuildMemberFitResultList(
        member_result_list,
        options.local_fit_model) };
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, options.local_fit_model) };
    result.member_results = DecodeMemberGaussianResults(raw_result, options.local_fit_model);
    return result;
}

} // namespace rhbm_gem::core
