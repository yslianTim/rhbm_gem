#include <rhbm_gem/utils/hrl/GaussianEstimator.hpp>

#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::gaussian_estimator {
namespace {
constexpr std::size_t kAlphaRSubsetSize{ 5 };
constexpr std::size_t kAlphaGSubsetSize{ 10 };

RHBMExecutionOptions MakeExecutionOptions(const CrossValidationOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = true;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

rhbm_trainer::AlphaTrainer MakeAlphaTrainer(const CrossValidationOptions & options)
{
    rhbm_trainer::AlphaTrainer trainer{
        options.alpha_min,
        options.alpha_max,
        options.alpha_step
    };
    return trainer;
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
    amplitude_values.reserve(alpha_list.size());
    width_values.reserve(alpha_list.size());
    for (Eigen::Index column = 0; column < gaus_bias_matrix.cols(); ++column)
    {
        amplitude_values.emplace_back(gaus_bias_matrix(0, column));
        width_values.emplace_back(gaus_bias_matrix(1, column));
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
    const CrossValidationOptions & options)
{
    if (!output_study_plot) return false;
    if (!options.study_plot_dir.empty()) return true;

    Logger::Log(LogLevel::Debug,
        "Alpha training report was skipped because no study plot directory was provided.");
    return false;
}

rhbm_trainer::AlphaTrainer::AlphaTrainingOptions MakeTrainingOptions(
    std::size_t subset_size,
    const CrossValidationOptions & options)
{
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions training_options;
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

rhbm_trainer::AlphaTrainer::AlphaRunOptions MakeStudyOptions(
    const CrossValidationOptions & options)
{
    rhbm_trainer::AlphaTrainer::AlphaRunOptions study_options;
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

RHBMMemberDataset BuildMemberDataset(
    const LocalPotentialSampleList & sample_entries,
    const CrossValidationOptions & options)
{
    return rhbm_helper::BuildMemberDataset(
        sample_entries, options.fit_range_min, options.fit_range_max);
}

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const CrossValidationOptions & options)
{
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(BuildMemberDataset(sample_entries, options));
    }
    return dataset_list;
}

RHBMBetaEstimateResult EstimateMemberBeta(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const CrossValidationOptions & options)
{
    return rhbm_helper::EstimateBetaMDPDE(
        alpha_r,
        BuildMemberDataset(sample_entries, options),
        MakeExecutionOptions(options));
}

std::vector<RHBMBetaEstimateResult> EstimateMemberBetaList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<double> & alpha_r_list,
    const CrossValidationOptions & options)
{
    if (sample_entries_list.size() != alpha_r_list.size())
    {
        throw std::invalid_argument("sample_entries_list and alpha_r_list sizes are inconsistent.");
    }

    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        fit_result_list.emplace_back(
            EstimateMemberBeta(sample_entries_list.at(i), alpha_r_list.at(i), options));
    }
    return fit_result_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result)
{
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{
            linearization_service::DecodeParameterVector(fit_result.beta_ols),
            GaussianModel3DUncertainty{}
        },
        GaussianModel3DWithUncertainty{
            linearization_service::DecodeParameterVector(fit_result.beta_mdpde),
            GaussianModel3DUncertainty{}
        },
        false,
        0.0
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean),
        linearization_service::DecodeParameterVector(result.mu_mdpde),
        linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda)
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (result.capital_sigma_posterior_list.size() != member_count ||
        result.outlier_flag_array.rows() != static_cast<Eigen::Index>(member_count) ||
        result.statistical_distance_array.rows() != static_cast<Eigen::Index>(member_count))
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }

    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(static_cast<std::size_t>(i)))
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

RHBMGroupEstimationResult EstimateGroupFromSamples(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<double> & alpha_r_list,
    double alpha_g,
    const CrossValidationOptions & options)
{
    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, options) };
    const auto fit_result_list{ EstimateMemberBetaList(sample_entries_list, alpha_r_list, options) };
    return rhbm_helper::EstimateGroup(
        alpha_g,
        rhbm_helper::BuildGroupInput(dataset_list, fit_result_list),
        MakeExecutionOptions(options));
}

} // namespace

double CrossValidationAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const CrossValidationOptions & options,
    bool output_study_plot)
{
    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, options) };
    const auto trainer{ MakeAlphaTrainer(options) };
    const auto training_result{
        trainer.TrainAlphaR(dataset_list, MakeTrainingOptions(kAlphaRSubsetSize, options))
    };

    if (ShouldOutputStudyPlot(output_study_plot, options))
    {
        const auto bias_matrix{
            trainer.StudyAlphaRBias(dataset_list, MakeStudyOptions(options))
        };
        (void)EmitTrainingReportIfRequested(
            bias_matrix, trainer.AlphaGrid(), "#alpha_{r}", "Deviation with OLS",
            options.study_plot_dir, "alpha_r_bias.pdf");
    }

    return training_result.best_alpha;
}

double CrossValidationAlphaG(
    const std::vector<std::vector<LocalPotentialSampleList>> & sample_group_list,
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const CrossValidationOptions & options,
    bool output_study_plot)
{
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
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(
                    member_result.mdpde.GetModel()));
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

    const auto trainer{ MakeAlphaTrainer(options) };
    const auto training_result{
        trainer.TrainAlphaG(beta_group_list, MakeTrainingOptions(kAlphaGSubsetSize, options))
    };

    if (ShouldOutputStudyPlot(output_study_plot, options))
    {
        const auto bias_matrix{
            trainer.StudyAlphaGBias(beta_group_list, MakeStudyOptions(options))
        };
        (void)EmitTrainingReportIfRequested(
            bias_matrix, trainer.AlphaGrid(), "#alpha_{g}", "Deviation with Mean",
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
    const CrossValidationOptions & options)
{
    return DecodeLocalGaussianResult(
        alpha_r,
        EstimateMemberBeta(sample_entries, alpha_r, options));
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<double> & alpha_r_list,
    double alpha_g,
    const CrossValidationOptions & options)
{
    const auto raw_result{
        EstimateGroupFromSamples(sample_entries_list, alpha_r_list, alpha_g, options)
    };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result) };
    result.member_results = DecodeMemberGaussianResults(raw_result);
    return result;
}

} // namespace rhbm_gem::gaussian_estimator
