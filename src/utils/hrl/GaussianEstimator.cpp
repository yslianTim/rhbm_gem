#include <rhbm_gem/utils/hrl/GaussianEstimator.hpp>

#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>

#include <filesystem>
#include <optional>
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
    const CrossValidationOptions & options,
    bool output_progress)
{
    rhbm_trainer::AlphaTrainer::AlphaTrainingOptions training_options;
    training_options.subset_size = subset_size;
    training_options.execution_options = MakeExecutionOptions(options);
    if (output_progress)
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
    study_options.progress_callback =
        [](std::size_t completed, std::size_t total)
        {
            Logger::ProgressPercent(completed, total);
        };
    return study_options;
}

} // namespace

double CrossValidationAlphaR(
    const LocalPotentialSampleList & sample_entries,
    bool output_study_plot)
{
    const CrossValidationOptions options;
    const std::vector<RHBMMemberDataset> dataset_list{
        rhbm_helper::BuildMemberDataset(sample_entries, options.fit_range_min, options.fit_range_max)
    };
    return CrossValidationAlphaR(dataset_list, options, output_study_plot);
}

double CrossValidationAlphaG(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    bool output_study_plot)
{
    const CrossValidationOptions options;
    return CrossValidationAlphaG(beta_group_list, options, output_study_plot);
}

double CrossValidationAlphaR(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const CrossValidationOptions & options,
    bool output_study_plot)
{
    const auto trainer{ MakeAlphaTrainer(options) };
    const auto training_result{
        trainer.TrainAlphaR(dataset_list, MakeTrainingOptions(kAlphaRSubsetSize, options, false))
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
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    const CrossValidationOptions & options,
    bool output_study_plot)
{
    Logger::Log(LogLevel::Info,
        "Run Alpha_G Training with " + std::to_string(beta_group_list.size()) + " groups.");

    if (beta_group_list.empty())
    {
        Logger::Log(LogLevel::Warning,
            "Skip Alpha_G Training because no eligible groups were selected.");
        Logger::Log(LogLevel::Info,
            "Alpha_G Training Results Summary: minimum mu error sum alpha_g = "
            + std::to_string(options.alpha_min));
        return options.alpha_min;
    }

    const auto trainer{ MakeAlphaTrainer(options) };
    const auto training_result{
        trainer.TrainAlphaG(beta_group_list, MakeTrainingOptions(kAlphaGSubsetSize, options, true))
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

    Logger::Log(LogLevel::Info,
        "Alpha_G Training Results Summary: minimum mu error sum alpha_g = "
        + std::to_string(training_result.best_alpha));
    return training_result.best_alpha;
}

} // namespace rhbm_gem::gaussian_estimator
