#include "command/PotentialAnalysisCommand.hpp"

#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem::detail {

std::vector<double> BuildOrderedAlphaRTrainingList()
{
    std::vector<double> alpha_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    std::sort(alpha_list.begin(), alpha_list.end());
    return alpha_list;
}

std::vector<double> BuildOrderedAlphaGTrainingList()
{
    std::vector<double> alpha_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    std::sort(alpha_list.begin(), alpha_list.end());
    return alpha_list;
}

std::filesystem::path ResolveTrainingReportPath(
    const PotentialAnalysisCommandOptions & options,
    std::string_view report_file_name)
{
    if (options.training_report_dir.empty())
    {
        return {};
    }
    return options.training_report_dir / std::filesystem::path{ report_file_name };
}

bool EmitTrainingReportIfRequested(
    const Eigen::MatrixXd & gaus_bias_matrix,
    const std::vector<double> & alpha_list,
    std::string_view x_label,
    std::string_view y_label,
    const PotentialAnalysisCommandOptions & options,
    std::string_view report_file_name)
{
    const auto report_path{ ResolveTrainingReportPath(options, report_file_name) };
    if (report_path.empty())
    {
        return false;
    }

    std::error_code ec;
    std::filesystem::create_directories(report_path.parent_path(), ec);
    if (ec)
    {
        Logger::Log(
            LogLevel::Warning,
            "Failed to create training report directory '" +
                report_path.parent_path().string() + "': " + ec.message());
        return false;
    }

    try
    {
        LocalPainter::PaintTemplate1(
            gaus_bias_matrix,
            alpha_list,
            std::string(x_label),
            std::string(y_label),
            report_path.string());
    }
    catch (const std::exception & ex)
    {
        Logger::Log(
            LogLevel::Warning,
            "Failed to emit training report '" + report_path.string() + "': " + ex.what());
        return false;
    }

    if (!std::filesystem::exists(report_path))
    {
        Logger::Log(
            LogLevel::Warning,
            "Training report output was requested but no file was produced: " +
                report_path.string());
        return false;
    }

    return true;
}

} // namespace rhbm_gem::detail
