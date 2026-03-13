#pragma once

#include <filesystem>
#include <string_view>
#include <vector>

#include <Eigen/Dense>

#include "command/PotentialAnalysisCommand.hpp"

namespace rhbm_gem::detail {

std::vector<double> BuildOrderedAlphaRTrainingList();
std::vector<double> BuildOrderedAlphaGTrainingList();

std::filesystem::path ResolveTrainingReportPath(
    const PotentialAnalysisCommandOptions & options,
    std::string_view report_file_name);

bool EmitTrainingReportIfRequested(
    const Eigen::MatrixXd & gaus_bias_matrix,
    const std::vector<double> & alpha_list,
    std::string_view x_label,
    std::string_view y_label,
    const PotentialAnalysisCommandOptions & options,
    std::string_view report_file_name);

} // namespace rhbm_gem::detail
