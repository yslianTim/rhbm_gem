#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

namespace HRLAlphaTrainer
{
Eigen::VectorXd EvaluateAlphaR(
    const HRLMemberDataset & dataset,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const HRLExecutionOptions & options = {}
);

Eigen::VectorXd EvaluateAlphaG(
    const std::vector<Eigen::VectorXd> & beta_list,
    std::size_t subset_size,
    const std::vector<double> & alpha_list,
    const HRLExecutionOptions & options = {}
);
} // namespace HRLAlphaTrainer
