#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::rhbm_trainer
{

struct RHBMTrainingOptions
{
    double alpha_min{ 0.0 };
    double alpha_max{ 2.0 };
    double alpha_step{ 0.1 };
    RHBMExecutionOptions execution_options{};
    std::size_t subset_size{ 5 };
};

struct RHBMTrainingResult
{
    double best_alpha{ 0.0 };
    Eigen::VectorXd error_sum_list;
    std::vector<double> alpha_grid;
};

RHBMTrainingResult CrossValidationAlphaR(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const RHBMTrainingOptions & options
);

RHBMTrainingResult CrossValidationAlphaG(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    const RHBMTrainingOptions & options
);

} // namespace rhbm_gem::rhbm_trainer
