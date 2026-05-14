#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::rhbm_trainer
{

using ProgressCallback = std::function<void(std::size_t completed, std::size_t total)>;

struct AlphaTrainingOptions
{
    RHBMExecutionOptions execution_options{};
    ProgressCallback progress_callback{};
    std::size_t subset_size{ 0 };
};

struct AlphaTrainingResult
{
    double best_alpha{ 0.0 };
    Eigen::VectorXd error_sum_list;
    std::vector<double> alpha_grid;
};

AlphaTrainingResult TrainAlphaR(
    const std::vector<RHBMMemberDataset> & dataset_list,
    double alpha_min,
    double alpha_max,
    double alpha_step,
    const AlphaTrainingOptions & options
);

AlphaTrainingResult TrainAlphaG(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    double alpha_min,
    double alpha_max,
    double alpha_step,
    const AlphaTrainingOptions & options
);

} // namespace rhbm_gem::rhbm_trainer
