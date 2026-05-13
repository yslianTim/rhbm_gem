#pragma once

#include <filesystem>
#include <vector>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::gaussian_estimator {

struct CrossValidationOptions
{
    double alpha_min{ 0.0 };
    double alpha_max{ 2.0 };
    double alpha_step{ 0.1 };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    int thread_size{ 1 };
    std::filesystem::path study_plot_dir{};
};

double CrossValidationAlphaR(
    const LocalPotentialSampleList & sample_entries,
    bool output_study_plot = false);

double CrossValidationAlphaG(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    bool output_study_plot = false);

double CrossValidationAlphaR(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const CrossValidationOptions & options,
    bool output_study_plot = false);

double CrossValidationAlphaG(
    const std::vector<std::vector<RHBMParameterVector>> & beta_group_list,
    const CrossValidationOptions & options,
    bool output_study_plot = false);

} // namespace rhbm_gem::gaussian_estimator
