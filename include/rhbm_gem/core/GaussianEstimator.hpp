#pragma once

#include <filesystem>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::core {

struct FitOptions
{
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
    double distance_min{ 0.0 };
    double distance_max{ 1.0 };
    int thread_size{ 1 };
    bool output_progress{ false };
    bool output_summary_log{ false };
    std::filesystem::path study_plot_dir{};
};

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options,
    bool output_study_plot = false);

double TrainAlphaG(
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options,
    bool output_study_plot = false);

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options);

LocalGaussianResult EstimateLocalGaussianWithIntercept(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options);

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options);

} // namespace rhbm_gem::core
