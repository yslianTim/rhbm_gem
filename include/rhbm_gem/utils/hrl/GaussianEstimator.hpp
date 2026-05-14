#pragma once

#include <filesystem>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::gaussian_estimator {

struct TrainingOptions
{
    double alpha_min{ 0.0 };
    double alpha_max{ 2.0 };
    double alpha_step{ 0.1 };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    int thread_size{ 1 };
    bool output_progress{ false };
    bool output_summary_log{ false };
    std::filesystem::path study_plot_dir{};
};

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const TrainingOptions & options,
    bool output_study_plot = false);

double TrainAlphaG(
    const std::vector<std::vector<LocalPotentialSampleList>> & sample_group_list,
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const TrainingOptions & options,
    bool output_study_plot = false);

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const TrainingOptions & options);

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const TrainingOptions & options);

} // namespace rhbm_gem::gaussian_estimator
