#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {
class ModelObject;

namespace core {

struct FitOptions
{
    double distance_min{ 0.0 };
    double distance_max{ 1.0 };
    int thread_size{ 1 };
    bool quiet_mode{ false };
    bool exclude_hydrogen{ false };
    std::optional<std::filesystem::path> local_fitting_result_csv_path{};
};

struct LocalGaussianDesignTemplate
{
    std::size_t source_sample_count{ 0 };
    std::vector<std::size_t> source_sample_index_list{};
    std::vector<double> distance_list{};
    RHBMDesignMatrix design_matrix{};
};

LocalGaussianDesignTemplate BuildLocalGaussianDesignTemplate(
    const LocalPotentialSampleList & sample_entries,
    double range_min,
    double range_max);

RHBMMemberDataset BuildLocalGaussianPreparedDataset(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    const GaussianModel3D & offset_model);

LocalGaussianResult EstimateLocalGaussianPrepared(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model = GaussianModel3D{ 0.0, 1.0, 0.0 });

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options);

double TrainAlphaG(
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options);

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model = GaussianModel3D{ 0.0, 1.0, 0.0 });

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options);

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);

void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);

void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);

bool RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options);

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options);

} // namespace core

} // namespace rhbm_gem
