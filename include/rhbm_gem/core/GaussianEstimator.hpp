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
    std::size_t second_stage_boundary_halo_depth{ 1 };
    bool enable_second_stage_dependency_polish{ true };
    std::size_t second_stage_dependency_polish_max_iterations{ 10 };
    std::optional<std::filesystem::path> result_csv_path{};
};

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options);

double TrainAlphaG(
    const std::vector<std::vector<GaussianModel3D>> & model_group_list,
    const FitOptions & options);

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model = GaussianModel3D{ 0.0, 1.0, 0.0 });

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<GroupGaussianMemberInput> & member_list,
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
    const FitOptions & options);

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options);

} // namespace core

} // namespace rhbm_gem
