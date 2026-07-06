#pragma once

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
};

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

LocalGaussianResult EstimateLocalGaussianWithOffset(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    double offset_initial = 0.0);

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options);

void RunLocalAlphaTraining(ModelObject & model_object, const FitOptions & options);
void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options);
void RunFirstStageLocalFitting(ModelObject & model_object, const FitOptions & options, bool fit_offset = false);
void RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options);
void RunThirdStageLocalFitting(ModelObject & model_object, const FitOptions & options);
void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options);
void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options);

} // namespace core

} // namespace rhbm_gem
