#pragma once

#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {
class ModelObject;

namespace core {

struct FitOptions
{
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
    double distance_min{ 0.0 };
    double distance_max{ 1.0 };
    int thread_size{ 1 };
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

void RunLocalAlphaTraining(ModelObject & model_object, const FitOptions & options);
void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options);
void RunLocalPotentialFitting(ModelObject & model_object, const FitOptions & options);
void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options);

} // namespace core

} // namespace rhbm_gem
