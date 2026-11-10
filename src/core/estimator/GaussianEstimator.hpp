#pragma once

#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {
class ModelObject;
class AtomObject;

struct GroupGaussianMemberInput
{
    LocalPotentialSampleList sample_entries{};
    double alpha_r{ 0.0 };
    GaussianModel3D local_model{};
};

namespace core {

struct FitOptions
{
    int thread_size{ 1 };
    bool quiet_mode{ false };
    int joint_initialization_thread_size{ 1 };
    SphereSamplingMethod first_stage_sampling_method{ SphereSamplingMethod::FibonacciDeterministic };
};

double TrainAlphaR(
    const LocalPotentialSampleList & sample_entries,
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

// Reconstructs a Joint group member's posterior signal with its fixed Joint charge profile.
double EvaluateGroupMemberResponseAtDistance(const GroupGaussianMemberResult &, double distance);

void RunFirstStageLocalAlphaTraining(ModelObject &, const FitOptions &);
void RunFirstStageLocalFitting(ModelObject &, const FitOptions &);

void RunFirstStageLocalAlphaTraining(ModelObject &, const FitOptions &,
    const std::vector<AtomObject *> & atoms);
void RunFirstStageLocalFitting(ModelObject &, const FitOptions &,
    const std::vector<AtomObject *> & atoms);

} // namespace core

} // namespace rhbm_gem
