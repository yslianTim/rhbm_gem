#pragma once

#include <cstddef>
#include <rhbm_gem/core/PotentialEstimator.hpp>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {
class ModelObject;
class MapObject;
class AtomObject;

namespace core {

struct FitOptions
{
    PotentialEstimator estimator{ PotentialEstimator::JOINT_COMPONENTS };
    int thread_size{ 1 };
    bool quiet_mode{ false };
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

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);

void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage);

void RunLocalAlphaTraining(ModelObject &, const FitOptions &, FittingStage,
    const std::vector<AtomObject *> & atoms);
void RunFixedOffsetLocalFitting(ModelObject &, const FitOptions &, FittingStage,
    const std::vector<AtomObject *> & atoms);
void RunPotentialFittingWorkflow(MapObject &, ModelObject &, const FitOptions &);

void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options);

} // namespace core

} // namespace rhbm_gem
