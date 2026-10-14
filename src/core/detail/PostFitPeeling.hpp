#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <map>
#include <span>
#include <string>
#include <vector>
namespace rhbm_gem::core::detail {
struct PostJointLocalSample
{
    std::size_t raw_sample_index{};
    double distance{};
    double raw_response{};
    std::optional<double> peeled_response;
    std::optional<double> fixed_charge_response;
    std::optional<double> gaussian_response;
    std::string reason;
};

struct PostJointLocalSamples
{
    EstimateSource source;
    JointFixedChargeProfile fixed_charge;
    std::vector<PostJointLocalSample> samples;

    LocalPotentialSampleList BuildGaussianSamples(const LocalPotentialSampleList &) const;
};

// Reads geometry, immutable support, raw samples and fixed Second states only.
std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const MapObject &, const ModelObject &, const JointProblem &,
    std::optional<std::span<const std::size_t>> outputs = std::nullopt,
    const JointAnalysisResult * result = nullptr);

PostJointLocalSamples BuildPostJointLocalSamples(
    const MapObject &, const ModelObject &, const JointProblem &, int atom_id,
    const PostFitPeelingResult &, const JointAnalysisResult &);
}
