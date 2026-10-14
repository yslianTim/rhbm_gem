#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <map>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>
namespace rhbm_gem::core::detail {
struct JointPostprocessingMetrics
{
    std::size_t observation_row_mapping_builds{};
    std::size_t contributor_mapping_builds{};
    std::size_t joint_stage_mapping_builds{};
    std::size_t target_preparation_count{};
    std::size_t total_prepared_sample_count{};
    std::size_t valid_local_sample_count{};
    std::size_t invalid_local_sample_count{};
};

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

// Owns compact, immutable lookups shared by all targets in one Joint run.
class JointPostprocessingContext
{
public:
    JointPostprocessingContext(const MapObject &, const ModelObject &, const JointProblem &,
        const JointAnalysisResult *, JointPostprocessingMetrics * = nullptr);

    std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
        std::optional<std::span<const std::size_t>> outputs = std::nullopt) const;
    PostJointLocalSamples BuildPostJointLocalSamples(
        int atom_id, const PostFitPeelingResult &) const;

private:
    const MapObject & m_geometry;
    const ModelObject & m_model;
    const JointProblem & m_problem;
    const JointAnalysisResult * m_snapshot;
    JointPostprocessingMetrics * m_metrics;
    std::unordered_map<std::size_t, std::size_t> m_row_index_by_grid;
    std::unordered_map<int, std::size_t> m_contributor_index_by_id;
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
