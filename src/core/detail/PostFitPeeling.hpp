#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <cstddef>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>
namespace rhbm_gem::core::detail {
struct JointPostprocessingMetrics
{
    std::size_t snapshot_identity_validation_count{};
    std::size_t full_row_identity_scan_count{};
    std::size_t observation_row_mapping_builds{};
    std::size_t contributor_mapping_builds{};
    std::size_t joint_stage_mapping_builds{};
    std::size_t target_preparation_count{};
    std::size_t total_prepared_sample_count{};
    std::size_t valid_local_sample_count{};
    std::size_t invalid_local_sample_count{};
    std::size_t estimated_context_lookup_bytes{};
    std::size_t estimated_peeling_result_bytes{};
    std::size_t estimated_peak_postprocessing_payload_bytes{};
    double context_construction_seconds{};
    double peeling_seconds{};
    double charge_correction_seconds{};
    double local_alpha_training_seconds{};
    double local_mdpde_fitting_seconds{};
    double group_alpha_training_seconds{};
    double group_web_seconds{};
    double total_postprocessing_seconds{};
    std::optional<std::size_t> rss_before_postprocessing_bytes;
    std::optional<std::size_t> rss_after_context_bytes;
    std::optional<std::size_t> rss_after_peeling_bytes;
    std::optional<std::size_t> rss_after_local_bytes;
    std::optional<std::size_t> rss_after_group_bytes;
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

// Borrows one immutable Joint run and owns compact lookups shared by its targets.
// The geometry, model endpoints, problem input, and snapshot must outlive this
// context and remain unchanged until postprocessing completes.
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
