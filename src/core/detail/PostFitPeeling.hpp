#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>
namespace rhbm_gem::core::detail {
struct HaloMedianBuildResult;

struct JointPeelingContributor
{
    int atom_id{};
    std::array<double,3> position{};
    std::optional<std::size_t> contributor_index;
    bool target{};
    bool reconstruction_eligible{};
};

struct JointPeelingContributorUniverse
{
    std::string joint_run_id;
    std::vector<JointPeelingContributor> contributors;
};

JointPeelingContributorUniverse BuildJointPeelingContributorUniverse(
    const ModelObject &,const JointProblem &,const JointAnalysisResult *);

struct JointPostprocessingMetrics
{
    std::size_t snapshot_identity_validation_count{};
    std::size_t full_row_identity_scan_count{};
    std::size_t observation_row_mapping_builds{};
    std::size_t contributor_mapping_builds{};
    std::size_t joint_stage_mapping_builds{};
    std::size_t halo_median_valid_target_count{};
    std::size_t halo_median_geometry_local_halos{};
    std::size_t halo_median_unavailable_halos{};
    std::size_t reconstruction_lookup_atom_count{};
    std::size_t reconstruction_only_node_evaluation_count{};
    std::size_t reconstruction_node_cache_hit_count{};
    std::size_t reconstruction_node_cache_peak_entries{};
    std::size_t estimated_reconstruction_workspace_bytes{};
    // Peeling availability diagnostics.  The historical sample reason
    // "outside-joint-domain" remains the serialized compatibility label;
    // these counters split it into map geometry and original-domain causes.
    std::size_t raw_peeling_sample_count{};
    std::size_t nonzero_stencil_node_count{};
    std::size_t available_peeling_sample_count{};
    std::size_t unavailable_peeling_sample_count{};
    std::size_t map_geometry_boundary_samples{};
    std::size_t map_geometry_boundary_nodes{};
    std::size_t original_support_boundary_samples{};
    std::size_t original_support_boundary_nodes{};
    std::size_t fitting_mask_only_intersection_samples{};
    std::size_t fitting_mask_only_intersection_nodes{};
    std::size_t mapping_mismatch_samples{};
    std::size_t missing_contributor_samples{};
    std::size_t other_peeling_failure_samples{};
    std::size_t tail_target_count{};
    std::size_t tail_full_coverage_targets{};
    std::size_t tail_partial_coverage_targets{};
    std::size_t tail_zero_coverage_targets{};
    std::size_t tail_no_sample_targets{};
    std::size_t tail_ratio_available_targets{};
    std::size_t tail_raw_sample_count{};
    std::size_t tail_valid_sample_count{};
    std::size_t tail_unavailable_sample_count{};
    double tail_raw_response_sum{};
    double tail_peeling_response_sum{};
    std::map<std::string, std::size_t> tail_coverage_distribution;
    std::map<std::string, std::size_t> tail_unavailable_reasons;
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
        const JointAnalysisResult *, JointPostprocessingMetrics * = nullptr,
        const HaloMedianBuildResult * = nullptr);
    JointPostprocessingContext(const MapObject &, const ModelObject &, const JointProblem &,
        const JointAnalysisResult *, JointPostprocessingMetrics *, const HaloMedianBuildResult *,
        JointPeelingContributorUniverse);

    std::vector<std::pair<int, PostFitPeelingResult>> BuildPostFitPeelingBatch(
        std::optional<std::span<const std::size_t>> outputs = std::nullopt) const;
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
    const HaloMedianBuildResult * m_halo_medians;
    JointPeelingContributorUniverse m_contributors;
    std::unordered_map<std::size_t, std::size_t> m_row_index_by_grid;
    std::unordered_map<std::size_t, std::size_t> m_fitting_row_index_by_grid;
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
