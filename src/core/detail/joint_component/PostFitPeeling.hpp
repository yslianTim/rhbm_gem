#pragma once
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <vector>
namespace rhbm_gem::core::detail {
struct HaloMedianBuildResult;

inline constexpr double PostJointCandidateRadius = 2.0;
std::vector<PostJointGridObservation> CollectPostJointGridObservations(
    const MapObject &, const JointProblemInput &, std::size_t target,
    const std::array<double,3> & target_position);

struct JointPeelingContributor
{
    int atom_id{};
    std::array<double,3> position{};
    std::size_t contributor_index{};
    bool target{};
    bool reconstruction_eligible{};
};

struct JointPeelingContributorUniverse
{
    std::string joint_run_id;
    // BuildJointPeelingContributorUniverse returns contributors sorted by atom_id.
    std::vector<JointPeelingContributor> contributors;
};

JointPeelingContributorUniverse BuildJointPeelingContributorUniverse(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &);

struct JointPostprocessingMetrics
{
    std::size_t halo_median_valid_target_count{};
    std::size_t halo_median_geometry_local_halos{};
    std::size_t halo_median_unavailable_halos{};
    std::size_t estimated_peeling_workspace_bytes{};
    std::size_t grid_peeling_sample_count{};
    std::size_t available_peeling_sample_count{};
    std::size_t unavailable_peeling_sample_count{};
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

struct PostJointLocalSample : PostJointGridObservation
{
    double distance{};
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

    LocalPotentialSampleList BuildGaussianSamples() const;
};

// Borrows one immutable Joint run and owns compact lookups shared by its targets.
// The geometry, model endpoints, problem input, snapshot, and optional Halo donor
// workspace must outlive this context and remain unchanged until postprocessing completes.
class JointPostprocessingContext
{
public:
    JointPostprocessingContext(const MapObject &, const ModelObject &, const JointProblem &,
        const JointAnalysisResult &, JointPostprocessingMetrics *, const HaloMedianBuildResult *,
        JointPeelingContributorUniverse);

    std::vector<std::pair<int, PostFitPeelingResult>> BuildPostFitPeelingBatch(
        std::span<const std::size_t> outputs) const;
    PostJointLocalSamples BuildPostJointLocalSamples(
        int atom_id, const PostFitPeelingResult &) const;

private:
    const MapObject & m_geometry;
    const ModelObject & m_model;
    const JointProblem & m_problem;
    const JointAnalysisResult & m_snapshot;
    JointPostprocessingMetrics * m_metrics;
    const HaloMedianBuildResult * m_halo_medians;
    JointPeelingContributorUniverse m_contributors;
};

}
