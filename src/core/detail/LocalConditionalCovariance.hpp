#pragma once

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem::core::detail {

struct LocalConditionalCovariance
{
    EstimateSource source;
    EvidenceStatus status{ EvidenceStatus::Unavailable };
    std::string method{ "local-conditional-iid-ls" };
    std::string reason{ "not-computed" };
    // Covariance in GroupParameterEvidence's (beta_1, beta_2) coordinates.
    std::optional<Eigen::Matrix2d> covariance;
    std::optional<double> residual_variance;
    std::size_t observation_count{};
    std::size_t degrees_of_freedom{};
    std::size_t temporary_workspace_estimate_bytes{};
};

std::optional<std::size_t> EstimateLocalConditionalCovarianceWorkspaceBytes(std::size_t observations);

LocalConditionalCovariance ComputeLocalConditionalCovariance(
    const LocalStageEstimate &, const LocalPotentialSampleList &,
    const PostFitPeelingResult &, bool sample_geometry_available);

GroupParameterEvidence BuildLocalGroupEvidence(
    const LocalStageEstimate &, const PostFitPeelingResult *,
    const LocalConditionalCovariance &);

} // namespace rhbm_gem::core::detail
