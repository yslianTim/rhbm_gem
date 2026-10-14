#include "LocalConditionalCovariance.hpp"

#include "joint_component/Numerics.hpp"

#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::detail {
namespace {

bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role;
}

bool PositiveDefinite(const Eigen::Matrix2d & covariance)
{
    if (!covariance.allFinite()) return false;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> spectrum(covariance);
    return spectrum.info() == Eigen::Success && spectrum.eigenvalues()(0) >
        std::numeric_limits<double>::epsilon() * 2.0 * spectrum.eigenvalues()(1);
}

} // namespace

std::optional<std::size_t> EstimateLocalConditionalCovarianceWorkspaceBytes(std::size_t observations)
{
    constexpr std::size_t bytes_per_observation = 4 * sizeof(double);
    if (observations > std::numeric_limits<std::size_t>::max() / bytes_per_observation) return std::nullopt;
    return observations * bytes_per_observation;
}

LocalConditionalCovariance ComputeLocalConditionalCovariance(
    const LocalStageEstimate & stage, const LocalPotentialSampleList & raw_samples,
    const PostFitPeelingResult & peeling, bool sample_geometry_available)
{
    LocalConditionalCovariance result;
    result.source = stage.source;
    if (stage.source.method != EstimateMethod::JointComponents)
    { result.reason = "not-a-joint-estimate"; return result; }
    if (!stage.point) { result.reason = "missing-point-state"; return result; }
    if (stage.convergence != JointCheckStatus::Passed)
    { result.reason = "joint-endpoint-not-converged"; return result; }
    const auto & point = *stage.point;
    const double amplitude = point.GetAmplitude(), width = point.GetWidth(), charge = point.GetOffset();
    if (!(amplitude > 0.0) || !std::isfinite(amplitude))
    { result.reason = "nonpositive-amplitude"; return result; }
    if (!(width > 0.0) || !std::isfinite(width))
    { result.reason = "nonpositive-width"; return result; }
    if (!std::isfinite(charge)) { result.reason = "nonfinite-charge"; return result; }
    if (!SameSource(stage.source, peeling.source))
    { result.reason = "peeling-source-mismatch"; return result; }
    if (!sample_geometry_available) { result.reason = "sample-geometry-unavailable"; return result; }
    if (raw_samples.size() != peeling.samples.size())
    { result.reason = "peeling-sample-count-mismatch"; return result; }

    Eigen::Index available = 0;
    for (const auto & sample : peeling.samples)
        if (sample.response && std::isfinite(*sample.response)) ++available;
    result.observation_count = static_cast<std::size_t>(available);
    if (available <= 2) { result.reason = "insufficient-observations"; return result; }
    result.degrees_of_freedom = static_cast<std::size_t>(available - 2);
    const auto workspace_bytes = EstimateLocalConditionalCovarianceWorkspaceBytes(result.observation_count);
    if (!workspace_bytes) { result.reason = "local-workspace-size-overflow"; return result; }
    result.temporary_workspace_estimate_bytes = *workspace_bytes;

    Eigen::MatrixXd jacobian(available, 2);
    long double residual_sum_squares = 0.0L;
    Eigen::Index row = 0;
    for (std::size_t i = 0; i < peeling.samples.size(); ++i)
    {
        const auto & sample = peeling.samples[i];
        if (!sample.response || !std::isfinite(*sample.response)) continue;
        const double distance = raw_samples[i].point.distance;
        if (!(distance >= 0.0) || !std::isfinite(distance))
        { result.reason = "invalid-sample-distance"; return result; }
        const auto basis = joint_component::EvaluateKernel(distance * distance, width, 2.5);
        const double prediction = amplitude * basis.gaussian + charge * basis.charge;
        jacobian(row, 0) = basis.gaussian;
        jacobian(row, 1) =
            (amplitude * basis.gaussian_log_width + charge * basis.charge_log_width) / width;
        const double residual = *sample.response - prediction;
        if (!std::isfinite(prediction) || !jacobian.row(row).allFinite() || !std::isfinite(residual))
        { result.reason = "nonfinite-local-design"; return result; }
        residual_sum_squares += static_cast<long double>(residual) * residual;
        ++row;
    }

    const Eigen::Vector2d column_norms = jacobian.colwise().norm();
    if (!column_norms.allFinite() || (column_norms.array() <= 0.0).any())
    { result.reason = "rank-deficient-local-design"; return result; }
    jacobian.col(0) /= column_norms(0);
    jacobian.col(1) /= column_norms(1);

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(jacobian);
    const double relative_threshold = std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max<Eigen::Index>(available, 2));
    qr.setThreshold(relative_threshold);
    if (qr.rank() != 2)
    { result.reason = "rank-deficient-local-design"; return result; }

    Eigen::Matrix2d r = qr.matrixR().topLeftCorner<2, 2>().template triangularView<Eigen::Upper>();
    r(1, 0) = 0.0;
    const Eigen::JacobiSVD<Eigen::Matrix2d> spectrum(r);
    if (spectrum.info() != Eigen::Success || !spectrum.singularValues().allFinite() ||
        spectrum.singularValues()(1) <= relative_threshold * spectrum.singularValues()(0))
    { result.reason = "ill-conditioned-local-design"; return result; }

    const long double variance = residual_sum_squares / static_cast<long double>(result.degrees_of_freedom);
    result.residual_variance = static_cast<double>(variance);
    if (!(variance > 0.0L) || !std::isfinite(variance) || !std::isfinite(*result.residual_variance))
    { result.reason = "residual-variance-unavailable"; return result; }

    const Eigen::Matrix2d inverse_r = r.template triangularView<Eigen::Upper>().solve(Eigen::Matrix2d::Identity());
    const Eigen::Matrix2d pivoted_inverse = inverse_r * inverse_r.transpose();
    const Eigen::Matrix2d inverse_normalized = qr.colsPermutation() * pivoted_inverse * qr.colsPermutation().transpose();
    const Eigen::Matrix2d inverse_norms = column_norms.cwiseInverse().asDiagonal();
    const Eigen::Matrix2d covariance_ab = static_cast<double>(variance) *
        inverse_norms * inverse_normalized * inverse_norms;
    Eigen::Matrix2d transform;
    transform << 1.0 / amplitude, -3.0 / width,
        0.0, -2.0 / (width * width * width);
    Eigen::Matrix2d covariance = transform * covariance_ab * transform.transpose();
    covariance = 0.5 * (covariance + covariance.transpose()).eval();
    if (!covariance.allFinite()) { result.reason = "nonfinite-parameter-covariance"; return result; }
    if (!PositiveDefinite(covariance)) { result.reason = "singular-parameter-covariance"; return result; }

    result.covariance = std::move(covariance);
    result.status = EvidenceStatus::Available;
    result.reason.clear();
    return result;
}

GroupParameterEvidence BuildLocalGroupEvidence(
    const LocalStageEstimate & stage, const PostFitPeelingResult * peeling,
    const LocalConditionalCovariance & covariance)
{
    GroupParameterEvidence evidence;
    evidence.atom_id = stage.source.atom_id;
    evidence.component_id = stage.source.component_id;
    evidence.source_id = stage.source.run_id;
    evidence.uncertainty_method = covariance.method;
    evidence.status = EvidenceStatus::Ineligible;
    if (stage.source.method != EstimateMethod::JointComponents)
    { evidence.reason = "not-a-joint-estimate"; return evidence; }
    if (stage.source.role != FittingRole::Target) { evidence.reason = "not-a-target"; return evidence; }
    if (!stage.point) { evidence.reason = "missing-point-state"; return evidence; }
    if (stage.convergence != JointCheckStatus::Passed)
    { evidence.reason = "joint-endpoint-not-converged"; return evidence; }
    if (!(stage.point->GetAmplitude() > 0.0)) { evidence.reason = "nonpositive-amplitude"; return evidence; }
    if (!(stage.point->GetWidth() > 0.0)) { evidence.reason = "nonpositive-width"; return evidence; }
    evidence.estimate = linearization_service::EncodeGaussianToParameterVector(*stage.point);
    if (!evidence.estimate->allFinite())
    { evidence.estimate.reset(); evidence.reason = "nonfinite-encoding"; return evidence; }
    if (!peeling) { evidence.status = EvidenceStatus::Unavailable; evidence.reason = "peeling-unavailable"; return evidence; }
    if (!SameSource(stage.source, peeling->source) || !SameSource(stage.source, covariance.source))
    { evidence.status = EvidenceStatus::Unavailable; evidence.reason = "evidence-source-mismatch"; return evidence; }
    if (covariance.status != EvidenceStatus::Available || !covariance.covariance)
    { evidence.status = EvidenceStatus::Unavailable; evidence.reason = covariance.reason; return evidence; }
    if (!PositiveDefinite(*covariance.covariance))
    { evidence.status = EvidenceStatus::Unavailable; evidence.reason = "singular-parameter-covariance"; return evidence; }
    evidence.covariance = covariance.covariance;
    evidence.status = EvidenceStatus::Available;
    return evidence;
}

} // namespace rhbm_gem::core::detail
