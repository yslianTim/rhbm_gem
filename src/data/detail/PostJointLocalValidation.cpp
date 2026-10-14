#include "PostJointLocalValidation.hpp"
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <set>
#include <utility>

namespace rhbm_gem::data_internal {
namespace {
bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role;
}

bool PositiveDefinite(const Eigen::MatrixXd & matrix)
{
    if (matrix.rows() != 2 || matrix.cols() != 2 || !matrix.allFinite() ||
        !matrix.isApprox(matrix.transpose(), 1.0e-12)) return false;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(matrix);
    return spectrum.info() == Eigen::Success && spectrum.eigenvalues()(0) >
        std::numeric_limits<double>::epsilon() * static_cast<double>(matrix.rows()) *
            spectrum.eigenvalues().maxCoeff();
}

bool FiniteRefinement(const std::optional<RHBMEndpointRefinementDiagnostics> & refinement)
{
    if (!refinement) return true;
    const auto finite = [](const std::optional<double> & value) {
        return !value || std::isfinite(*value);
    };
    if (refinement->candidate_equation_evaluations < 0 || refinement->reference_updates < 0 ||
        !finite(refinement->original_residual) || !finite(refinement->candidate_residual) ||
        !finite(refinement->reference_residual) || !finite(refinement->weight_max_difference))
        return false;
    return !refinement->relative_coordinate_difference ||
        std::all_of(refinement->relative_coordinate_difference->begin(),
            refinement->relative_coordinate_difference->end(), [](double value) { return std::isfinite(value); });
}

PostJointLocalEvidenceAssessment Invalid(std::string reason)
{
    return {false, false, std::move(reason)};
}
}

PostJointLocalEvidenceAssessment AssessPostJointLocalResult(
    const PostJointLocalResult & result, int atom_id, const LocalStageEstimate & joint,
    const PostFitPeelingResult * peeling, const LocalPotentialSampleList & raw_samples)
{
    const auto atom_identity = std::to_string(atom_id);
    if (joint.source.method != EstimateMethod::JointComponents || joint.source.role != FittingRole::Target ||
        joint.source.atom_id != atom_identity || joint.source.component_id.empty() || joint.source.run_id.empty() ||
        !SameSource(result.source, joint.source))
        return Invalid("Post-Joint local result source mismatch.");
    if (result.status == EvidenceStatus::NotRun)
        return Invalid("Post-Joint local result status was not recorded.");
    if ((result.status == EvidenceStatus::Available && !result.reason.empty()) ||
        (result.status != EvidenceStatus::Available && result.reason.empty()))
        return Invalid("Post-Joint local result status/reason mismatch.");
    if (result.sample_count != raw_samples.size() || result.effective_sample_count > result.sample_count)
        return Invalid("Post-Joint local sample count mismatch.");
    if (peeling)
    {
        if (!SameSource(peeling->source, result.source) || peeling->samples.size() != raw_samples.size())
            return Invalid("Post-Joint peeling/raw sample pairing mismatch.");
        for (const auto & sample : peeling->samples)
            if (sample.response && !std::isfinite(*sample.response))
                return Invalid("Post-Joint peeling contains a nonfinite response.");
    }
    if (result.status == EvidenceStatus::Available && !peeling)
        return Invalid("Available post-Joint local evidence requires current peeling.");

    std::set<std::size_t> failure_indices;
    for (const auto & failure : result.sample_failures)
        if (failure.raw_sample_index >= raw_samples.size() ||
            !failure_indices.insert(failure.raw_sample_index).second ||
            !std::isfinite(failure.distance) || failure.distance != raw_samples[failure.raw_sample_index].point.distance ||
            failure.reason.empty())
            return Invalid("Invalid post-Joint local sample failure pairing.");

    if (joint.point)
    {
        if (!std::isfinite(joint.point->GetOffset()) || !std::isfinite(joint.point->GetWidth()) ||
            joint.point->GetWidth() <= 0.0 ||
            result.fixed_charge.coefficient != joint.point->GetOffset() ||
            result.fixed_charge.width != joint.point->GetWidth())
            return Invalid("Post-Joint fixed charge profile does not match its Joint endpoint.");
    }

    if (result.status != EvidenceStatus::Available)
    {
        if (result.information)
            return Invalid("Unavailable post-Joint local evidence contains group measurement information.");
        return {true, false, result.reason};
    }

    const auto & point = result.fit.mdpde.GetModel();
    const auto & ols = result.fit.ols.GetModel();
    const auto & diagnostics = result.fit.diagnostics;
    const auto & information = result.information;
    if (!joint.point || joint.convergence != JointCheckStatus::Passed)
        return Invalid("Available post-Joint local evidence requires a qualified Joint endpoint.");
    if (!std::isfinite(result.fixed_charge.coefficient) || !std::isfinite(result.fixed_charge.width) ||
        result.fixed_charge.width <= 0.0)
        return Invalid("Available post-Joint local evidence has an incomplete fixed charge profile.");
    if (!std::isfinite(result.fit.alpha_r) || result.fit.alpha_r < 0.0)
        return Invalid("Available post-Joint local evidence has invalid alpha.");
    for (const auto * model : {&point, &ols})
        if (!std::isfinite(model->GetAmplitude()) || model->GetAmplitude() <= 0.0 ||
            !std::isfinite(model->GetWidth()) || model->GetWidth() <= 0.0 ||
            !std::isfinite(model->GetOffset()))
            return Invalid("Available post-Joint local evidence has invalid Gaussian parameters.");
    if (!diagnostics || diagnostics->status != RHBMEstimationStatus::SUCCESS ||
        GetSolveQualification(diagnostics->status, diagnostics->refinement) == RHBMSolveQualification::Unqualified ||
        !FiniteRefinement(diagnostics->refinement) ||
        !std::isfinite(diagnostics->sigma_square) || diagnostics->sigma_square <= 0.0 ||
        diagnostics->sigma_square == std::numeric_limits<double>::max() || diagnostics->iterations.iterations < 0 ||
        (diagnostics->iterations.squared_beta_change &&
            (!std::isfinite(*diagnostics->iterations.squared_beta_change) ||
             *diagnostics->iterations.squared_beta_change < 0.0)) ||
        (diagnostics->iterations.relative_variance_change &&
            (!std::isfinite(*diagnostics->iterations.relative_variance_change) ||
             *diagnostics->iterations.relative_variance_change < 0.0)))
        return Invalid("Available post-Joint local evidence lacks qualified diagnostics.");
    if (result.effective_sample_count <= 2 || !information ||
        information->moment.size() != 2 || !information->moment.allFinite() ||
        !PositiveDefinite(information->precision))
        return Invalid("Available post-Joint local evidence lacks positive-definite measurement information.");

    try
    {
        if (!linearization_service::EncodeGaussianToParameterVector(point).allFinite())
            return Invalid("Available post-Joint local beta cannot be encoded for group inference.");
    }
    catch (const std::exception &)
    {
        return Invalid("Available post-Joint local beta cannot be encoded for group inference.");
    }
    return {true, true, {}};
}
} // namespace rhbm_gem::data_internal
