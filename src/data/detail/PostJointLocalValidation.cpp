#include "PostJointLocalValidation.hpp"
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <numbers>
#include <set>
#include <utility>

namespace rhbm_gem::data_internal {
namespace {
bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role && a.fixed_from_first_stage==b.fixed_from_first_stage;
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

std::string ValidatePostFitPeelingResult(const PostFitPeelingResult & peeling,
    const LocalStageEstimate & joint, const JointAnalysisResult * snapshot,
    const std::array<double,3> * target_position)
{
    if (peeling.mode!="original-grid-v1" || !SameSource(peeling.source,joint.source) ||
        joint.source.method!=EstimateMethod::JointComponents || joint.source.role!=FittingRole::Target)
        return "Post-Joint grid source/contract mismatch.";
    std::set<std::size_t> grids;
    std::optional<std::size_t> previous;
    for (const auto & sample:peeling.samples)
    {
        const auto & observation=sample.observation;
        if ((previous && observation.original_row_index<=*previous) ||
            !grids.insert(observation.global_grid_index).second ||
            !std::isfinite(observation.squared_distance) || observation.squared_distance<0 ||
            observation.squared_distance>4.0 || !std::isfinite(observation.raw_response) ||
            !std::all_of(observation.position.begin(),observation.position.end(),
                [](double value){return std::isfinite(value);}))
            return "Invalid Post-Joint grid identity/order/observation.";
        previous=observation.original_row_index;
        if ((sample.response && (!std::isfinite(*sample.response) || !sample.reason.empty())) ||
            (!sample.response && sample.reason.empty()))
            return "Invalid Post-Joint available/unavailable response pairing.";
        if (snapshot && (observation.original_row_index>=snapshot->row_ids.size() ||
            snapshot->row_ids[observation.original_row_index]!=std::to_string(observation.global_grid_index)))
            return "Post-Joint Original row mapping mismatch.";
        if (snapshot && std::all_of(snapshot->metadata.grid_size.begin(),snapshot->metadata.grid_size.end(),
                [](int size){return size>0;}))
        {
            auto grid=observation.global_grid_index;
            for(std::size_t axis=0;axis<3;++axis)
            {
                const auto size=static_cast<std::size_t>(snapshot->metadata.grid_size[axis]);
                const auto coordinate=grid%size; grid/=size;
                const double position=snapshot->metadata.origin[axis]+
                    snapshot->metadata.grid_spacing[axis]*static_cast<double>(coordinate);
                if (!std::isfinite(position) || std::abs(position-observation.position[axis])>
                    1.0e-12*std::max(1.0,std::abs(position)))
                    return "Post-Joint grid index/position mismatch.";
            }
            if(grid) return "Post-Joint grid index outside map geometry.";
        }
        if (target_position)
        {
            double square{};
            for (std::size_t axis=0;axis<3;++axis)
            {
                const double delta=observation.position[axis]-(*target_position)[axis];
                square+=delta*delta;
            }
            if (std::abs(square-observation.squared_distance)>1.0e-12*std::max(1.0,square))
                return "Post-Joint target/grid geometry mismatch.";
        }
    }
    return {};
}

PostJointLocalEvidenceAssessment AssessPostJointLocalResult(
    const PostJointLocalResult & result, int atom_id, const LocalStageEstimate & joint,
    const PostFitPeelingResult * peeling)
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
    if (result.effective_sample_count > result.positive_sample_count ||
        result.positive_sample_count > result.corrected_sample_count ||
        result.corrected_sample_count > result.sample_count ||
        !std::isfinite(result.fitting_radius) || result.fitting_radius<=0 || result.fitting_radius>2.0)
        return Invalid("Post-Joint local sample counts/range mismatch.");
    if (peeling)
    {
        const auto reason=ValidatePostFitPeelingResult(*peeling,joint);
        if (!reason.empty()) return Invalid(reason);
        if (result.sample_count!=peeling->samples.size())
            return Invalid("Post-Joint local/grid sample count mismatch.");
    }
    else if (result.sample_count || result.status==EvidenceStatus::Available)
        return Invalid("Post-Joint local evidence requires current grid peeling.");

    std::set<std::size_t> failure_indices;
    for (const auto & failure : result.sample_failures)
    {
        if (!peeling || !failure_indices.insert(failure.original_row_index).second ||
            !std::isfinite(failure.distance) || failure.reason.empty())
            return Invalid("Invalid post-Joint local sample failure pairing.");
        const auto found=std::lower_bound(peeling->samples.begin(),peeling->samples.end(),failure.original_row_index,
            [](const auto & sample,std::size_t row){return sample.observation.original_row_index<row;});
        if (found==peeling->samples.end() || found->observation.original_row_index!=failure.original_row_index ||
            found->observation.global_grid_index!=failure.global_grid_index ||
            std::sqrt(found->observation.squared_distance)!=failure.distance)
            return Invalid("Invalid post-Joint local sample failure identity.");
    }
    if (result.status==EvidenceStatus::Available)
    {
        std::size_t corrected{},positive{},retained{};
        for(const auto & sample:peeling->samples)
        {
            const auto & observation=sample.observation;
            const double distance=std::sqrt(observation.squared_distance);
            // Independent scalar replay of the fixed charge; C is not an offset.
            const double width=result.fixed_charge.width;
            if(!std::isfinite(width) || width<=0 || !std::isfinite(result.fixed_charge.coefficient))
                return Invalid("Invalid Post-Joint fixed charge profile.");
            const double kernel=distance<1e-5 ? std::sqrt(2.0/std::numbers::pi)/width :
                std::erf(distance/width/std::sqrt(2.0))/distance;
            const double charge=result.fixed_charge.coefficient*kernel;
            const double response=sample.response ? *sample.response-charge :
                std::numeric_limits<double>::quiet_NaN();
            const bool failed=!std::isfinite(response);
            if(failure_indices.contains(observation.original_row_index)!=failed)
                return Invalid("Post-Joint correction/failure identity mismatch.");
            if(failed) continue;
            ++corrected;
            if(response>0)
            {
                ++positive;
                if(distance<=result.fitting_radius) ++retained;
            }
        }
        if(corrected!=result.corrected_sample_count || positive!=result.positive_sample_count ||
            retained!=result.effective_sample_count)
            return Invalid("Post-Joint corrected/positive/retained dataset mismatch.");
    }

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
            (!std::isfinite(model->GetOffset()) || model->GetOffset()!=0.0))
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
