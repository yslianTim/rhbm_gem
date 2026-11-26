#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rhbm_gem {

namespace {

const LocalPotentialEntry & RequireLocalEntry(
    const LocalPotentialEntry * entry,
    const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
}

const GaussianModel3D & RequireStageModel(const LocalStageEstimate & estimate)
{
    if (!estimate.point) throw std::runtime_error("Stage estimate unavailable: " + estimate.reason);
    return *estimate.point;
}

// Accumulate in units of a power of two, independently for each component.
// The significand sum can grow only with the number of samples.
struct ScaledMagnitudeSum
{
    double sum{}, correction{}, largest{};
    int exponent{};

    void Add(double magnitude)
    {
        if (magnitude == 0.0) return;
        if (magnitude > largest)
        {
            int next_exponent{};
            std::frexp(magnitude, &next_exponent);
            sum = std::scalbn(sum, exponent - next_exponent);
            correction = std::scalbn(correction, exponent - next_exponent);
            exponent = next_exponent;
            largest = magnitude;
        }
        const double adjusted = std::scalbn(magnitude, -exponent) - correction;
        const double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
};

} // namespace

AtomLocalPotentialView::AtomLocalPotentialView(const AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

AtomLocalPotentialView AtomLocalPotentialView::For(const AtomObject & atom_object)
{
    return AtomLocalPotentialView(&atom_object);
}

bool AtomLocalPotentialView::IsAvailable() const
{
    return FindEntry() != nullptr;
}

const LocalPotentialEntry * AtomLocalPotentialView::FindEntry() const
{
    if (m_atom_object != nullptr && m_atom_object->m_owner_model != nullptr)
    {
        return ModelAnalysisData::Of(*m_atom_object->m_owner_model).FindAtomLocalEntry(*m_atom_object);
    }
    return nullptr;
}

const LocalPotentialEntry & AtomLocalPotentialView::RequireEntry(const char * context) const
{
    return RequireLocalEntry(FindEntry(), context);
}

const LocalStageEstimate & AtomLocalPotentialView::GetInitialStageEstimate() const
{
    return RequireEntry("Stage estimate").InitialStageEstimate();
}

bool AtomLocalPotentialView::HasInitialStageModel() const
{
    return IsAvailable() && GetInitialStageEstimate().point.has_value();
}

const GaussianModel3D & AtomLocalPotentialView::GetInitialStageModel() const
{
    return RequireStageModel(GetInitialStageEstimate());
}

const LocalStageEstimate & AtomLocalPotentialView::GetJointStageEstimate() const
{
    return RequireEntry("Stage estimate").JointStageEstimate();
}

bool AtomLocalPotentialView::HasJointStageModel() const
{
    return IsAvailable() && GetJointStageEstimate().point.has_value();
}

const GaussianModel3D & AtomLocalPotentialView::GetJointStageModel() const
{
    return RequireStageModel(GetJointStageEstimate());
}

const std::optional<PostFitPeelingResult> & AtomLocalPotentialView::GetPostFitPeeling() const
{
    return RequireEntry("Post-fit peeling").PostFitPeeling();
}

const std::optional<PostJointLocalResult> & AtomLocalPotentialView::GetPostJointLocalResult() const
{
    return RequireEntry("Post-Joint local result").PostJointLocalResult();
}

const std::optional<PostJointQScoreResult> & AtomLocalPotentialView::GetPostJointQScore() const
{
    return RequireEntry("Post-Joint QScore").PostJointQScore();
}

bool AtomLocalPotentialView::HasSampleGeometry() const
{
    return IsAvailable() && RequireEntry("Sample geometry").SampleGeometryAvailable();
}


const LocalPotentialEntry & AtomLocalPotentialView::RequireInitialStageEntry(const char * diagnostic_error) const
{
    const auto & entry = RequireEntry("Stage estimate");
    if (entry.InitialStageEstimate().source.method == EstimateMethod::JointComponents)
        throw std::runtime_error(diagnostic_error);
    return entry;
}

LocalGaussianResult AtomLocalPotentialView::GetInitialStageGaussianResult() const
{
    return RequireInitialStageEntry("OLS/MDPDE diagnostics unavailable for Joint stage.").InitialStageGaussianResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetInitialStageOLS() const
{
    return RequireInitialStageEntry("OLS diagnostics unavailable for Joint stage.").InitialStageOLSModel();
}

const std::optional<GroupGaussianMemberResult> & AtomLocalPotentialView::GetGroupMemberResult() const
{
    return RequireEntry("Group Gaussian member result").GroupMemberResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetInitialStageMDPDE() const
{
    return RequireInitialStageEntry("MDPDE diagnostics unavailable for Joint stage.").InitialStageMDPDEModel();
}

LocalPotentialSampleList AtomLocalPotentialView::GetSelectedSamplingEntries() const
{
    const auto & sampling_entries{ RequireEntry("Local sampling entries").SamplingEntries() };
    LocalPotentialSampleList selected_entries;
    selected_entries.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        if (sample.point.is_selected)
        {
            selected_entries.emplace_back(sample);
        }
    }
    return selected_entries;
}

const LocalPotentialSampleList & AtomLocalPotentialView::GetSamplingEntriesView() const
{
    return RequireEntry("Local sampling entries").SamplingEntries();
}

LocalPotentialSampleList AtomLocalPotentialView::GetPeeledObservationEntries() const
{
    return RequireEntry("Local peeled observation entries").PeeledObservationEntries();
}

PeelingRatioAssessment AtomLocalPotentialView::AssessLocalFittingPeelingRatio(PeelingRegion region) const
{
    PeelingRatioAssessment assessment;
    const auto & entry = RequireEntry("Post-fit peeling");
    const auto & peeling = entry.PostFitPeeling();
    const auto & stage = entry.JointStageEstimate();
    const bool joint_unqualified = stage.source.method == EstimateMethod::JointComponents &&
        (!stage.point || stage.convergence != JointCheckStatus::Passed);
    if (!peeling)
    {
        const auto & snapshot = ModelAnalysisData::Of(*m_atom_object->m_owner_model).joint_result;
        if (snapshot && !IsQualifiedJointAnalysisResult(*snapshot))
            assessment.reason = "joint-qualification-failed";
        else if (snapshot && snapshot->selection_domain &&
            snapshot->selection_domain->FittingCount() < snapshot->atom_ids.size() &&
            !snapshot->halo_median_reconstruction)
            assessment.reason = "halo-median-reconstruction-unavailable";
        else
            assessment.reason = joint_unqualified ? "joint-qualification-failed" : "missing-peeling-result";
        return assessment;
    }

    const auto bounds = GetPeelingRegionBounds(region);
    ScaledMagnitudeSum neighbor, residual;
    std::string missing_reason;
    for (const auto & sample : peeling->samples)
    {
        const double squared_distance = sample.observation.squared_distance;
        if (!std::isfinite(squared_distance) || squared_distance < 0.0)
        {
            assessment.reason = "invalid-squared-distance";
            return assessment;
        }
        if (!bounds.ContainsSquaredDistance(squared_distance)) continue;
        ++assessment.total_samples;
        const double raw = sample.observation.map_value;
        if (!std::isfinite(raw) || (sample.response && !std::isfinite(*sample.response)))
        {
            if (missing_reason.empty()) missing_reason = "nonfinite-response";
            continue;
        }
        if (!sample.response)
        {
            if (missing_reason.empty()) missing_reason = sample.reason;
            continue;
        }
        const double peeled = *sample.response;
        ++assessment.valid_samples;
        // Subtract like-signed values directly to preserve close differences.
        // With opposite signs, |raw - peeled| is the sum of the two magnitudes.
        if (std::signbit(raw) == std::signbit(peeled)) neighbor.Add(std::abs(raw - peeled));
        else
        {
            neighbor.Add(std::abs(raw));
            neighbor.Add(std::abs(peeled));
        }
        residual.Add(std::abs(peeled));
    }
    if (!assessment.total_samples)
        assessment.reason = "no-region-samples";
    else if (assessment.valid_samples != assessment.total_samples)
        assessment.reason = "peeling-coverage-incomplete:" +
            (missing_reason.empty() ? std::string("peeling-unavailable") : missing_reason);
    else if (joint_unqualified)
        assessment.reason = "joint-qualification-failed";
    else if (neighbor.sum == 0.0 && residual.sum == 0.0)
        assessment.reason = "zero-component-magnitude";
    if (!assessment.reason.empty()) return assessment;

    if (neighbor.sum == 0.0 || residual.sum == 0.0)
    {
        assessment.ratio = neighbor.sum == 0.0 ? 0.0 : 1.0;
        return assessment;
    }
    int neighbor_sum_exponent{}, residual_sum_exponent{};
    const double neighbor_fraction = std::frexp(neighbor.sum, &neighbor_sum_exponent);
    const double residual_fraction = std::frexp(residual.sum, &residual_sum_exponent);
    neighbor_sum_exponent += neighbor.exponent;
    residual_sum_exponent += residual.exponent;
    const int common_exponent = std::max(neighbor_sum_exponent, residual_sum_exponent);
    const double total = std::scalbn(neighbor_fraction, neighbor_sum_exponent - common_exponent) +
        std::scalbn(residual_fraction, residual_sum_exponent - common_exponent);
    // Scale the quotient last so a subnormal Ratio rounds only at the output scale.
    const double ratio = std::scalbn(neighbor_fraction / total, neighbor_sum_exponent - common_exponent);
    if (!std::isfinite(ratio) || ratio < 0.0 || ratio > 1.0)
    {
        assessment.reason = "peeling-ratio-numerical-failure";
        return assessment;
    }
    assessment.ratio = ratio;
    return assessment;
}

double AtomLocalPotentialView::GetInitialStageAlphaR() const
{
    return RequireInitialStageEntry("MDPDE alpha unavailable for Joint stage.").InitialStageAlphaR();
}

} // namespace rhbm_gem
