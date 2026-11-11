#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
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

LocalPotentialSampleList ApplySamplingEntrySelection(
    const LocalPotentialSampleList & sampling_entries,
    bool apply_selection)
{
    if (!apply_selection)
    {
        return sampling_entries;
    }

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

bool IsWithinRange(const SamplingPoint & point, double distance_min, double distance_max)
{
    return std::isfinite(point.distance) && point.distance >= distance_min &&
        point.distance <= distance_max;
}

bool StoreResponseSum(long double sum, double & value)
{
    if (!std::isfinite(sum) || std::abs(sum) > std::numeric_limits<double>::max()) return false;
    value = static_cast<double>(sum);
    return std::isfinite(value);
}

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

LocalPotentialSampleList AtomLocalPotentialView::GetRawSamplingEntries(bool apply_selection) const
{
    const auto & entry{ RequireEntry("Local raw sampling entries") };
    return ApplySamplingEntrySelection(entry.RawSamplingEntries(), apply_selection);
}

const LocalPotentialSampleList & AtomLocalPotentialView::GetRawSamplingEntriesView() const
{
    return RequireEntry("Local raw sampling entries").RawSamplingEntries();
}

LocalPotentialSampleList AtomLocalPotentialView::GetPeelingSamplingEntries(bool apply_selection) const
{
    const auto & entry{ RequireEntry("Local peeling sampling entries") };
    return ApplySamplingEntrySelection(entry.PeelingSamplingEntries(), apply_selection);
}

PeelingRatioAssessment AtomLocalPotentialView::AssessLocalFittingPeelingRatio(
    double distance_min,
    double distance_max) const
{
    numeric_validation::RequireFiniteNonNegativeRange(
        distance_min,
        distance_max,
        "peeling ratio distance range");

    PeelingRatioAssessment assessment;
    const auto & raw = GetRawSamplingEntriesView();
    if (const auto & peeling = GetPostFitPeeling(); peeling)
    {
        long double raw_sum{};
        long double valid_raw_sum{};
        long double peeled_sum{};
        std::string missing_reason;
        for (const auto & sample:peeling->samples)
        {
            const double distance=std::sqrt(sample.observation.squared_distance);
            if (distance<distance_min || distance>distance_max) continue;
            ++assessment.total_samples;
            if (std::isfinite(sample.observation.raw_response)) raw_sum+=sample.observation.raw_response;
            if (!sample.response)
            {
                if (missing_reason.empty()) missing_reason=sample.reason;
                continue;
            }
            if (!std::isfinite(sample.observation.raw_response) || !std::isfinite(*sample.response))
            {
                if (missing_reason.empty()) missing_reason="nonfinite-response";
                continue;
            }
            valid_raw_sum+=sample.observation.raw_response;
            peeled_sum+=*sample.response;
            ++assessment.valid_samples;
        }
        if (!StoreResponseSum(raw_sum, assessment.raw_response_sum) ||
            !StoreResponseSum(valid_raw_sum, assessment.valid_raw_response_sum) ||
            !StoreResponseSum(peeled_sum, assessment.peeled_response_sum))
        {
            assessment.reason = "nonfinite-response-sum";
            return assessment;
        }
        assessment.peeling_samples = assessment.valid_samples;
        if (!assessment.total_samples)
        {
            assessment.reason = "no-tail-samples";
            return assessment;
        }
        if (assessment.valid_samples != assessment.total_samples)
        {
            assessment.reason = "tail-coverage-incomplete:" +
                (missing_reason.empty() ? std::string("peeling-unavailable") : missing_reason);
            return assessment;
        }
        if (assessment.valid_raw_response_sum == 0.0)
        {
            assessment.reason = "invalid-zero-denominator";
            return assessment;
        }
        const long double ratio = 1.0L - peeled_sum / valid_raw_sum;
        if (!std::isfinite(ratio) || std::abs(ratio) > std::numeric_limits<double>::max())
        {
            assessment.reason = "nonfinite-ratio";
            return assessment;
        }
        assessment.ratio = static_cast<double>(ratio);
        return assessment;
    }

    const auto & stage = GetJointStageEstimate();
    assessment.reason = stage.source.method == EstimateMethod::JointComponents &&
        (!stage.point || stage.convergence != JointCheckStatus::Passed) ?
        "joint-qualification-failed" : "missing-peeling-result";
    long double raw_sum{};
    for (const auto & sample : raw)
    {
        if (!IsWithinRange(sample.point, distance_min, distance_max)) continue;
        ++assessment.total_samples;
        if (std::isfinite(sample.response)) raw_sum += sample.response;
    }
    if (!StoreResponseSum(raw_sum, assessment.raw_response_sum))
        assessment.reason = "nonfinite-response-sum";
    return assessment;
}

std::optional<double> AtomLocalPotentialView::GetLocalFittingPeelingRatio(
    double distance_min,
    double distance_max) const
{
    return AssessLocalFittingPeelingRatio(distance_min, distance_max).ratio;
}

int AtomLocalPotentialView::GetNeighborCountForPeeling() const
{
    return RequireEntry("Local peeling neighbor count").NeighborCountForPeeling();
}

double AtomLocalPotentialView::GetInitialStageAlphaR() const
{
    return RequireInitialStageEntry("MDPDE alpha unavailable for Joint stage.").InitialStageAlphaR();
}

} // namespace rhbm_gem
