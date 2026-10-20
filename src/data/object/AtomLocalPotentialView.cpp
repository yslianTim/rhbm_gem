#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>

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

using SampleIdentity = std::tuple<double, double, double, double, bool>;

struct ResponseAggregate
{
    std::size_t count{};
    long double response_sum{};
    bool responses_finite{ true };
};

bool IsWithinRange(const SamplingPoint & point, double distance_min, double distance_max)
{
    return std::isfinite(point.distance) && point.distance >= distance_min &&
        point.distance <= distance_max;
}

bool HasFinitePosition(const SamplingPoint & point)
{
    return std::all_of(point.position.begin(), point.position.end(),
        [](double coordinate) { return std::isfinite(coordinate); });
}

SampleIdentity MakeSampleIdentity(const SamplingPoint & point)
{
    return { point.distance, point.position[0], point.position[1], point.position[2],
        point.is_selected };
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

const LocalStageEstimate & AtomLocalPotentialView::GetStageEstimate(FittingStage stage) const
{
    return RequireEntry("Stage estimate").StageEstimate(stage);
}

bool AtomLocalPotentialView::HasFinalModel(FittingStage stage) const
{
    return IsAvailable() && GetStageEstimate(stage).point.has_value();
}

const GaussianModel3D & AtomLocalPotentialView::GetFinalModel(FittingStage stage) const
{
    const auto & estimate = GetStageEstimate(stage);
    if (!estimate.point) throw std::runtime_error("Stage estimate unavailable: " + estimate.reason);
    return *estimate.point;
}

const std::optional<PostFitPeelingResult> & AtomLocalPotentialView::GetPostFitPeeling() const
{
    return RequireEntry("Post-fit peeling").PostFitPeeling();
}

const std::optional<PostJointLocalResult> & AtomLocalPotentialView::GetPostJointLocalResult() const
{
    return RequireEntry("Post-Joint local result").PostJointLocalResult();
}

bool AtomLocalPotentialView::HasSampleGeometry() const
{
    return IsAvailable() && RequireEntry("Sample geometry").SampleGeometryAvailable();
}


LocalGaussianResult AtomLocalPotentialView::GetGaussianResult(FittingStage stage) const
{
    if (GetStageEstimate(stage).source.method == EstimateMethod::JointComponents)
        throw std::runtime_error("OLS/MDPDE diagnostics unavailable for Joint stage.");
    return RequireEntry("Local Gaussian result").GaussianResult(stage);
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS(FittingStage stage) const
{
    if (GetStageEstimate(stage).source.method == EstimateMethod::JointComponents)
        throw std::runtime_error("OLS diagnostics unavailable for Joint stage.");
    return RequireEntry("Local Gaussian result").OLSModel(stage);
}

const std::optional<GroupGaussianMemberResult> & AtomLocalPotentialView::GetGroupMemberResult() const
{
    return RequireEntry("Group Gaussian member result").GroupMemberResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE(FittingStage stage) const
{
    if (GetStageEstimate(stage).source.method == EstimateMethod::JointComponents)
        throw std::runtime_error("MDPDE diagnostics unavailable for Joint stage.");
    return RequireEntry("Local Gaussian result").MDPDEModel(stage);
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

LocalPotentialSampleList AtomLocalPotentialView::GetSamplingEntries(FittingStage stage) const
{
    switch (stage)
    {
        case FittingStage::First:
            return GetRawSamplingEntries(true);
        case FittingStage::Second:
            return GetPeelingSamplingEntries(false);
    }
    throw std::invalid_argument("Unknown local fitting stage.");
}

bool AtomLocalPotentialView::HasEnoughSamplingEntriesInRange(
    FittingStage stage,
    double distance_min,
    double distance_max,
    std::size_t minimum_sample_count) const
{
    std::size_t count{ 0 };
    for (const auto & sample : GetSamplingEntries(stage))
    {
        if (sample.point.distance < distance_min || sample.point.distance > distance_max) continue;
        count++;
        if (count >= minimum_sample_count) return true;
    }
    return false;
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
        if (raw.size() != peeling->samples.size())
        {
            for (const auto & sample : raw)
                if (!(sample.point.distance < distance_min || sample.point.distance > distance_max))
                    ++assessment.total_samples;
            assessment.reason = "raw-peeling-sample-count-mismatch";
            return assessment;
        }
        long double raw_sum{};
        long double valid_raw_sum{};
        long double peeled_sum{};
        std::string missing_reason;
        for (std::size_t i = 0; i < raw.size(); ++i)
        {
            if (!IsWithinRange(raw[i].point, distance_min, distance_max)) continue;
            ++assessment.total_samples;
            if (std::isfinite(raw[i].response)) raw_sum += raw[i].response;
            if (!peeling->samples[i].response)
            {
                if (missing_reason.empty()) missing_reason = peeling->samples[i].reason;
                continue;
            }
            if (!std::isfinite(raw[i].response) || !std::isfinite(*peeling->samples[i].response))
            {
                if (missing_reason.empty()) missing_reason = "nonfinite-response";
                continue;
            }
            valid_raw_sum += raw[i].response;
            peeled_sum += *peeling->samples[i].response;
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

    const auto & stage = GetStageEstimate(FittingStage::Second);
    if (stage.source.method == EstimateMethod::JointComponents)
    {
        assessment.reason = !stage.point || stage.convergence != JointCheckStatus::Passed ?
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

    std::map<SampleIdentity, ResponseAggregate> raw_by_identity;
    std::map<SampleIdentity, ResponseAggregate> peeling_by_identity;
    long double raw_sum{};
    long double peeled_sum{};
    std::size_t peeling_tail_count{};
    bool invalid_sample_geometry{};
    bool nonfinite_response{};
    for (const auto & sample : raw)
    {
        if (!IsWithinRange(sample.point, distance_min, distance_max)) continue;
        ++assessment.total_samples;
        if (std::isfinite(sample.response)) raw_sum += sample.response;
        else nonfinite_response = true;
        if (!HasFinitePosition(sample.point))
        {
            invalid_sample_geometry = true;
            continue;
        }
        auto & aggregate = raw_by_identity[MakeSampleIdentity(sample.point)];
        ++aggregate.count;
        if (std::isfinite(sample.response)) aggregate.response_sum += sample.response;
        else aggregate.responses_finite = false;
    }
    const auto peeling = GetPeelingSamplingEntries(false);
    for (const auto & sample : peeling)
    {
        if (!IsWithinRange(sample.point, distance_min, distance_max)) continue;
        ++peeling_tail_count;
        if (std::isfinite(sample.response)) peeled_sum += sample.response;
        else nonfinite_response = true;
        if (std::isfinite(sample.response)) ++assessment.peeling_samples;
        if (!HasFinitePosition(sample.point))
        {
            invalid_sample_geometry = true;
            continue;
        }
        auto & aggregate = peeling_by_identity[MakeSampleIdentity(sample.point)];
        ++aggregate.count;
        if (std::isfinite(sample.response)) aggregate.response_sum += sample.response;
        else aggregate.responses_finite = false;
    }
    if (!StoreResponseSum(raw_sum, assessment.raw_response_sum) ||
        !StoreResponseSum(peeled_sum, assessment.peeled_response_sum))
    {
        assessment.reason = "nonfinite-response-sum";
        return assessment;
    }
    if (!assessment.total_samples)
    {
        assessment.reason = "no-tail-samples";
        return assessment;
    }
    if (!peeling_tail_count)
    {
        assessment.reason = "no-peeling-samples";
        return assessment;
    }
    if (!HasSampleGeometry())
    {
        assessment.reason = "sample-geometry-unavailable";
        return assessment;
    }
    if (invalid_sample_geometry)
    {
        assessment.reason = "nonfinite-sample-geometry";
        return assessment;
    }

    bool identity_mismatch = raw_by_identity.size() != peeling_by_identity.size();
    long double valid_raw_sum{};
    long double valid_peeled_sum{};
    for (const auto & [identity, raw_aggregate] : raw_by_identity)
    {
        const auto peeled_aggregate = peeling_by_identity.find(identity);
        if (peeled_aggregate == peeling_by_identity.end() ||
            peeled_aggregate->second.count != raw_aggregate.count)
        {
            identity_mismatch = true;
            continue;
        }
        if (!raw_aggregate.responses_finite || !peeled_aggregate->second.responses_finite)
        {
            nonfinite_response = true;
            continue;
        }
        assessment.valid_samples += raw_aggregate.count;
        valid_raw_sum += raw_aggregate.response_sum;
        valid_peeled_sum += peeled_aggregate->second.response_sum;
    }
    if (!StoreResponseSum(valid_raw_sum, assessment.valid_raw_response_sum))
    {
        assessment.reason = "nonfinite-response-sum";
        return assessment;
    }
    if (identity_mismatch)
    {
        assessment.reason = "peeling-sample-identity-mismatch";
        return assessment;
    }
    if (nonfinite_response || assessment.valid_samples != assessment.total_samples)
    {
        assessment.reason = "nonfinite-response";
        return assessment;
    }
    if (assessment.valid_raw_response_sum == 0.0)
    {
        assessment.reason = "invalid-zero-denominator";
        return assessment;
    }

    const long double ratio{ 1.0L - valid_peeled_sum / valid_raw_sum };
    if (!std::isfinite(ratio) || std::abs(ratio) > std::numeric_limits<double>::max())
    {
        assessment.reason = "nonfinite-ratio";
        return assessment;
    }
    assessment.ratio = static_cast<double>(ratio);
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

double AtomLocalPotentialView::GetAlphaR(FittingStage stage) const
{
    if (GetStageEstimate(stage).source.method == EstimateMethod::JointComponents)
        throw std::runtime_error("MDPDE alpha unavailable for Joint stage.");
    return RequireEntry("Local Gaussian result").AlphaR(stage);
}

} // namespace rhbm_gem
