#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

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

const LocalGaussianResult & AtomLocalPotentialView::GetGaussianResult(
    FittingStage stage) const
{
    return RequireEntry("Local Gaussian result").GaussianResult(stage);
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS(
    FittingStage stage) const
{
    return RequireEntry("Local estimate OLS").GaussianResult(stage).ols.GetModel();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE(
    FittingStage stage) const
{
    return RequireEntry("Local estimate MDPDE").GaussianResult(stage).mdpde.GetModel();
}

LocalPotentialSampleList AtomLocalPotentialView::GetRawSamplingEntries(bool apply_selection) const
{
    const auto & entry{ RequireEntry("Local raw sampling entries") };
    return ApplySamplingEntrySelection(entry.RawSamplingEntries(), apply_selection);
}

LocalPotentialSampleList AtomLocalPotentialView::GetPeelingSamplingEntries(bool apply_selection) const
{
    const auto & entry{ RequireEntry("Local peeling sampling entries") };
    return ApplySamplingEntrySelection(entry.PeelingSamplingEntries(), apply_selection);
}

LocalPotentialSampleList AtomLocalPotentialView::GetSamplingEntries(
    FittingStage stage) const
{
    switch (stage)
    {
        case FittingStage::First:
            return GetRawSamplingEntries();
        case FittingStage::Second:
        case FittingStage::Third:
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
    const auto & entry{ RequireEntry("Local fitting samples") };
    const LocalPotentialSampleList * sample_entries{ nullptr };
    bool apply_selection{ false };
    switch (stage)
    {
        case FittingStage::First:
            sample_entries = &entry.RawSamplingEntries();
            apply_selection = true;
            break;
        case FittingStage::Second:
        case FittingStage::Third:
            sample_entries = &entry.PeelingSamplingEntries();
            break;
        default:
            throw std::invalid_argument("Unknown local fitting stage.");
    }

    std::size_t count{ 0 };
    for (const auto & sample : *sample_entries)
    {
        if (apply_selection && !sample.point.is_selected) continue;
        if (sample.point.distance < distance_min || sample.point.distance > distance_max)
        {
            continue;
        }
        count++;
        if (count >= minimum_sample_count)
        {
            return true;
        }
    }
    return false;
}

std::optional<double> AtomLocalPotentialView::GetLocalFittingPeelingRatio(
    bool peeling_applied,
    double distance_min,
    double distance_max) const
{
    numeric_validation::RequireFiniteNonNegativeRange(
        distance_min,
        distance_max,
        "peeling ratio distance range");
    if (!peeling_applied)
    {
        return std::nullopt;
    }

    const auto raw_sampling_entries{ GetRawSamplingEntries(false) };
    const auto peeling_sampling_entries{ GetPeelingSamplingEntries(false) };
    double raw_sum{ 0.0 };
    std::size_t raw_sample_count{ 0 };
    for (const auto & sample : raw_sampling_entries)
    {
        if (sample.point.distance < distance_min
            || sample.point.distance > distance_max
            || !std::isfinite(sample.point.distance))
        {
            continue;
        }
        raw_sum += sample.response;
        ++raw_sample_count;
    }
    double peeling_sum{ 0.0 };
    std::size_t peeling_sample_count{ 0 };
    for (const auto & sample : peeling_sampling_entries)
    {
        if (sample.point.distance < distance_min
            || sample.point.distance > distance_max
            || !std::isfinite(sample.point.distance))
        {
            continue;
        }
        peeling_sum += sample.response;
        ++peeling_sample_count;
    }
    if (raw_sample_count == 0
        || peeling_sample_count == 0
        || !std::isfinite(raw_sum)
        || !std::isfinite(peeling_sum)
        || raw_sum == 0.0)
    {
        return std::nullopt;
    }

    const auto ratio{ (raw_sum - peeling_sum) / raw_sum };
    return std::isfinite(ratio) ? std::optional<double>{ ratio } : std::nullopt;
}

int AtomLocalPotentialView::GetNeighborCountForPeeling() const
{
    return RequireEntry("Local peeling neighbor count").NeighborCountForPeeling();
}

double AtomLocalPotentialView::GetAlphaR(FittingStage stage) const
{
    return RequireEntry("Local alpha-r").GaussianResult(stage).alpha_r;
}

} // namespace rhbm_gem
