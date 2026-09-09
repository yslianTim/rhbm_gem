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

const LocalGaussianResult & AtomLocalPotentialView::GetGaussianResult(FittingStage stage) const
{
    return RequireEntry("Local Gaussian result").GaussianResult(stage);
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS(FittingStage stage) const
{
    return RequireEntry("Local estimate OLS").GaussianResult(stage).ols.GetModel();
}

const std::optional<GroupGaussianMemberResult> & AtomLocalPotentialView::GetGroupMemberResult() const
{
    return RequireEntry("Group Gaussian member result").GroupMemberResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE(FittingStage stage) const
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

std::optional<double> AtomLocalPotentialView::GetLocalFittingPeelingRatio(
    bool peeling_applied,
    double distance_min,
    double distance_max) const
{
    numeric_validation::RequireFiniteNonNegativeRange(
        distance_min,
        distance_max,
        "peeling ratio distance range");
    if (!peeling_applied) return std::nullopt;

    double raw_sum{ 0.0 };
    std::size_t raw_sample_count{ 0 };
    for (const auto & sample : GetRawSamplingEntries(false))
    {
        if (sample.point.distance < distance_min || sample.point.distance > distance_max) continue;
        raw_sum += sample.response;
        ++raw_sample_count;
    }
    double peeling_sum{ 0.0 };
    std::size_t peeling_sample_count{ 0 };
    for (const auto & sample : GetPeelingSamplingEntries(false))
    {
        if (sample.point.distance < distance_min || sample.point.distance > distance_max) continue;
        peeling_sum += sample.response;
        ++peeling_sample_count;
    }
    if (raw_sample_count == 0 || peeling_sample_count == 0 || raw_sum == 0.0) return std::nullopt;

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
