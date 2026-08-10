#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>

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

AtomLocalPotentialView AtomLocalPotentialView::RequireFor(const AtomObject & atom_object)
{
    auto view{ AtomLocalPotentialView::For(atom_object) };
    (void)view.RequireEntry("Atom local analysis");
    return view;
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
    LocalFittingStage stage) const
{
    return RequireEntry("Local Gaussian result").GaussianResult(stage);
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS(
    LocalFittingStage stage) const
{
    return RequireEntry("Local estimate OLS").GaussianResult(stage).ols.GetModel();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE(
    LocalFittingStage stage) const
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

int AtomLocalPotentialView::GetNeighborCountForPeeling() const
{
    return RequireEntry("Local peeling neighbor count").NeighborCountForPeeling();
}

double AtomLocalPotentialView::GetAlphaR(LocalFittingStage stage) const
{
    return RequireEntry("Local alpha-r").GaussianResult(stage).alpha_r;
}

} // namespace rhbm_gem
