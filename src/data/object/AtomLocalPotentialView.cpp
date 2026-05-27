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

const LocalGaussianResult & AtomLocalPotentialView::GetGaussianResult() const
{
    return RequireEntry("Local Gaussian result").GaussianResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS() const
{
    return RequireEntry("Local estimate OLS").GaussianResult().ols.GetModel();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE() const
{
    return RequireEntry("Local estimate MDPDE").GaussianResult().mdpde.GetModel();
}

LocalPotentialSampleList AtomLocalPotentialView::GetSamplingEntries(bool apply_selection) const
{
    const auto & sampling_entries{ RequireEntry("Local sampling entries").SamplingEntries() };
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

double AtomLocalPotentialView::GetAlphaR() const
{
    return RequireEntry("Local alpha-r").GaussianResult().alpha_r;
}

std::optional<LocalPotentialAnnotation> AtomLocalPotentialView::FindAnnotation(
    const std::string & key) const
{
    const auto * annotation{ RequireEntry("Local annotation").FindAnnotation(key) };
    if (annotation == nullptr)
    {
        return std::nullopt;
    }
    return LocalPotentialAnnotation{
        annotation->gaussian,
        annotation->is_outlier,
        annotation->statistical_distance
    };
}

} // namespace rhbm_gem
