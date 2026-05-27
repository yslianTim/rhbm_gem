#include <rhbm_gem/data/object/AtomLocalPotentialEditor.hpp>

#include "data/detail/LocalPotentialEntry.hpp"

#include <utility>

namespace rhbm_gem {

AtomLocalPotentialEditor::AtomLocalPotentialEditor(LocalPotentialEntry & entry) :
    m_entry{ &entry }
{
}

void AtomLocalPotentialEditor::SetSamplingEntries(LocalPotentialSampleList value)
{
    m_entry->SetSamplingEntries(std::move(value));
}

void AtomLocalPotentialEditor::SetGaussianResult(LocalGaussianResult value)
{
    m_entry->SetGaussianResult(std::move(value));
}

void AtomLocalPotentialEditor::SetAlphaR(double value)
{
    m_entry->SetAlphaR(value);
}

} // namespace rhbm_gem
