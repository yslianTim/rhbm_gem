#include <rhbm_gem/data/object/AtomLocalPotentialEditor.hpp>

#include "data/detail/LocalPotentialEntry.hpp"

#include <utility>

namespace rhbm_gem {

AtomLocalPotentialEditor::AtomLocalPotentialEditor(LocalPotentialEntry & entry) :
    m_entry{ &entry }
{
}

void AtomLocalPotentialEditor::SetRawSamplingEntries(LocalPotentialSampleList value)
{
    m_entry->SetRawSamplingEntries(std::move(value));
}

void AtomLocalPotentialEditor::SetPeelingSamplingEntries(LocalPotentialSampleList value)
{
    m_entry->SetPeelingSamplingEntries(std::move(value));
}

void AtomLocalPotentialEditor::SetNeighborCountForPeeling(int value)
{
    m_entry->SetNeighborCountForPeeling(value);
}

void AtomLocalPotentialEditor::SetGaussianResult(
    FittingStage stage,
    LocalGaussianResult value)
{
    m_entry->SetGaussianResult(stage, std::move(value));
}

void AtomLocalPotentialEditor::SetAlphaR(FittingStage stage, double value)
{
    m_entry->SetAlphaR(stage, value);
}

} // namespace rhbm_gem
