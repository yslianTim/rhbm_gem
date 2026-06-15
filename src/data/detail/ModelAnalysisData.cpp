#include "data/detail/ModelAnalysisData.hpp"

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

ModelAnalysisData::ModelAnalysisData() = default;
ModelAnalysisData::~ModelAnalysisData() = default;

ModelAnalysisData & ModelAnalysisData::Of(ModelObject & model_object)
{
    return *model_object.m_analysis_data;
}

const ModelAnalysisData & ModelAnalysisData::Of(const ModelObject & model_object)
{
    return *model_object.m_analysis_data;
}

void ModelAnalysisData::Clear()
{
    m_atom_group_entry = AtomGroupPotentialEntry{};
    m_atom_local_entry_map.clear();
}

AtomGroupPotentialEntry & ModelAnalysisData::AtomGroupEntry()
{
    return m_atom_group_entry;
}

const AtomGroupPotentialEntry & ModelAnalysisData::AtomGroupEntry() const
{
    return m_atom_group_entry;
}

LocalPotentialEntry & ModelAnalysisData::EnsureAtomLocalEntry(const AtomObject & atom_object)
{
    auto & entry{ m_atom_local_entry_map[BuildAtomFitStateKey(atom_object)] };
    if (entry == nullptr)
    {
        entry = std::make_unique<LocalPotentialEntry>();
    }
    return *entry;
}

void ModelAnalysisData::SetAtomLocalEntry(
    const AtomObject & atom_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    m_atom_local_entry_map[BuildAtomFitStateKey(atom_object)] = std::move(entry);
}

LocalPotentialEntry * ModelAnalysisData::FindAtomLocalEntry(const AtomObject & atom_object)
{
    const auto iter{ m_atom_local_entry_map.find(BuildAtomFitStateKey(atom_object)) };
    return iter == m_atom_local_entry_map.end() || iter->second == nullptr ? nullptr : iter->second.get();
}

const LocalPotentialEntry * ModelAnalysisData::FindAtomLocalEntry(const AtomObject & atom_object) const
{
    const auto iter{ m_atom_local_entry_map.find(BuildAtomFitStateKey(atom_object)) };
    return iter == m_atom_local_entry_map.end() || iter->second == nullptr ? nullptr : iter->second.get();
}

ModelAnalysisData::AtomLocalEntryMap & ModelAnalysisData::AtomLocalEntries()
{
    return m_atom_local_entry_map;
}

const ModelAnalysisData::AtomLocalEntryMap & ModelAnalysisData::AtomLocalEntries() const
{
    return m_atom_local_entry_map;
}

int ModelAnalysisData::BuildAtomFitStateKey(const AtomObject & atom_object)
{
    return atom_object.GetSerialID();
}

} // namespace rhbm_gem
