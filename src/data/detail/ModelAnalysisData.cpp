#include "data/detail/ModelAnalysisData.hpp"

#include "data/detail/AtomClassifier.hpp"
#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#include <string>

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
    ClearTransientFitStates();
    m_atom_group_entry_map.clear();
    m_atom_local_entry_map.clear();
}

void ModelAnalysisData::ClearTransientFitStates()
{
    for (auto & [serial_id, entry] : m_atom_local_entry_map)
    {
        (void)serial_id;
        if (entry != nullptr)
        {
            entry->ClearTransientFitState();
        }
    }
}

void ModelAnalysisData::RebuildAtomGroupEntriesFromSelection(const ModelObject & model_object)
{
    m_atom_group_entry_map.clear();
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto & group_entry{ EnsureAtomGroupEntry(class_key) };
        for (auto * atom : model_object.GetSelectedAtoms())
        {
            const auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_entry.AddMember(group_key, *atom);
        }
    }
}

AtomGroupPotentialEntry & ModelAnalysisData::EnsureAtomGroupEntry(const std::string & class_key)
{
    auto [iter, inserted]{ m_atom_group_entry_map.try_emplace(class_key) };
    (void)inserted;
    return iter->second;
}

AtomGroupPotentialEntry * ModelAnalysisData::FindAtomGroupEntry(const std::string & class_key)
{
    const auto iter{ m_atom_group_entry_map.find(class_key) };
    return iter == m_atom_group_entry_map.end() ? nullptr : &iter->second;
}

const AtomGroupPotentialEntry * ModelAnalysisData::FindAtomGroupEntry(
    const std::string & class_key) const
{
    const auto iter{ m_atom_group_entry_map.find(class_key) };
    return iter == m_atom_group_entry_map.end() ? nullptr : &iter->second;
}

const ModelAnalysisData::AtomGroupEntryMap & ModelAnalysisData::AtomGroupEntries() const
{
    return m_atom_group_entry_map;
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

int ModelAnalysisData::BuildAtomFitStateKey(const AtomObject & atom_object)
{
    return atom_object.GetSerialID();
}

} // namespace rhbm_gem
