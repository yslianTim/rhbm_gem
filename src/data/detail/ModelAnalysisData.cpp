#include "data/detail/ModelAnalysisData.hpp"

#include "data/detail/AtomClassifier.hpp"
#include "data/detail/BondClassifier.hpp"
#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem {

namespace {

template <typename EntryT>
std::vector<GroupKey> CollectGroupKeys(const EntryT * entry)
{
    return entry == nullptr ? std::vector<GroupKey>{} : entry->CollectGroupKeys();
}

} // namespace

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

ModelObject * ModelAnalysisData::OwnerOf(const AtomObject & atom_object)
{
    return atom_object.m_owner_model;
}

ModelObject * ModelAnalysisData::OwnerOf(const BondObject & bond_object)
{
    return bond_object.m_owner_model;
}

const LocalPotentialEntry * ModelAnalysisData::FindLocalEntry(const AtomObject & atom_object)
{
    auto * owner{ OwnerOf(atom_object) };
    return owner == nullptr ? nullptr : Of(*owner).FindAtomLocalEntry(atom_object);
}

const LocalPotentialEntry * ModelAnalysisData::FindLocalEntry(const BondObject & bond_object)
{
    auto * owner{ OwnerOf(bond_object) };
    return owner == nullptr ? nullptr : Of(*owner).FindBondLocalEntry(bond_object);
}

const LocalPotentialEntry & ModelAnalysisData::RequireLocalEntry(
    const LocalPotentialEntry * entry,
    const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
}

const LocalPotentialEntry & ModelAnalysisData::RequireLocalEntry(const AtomObject & atom_object)
{
    return RequireLocalEntry(FindLocalEntry(atom_object), "Atom local entry");
}

const LocalPotentialEntry & ModelAnalysisData::RequireLocalEntry(const BondObject & bond_object)
{
    return RequireLocalEntry(FindLocalEntry(bond_object), "Bond local entry");
}

void ModelAnalysisData::Clear()
{
    ClearTransientFitStates();
    m_atom_group_entry_map.clear();
    m_atom_local_entry_map.clear();
    m_bond_group_entry_map.clear();
    m_bond_local_entry_map.clear();
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
    for (auto & [key, entry] : m_bond_local_entry_map)
    {
        (void)key;
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

void ModelAnalysisData::RebuildBondGroupEntriesFromSelection(const ModelObject & model_object)
{
    m_bond_group_entry_map.clear();
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto & group_entry{ EnsureBondGroupEntry(class_key) };
        for (auto * bond : model_object.GetSelectedBonds())
        {
            const auto group_key{ BondClassifier::GetGroupKeyInClass(bond, class_key) };
            group_entry.AddMember(group_key, *bond);
        }
    }
}

std::vector<GroupKey> ModelAnalysisData::CollectAtomGroupKeys(const std::string & class_key) const
{
    return CollectGroupKeys(FindAtomGroupEntry(class_key));
}

std::vector<GroupKey> ModelAnalysisData::CollectBondGroupKeys(const std::string & class_key) const
{
    return CollectGroupKeys(FindBondGroupEntry(class_key));
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

BondGroupPotentialEntry & ModelAnalysisData::EnsureBondGroupEntry(const std::string & class_key)
{
    auto [iter, inserted]{ m_bond_group_entry_map.try_emplace(class_key) };
    (void)inserted;
    return iter->second;
}

BondGroupPotentialEntry * ModelAnalysisData::FindBondGroupEntry(const std::string & class_key)
{
    const auto iter{ m_bond_group_entry_map.find(class_key) };
    return iter == m_bond_group_entry_map.end() ? nullptr : &iter->second;
}

const BondGroupPotentialEntry * ModelAnalysisData::FindBondGroupEntry(
    const std::string & class_key) const
{
    const auto iter{ m_bond_group_entry_map.find(class_key) };
    return iter == m_bond_group_entry_map.end() ? nullptr : &iter->second;
}

const ModelAnalysisData::BondGroupEntryMap & ModelAnalysisData::BondGroupEntries() const
{
    return m_bond_group_entry_map;
}

LocalPotentialEntry & ModelAnalysisData::EnsureBondLocalEntry(const BondObject & bond_object)
{
    auto & entry{ m_bond_local_entry_map[BuildBondFitStateKey(bond_object)] };
    if (entry == nullptr)
    {
        entry = std::make_unique<LocalPotentialEntry>();
    }
    return *entry;
}

void ModelAnalysisData::SetBondLocalEntry(
    const BondObject & bond_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    m_bond_local_entry_map[BuildBondFitStateKey(bond_object)] = std::move(entry);
}

LocalPotentialEntry * ModelAnalysisData::FindBondLocalEntry(const BondObject & bond_object)
{
    const auto iter{ m_bond_local_entry_map.find(BuildBondFitStateKey(bond_object)) };
    return iter == m_bond_local_entry_map.end() || iter->second == nullptr ? nullptr : iter->second.get();
}

const LocalPotentialEntry * ModelAnalysisData::FindBondLocalEntry(const BondObject & bond_object) const
{
    const auto iter{ m_bond_local_entry_map.find(BuildBondFitStateKey(bond_object)) };
    return iter == m_bond_local_entry_map.end() || iter->second == nullptr ? nullptr : iter->second.get();
}

int ModelAnalysisData::BuildAtomFitStateKey(const AtomObject & atom_object)
{
    return atom_object.GetSerialID();
}

ModelAnalysisData::BondFitStateKey ModelAnalysisData::BuildBondFitStateKey(
    const BondObject & bond_object)
{
    return { bond_object.GetAtomSerialID1(), bond_object.GetAtomSerialID2() };
}

} // namespace rhbm_gem
