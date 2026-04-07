#include "data/detail/ModelAnalysisData.hpp"

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/LocalPotentialFitState.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <stdexcept>
#include <string>

namespace rhbm_gem {

void AtomAnalysisStore::ClearTransientFitStates()
{
    m_fit_state_map.clear();
}

void AtomAnalysisStore::Clear()
{
    m_group_entry_map.clear();
    m_local_entry_map.clear();
    ClearTransientFitStates();
}

GroupPotentialEntry & AtomAnalysisStore::EnsureGroupEntry(const std::string & class_key)
{
    return m_group_entry_map[class_key];
}

GroupPotentialEntry * AtomAnalysisStore::FindGroupEntry(const std::string & class_key)
{
    const auto iter{ m_group_entry_map.find(class_key) };
    return (iter == m_group_entry_map.end()) ? nullptr : &iter->second;
}

const GroupPotentialEntry * AtomAnalysisStore::FindGroupEntry(const std::string & class_key) const
{
    const auto iter{ m_group_entry_map.find(class_key) };
    return (iter == m_group_entry_map.end()) ? nullptr : &iter->second;
}

const AtomAnalysisStore::GroupEntryMap & AtomAnalysisStore::Entries() const
{
    return m_group_entry_map;
}

LocalPotentialEntry & AtomAnalysisStore::EnsureLocalEntry(const AtomObject & atom_object)
{
    auto & entry{ m_local_entry_map[BuildFitStateKey(atom_object)] };
    if (entry == nullptr)
    {
        entry = std::make_unique<LocalPotentialEntry>();
    }
    return *entry;
}

void AtomAnalysisStore::SetLocalEntry(
    const AtomObject & atom_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    m_local_entry_map[BuildFitStateKey(atom_object)] = std::move(entry);
}

LocalPotentialEntry * AtomAnalysisStore::FindLocalEntry(const AtomObject & atom_object)
{
    const auto iter{ m_local_entry_map.find(BuildFitStateKey(atom_object)) };
    return (iter == m_local_entry_map.end() || iter->second == nullptr) ? nullptr : iter->second.get();
}

const LocalPotentialEntry * AtomAnalysisStore::FindLocalEntry(const AtomObject & atom_object) const
{
    const auto iter{ m_local_entry_map.find(BuildFitStateKey(atom_object)) };
    return (iter == m_local_entry_map.end() || iter->second == nullptr) ? nullptr : iter->second.get();
}

const AtomAnalysisStore::LocalEntryMap & AtomAnalysisStore::LocalEntries() const
{
    return m_local_entry_map;
}

LocalPotentialFitState & AtomAnalysisStore::EnsureFitState(const AtomObject & atom_object)
{
    return m_fit_state_map[BuildFitStateKey(atom_object)];
}

LocalPotentialFitState * AtomAnalysisStore::FindFitState(const AtomObject & atom_object)
{
    const auto iter{ m_fit_state_map.find(BuildFitStateKey(atom_object)) };
    return (iter == m_fit_state_map.end()) ? nullptr : &iter->second;
}

const LocalPotentialFitState * AtomAnalysisStore::FindFitState(const AtomObject & atom_object) const
{
    const auto iter{ m_fit_state_map.find(BuildFitStateKey(atom_object)) };
    return (iter == m_fit_state_map.end()) ? nullptr : &iter->second;
}

int AtomAnalysisStore::BuildFitStateKey(const AtomObject & atom_object)
{
    return atom_object.GetSerialID();
}

void BondAnalysisStore::ClearTransientFitStates()
{
    m_fit_state_map.clear();
}

void BondAnalysisStore::Clear()
{
    m_group_entry_map.clear();
    m_local_entry_map.clear();
    ClearTransientFitStates();
}

GroupPotentialEntry & BondAnalysisStore::EnsureGroupEntry(const std::string & class_key)
{
    return m_group_entry_map[class_key];
}

GroupPotentialEntry * BondAnalysisStore::FindGroupEntry(const std::string & class_key)
{
    const auto iter{ m_group_entry_map.find(class_key) };
    return (iter == m_group_entry_map.end()) ? nullptr : &iter->second;
}

const GroupPotentialEntry * BondAnalysisStore::FindGroupEntry(const std::string & class_key) const
{
    const auto iter{ m_group_entry_map.find(class_key) };
    return (iter == m_group_entry_map.end()) ? nullptr : &iter->second;
}

const BondAnalysisStore::GroupEntryMap & BondAnalysisStore::Entries() const
{
    return m_group_entry_map;
}

LocalPotentialEntry & BondAnalysisStore::EnsureLocalEntry(const BondObject & bond_object)
{
    auto & entry{ m_local_entry_map[BuildFitStateKey(bond_object)] };
    if (entry == nullptr)
    {
        entry = std::make_unique<LocalPotentialEntry>();
    }
    return *entry;
}

void BondAnalysisStore::SetLocalEntry(
    const BondObject & bond_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    m_local_entry_map[BuildFitStateKey(bond_object)] = std::move(entry);
}

LocalPotentialEntry * BondAnalysisStore::FindLocalEntry(const BondObject & bond_object)
{
    const auto iter{ m_local_entry_map.find(BuildFitStateKey(bond_object)) };
    return (iter == m_local_entry_map.end() || iter->second == nullptr) ? nullptr : iter->second.get();
}

const LocalPotentialEntry * BondAnalysisStore::FindLocalEntry(const BondObject & bond_object) const
{
    const auto iter{ m_local_entry_map.find(BuildFitStateKey(bond_object)) };
    return (iter == m_local_entry_map.end() || iter->second == nullptr) ? nullptr : iter->second.get();
}

const BondAnalysisStore::LocalEntryMap & BondAnalysisStore::LocalEntries() const
{
    return m_local_entry_map;
}

LocalPotentialFitState & BondAnalysisStore::EnsureFitState(const BondObject & bond_object)
{
    return m_fit_state_map[BuildFitStateKey(bond_object)];
}

LocalPotentialFitState * BondAnalysisStore::FindFitState(const BondObject & bond_object)
{
    const auto iter{ m_fit_state_map.find(BuildFitStateKey(bond_object)) };
    return (iter == m_fit_state_map.end()) ? nullptr : &iter->second;
}

const LocalPotentialFitState * BondAnalysisStore::FindFitState(const BondObject & bond_object) const
{
    const auto iter{ m_fit_state_map.find(BuildFitStateKey(bond_object)) };
    return (iter == m_fit_state_map.end()) ? nullptr : &iter->second;
}

BondAnalysisStore::FitStateKey BondAnalysisStore::BuildFitStateKey(const BondObject & bond_object)
{
    return { bond_object.GetAtomSerialID1(), bond_object.GetAtomSerialID2() };
}

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
    return owner == nullptr ? nullptr : Of(*owner).Atoms().FindLocalEntry(atom_object);
}

const LocalPotentialEntry * ModelAnalysisData::FindLocalEntry(const BondObject & bond_object)
{
    auto * owner{ OwnerOf(bond_object) };
    return owner == nullptr ? nullptr : Of(*owner).Bonds().FindLocalEntry(bond_object);
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
    m_atoms.Clear();
    m_bonds.Clear();
}

void ModelAnalysisData::ClearTransientFitStates()
{
    m_atoms.ClearTransientFitStates();
    m_bonds.ClearTransientFitStates();
}

} // namespace rhbm_gem
