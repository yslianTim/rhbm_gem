#include "data/object/ModelAnalysisState.hpp"

#include "data/object/LocalPotentialFitState.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>

#include <stdexcept>

namespace rhbm_gem {

ModelAnalysisState::ModelAnalysisState() = default;

ModelAnalysisState::~ModelAnalysisState() = default;

void ModelAnalysisState::Clear()
{
    ClearGroupEntries();
    ClearFitStates();
}

void ModelAnalysisState::ClearGroupEntries()
{
    m_atom_group_potential_entry_map.clear();
    m_bond_group_potential_entry_map.clear();
}

void ModelAnalysisState::ClearFitStates()
{
    m_atom_fit_state_map.clear();
    m_bond_fit_state_map.clear();
}

void ModelAnalysisState::SetAtomGroupPotentialEntry(
    const std::string & class_key,
    std::unique_ptr<GroupPotentialEntry> entry)
{
    m_atom_group_potential_entry_map[class_key] = std::move(entry);
}

void ModelAnalysisState::SetBondGroupPotentialEntry(
    const std::string & class_key,
    std::unique_ptr<GroupPotentialEntry> entry)
{
    m_bond_group_potential_entry_map[class_key] = std::move(entry);
}

GroupPotentialEntry * ModelAnalysisState::GetAtomGroupPotentialEntry(const std::string & class_key) const
{
    return m_atom_group_potential_entry_map.at(class_key).get();
}

GroupPotentialEntry * ModelAnalysisState::GetBondGroupPotentialEntry(const std::string & class_key) const
{
    return m_bond_group_potential_entry_map.at(class_key).get();
}

const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
ModelAnalysisState::GetAtomGroupPotentialEntryMap() const
{
    return m_atom_group_potential_entry_map;
}

const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
ModelAnalysisState::GetBondGroupPotentialEntryMap() const
{
    return m_bond_group_potential_entry_map;
}

LocalPotentialFitState & ModelAnalysisState::EnsureAtomFitState(const AtomObject & atom_object)
{
    const auto key{ BuildAtomFitStateKey(atom_object) };
    auto & slot{ m_atom_fit_state_map[key] };
    if (slot == nullptr)
    {
        slot = std::make_unique<LocalPotentialFitState>();
    }
    return *slot;
}

LocalPotentialFitState & ModelAnalysisState::EnsureBondFitState(const BondObject & bond_object)
{
    const auto key{ BuildBondFitStateKey(bond_object) };
    auto & slot{ m_bond_fit_state_map[key] };
    if (slot == nullptr)
    {
        slot = std::make_unique<LocalPotentialFitState>();
    }
    return *slot;
}

LocalPotentialFitState * ModelAnalysisState::FindAtomFitState(const AtomObject & atom_object)
{
    const auto key{ BuildAtomFitStateKey(atom_object) };
    const auto iter{ m_atom_fit_state_map.find(key) };
    if (iter == m_atom_fit_state_map.end())
    {
        return nullptr;
    }
    return iter->second.get();
}

LocalPotentialFitState * ModelAnalysisState::FindBondFitState(const BondObject & bond_object)
{
    const auto key{ BuildBondFitStateKey(bond_object) };
    const auto iter{ m_bond_fit_state_map.find(key) };
    if (iter == m_bond_fit_state_map.end())
    {
        return nullptr;
    }
    return iter->second.get();
}

const LocalPotentialFitState * ModelAnalysisState::FindAtomFitState(const AtomObject & atom_object) const
{
    const auto key{ BuildAtomFitStateKey(atom_object) };
    const auto iter{ m_atom_fit_state_map.find(key) };
    if (iter == m_atom_fit_state_map.end())
    {
        return nullptr;
    }
    return iter->second.get();
}

const LocalPotentialFitState * ModelAnalysisState::FindBondFitState(const BondObject & bond_object) const
{
    const auto key{ BuildBondFitStateKey(bond_object) };
    const auto iter{ m_bond_fit_state_map.find(key) };
    if (iter == m_bond_fit_state_map.end())
    {
        return nullptr;
    }
    return iter->second.get();
}

int ModelAnalysisState::BuildAtomFitStateKey(const AtomObject & atom_object)
{
    return atom_object.GetSerialID();
}

std::pair<int, int> ModelAnalysisState::BuildBondFitStateKey(const BondObject & bond_object)
{
    return { bond_object.GetAtomSerialID1(), bond_object.GetAtomSerialID2() };
}

} // namespace rhbm_gem
