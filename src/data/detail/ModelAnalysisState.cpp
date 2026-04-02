#include "data/detail/ModelAnalysisState.hpp"

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialFitState.hpp"
#include "data/detail/ModelObjectAccess.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

namespace rhbm_gem {

void AtomAnalysisStore::ClearGroupEntries()
{
    m_group_entry_map.clear();
}

void AtomAnalysisStore::ClearFitStates()
{
    m_fit_state_map.clear();
}

void AtomAnalysisStore::Clear()
{
    ClearGroupEntries();
    ClearFitStates();
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

void BondAnalysisStore::ClearGroupEntries()
{
    m_group_entry_map.clear();
}

void BondAnalysisStore::ClearFitStates()
{
    m_fit_state_map.clear();
}

void BondAnalysisStore::Clear()
{
    ClearGroupEntries();
    ClearFitStates();
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

ModelAnalysisState::ModelAnalysisState() = default;

ModelAnalysisState::~ModelAnalysisState() = default;

void ModelAnalysisState::Clear()
{
    m_atoms.Clear();
    m_bonds.Clear();
}

void ModelAnalysisState::ClearFitStates()
{
    m_atoms.ClearFitStates();
    m_bonds.ClearFitStates();
}

void ModelObjectAccess::ClearAnalysisFitStates(ModelObject & model_object)
{
    model_object.m_analysis_state->ClearFitStates();
}

} // namespace rhbm_gem
