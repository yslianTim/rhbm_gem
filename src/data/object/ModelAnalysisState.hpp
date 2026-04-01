#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class GroupPotentialEntry;
class LocalPotentialFitState;

class ModelAnalysisState
{
    std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> m_atom_group_potential_entry_map;
    std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> m_bond_group_potential_entry_map;
    std::unordered_map<int, std::unique_ptr<LocalPotentialFitState>> m_atom_fit_state_map;
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialFitState>> m_bond_fit_state_map;

public:
    ModelAnalysisState();
    ~ModelAnalysisState();

    void Clear();
    void ClearGroupEntries();
    void ClearFitStates();

    void SetAtomGroupPotentialEntry(
        const std::string & class_key,
        std::unique_ptr<GroupPotentialEntry> entry);
    void SetBondGroupPotentialEntry(
        const std::string & class_key,
        std::unique_ptr<GroupPotentialEntry> entry);

    GroupPotentialEntry * GetAtomGroupPotentialEntry(const std::string & class_key) const;
    GroupPotentialEntry * GetBondGroupPotentialEntry(const std::string & class_key) const;

    const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
    GetAtomGroupPotentialEntryMap() const;
    const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
    GetBondGroupPotentialEntryMap() const;

    LocalPotentialFitState & EnsureAtomFitState(const AtomObject & atom_object);
    LocalPotentialFitState & EnsureBondFitState(const BondObject & bond_object);
    LocalPotentialFitState * FindAtomFitState(const AtomObject & atom_object);
    LocalPotentialFitState * FindBondFitState(const BondObject & bond_object);
    const LocalPotentialFitState * FindAtomFitState(const AtomObject & atom_object) const;
    const LocalPotentialFitState * FindBondFitState(const BondObject & bond_object) const;

private:
    static int BuildAtomFitStateKey(const AtomObject & atom_object);
    static std::pair<int, int> BuildBondFitStateKey(const BondObject & bond_object);
};

} // namespace rhbm_gem
