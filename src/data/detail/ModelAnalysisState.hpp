#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <utility>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class GroupPotentialEntry;
class LocalPotentialFitState;
class ModelObjectAccess;

class AtomAnalysisStore
{
public:
    using GroupEntryMap = std::unordered_map<std::string, GroupPotentialEntry>;
    using FitStateMap = std::unordered_map<int, LocalPotentialFitState>;

private:
    GroupEntryMap m_group_entry_map;
    FitStateMap m_fit_state_map;

public:
    GroupPotentialEntry & EnsureGroupEntry(const std::string & class_key);
    GroupPotentialEntry * FindGroupEntry(const std::string & class_key);
    const GroupPotentialEntry * FindGroupEntry(const std::string & class_key) const;
    const GroupEntryMap & Entries() const;

    LocalPotentialFitState & EnsureFitState(const AtomObject & atom_object);
    LocalPotentialFitState * FindFitState(const AtomObject & atom_object);
    const LocalPotentialFitState * FindFitState(const AtomObject & atom_object) const;

private:
    friend class ModelAnalysisState;

    void ClearGroupEntries();
    void ClearFitStates();
    void Clear();

    static int BuildFitStateKey(const AtomObject & atom_object);
};

class BondAnalysisStore
{
public:
    using GroupEntryMap = std::unordered_map<std::string, GroupPotentialEntry>;
    using FitStateKey = std::pair<int, int>;
    using FitStateMap = std::map<FitStateKey, LocalPotentialFitState>;

private:
    GroupEntryMap m_group_entry_map;
    FitStateMap m_fit_state_map;

public:
    GroupPotentialEntry & EnsureGroupEntry(const std::string & class_key);
    GroupPotentialEntry * FindGroupEntry(const std::string & class_key);
    const GroupPotentialEntry * FindGroupEntry(const std::string & class_key) const;
    const GroupEntryMap & Entries() const;

    LocalPotentialFitState & EnsureFitState(const BondObject & bond_object);
    LocalPotentialFitState * FindFitState(const BondObject & bond_object);
    const LocalPotentialFitState * FindFitState(const BondObject & bond_object) const;

private:
    friend class ModelAnalysisState;

    void ClearGroupEntries();
    void ClearFitStates();
    void Clear();

    static FitStateKey BuildFitStateKey(const BondObject & bond_object);
};

class ModelAnalysisState
{
    AtomAnalysisStore m_atoms;
    BondAnalysisStore m_bonds;

public:
    ModelAnalysisState();
    ~ModelAnalysisState();

    void Clear();

    AtomAnalysisStore & Atoms() { return m_atoms; }
    const AtomAnalysisStore & Atoms() const { return m_atoms; }
    BondAnalysisStore & Bonds() { return m_bonds; }
    const BondAnalysisStore & Bonds() const { return m_bonds; }

private:
    friend class ModelObjectAccess;

    void ClearFitStates();
};

} // namespace rhbm_gem
