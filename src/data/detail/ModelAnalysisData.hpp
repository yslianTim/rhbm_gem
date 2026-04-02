#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class GroupPotentialEntry;
class LocalPotentialEntry;
class LocalPotentialFitState;
class ModelAnalysisAccess;
class ModelObject;

class AtomAnalysisStore
{
public:
    using GroupEntryMap = std::unordered_map<std::string, GroupPotentialEntry>;
    using FitStateMap = std::unordered_map<int, LocalPotentialFitState>;
    using LocalEntryMap = std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>;

private:
    GroupEntryMap m_group_entry_map;
    FitStateMap m_fit_state_map;
    LocalEntryMap m_local_entry_map;

public:
    GroupPotentialEntry & EnsureGroupEntry(const std::string & class_key);
    GroupPotentialEntry * FindGroupEntry(const std::string & class_key);
    const GroupPotentialEntry * FindGroupEntry(const std::string & class_key) const;
    const GroupEntryMap & Entries() const;

    LocalPotentialEntry & EnsureLocalEntry(const AtomObject & atom_object);
    void SetLocalEntry(const AtomObject & atom_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindLocalEntry(const AtomObject & atom_object);
    const LocalPotentialEntry * FindLocalEntry(const AtomObject & atom_object) const;
    const LocalEntryMap & LocalEntries() const;

    LocalPotentialFitState & EnsureFitState(const AtomObject & atom_object);
    LocalPotentialFitState * FindFitState(const AtomObject & atom_object);
    const LocalPotentialFitState * FindFitState(const AtomObject & atom_object) const;

private:
    friend class ModelAnalysisData;

    void ClearGroupEntries();
    void ClearLocalEntries();
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
    using LocalEntryMap = std::map<FitStateKey, std::unique_ptr<LocalPotentialEntry>>;

private:
    GroupEntryMap m_group_entry_map;
    FitStateMap m_fit_state_map;
    LocalEntryMap m_local_entry_map;

public:
    GroupPotentialEntry & EnsureGroupEntry(const std::string & class_key);
    GroupPotentialEntry * FindGroupEntry(const std::string & class_key);
    const GroupPotentialEntry * FindGroupEntry(const std::string & class_key) const;
    const GroupEntryMap & Entries() const;

    LocalPotentialEntry & EnsureLocalEntry(const BondObject & bond_object);
    void SetLocalEntry(const BondObject & bond_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindLocalEntry(const BondObject & bond_object);
    const LocalPotentialEntry * FindLocalEntry(const BondObject & bond_object) const;
    const LocalEntryMap & LocalEntries() const;

    LocalPotentialFitState & EnsureFitState(const BondObject & bond_object);
    LocalPotentialFitState * FindFitState(const BondObject & bond_object);
    const LocalPotentialFitState * FindFitState(const BondObject & bond_object) const;

private:
    friend class ModelAnalysisData;

    void ClearGroupEntries();
    void ClearLocalEntries();
    void ClearFitStates();
    void Clear();

    static FitStateKey BuildFitStateKey(const BondObject & bond_object);
};

class ModelAnalysisData
{
    AtomAnalysisStore m_atoms;
    BondAnalysisStore m_bonds;

public:
    ModelAnalysisData();
    ~ModelAnalysisData();

    void Clear();

    AtomAnalysisStore & Atoms() { return m_atoms; }
    const AtomAnalysisStore & Atoms() const { return m_atoms; }
    BondAnalysisStore & Bonds() { return m_bonds; }
    const BondAnalysisStore & Bonds() const { return m_bonds; }

private:
    friend class ModelAnalysisAccess;

    void ClearFitStates();
};

} // namespace rhbm_gem
