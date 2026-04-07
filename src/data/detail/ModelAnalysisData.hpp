#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class GroupPotentialEntry;
class LocalPotentialEntry;
class LocalPotentialFitState;
class ModelObject;

class ModelAnalysisData
{
public:
    using AtomGroupEntryMap = std::unordered_map<std::string, GroupPotentialEntry>;
    using AtomLocalEntryMap = std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>;
    using AtomFitStateMap = std::unordered_map<int, LocalPotentialFitState>;
    using BondGroupEntryMap = std::unordered_map<std::string, GroupPotentialEntry>;
    using BondFitStateKey = std::pair<int, int>;
    using BondLocalEntryMap = std::map<BondFitStateKey, std::unique_ptr<LocalPotentialEntry>>;
    using BondFitStateMap = std::map<BondFitStateKey, LocalPotentialFitState>;

private:
    AtomGroupEntryMap m_atom_group_entry_map;
    AtomLocalEntryMap m_atom_local_entry_map;
    AtomFitStateMap m_atom_fit_state_map;
    BondGroupEntryMap m_bond_group_entry_map;
    BondLocalEntryMap m_bond_local_entry_map;
    BondFitStateMap m_bond_fit_state_map;

public:
    ModelAnalysisData();
    ~ModelAnalysisData();

    static ModelAnalysisData & Of(ModelObject & model_object);
    static const ModelAnalysisData & Of(const ModelObject & model_object);

    static ModelObject * OwnerOf(const AtomObject & atom_object);
    static ModelObject * OwnerOf(const BondObject & bond_object);

    static const LocalPotentialEntry * FindLocalEntry(const AtomObject & atom_object);
    static const LocalPotentialEntry * FindLocalEntry(const BondObject & bond_object);
    static const LocalPotentialEntry & RequireLocalEntry(
        const LocalPotentialEntry * entry,
        const char * context);
    static const LocalPotentialEntry & RequireLocalEntry(const AtomObject & atom_object);
    static const LocalPotentialEntry & RequireLocalEntry(const BondObject & bond_object);

    void Clear();
    void ClearTransientFitStates();
    void RebuildAtomGroupEntriesFromSelection(const ModelObject & model_object);
    void RebuildBondGroupEntriesFromSelection(const ModelObject & model_object);
    std::vector<GroupKey> CollectAtomGroupKeys(const std::string & class_key) const;
    std::vector<GroupKey> CollectBondGroupKeys(const std::string & class_key) const;

    GroupPotentialEntry & EnsureAtomGroupEntry(const std::string & class_key);
    GroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key);
    const GroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key) const;
    const AtomGroupEntryMap & AtomGroupEntries() const;

    LocalPotentialEntry & EnsureAtomLocalEntry(const AtomObject & atom_object);
    void SetAtomLocalEntry(const AtomObject & atom_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindAtomLocalEntry(const AtomObject & atom_object);
    const LocalPotentialEntry * FindAtomLocalEntry(const AtomObject & atom_object) const;

    LocalPotentialFitState & EnsureAtomFitState(const AtomObject & atom_object);
    LocalPotentialFitState * FindAtomFitState(const AtomObject & atom_object);
    const LocalPotentialFitState * FindAtomFitState(const AtomObject & atom_object) const;

    GroupPotentialEntry & EnsureBondGroupEntry(const std::string & class_key);
    GroupPotentialEntry * FindBondGroupEntry(const std::string & class_key);
    const GroupPotentialEntry * FindBondGroupEntry(const std::string & class_key) const;
    const BondGroupEntryMap & BondGroupEntries() const;

    LocalPotentialEntry & EnsureBondLocalEntry(const BondObject & bond_object);
    void SetBondLocalEntry(const BondObject & bond_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindBondLocalEntry(const BondObject & bond_object);
    const LocalPotentialEntry * FindBondLocalEntry(const BondObject & bond_object) const;

    LocalPotentialFitState & EnsureBondFitState(const BondObject & bond_object);
    LocalPotentialFitState * FindBondFitState(const BondObject & bond_object);
    const LocalPotentialFitState * FindBondFitState(const BondObject & bond_object) const;

private:
    static int BuildAtomFitStateKey(const AtomObject & atom_object);
    static BondFitStateKey BuildBondFitStateKey(const BondObject & bond_object);
};

} // namespace rhbm_gem
