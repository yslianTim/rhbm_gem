#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "data/detail/GroupPotentialEntry.hpp"
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class LocalPotentialEntry;
class ModelObject;

class ModelAnalysisData
{
public:
    using AtomGroupEntryMap = std::unordered_map<std::string, AtomGroupPotentialEntry>;
    using AtomLocalEntryMap = std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>;
    using BondGroupEntryMap = std::unordered_map<std::string, BondGroupPotentialEntry>;
    using BondFitStateKey = std::pair<int, int>;
    using BondLocalEntryMap = std::map<BondFitStateKey, std::unique_ptr<LocalPotentialEntry>>;

private:
    AtomGroupEntryMap m_atom_group_entry_map;
    AtomLocalEntryMap m_atom_local_entry_map;
    BondGroupEntryMap m_bond_group_entry_map;
    BondLocalEntryMap m_bond_local_entry_map;

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

    AtomGroupPotentialEntry & EnsureAtomGroupEntry(const std::string & class_key);
    AtomGroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key);
    const AtomGroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key) const;
    const AtomGroupEntryMap & AtomGroupEntries() const;

    LocalPotentialEntry & EnsureAtomLocalEntry(const AtomObject & atom_object);
    void SetAtomLocalEntry(const AtomObject & atom_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindAtomLocalEntry(const AtomObject & atom_object);
    const LocalPotentialEntry * FindAtomLocalEntry(const AtomObject & atom_object) const;

    BondGroupPotentialEntry & EnsureBondGroupEntry(const std::string & class_key);
    BondGroupPotentialEntry * FindBondGroupEntry(const std::string & class_key);
    const BondGroupPotentialEntry * FindBondGroupEntry(const std::string & class_key) const;
    const BondGroupEntryMap & BondGroupEntries() const;

    LocalPotentialEntry & EnsureBondLocalEntry(const BondObject & bond_object);
    void SetBondLocalEntry(const BondObject & bond_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindBondLocalEntry(const BondObject & bond_object);
    const LocalPotentialEntry * FindBondLocalEntry(const BondObject & bond_object) const;

private:
    static int BuildAtomFitStateKey(const AtomObject & atom_object);
    static BondFitStateKey BuildBondFitStateKey(const BondObject & bond_object);
};

} // namespace rhbm_gem
