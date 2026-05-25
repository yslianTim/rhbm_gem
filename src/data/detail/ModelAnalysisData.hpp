#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "data/detail/GroupPotentialEntry.hpp"
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class LocalPotentialEntry;
class ModelObject;

class ModelAnalysisData
{
public:
    using AtomGroupEntryMap = std::unordered_map<std::string, AtomGroupPotentialEntry>;
    using AtomLocalEntryMap = std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>;

private:
    AtomGroupEntryMap m_atom_group_entry_map;
    AtomLocalEntryMap m_atom_local_entry_map;

public:
    ModelAnalysisData();
    ~ModelAnalysisData();

    static ModelAnalysisData & Of(ModelObject & model_object);
    static const ModelAnalysisData & Of(const ModelObject & model_object);

    void Clear();
    AtomGroupPotentialEntry & EnsureAtomGroupEntry(const std::string & class_key);
    AtomGroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key);
    const AtomGroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key) const;
    AtomGroupEntryMap & AtomGroupEntries();
    const AtomGroupEntryMap & AtomGroupEntries() const;

    LocalPotentialEntry & EnsureAtomLocalEntry(const AtomObject & atom_object);
    void SetAtomLocalEntry(const AtomObject & atom_object, std::unique_ptr<LocalPotentialEntry> entry);
    LocalPotentialEntry * FindAtomLocalEntry(const AtomObject & atom_object);
    const LocalPotentialEntry * FindAtomLocalEntry(const AtomObject & atom_object) const;
    AtomLocalEntryMap & AtomLocalEntries();
    const AtomLocalEntryMap & AtomLocalEntries() const;

private:
    static int BuildAtomFitStateKey(const AtomObject & atom_object);
};

} // namespace rhbm_gem
