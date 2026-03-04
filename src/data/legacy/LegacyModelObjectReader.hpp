#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;
class AtomObject;
class BondObject;
class LocalPotentialEntry;

class LegacyModelObjectReader
{
    SQLiteWrapper * m_database;
    mutable std::unordered_set<std::string> m_table_cache;

public:
    explicit LegacyModelObjectReader(SQLiteWrapper * database);

    std::unique_ptr<ModelObject> Load(const std::string & key_tag);

private:
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & key_tag);
    std::vector<std::unique_ptr<BondObject>> LoadBondObjectList(
        const std::string & key_tag, const ModelObject * model_object);

    void LoadChemicalComponentEntryList(ModelObject * model_obj, const std::string & table_name);
    void LoadComponentAtomEntryList(ModelObject * model_obj, const std::string & table_name);
    void LoadComponentBondEntryList(ModelObject * model_obj, const std::string & table_name);

    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>
    LoadAtomLocalPotentialEntryMap(const std::string & table_name);
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>>
    LoadBondLocalPotentialEntryMap(const std::string & table_name);
    void LoadAtomLocalPotentialEntrySubList(
        const std::string & table_name,
        const std::string & class_key,
        std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map);
    void LoadBondLocalPotentialEntrySubList(
        const std::string & table_name,
        const std::string & class_key,
        std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map);
    void LoadAtomGroupPotentialEntryList(
        ModelObject * model_obj, const std::string & class_key, const std::string & table_name);
    void LoadBondGroupPotentialEntryList(
        ModelObject * model_obj, const std::string & class_key, const std::string & table_name);

    bool TableExists(const std::string & table_name) const;
    std::string SanitizeTableName(const std::string & key_tag) const;
};

} // namespace rhbm_gem
