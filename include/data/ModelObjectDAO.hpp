#pragma once

#include <memory>
#include <vector>
#include <map>
#include <utility>
#include <unordered_map>
#include <unordered_set>

#include "DataObjectDAOBase.hpp"

class SQLiteWrapper;
class ModelObject;
class AtomObject;
class BondObject;
class LocalPotentialEntry;
class GroupPotentialEntry;

class ModelObjectDAO : public DataObjectDAOBase
{
    SQLiteWrapper * m_database;
    mutable std::unordered_set<std::string> m_table_cache;

public:
    ModelObjectDAO(SQLiteWrapper * db_manager);
    ~ModelObjectDAO();

    void Save(const DataObjectBase * obj, const std::string & key_tag) override;
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) override;
    
private:
    void SaveAtomObjectList(const ModelObject * model_obj, const std::string & table_name);
    void SaveBondObjectList(const ModelObject * model_obj, const std::string & table_name);
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & key_tag);
    std::vector<std::unique_ptr<BondObject>> LoadBondObjectList(const std::string & key_tag, const ModelObject * model_obj);

    void SaveChemicalComponentEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveComponentAtomEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveComponentBondEntryList(const ModelObject * model_obj, const std::string & table_name);
    void LoadChemicalComponentEntryList(ModelObject * model_obj, const std::string & table_name);
    void LoadComponentAtomEntryList(ModelObject * model_obj, const std::string & table_name);
    void LoadComponentBondEntryList(ModelObject * model_obj, const std::string & table_name);

    void SaveAtomLocalPotentialEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveBondLocalPotentialEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveAtomLocalPotentialEntrySubList(const ModelObject * model_obj, const std::string & table_name, const std::string & class_key);
    void SaveBondLocalPotentialEntrySubList(const ModelObject * model_obj, const std::string & table_name, const std::string & class_key);
    void SaveAtomGroupPotentialEntryList(const GroupPotentialEntry * group_entry, const std::string & table_name);
    void SaveBondGroupPotentialEntryList(const GroupPotentialEntry * group_entry, const std::string & table_name);

    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> LoadAtomLocalPotentialEntryMap(const std::string & table_name);
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> LoadBondLocalPotentialEntryMap(const std::string & table_name);
    void LoadAtomLocalPotentialEntrySubList(const std::string & table_name, const std::string & class_key, std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map);
    void LoadBondLocalPotentialEntrySubList(const std::string & table_name, const std::string & class_key, std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map);
    void LoadAtomGroupPotentialEntryList(ModelObject * model_obj, const std::string & class_key, const std::string & table_name);
    void LoadBondGroupPotentialEntryList(ModelObject * model_obj, const std::string & class_key, const std::string & table_name);


    bool TableExists(const std::string & table_name) const;
    std::string SanitizeTableName(const std::string & key_tag) const;
    
};
