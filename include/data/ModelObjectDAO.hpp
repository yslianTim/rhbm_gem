#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "DataObjectDAOBase.hpp"

class SQLiteWrapper;
class ModelObject;
class AtomObject;
class AtomicPotentialEntry;
class GroupPotentialEntry;
class ModelObjectDAO : public DataObjectDAOBase
{
    SQLiteWrapper * m_database;
    std::unordered_set<std::string> m_table_cache;

public:
    ModelObjectDAO(SQLiteWrapper * db_manager);
    ~ModelObjectDAO();

    void Save(const DataObjectBase * obj) override;
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) override;
    
private:
    void CreateModelObjectListTable(const std::string & table_name);
    void CreateAtomObjectListTable(const std::string & table_name);
    void CreateAtomicPotentialEntryListTable(const std::string & table_name);
    void CreateAtomicPotentialEntrySubListTable(const std::string & table_name);
    void CreateGroupPotentialEntryListTable(const std::string & table_name);

    void SaveAtomObjectList(const ModelObject * model_obj, const std::string & table_name);
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & key_tag);

    void SaveAtomicPotentialEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveAtomicPotentialEntrySubList(const ModelObject * model_obj, const std::string & table_name, const std::string & class_key);
    void SaveGroupPotentialEntryList(const GroupPotentialEntry * group_entry, const std::string & table_name);

    std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>> LoadAtomicPotentialEntryMap(const std::string & table_name);
    void LoadAtomicPotentialEntrySubList(const std::string & table_name, const std::string & class_key, std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>> & entry_map);
    void LoadGroupPotentialEntryList(ModelObject * model_obj, const std::string & class_key, const std::string & table_name);
    void LoadGroupPotentialEntrySubList(ModelObject * model_obj, const std::string & class_key);

    bool TableExists(const std::string & table_name) const;
    std::string SanitizeTableName(const std::string & key_tag) const;
    
};
