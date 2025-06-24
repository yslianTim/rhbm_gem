#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include "DataObjectDAOBase.hpp"

class SQLiteWrapper;
class ModelObject;
class AtomObject;
class AtomicPotentialEntry;
class GroupPotentialEntry;
class AtomicChargeEntry;
class GroupChargeEntry;
class ModelObjectDAO : public DataObjectDAOBase
{
    SQLiteWrapper * m_database;

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
    void CreateAtomicChargeEntryListTable(const std::string & table_name);
    void CreateAtomicChargeEntrySubListTable(const std::string & table_name);
    void CreateGroupChargeEntryListTable(const std::string & table_name);

    void SaveAtomObjectList(const ModelObject * model_obj, const std::string & table_name);
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & key_tag);

    void SaveAtomicPotentialEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveAtomicPotentialEntrySubList(const ModelObject * model_obj, const std::string & table_name, const std::string & class_key);
    void SaveGroupPotentialEntryList(const GroupPotentialEntry * group_entry, const std::string & table_name);

    void SaveAtomicChargeEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveAtomicChargeEntrySubList(const ModelObject * model_obj, const std::string & table_name, const std::string & class_key);
    void SaveGroupChargeEntryList(const GroupChargeEntry * group_entry, const std::string & table_name);

    std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>> LoadAtomicPotentialEntryMap(const std::string & table_name);
    void LoadAtomicPotentialEntrySubList(const std::string & table_name, const std::string & class_key, std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>> & entry_map);
    void LoadGroupPotentialEntryList(ModelObject * model_obj, const std::string & class_key, const std::string & table_name);
    void LoadGroupPotentialEntrySubList(ModelObject * model_obj, const std::string & class_key);

    std::unordered_map<int, std::unique_ptr<AtomicChargeEntry>> LoadAtomicChargeEntryMap(const std::string & table_name);
    void LoadAtomicChargeEntrySubList(const std::string & table_name, const std::string & class_key, std::unordered_map<int, std::unique_ptr<AtomicChargeEntry>> & entry_map);
    void LoadGroupChargeEntryList(ModelObject * model_obj, const std::string & class_key, const std::string & table_name);
    void LoadGroupChargeEntrySubList(ModelObject * model_obj, const std::string & class_key);

    bool TableExists(const std::string & table_name) const;
    
};
