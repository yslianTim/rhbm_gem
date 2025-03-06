#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include "DataObjectDAOBase.hpp"

class SQLiteWrapper;
class ModelObject;
class AtomObject;
class GroupPotentialEntryBase;
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
    void CreateAtomicPotentialEntrySubListTable(const std::string & table_name, const std::string & class_key);
    void CreateGroupPotentialEntryElementClassListTable(const std::string & table_name);
    void CreateGroupPotentialEntryResidueClassListTable(const std::string & table_name);
    void SaveAtomObjectList(const ModelObject * model_obj, const std::string & table_name);
    void SaveAtomicPotentialEntryList(const ModelObject * model_obj, const std::string & table_name);
    void SaveAtomicPotentialEntrySubList(const ModelObject * model_obj, const std::string & table_name, const std::string & class_key);
    void SaveGroupPotentialEntryElementClassList(const GroupPotentialEntryBase * group_entry, const std::string & table_name);
    void SaveGroupPotentialEntryResidueClassList(const GroupPotentialEntryBase * group_entry, const std::string & table_name);
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & table_name);
    
};