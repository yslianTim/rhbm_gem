#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include "DataObjectDAOBase.hpp"

class DatabaseManager;
class SQLiteWrapper;
class ModelObject;
class AtomObject;

class ModelObjectDAO : public DataObjectDAOBase
{
    DatabaseManager & m_db_manager;

public:
    ModelObjectDAO(DatabaseManager & db_manager);
    ~ModelObjectDAO();

    void Save(const DataObjectBase * obj, const std::string & key_tag="") override;
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) override;
    
private:
    void CreateModelObjectListTable(SQLiteWrapper * db);
    void CreateAtomObjectListTable(SQLiteWrapper * db, const std::string & table_name);
    void SaveAtomObjectList(const ModelObject * model_obj);
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & key_tag);
    
};