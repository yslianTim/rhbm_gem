#pragma once

#include <string>
#include <memory>

class SQLiteWrapper;
class DataObjectBase;
class DataObjectDAOBase;

class DatabaseManager
{
    std::unique_ptr<SQLiteWrapper> m_database;

public:
    explicit DatabaseManager(const std::string & database_path);
    ~DatabaseManager();

    void SaveDataObject(const DataObjectBase * obj);
    std::unique_ptr<DataObjectBase> LoadDataObject(const std::string & key_tag);

    SQLiteWrapper * GetDatabase(void) { return m_database.get(); }
    std::unique_ptr<DataObjectDAOBase> CreateDataObjectDAO(const DataObjectBase * obj);
    
private:
    
    
};