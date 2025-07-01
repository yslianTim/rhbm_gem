#pragma once

#include <string>
#include <memory>
#include <filesystem>

class SQLiteWrapper;
class DataObjectBase;
class DataObjectDAOBase;

class DatabaseManager
{
    std::filesystem::path m_database_path;
    std::unique_ptr<SQLiteWrapper> m_database;

public:
    explicit DatabaseManager(const std::filesystem::path & database_path);
    ~DatabaseManager();

    void SaveDataObject(const DataObjectBase * obj);
    std::unique_ptr<DataObjectBase> LoadDataObject(const std::string & key_tag);

    const std::filesystem::path & GetDatabasePath(void) const { return m_database_path; }
    SQLiteWrapper * GetDatabase(void) { return m_database.get(); }
    std::unique_ptr<DataObjectDAOBase> CreateDataObjectDAO(const DataObjectBase * obj);
    
};
