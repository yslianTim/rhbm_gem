#pragma once

#include <string>
#include <memory>
#include <filesystem>
#include <unordered_map>
#include <typeindex>
#include <mutex>

namespace rhbm_gem {

class SQLiteWrapper;
class DataObjectBase;
class DataObjectDAOBase;

class DatabaseManager
{
    std::filesystem::path m_database_path;
    std::unique_ptr<SQLiteWrapper> m_database;
    std::unordered_map<std::type_index, std::shared_ptr<DataObjectDAOBase>> m_dao_cache;
    mutable std::mutex m_mutex;     // Protects m_dao_cache
    mutable std::mutex m_db_mutex;  // Protects database operations

public:
    explicit DatabaseManager(const std::filesystem::path & database_path);
    ~DatabaseManager();
    void SaveDataObject(const DataObjectBase * data_object, const std::string & key_tag);
    std::unique_ptr<DataObjectBase> LoadDataObject(const std::string & key_tag);
    const std::filesystem::path & GetDatabasePath(void) const { return m_database_path; }
    SQLiteWrapper * GetDatabase(void) { return m_database.get(); }
    std::shared_ptr<DataObjectDAOBase> CreateDataObjectDAO(const DataObjectBase * data_object);
    std::shared_ptr<DataObjectDAOBase> CreateDataObjectDAO(const std::string & object_type);
    
};

} // namespace rhbm_gem
