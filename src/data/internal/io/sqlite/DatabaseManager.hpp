#pragma once

#include <string>
#include <memory>
#include <filesystem>
#include <unordered_map>
#include <typeindex>
#include <mutex>

#include "internal/migration/DatabaseSchemaManager.hpp"

namespace rhbm_gem {

class SQLiteWrapper;
class DataObjectBase;
class DataObjectDAOBase;
class DataObjectDAOFactoryRegistry;

class DatabaseManager
{
    std::filesystem::path m_database_path;
    std::unique_ptr<SQLiteWrapper> m_database;
    std::shared_ptr<const DataObjectDAOFactoryRegistry> m_dao_factory_registry;
    std::unordered_map<std::type_index, std::shared_ptr<DataObjectDAOBase>> m_dao_cache;
    DatabaseSchemaVersion m_schema_version;
    mutable std::mutex m_mutex;     // Protects m_dao_cache
    mutable std::mutex m_db_mutex;  // Protects database operations

public:
    DatabaseManager(
        const std::filesystem::path & database_path,
        std::shared_ptr<const DataObjectDAOFactoryRegistry> dao_factory_registry);
    ~DatabaseManager();
    void SaveDataObject(const DataObjectBase * data_object, const std::string & key_tag);
    std::unique_ptr<DataObjectBase> LoadDataObject(const std::string & key_tag);
    const std::filesystem::path & GetDatabasePath() const { return m_database_path; }
    DatabaseSchemaVersion GetSchemaVersion() const { return m_schema_version; }

private:
    std::shared_ptr<DataObjectDAOBase> CreateDataObjectDAO(const std::string & object_type);
};

} // namespace rhbm_gem
