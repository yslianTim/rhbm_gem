#pragma once

#include <string>
#include <memory>
#include <filesystem>
#include <mutex>

#include "internal/sqlite/DatabaseSchemaManager.hpp"

namespace rhbm_gem {

class DataObjectBase;
class SQLiteWrapper;
class ModelObjectStorage;
class MapObjectStorage;

class DatabaseManager
{
    std::filesystem::path m_database_path;
    std::unique_ptr<SQLiteWrapper> m_database;
    std::shared_ptr<ModelObjectStorage> m_model_store;
    std::shared_ptr<MapObjectStorage> m_map_store;
    DatabaseSchemaVersion m_schema_version;
    mutable std::mutex m_db_mutex;  // Protects database operations

public:
    explicit DatabaseManager(const std::filesystem::path & database_path);
    ~DatabaseManager();
    void SaveDataObject(const DataObjectBase * data_object, const std::string & key_tag);
    std::unique_ptr<DataObjectBase> LoadDataObject(const std::string & key_tag);
    const std::filesystem::path & GetDatabasePath() const { return m_database_path; }
    DatabaseSchemaVersion GetSchemaVersion() const { return m_schema_version; }
};

} // namespace rhbm_gem
