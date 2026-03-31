#pragma once

#include <filesystem>
#include <memory>
#include <mutex>
#include <string>

namespace rhbm_gem {

class DataObjectBase;
class SQLiteWrapper;
class ModelObjectStorage;

class SQLitePersistence
{
    std::filesystem::path m_database_path;
    std::unique_ptr<SQLiteWrapper> m_database;
    std::unique_ptr<ModelObjectStorage> m_model_store;
    mutable std::mutex m_db_mutex;

public:
    explicit SQLitePersistence(const std::filesystem::path & database_path);
    ~SQLitePersistence();

    SQLitePersistence(const SQLitePersistence &) = delete;
    SQLitePersistence & operator=(const SQLitePersistence &) = delete;

    void Save(const DataObjectBase & data_object, const std::string & key_tag);
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag);

    const std::filesystem::path & GetDatabasePath() const noexcept { return m_database_path; }
};

} // namespace rhbm_gem
