#pragma once

#include <filesystem>
#include <memory>
#include <mutex>
#include <string>

namespace rhbm_gem {

class MapObject;
class SQLiteWrapper;
class ModelObject;

class SQLitePersistence
{
    std::filesystem::path m_database_path;
    std::unique_ptr<SQLiteWrapper> m_database;
    mutable std::mutex m_db_mutex;

public:
    explicit SQLitePersistence(const std::filesystem::path & database_path);
    ~SQLitePersistence();

    SQLitePersistence(const SQLitePersistence &) = delete;
    SQLitePersistence & operator=(const SQLitePersistence &) = delete;

    void SaveModel(const ModelObject & model_object, const std::string & key_tag);
    void SaveMap(const MapObject & map_object, const std::string & key_tag);
    std::unique_ptr<ModelObject> LoadModel(const std::string & key_tag);
    std::unique_ptr<MapObject> LoadMap(const std::string & key_tag);

    const std::filesystem::path & GetDatabasePath() const noexcept { return m_database_path; }
};

} // namespace rhbm_gem
