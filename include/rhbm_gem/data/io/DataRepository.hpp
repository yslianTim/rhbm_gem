#pragma once

#include <filesystem>
#include <memory>
#include <string>

namespace rhbm_gem {

class MapObject;
class ModelObject;
class SQLitePersistence;

class DataRepository
{
    std::unique_ptr<SQLitePersistence> m_database;

public:
    explicit DataRepository(const std::filesystem::path & database_path);
    ~DataRepository();
    DataRepository(const DataRepository &) = delete;
    DataRepository & operator=(const DataRepository &) = delete;

    std::unique_ptr<ModelObject> LoadModel(const std::string & key_tag) const;
    std::unique_ptr<MapObject> LoadMap(const std::string & key_tag) const;
    void SaveModel(const ModelObject & model_object, const std::string & key_tag) const;
    void SaveMap(const MapObject & map_object, const std::string & key_tag) const;
    
};

} // namespace rhbm_gem
