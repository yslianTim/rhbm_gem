#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class SQLiteWrapper;
class MapObject;

class MapObjectStorage {
    SQLiteWrapper* m_database;

  public:
    explicit MapObjectStorage(SQLiteWrapper* database);
    ~MapObjectStorage();
    static void CreateTables(SQLiteWrapper& database);
    void Save(const MapObject& map_object, const std::string& key_tag);
    std::unique_ptr<MapObject> Load(const std::string& key_tag);
};

} // namespace rhbm_gem
