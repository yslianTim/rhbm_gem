#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;

class ModelObjectStorage {
    SQLiteWrapper* m_database;

  public:
    explicit ModelObjectStorage(SQLiteWrapper* db_manager);
    ~ModelObjectStorage();

    static void CreateTables(SQLiteWrapper& database);

    void Save(const ModelObject& obj, const std::string& key_tag);
    std::unique_ptr<ModelObject> Load(const std::string& key_tag);
};

} // namespace rhbm_gem
