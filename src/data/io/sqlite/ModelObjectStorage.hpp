#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;

class ModelObjectStorage
{
public:
    static void CreateTables(SQLiteWrapper& database);
    static void Save(SQLiteWrapper& database, const ModelObject& obj, const std::string& key_tag);
    static std::unique_ptr<ModelObject> Load(SQLiteWrapper& database, const std::string& key_tag);
    
};

} // namespace rhbm_gem
