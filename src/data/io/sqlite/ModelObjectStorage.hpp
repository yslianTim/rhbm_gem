#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;

namespace model_storage {

void CreateTables(SQLiteWrapper & database);
void Save(SQLiteWrapper & database, const ModelObject & obj, const std::string & key_tag);
std::unique_ptr<ModelObject> Load(SQLiteWrapper & database, const std::string & key_tag);

} // namespace model_storage

} // namespace rhbm_gem
