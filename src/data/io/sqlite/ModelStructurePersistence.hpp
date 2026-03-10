#pragma once

#include <string>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;

namespace model_io {

void SaveStructure(SQLiteWrapper & database, const ModelObject & model_obj, const std::string & key_tag);
void LoadStructure(SQLiteWrapper & database, ModelObject & model_obj, const std::string & key_tag);

} // namespace model_io

} // namespace rhbm_gem
