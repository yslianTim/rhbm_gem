#pragma once

#include <string>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;

namespace model_io {

void SaveAnalysis(SQLiteWrapper & database, const ModelObject & model_obj, const std::string & key_tag);
void LoadAnalysis(SQLiteWrapper & database, ModelObject & model_obj, const std::string & key_tag);

} // namespace model_io

} // namespace rhbm_gem
