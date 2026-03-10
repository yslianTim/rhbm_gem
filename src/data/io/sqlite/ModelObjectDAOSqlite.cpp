#include "internal/io/sqlite/ModelObjectDAOSqlite.hpp"

#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/dispatch/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "internal/io/sqlite/SQLiteWrapper.hpp"

#include "io/sqlite/ModelAnalysisPersistence.hpp"
#include "io/sqlite/ModelSchemaSql.hpp"
#include "io/sqlite/ModelStructurePersistence.hpp"

#include <stdexcept>
#include <string>

namespace {

void DeleteRowsForKey(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & table_name,
    const std::string & key_tag)
{
    database.Prepare(
        std::string(rhbm_gem::model_io::kDeleteRowsForKeySqlPrefix) +
        table_name +
        std::string(rhbm_gem::model_io::kDeleteRowsForKeySqlSuffix));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.StepOnce();
}

} // namespace

namespace rhbm_gem {

ModelObjectDAOSqlite::ModelObjectDAOSqlite(SQLiteWrapper * db_manager) :
    m_database{ db_manager }
{
}

ModelObjectDAOSqlite::~ModelObjectDAOSqlite()
{
}

void ModelObjectDAOSqlite::EnsureSchema(SQLiteWrapper & database)
{
    for (const auto create_sql : model_io::kCreateModelTableSqlList)
    {
        database.Execute(std::string(create_sql));
    }
}

void ModelObjectDAOSqlite::Save(const DataObjectBase & obj, const std::string & key_tag)
{
    const auto & model_obj{ ExpectModelObject(obj, "ModelObjectDAOSqlite::Save()") };

    for (const auto table_name : model_io::kModelTablesScopedByKey)
    {
        DeleteRowsForKey(*m_database, std::string(table_name), key_tag);
    }

    model_io::SaveStructure(*m_database, model_obj, key_tag);
    model_io::SaveAnalysis(*m_database, model_obj, key_tag);
}

std::unique_ptr<DataObjectBase> ModelObjectDAOSqlite::Load(const std::string & key_tag)
{
    auto model_object{ std::make_unique<ModelObject>() };
    model_io::LoadStructure(*m_database, *model_object, key_tag);
    model_io::LoadAnalysis(*m_database, *model_object, key_tag);
    return model_object;
}

} // namespace rhbm_gem
