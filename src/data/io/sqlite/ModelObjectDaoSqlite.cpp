#include "internal/ModelObjectDAOv2.hpp"

#include <rhbm_gem/data/DataObjectBase.hpp>
#include <rhbm_gem/data/DataObjectDispatch.hpp>
#include <rhbm_gem/data/ModelObject.hpp>
#include "internal/SQLiteWrapper.hpp"

#include "model_io/ModelAnalysisPersistence.hpp"
#include "model_io/ModelSchemaSql.hpp"
#include "model_io/ModelStructurePersistence.hpp"

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

ModelObjectDAOv2::ModelObjectDAOv2(SQLiteWrapper * db_manager) :
    m_database{ db_manager }
{
}

ModelObjectDAOv2::~ModelObjectDAOv2()
{
}

void ModelObjectDAOv2::EnsureSchema(SQLiteWrapper & database)
{
    for (const auto create_sql : model_io::kCreateModelTableSqlList)
    {
        database.Execute(std::string(create_sql));
    }
}

void ModelObjectDAOv2::Save(const DataObjectBase & obj, const std::string & key_tag)
{
    const auto & model_obj{ ExpectModelObject(obj, "ModelObjectDAOv2::Save()") };

    for (const auto table_name : model_io::kModelTablesScopedByKey)
    {
        DeleteRowsForKey(*m_database, std::string(table_name), key_tag);
    }

    model_io::SaveStructure(*m_database, model_obj, key_tag);
    model_io::SaveAnalysis(*m_database, model_obj, key_tag);
}

std::unique_ptr<DataObjectBase> ModelObjectDAOv2::Load(const std::string & key_tag)
{
    auto model_object{ std::make_unique<ModelObject>() };
    model_io::LoadStructure(*m_database, *model_object, key_tag);
    model_io::LoadAnalysis(*m_database, *model_object, key_tag);
    return model_object;
}

} // namespace rhbm_gem
