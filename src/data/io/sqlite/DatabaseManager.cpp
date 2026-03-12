#include "internal/io/sqlite/DatabaseManager.hpp"
#include "internal/migration/DatabaseSchemaManager.hpp"
#include "internal/io/sqlite/MapObjectDAO.hpp"
#include "internal/io/sqlite/ModelObjectDAOSqlite.hpp"
#include "internal/io/sqlite/SQLiteWrapper.hpp"
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/dispatch/DataObjectDispatch.hpp>

namespace rhbm_gem {

DatabaseManager::DatabaseManager(const std::filesystem::path & database_path) :
    m_database_path{ database_path },
    m_database{ nullptr },
    m_model_store{ nullptr },
    m_map_store{ nullptr },
    m_schema_version{ DatabaseSchemaVersion::NormalizedV2 }
{
    if (m_database_path.empty()) m_database_path = "database.sqlite";
    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
    m_model_store = std::make_shared<ModelObjectDAOSqlite>(m_database.get());
    m_map_store = std::make_shared<MapObjectDAO>(m_database.get());
    DatabaseSchemaManager schema_manager{ m_database.get() };
    m_schema_version = schema_manager.EnsureSchema();
}

DatabaseManager::~DatabaseManager()
{
}

void DatabaseManager::SaveDataObject(
    const DataObjectBase * data_object, const std::string & key_tag)
{
    if (data_object == nullptr)
    {
        throw std::runtime_error("Null data object pointer provided.");
    }
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    auto type_name{ GetCatalogTypeName(*data_object) };
    m_database->Prepare(
        "INSERT INTO object_catalog(key_tag, object_type) VALUES (?, ?) "
        "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;");
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, type_name);
    m_database->StepOnce();

    if (type_name == "model")
    {
        m_model_store->Save(*data_object, key_tag);
        return;
    }
    if (type_name == "map")
    {
        m_map_store->Save(*data_object, key_tag);
        return;
    }
    throw std::runtime_error("Unsupported data object type: " + type_name);
}

std::unique_ptr<DataObjectBase> DatabaseManager::LoadDataObject(
    const std::string & key_tag)
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    m_database->Prepare(
        "SELECT object_type FROM object_catalog WHERE key_tag = ? LIMIT 1;");
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    auto rc{ m_database->StepNext() };
    if (rc == SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    else if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    auto type_name{ m_database->GetColumn<std::string>(0) };
    if (type_name == "model")
    {
        return m_model_store->Load(key_tag);
    }
    if (type_name == "map")
    {
        return m_map_store->Load(key_tag);
    }
    throw std::runtime_error("Unsupported data object type: " + type_name);
}

} // namespace rhbm_gem
