#include "internal/DatabaseManager.hpp"
#include "internal/DatabaseSchemaManager.hpp"
#include "internal/SQLiteWrapper.hpp"
#include "internal/DataObjectDAOBase.hpp"
#include "internal/DataObjectDAOFactoryRegistry.hpp"
#include <rhbm_gem/data/DataObjectBase.hpp>
#include <rhbm_gem/data/DataObjectDispatch.hpp>
#include <rhbm_gem/utils/Logger.hpp>

namespace rhbm_gem {

DatabaseManager::DatabaseManager(const std::filesystem::path & database_path) :
    m_database_path{ database_path },
    m_database{ nullptr },
    m_schema_version{ DatabaseSchemaVersion::NormalizedV2 }
{
    if (m_database_path.empty()) m_database_path = "database.sqlite";
    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
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

    auto dao{ CreateDataObjectDAO(type_name) };
    dao->Save(*data_object, key_tag);
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
    auto dao{ CreateDataObjectDAO(type_name) };
    return dao->Load(key_tag);
}

std::shared_ptr<DataObjectDAOBase> DatabaseManager::CreateDataObjectDAO(
    const std::string & object_type)
{
    auto type{ DataObjectDAOFactoryRegistry::Instance().GetTypeIndex(object_type) };
    std::lock_guard<std::mutex> lock(m_mutex);
    auto iter{ m_dao_cache.find(type) };
    if (iter != m_dao_cache.end())
    {
        return iter->second;
    }
    auto dao_unique{ DataObjectDAOFactoryRegistry::Instance().CreateDAO(type, m_database.get()) };
    std::shared_ptr<DataObjectDAOBase> dao{ std::move(dao_unique) };
    m_dao_cache.emplace(type, dao);
    return dao;
}

} // namespace rhbm_gem
