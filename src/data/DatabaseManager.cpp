#include "DatabaseManager.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectDAOBase.hpp"
#include "DataObjectDAOFactoryRegistry.hpp"
#include "DataObjectBase.hpp"
#include "Logger.hpp"

namespace rhbm_gem {

DatabaseManager::DatabaseManager(const std::filesystem::path & database_path) :
    m_database_path{ database_path },
    m_database{ nullptr }
{
    if (m_database_path.empty()) m_database_path = "database.sqlite";
    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
    m_database->Execute(
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
}

DatabaseManager::~DatabaseManager()
{
}

void DatabaseManager::SaveDataObject(
    const DataObjectBase * data_object, const std::string & key_tag)
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    auto dao{ CreateDataObjectDAO(data_object) };
    dao->Save(data_object, key_tag);

    SQLiteWrapper::TransactionGuard transaction(*m_database);
    std::string type_name{
        DataObjectDAOFactoryRegistry::Instance().GetTypeName(std::type_index(typeid(*data_object)))
    };
    m_database->Prepare(
        "INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?) "
        "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;");
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, type_name);
    m_database->StepOnce();
}

std::unique_ptr<DataObjectBase> DatabaseManager::LoadDataObject(
    const std::string & key_tag)
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    m_database->Prepare(
        "SELECT object_type FROM object_metadata WHERE key_tag = ? LIMIT 1;");
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
    const DataObjectBase * data_object)
{
    if (data_object == nullptr)
    {
        throw std::runtime_error("Null data object pointer provided.");
    }
    auto type{ std::type_index(typeid(*data_object)) };
    return CreateDataObjectDAO(DataObjectDAOFactoryRegistry::Instance().GetTypeName(type));
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
