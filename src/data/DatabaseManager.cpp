#include "DatabaseManager.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectDAOBase.hpp"
#include "DataObjectDAOFactoryRegistry.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "Logger.hpp"

DatabaseManager::DatabaseManager(const std::filesystem::path & database_path) :
    m_database_path{ database_path },
    m_database{ nullptr }
{
    Logger::Log(LogLevel::Debug, "DatabaseManager::DatabaseManager() called");
    if (m_database_path.empty()) m_database_path = "database.sqlite";
    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
}

DatabaseManager::~DatabaseManager()
{
    Logger::Log(LogLevel::Debug, "DatabaseManager::~DatabaseManager() called");
}

void DatabaseManager::SaveDataObject(
    const DataObjectBase * data_object, const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DatabaseManager::SaveDataObject() called");
    auto dao{ CreateDataObjectDAO(data_object) };
    dao->Save(data_object, key_tag);
}

std::unique_ptr<DataObjectBase> DatabaseManager::LoadDataObject(
    const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DatabaseManager::LoadDataObject() called");
    auto table_exists = [this](const std::string & table_name) {
        m_database->Prepare(
            "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;");
        SQLiteWrapper::StatementGuard guard(*m_database);
        m_database->Bind<std::string>(1, table_name);
        auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepRow())
        {
            return true;
        }
        else if (rc == SQLiteWrapper::StepDone())
        {
            return false;
        }
        else
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }
    };

    auto has_entry = [this](const std::string & table, const std::string & key) {
        m_database->Prepare("SELECT key_tag FROM " + table + " WHERE key_tag = ? LIMIT 1;");
        SQLiteWrapper::StatementGuard guard(*m_database);
        m_database->Bind<std::string>(1, key);
        auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepRow()) return true;
        if (rc == SQLiteWrapper::StepDone()) return false;
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    };

    bool is_map{ false };
    bool is_model{ false };
    if (table_exists("map_list") && has_entry("map_list", key_tag))
    {
        is_map = true;
    }
    else if (table_exists("model_list") && has_entry("model_list", key_tag))
    {
        is_model = true;
    }
    else
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }

    if (is_map)
    {
        auto map_object{ std::make_unique<MapObject>() };
        auto dao{ CreateDataObjectDAO(map_object.get()) };
        return dao->Load(key_tag);
    }
    else if (is_model)
    {
        auto model_object{ std::make_unique<ModelObject>() };
        auto dao{ CreateDataObjectDAO(model_object.get()) };
        return dao->Load(key_tag);
    }
    else
    {
        return nullptr;
    }
}

DataObjectDAOBase * DatabaseManager::CreateDataObjectDAO(
    const DataObjectBase * data_object)
{
    Logger::Log(LogLevel::Debug, "DatabaseManager::CreateDataObjectDAO() called");
    if (data_object == nullptr)
    {
        throw std::runtime_error("Null data object pointer provided.");
    }
    auto type{ std::type_index(typeid(*data_object)) };
    std::lock_guard<std::mutex> lock(m_mutex);
    auto iter{ m_dao_cache.find(type) };
    if (iter != m_dao_cache.end())
    {
        return iter->second.get();
    }
    auto dao{ DataObjectDAOFactoryRegistry::Instance().CreateDAO(type, m_database.get()) };
    auto dao_ptr{ dao.get() };
    m_dao_cache.emplace(type, std::move(dao));
    return dao_ptr;
}
