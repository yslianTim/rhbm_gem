#include "DatabaseManager.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectDAOBase.hpp"
#include "DataObjectDAOFactoryRegistry.hpp"
#include "ModelObject.hpp"
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
    // To-Do : Implement case for MapObject type case
    auto model_object{ std::make_unique<ModelObject>() };
    auto dao{ CreateDataObjectDAO(model_object.get()) };
    return dao->Load(key_tag);
}

std::unique_ptr<DataObjectDAOBase> DatabaseManager::CreateDataObjectDAO(
    const DataObjectBase * data_object)
{
    Logger::Log(LogLevel::Debug, "DatabaseManager::CreateDataObjectDAO() called");
    if (data_object == nullptr)
    {
        throw std::runtime_error("Null data object pointer provided.");
    }
    auto type{ std::type_index(typeid(*data_object)) };
    return DataObjectDAOFactoryRegistry::Instance().CreateDAO(type, m_database.get());
}
