#include "DataObjectManager.hpp"
#include "FileProcessFactoryBase.hpp"
#include "FileProcessFactoryRegistry.hpp"
#include "FilePathHelper.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectVisitorBase.hpp"
#include "DatabaseManager.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

DataObjectManager::DataObjectManager(void) :
    m_db_manager{ nullptr }
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::DataObjectManager() called");
}

DataObjectManager::DataObjectManager(const std::filesystem::path & dbname) :
    m_db_manager{ std::make_unique<DatabaseManager>(dbname) }
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::DataObjectManager() with database path called");
    if (!m_db_manager->GetDatabase())
    {
        throw std::runtime_error("Failed to initialize database manager with path: " + dbname.string());
    }
}

DataObjectManager::~DataObjectManager()
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::~DataObjectManager() called");
}

void DataObjectManager::ProcessFile(
    const std::filesystem::path & filename, const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::ProcessFile() called");
    ScopeTimer timer("DataObjectManager::ProcessFile");
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(file_extension) };
    auto data_object{ factory->CreateDataObject(filename) };
    if (data_object == nullptr)
    {
        throw std::runtime_error("Failed to create data object");
    }
    data_object->SetKeyTag(key_tag);
    data_object->Display();

    AddDataObject(key_tag, std::move(data_object));
}

void DataObjectManager::ProduceFile(
    const std::filesystem::path & filename, const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::ProduceFile() called");
    ScopeTimer timer("DataObjectManager::ProduceFile");
    if (m_data_object_map.find(key_tag) == m_data_object_map.end())
    {
        Logger::Log(LogLevel::Warning,
                    "The data object with key tag: [" + key_tag + "] isn't presented, "
                    "no file will be produced.");
        return;
    }
    auto data_object{ m_data_object_map.at(key_tag).get() };
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(file_extension) };
    factory->OutputDataObject(filename, data_object);
};

void DataObjectManager::AddDataObject(
    const std::string & key_tag, std::unique_ptr<DataObjectBase> data_object)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::AddDataObject() called");
    if (!data_object)
    {
        Logger::Log(LogLevel::Error, "AddDataObject(): nullptr provided for key tag: " + key_tag);
        return;
    }
    if (m_data_object_map.find(key_tag) != m_data_object_map.end())
    {
        Logger::Log(LogLevel::Warning,
                    "The key tag: [" + key_tag + "] already presented in the data object map, "
                    "this data object will be replaced.");
    }
    m_data_object_map.insert_or_assign(key_tag, std::move(data_object));
}

void DataObjectManager::LoadDataObject(const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::LoadDataObject() called");
    ScopeTimer timer("DataObjectManager::LoadDataObject");
    if (m_db_manager == nullptr)
    {
        throw std::runtime_error("Database manager is not initialized.");
    }
    auto data_object{ m_db_manager->LoadDataObject(key_tag) };
    AddDataObject(key_tag, std::move(data_object));
}

void DataObjectManager::SaveDataObject(
    const std::string & key_tag, const std::string & renamed_key_tag) const
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::SaveDataObject() called");
    ScopeTimer timer("DataObjectManager::SaveDataObject");
    if (m_db_manager == nullptr)
    {
        throw std::runtime_error("Database manager is not initialized.");
    }
    if (m_data_object_map.find(key_tag) == m_data_object_map.end())
    {
        Logger::Log(LogLevel::Warning,
                    "The data object with key tag: [" + key_tag + "] isn't presented, "
                    "skip saving data object.");
        return;
    }

    if (renamed_key_tag != "")
    {
        Logger::Log(LogLevel::Info,
                    "The data object with key tag: [" + key_tag + "] will be renamed to: [" +
                    renamed_key_tag + "] and saved into database: " +
                    m_db_manager->GetDatabasePath().string());
        m_data_object_map.at(key_tag)->SetKeyTag(renamed_key_tag);
    }
    else
    {
        Logger::Log(LogLevel::Info, 
                    "The data object with key tag: [" + key_tag + "] will be saved into database: " +
                    m_db_manager->GetDatabasePath().string());
    }

    m_db_manager->SaveDataObject(m_data_object_map.at(key_tag).get());
}

void DataObjectManager::Accept(DataObjectVisitorBase * visitor)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::Accept() called");
    ScopeTimer timer("DataObjectManager::Accept");
    for (auto & [key, data_object] : m_data_object_map)
    {
        if (data_object)
        {
            data_object->Accept(visitor);
        }
    }
}

void DataObjectManager::Accept(DataObjectVisitorBase * visitor, const std::vector<std::string> & key_list)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::Accept() with key list called");
    ScopeTimer timer("DataObjectManager::Accept");
    for (const auto & key : key_list)
    {
        auto iter{ m_data_object_map.find(key) };
        if (iter == m_data_object_map.end())
        {
            Logger::Log(LogLevel::Warning, "Cannot find the data object with key tag: " + key);
            continue;
        }
        auto & data_object{ iter->second };
        if (data_object)
        {
            data_object->Accept(visitor);
        }
    }
}

void DataObjectManager::PrintDataObjectInfo(const std::string & key_tag) const
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::PrintDataObjectInfo() called");
    try
    {
        m_data_object_map.at(key_tag)->Display();
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Warning, ex.what());
    }
}

DataObjectBase * DataObjectManager::GetDataObjectPtr(const std::string & key_tag) const
{
    if (m_data_object_map.find(key_tag) == m_data_object_map.end())
    {
        throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
    }
    return m_data_object_map.at(key_tag).get();
}
