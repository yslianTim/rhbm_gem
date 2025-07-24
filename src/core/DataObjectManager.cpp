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

DataObjectManager::~DataObjectManager()
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::~DataObjectManager() called");
}

void DataObjectManager::SetDatabaseManager(const std::filesystem::path & dbname)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::SetDatabaseManager() called");
    m_db_manager = std::make_unique<DatabaseManager>(dbname);
    if (!m_db_manager->GetDatabase())
    {
        throw std::runtime_error("Failed to initialize database manager with path: " + dbname.string());
    }
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

    bool inserted{ AddDataObject(key_tag, std::move(data_object)) };
    if (inserted == false)
    {
        Logger::Log(LogLevel::Warning,
                    "Data object with key tag: [" + key_tag + "] overwritten or insertion failed.");
    }
}

void DataObjectManager::ProduceFile(
    const std::filesystem::path & filename, const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::ProduceFile() called");
    ScopeTimer timer("DataObjectManager::ProduceFile");
    if (HasDataObject(key_tag) == false)
    {
        Logger::Log(LogLevel::Warning,
                    "The data object with key tag: [" + key_tag + "] isn't presented, "
                    "no file will be produced.");
        return;
    }
    auto data_object{ GetDataObjectPtr(key_tag) };
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(file_extension) };
    factory->OutputDataObject(filename, data_object);
};

bool DataObjectManager::AddDataObject(
    const std::string & key_tag, std::unique_ptr<DataObjectBase> data_object)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::AddDataObject() called");
    if (!data_object)
    {
        Logger::Log(LogLevel::Error, "AddDataObject(): nullptr provided for key tag: " + key_tag);
        return false;
    }
    auto result{ m_data_object_map.insert_or_assign(key_tag, std::move(data_object)) };
    if (!result.second)
    {
        Logger::Log(LogLevel::Warning,
                    "The key tag: [" + key_tag + "] already presented in the data object map, "
                    "the data object has been replaced.");
    }
    return result.second;
}

bool DataObjectManager::HasDataObject(const std::string & key_tag) const
{
    return m_data_object_map.find(key_tag) != m_data_object_map.end();
}

void DataObjectManager::RemoveDataObject(const std::string & key_tag)
{
    m_data_object_map.erase(key_tag);
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
    bool inserted{ AddDataObject(key_tag, std::move(data_object)) };
    if (inserted == false)
    {
        Logger::Log(LogLevel::Warning,
                    "Data object with key tag: [" + key_tag + "] overwritten or insertion failed.");
    }
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
    if (HasDataObject(key_tag) == false)
    {
        Logger::Log(LogLevel::Warning,
                    "The data object with key tag: [" + key_tag + "] isn't presented, "
                    "skip saving data object.");
        return;
    }

    auto saved_key_tag{ renamed_key_tag.empty() ? key_tag : renamed_key_tag };
    if (renamed_key_tag != "")
    {
        Logger::Log(LogLevel::Info,
                    "The data object with key tag: [" + key_tag + "] will be renamed to: [" +
                    renamed_key_tag + "] and saved into database: " +
                    m_db_manager->GetDatabasePath().string());
    }
    else
    {
        Logger::Log(LogLevel::Info,
                    "The data object with key tag: [" + key_tag + "] will be saved into database: " +
                    m_db_manager->GetDatabasePath().string());
    }

    m_db_manager->SaveDataObject(GetDataObjectPtr(key_tag), saved_key_tag);
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
        if (HasDataObject(key) == false)
        {
            Logger::Log(LogLevel::Warning, "Cannot find the data object with key tag: " + key);
            continue;
        }
        auto data_object{ GetDataObjectPtr(key) };
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
        GetDataObjectPtr(key_tag)->Display();
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Warning, ex.what());
    }
}

DataObjectBase * DataObjectManager::GetDataObjectPtr(const std::string & key_tag)
{
    if (HasDataObject(key_tag) == false)
    {
        throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
    }
    return m_data_object_map.at(key_tag).get();
}

const DataObjectBase * DataObjectManager::GetDataObjectPtr(const std::string & key_tag) const
{
    if (HasDataObject(key_tag) == false)
    {
        throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
    }
    return m_data_object_map.at(key_tag).get();
}
