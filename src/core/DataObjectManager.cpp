#include "DataObjectManager.hpp"
#include "FileProcessFactoryBase.hpp"
#include "FileProcessFactoryRegistry.hpp"
#include "FilePathHelper.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectVisitorBase.hpp"
#include "DatabaseManager.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <utility>

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
    std::lock_guard<std::mutex> lock(m_db_mutex);
    if (m_db_manager && m_db_manager->GetDatabasePath() == dbname)
    {
        Logger::Log(LogLevel::Warning,
                    "Database already existed in the path: " + dbname.string()
                    + ", skip re-allocation of DatabaseManager.");
        return;
    }
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
    try
    {
        auto extension{ FilePathHelper::GetExtension(filename) };
        auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(extension) };
        auto data_object{ factory->CreateDataObject(filename) };
        if (!data_object)
        {
            throw std::runtime_error("Failed to create data object");
        }
        data_object->SetKeyTag(key_tag);
        data_object->Display();
        std::shared_ptr<DataObjectBase> shared_object{ std::move(data_object) };
        AddDataObject(key_tag, std::move(shared_object));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error("Failed to process file '" + filename.string() +
                                 "' for key tag '" + key_tag + "': " + ex.what());
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
    auto data_object{ GetDataObject(key_tag) };
    auto extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(extension) };
    factory->OutputDataObject(filename, data_object.get());
}

bool DataObjectManager::AddDataObject(
    const std::string & key_tag, std::shared_ptr<DataObjectBase> data_object)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::AddDataObject() called");
    if (!data_object)
    {
        Logger::Log(LogLevel::Error, "AddDataObject(): nullptr provided for key tag: " + key_tag);
        return false;
    }
    std::lock_guard<std::mutex> lock(m_map_mutex);
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
    std::lock_guard<std::mutex> lock(m_map_mutex);
    return m_data_object_map.find(key_tag) != m_data_object_map.end();
}

void DataObjectManager::LoadDataObject(const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::LoadDataObject() called");
    ScopeTimer timer("DataObjectManager::LoadDataObject");
    std::unique_ptr<DataObjectBase> data_object;
    {
        std::lock_guard<std::mutex> lock(m_db_mutex);
        if (m_db_manager == nullptr)
        {
            throw std::runtime_error("Database manager is not initialized.");
        }
        data_object = m_db_manager->LoadDataObject(key_tag);
    }
    try
    {
        std::shared_ptr<DataObjectBase> shared_object{ std::move(data_object) };
        bool inserted{ AddDataObject(key_tag, std::move(shared_object)) };
        if (inserted == false)
        {
            Logger::Log(LogLevel::Warning,
                        "Data object with key tag: [" + key_tag + "] overwritten or insertion failed.");
        }
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error("Failed to load data object with key tag '" + key_tag + "': "
                                 + ex.what());
    }
}

void DataObjectManager::SaveDataObject(
    const std::string & key_tag, const std::string & renamed_key_tag) const
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::SaveDataObject() called");
    ScopeTimer timer("DataObjectManager::SaveDataObject");
    std::lock_guard<std::mutex> db_lock(m_db_mutex);
    if (m_db_manager == nullptr)
    {
        throw std::runtime_error("Database manager is not initialized.");
    }

    std::shared_ptr<DataObjectBase> data_object;
    {
        std::lock_guard<std::mutex> lock(m_map_mutex);
        auto iter{ m_data_object_map.find(key_tag) };
        if (iter == m_data_object_map.end())
        {
            Logger::Log(LogLevel::Warning,
                        "The data object with key tag: [" + key_tag + "] isn't presented, "
                        "skip saving data object.");
            return;
        }
        data_object = iter->second;
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

    m_db_manager->SaveDataObject(data_object.get(), saved_key_tag);
}

void DataObjectManager::Accept(
    DataObjectVisitorBase * visitor, const std::vector<std::string> & key_tag_list)
{
    Logger::Log(LogLevel::Debug, "DataObjectManager::Accept() called");
    ScopeTimer timer("DataObjectManager::Accept");
    std::vector<DataObjectBase *> data_object_list;
    {
        std::lock_guard<std::mutex> lock(m_map_mutex);
        if (key_tag_list.empty())
        {
            for (auto & [key, data_object] : m_data_object_map)
            {
                if (data_object)
                {
                    data_object_list.push_back(data_object.get());
                }
            }
        }
        else
        {
            for (const auto & key : key_tag_list)
            {
                auto iter{ m_data_object_map.find(key) };
                if (iter == m_data_object_map.end())
                {
                    Logger::Log(LogLevel::Warning, "Cannot find the data object with key tag: " + key);
                    continue;
                }
                if (iter->second)
                {
                    data_object_list.push_back(iter->second.get());
                }
            }
        }
    }
    for (auto * data_object : data_object_list)
    {
        data_object->Accept(visitor);
    }
}

std::shared_ptr<DataObjectBase> DataObjectManager::GetDataObject(
    const std::string & key_tag)
{
    auto const_ptr{ std::as_const(*this).GetDataObject(key_tag) };
    return std::const_pointer_cast<DataObjectBase>(const_ptr);
}

std::shared_ptr<const DataObjectBase> DataObjectManager::GetDataObject(
    const std::string & key_tag) const
{
    std::lock_guard<std::mutex> lock(m_map_mutex);
    auto iter{ m_data_object_map.find(key_tag) };
    if (iter == m_data_object_map.end())
    {
        throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
    }
    return iter->second;
}

DatabaseManager * DataObjectManager::GetDatabaseManager(void) const
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    return m_db_manager.get();
}

const std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> &
DataObjectManager::GetDataObjectMap(void) const
{
    return m_data_object_map;
}
