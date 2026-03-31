#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include "sqlite/SQLitePersistence.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <utility>

namespace rhbm_gem {

namespace {

template <typename SharedPtrType>
std::vector<SharedPtrType> BuildDataObjectSnapshot(
    const std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> & data_object_map,
    std::mutex & map_mutex,
    const std::vector<std::string> & key_tag_list,
    bool deterministic_order)
{
    std::vector<SharedPtrType> data_object_list;
    std::lock_guard<std::mutex> lock(map_mutex);
    if (key_tag_list.empty())
    {
        if (deterministic_order)
        {
            std::vector<std::string> key_list;
            key_list.reserve(data_object_map.size());
            for (const auto & [key, data_object] : data_object_map)
            {
                (void)data_object;
                key_list.push_back(key);
            }
            std::sort(key_list.begin(), key_list.end());
            for (const auto & key : key_list)
            {
                auto iter{ data_object_map.find(key) };
                if (iter != data_object_map.end() && iter->second)
                {
                    data_object_list.push_back(iter->second);
                }
            }
        }
        else
        {
            for (const auto & [key, data_object] : data_object_map)
            {
                (void)key;
                if (data_object)
                {
                    data_object_list.push_back(data_object);
                }
            }
        }
        return data_object_list;
    }

    for (const auto & key : key_tag_list)
    {
        auto iter{ data_object_map.find(key) };
        if (iter == data_object_map.end())
        {
            Logger::Log(LogLevel::Warning, "Cannot find the data object with key tag: " + key);
            continue;
        }
        if (iter->second)
        {
            data_object_list.push_back(iter->second);
        }
    }
    return data_object_list;
}

std::shared_ptr<DataObjectBase> LookupDataObject(
    const std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> & data_object_map,
    std::mutex & map_mutex,
    const std::string & key_tag)
{
    std::lock_guard<std::mutex> lock(map_mutex);
    const auto iter{ data_object_map.find(key_tag) };
    if (iter == data_object_map.end())
    {
        throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
    }
    return iter->second;
}

} // namespace

DataObjectManager::DataObjectManager() = default;

DataObjectManager::~DataObjectManager() = default;

void DataObjectManager::ClearDataObjects()
{
    std::lock_guard<std::mutex> lock(m_map_mutex);
    m_data_object_map.clear();
}

void DataObjectManager::OpenDatabase(const std::filesystem::path & dbname)
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    if (m_db_manager && m_db_manager->GetDatabasePath() == dbname)
    {
        Logger::Log(LogLevel::Warning,
                    "Database already existed in the path: " + dbname.string()
                    + ", skip re-allocation of SQLitePersistence.");
        return;
    }
    m_db_manager = std::make_unique<SQLitePersistence>(dbname);
}

void DataObjectManager::LoadFileIntoMemory(
    const std::filesystem::path & filename, const std::string & key_tag)
{
    try
    {
        auto data_object{ ReadDataObject(filename) };
        data_object->SetKeyTag(key_tag);
        std::shared_ptr<DataObjectBase> shared_object{ std::move(data_object) };
        AddDataObject(key_tag, std::move(shared_object));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error("Failed to load file '" + filename.string() +
                                 "' for key tag '" + key_tag + "': " + ex.what());
    }
}

void DataObjectManager::WriteMemoryObjectToFile(
    const std::filesystem::path & filename, const std::string & key_tag) const
{
    try
    {
        auto data_object{ GetDataObject(key_tag) };
        WriteDataObject(filename, *data_object);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error("Failed to write stored data object to file '"
                                 + filename.string() + "' for key tag '"
                                 + key_tag + "': " + ex.what());
    }
}

bool DataObjectManager::AddDataObject(const std::string & key_tag, std::shared_ptr<DataObjectBase> data_object)
{
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
    std::lock_guard<std::mutex> db_lock(m_db_mutex);
    if (m_db_manager == nullptr)
    {
        throw std::runtime_error("Database manager is not initialized.");
    }

    std::shared_ptr<DataObjectBase> data_object;
    try
    {
        data_object = LookupDataObject(m_data_object_map, m_map_mutex, key_tag);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error("Failed to save data object with key tag '" + key_tag
                                 + "': " + ex.what());
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

void DataObjectManager::ForEachDataObject(
    const std::function<void(DataObjectBase &)> & callback,
    const std::vector<std::string> & key_tag_list)
{
    ForEachDataObject(callback, key_tag_list, IterateOptions{});
}

void DataObjectManager::ForEachDataObject(
    const std::function<void(DataObjectBase &)> & callback,
    const std::vector<std::string> & key_tag_list,
    const IterateOptions & options)
{
    if (!callback)
    {
        throw std::runtime_error("ForEachDataObject(): callback is empty.");
    }

    auto data_object_list{ BuildDataObjectSnapshot<std::shared_ptr<DataObjectBase>>(
        m_data_object_map, m_map_mutex, key_tag_list, options.deterministic_order) };
    for (auto & data_object : data_object_list)
    {
        callback(*data_object);
    }
}

void DataObjectManager::ForEachDataObject(
    const std::function<void(const DataObjectBase &)> & callback,
    const std::vector<std::string> & key_tag_list) const
{
    ForEachDataObject(callback, key_tag_list, IterateOptions{});
}

void DataObjectManager::ForEachDataObject(
    const std::function<void(const DataObjectBase &)> & callback,
    const std::vector<std::string> & key_tag_list,
    const IterateOptions & options) const
{
    if (!callback)
    {
        throw std::runtime_error("ForEachDataObject() const: callback is empty.");
    }

    auto data_object_list{ BuildDataObjectSnapshot<std::shared_ptr<const DataObjectBase>>(
        m_data_object_map, m_map_mutex, key_tag_list, options.deterministic_order) };
    for (const auto & data_object : data_object_list)
    {
        callback(*data_object);
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
    return LookupDataObject(m_data_object_map, m_map_mutex, key_tag);
}

} // namespace rhbm_gem
