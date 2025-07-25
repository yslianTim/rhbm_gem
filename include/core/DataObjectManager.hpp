#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <shared_mutex>

class DatabaseManager;
class FileProcessFactoryBase;
class DataObjectBase;
class DataObjectVisitorBase;

class DataObjectManager
{
    std::shared_ptr<DatabaseManager> m_db_manager;
    std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> m_data_object_map;
    mutable std::shared_mutex m_mutex; // Protects m_data_object_map and m_db_manager

public:
    DataObjectManager(void);
    ~DataObjectManager();
    DataObjectManager(const DataObjectManager &) = delete;
    DataObjectManager & operator=(const DataObjectManager &) = delete;
    DataObjectManager(DataObjectManager &&) noexcept = default;
    DataObjectManager & operator=(DataObjectManager &&) noexcept = default;
    void SetDatabaseManager(const std::filesystem::path & dbname);
    void SetDatabaseManager(std::shared_ptr<DatabaseManager> manager);
    void ProcessFile(const std::filesystem::path & filename, const std::string & key_tag);
    void ProduceFile(const std::filesystem::path & filename, const std::string & key_tag);
    bool AddDataObject(const std::string & key_tag, std::shared_ptr<DataObjectBase> data_object);
    bool HasDataObject(const std::string & key_tag) const;
    void RemoveDataObject(const std::string & key_tag);
    void LoadDataObject(const std::string & key_tag);
    void SaveDataObject(const std::string & key_tag, const std::string & renamed_key_tag="") const;
    void Accept(DataObjectVisitorBase * visitor);
    void Accept(DataObjectVisitorBase * visitor, const std::vector<std::string> & key_list);
    void PrintDataObjectInfo(const std::string & key_tag) const;
    std::shared_ptr<DataObjectBase> GetDataObject(const std::string & key_tag);
    std::shared_ptr<const DataObjectBase> GetDataObject(const std::string & key_tag) const;
    DatabaseManager * GetDatabaseManagerPtr(void);
    const DatabaseManager * GetDatabaseManagerPtr(void) const;
    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> GetTypedDataObject(const std::string & key_tag)
    {
        auto base_object{ GetDataObject(key_tag) };
        auto typed_object{ std::dynamic_pointer_cast<TypedDataObject>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }
    template <typename TypedDataObject>
    std::shared_ptr<const TypedDataObject> GetTypedDataObject(const std::string & key_tag) const
    {
        auto base_object{ GetDataObject(key_tag) };
        auto typed_object{ std::dynamic_pointer_cast<const TypedDataObject>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }
    std::unordered_map<std::string, std::shared_ptr<const DataObjectBase>> GetDataObjectMap(void) const;

};
