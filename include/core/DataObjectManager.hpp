#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>
#include <vector>

class DatabaseManager;
class FileProcessFactoryBase;
class DataObjectBase;
class DataObjectVisitorBase;

class DataObjectManager
{
    std::unique_ptr<DatabaseManager> m_db_manager;
    std::unordered_map<std::string, std::unique_ptr<DataObjectBase>> m_data_object_map;

public:
    DataObjectManager(void);
    ~DataObjectManager();
    DataObjectManager(const DataObjectManager &) = delete;
    DataObjectManager & operator=(const DataObjectManager &) = delete;
    DataObjectManager(DataObjectManager &&) noexcept = default;
    DataObjectManager & operator=(DataObjectManager &&) noexcept = default;
    void SetDatabaseManager(const std::filesystem::path & dbname);
    void ProcessFile(const std::filesystem::path & filename, const std::string & key_tag);
    void ProduceFile(const std::filesystem::path & filename, const std::string & key_tag);
    bool AddDataObject(const std::string & key_tag, std::unique_ptr<DataObjectBase> data_object);
    bool HasDataObject(const std::string & key_tag) const;
    void RemoveDataObject(const std::string & key_tag);
    void LoadDataObject(const std::string & key_tag);
    void SaveDataObject(const std::string & key_tag, const std::string & renamed_key_tag="") const;
    void Accept(DataObjectVisitorBase * visitor);
    void Accept(DataObjectVisitorBase * visitor, const std::vector<std::string> & key_list);
    void PrintDataObjectInfo(const std::string & key_tag) const;
    DataObjectBase * GetDataObjectPtr(const std::string & key_tag);
    const DataObjectBase * GetDataObjectPtr(const std::string & key_tag) const;
    template <typename TypedDataObject>
    TypedDataObject * GetTypedDataObjectPtr(const std::string & key_tag)
    {
        auto base_object{ GetDataObjectPtr(key_tag) };
        auto typed_object{ dynamic_cast<TypedDataObject *>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }
    template <typename TypedDataObject>
    const TypedDataObject * GetTypedDataObjectPtr(const std::string & key_tag) const
    {
        auto base_object{ GetDataObjectPtr(key_tag) };
        auto typed_object{ dynamic_cast<const TypedDataObject *>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }
    const std::unordered_map<std::string, std::unique_ptr<DataObjectBase>> & GetDataObjectMap(void) const { return m_data_object_map; }

};
