#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>

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
    DataObjectManager(const std::filesystem::path & dbname);
    ~DataObjectManager();
    void ProcessFile(const std::string & filename, const std::string & key_tag);
    void ProduceFile(const std::string & filename, const std::string & key_tag);
    void AddDataObject(const std::string & key_tag, std::unique_ptr<DataObjectBase> data_object);
    void LoadDataObject(const std::string & key_tag);
    void SaveDataObject(const std::string & key_tag, const std::string & renamed_key_tag="") const;
    void Accept(DataObjectVisitorBase * visitor);
    void PrintDataObjectInfo(const std::string & key_tag) const;
    DataObjectBase * GetDataObjectPtr(const std::string & key_tag);
    template <typename TypedDataObject>
    TypedDataObject * GetTypedDataObjectPtr(const std::string & key_tag)
    {
        auto base_object{ GetDataObjectPtr(key_tag) };
        auto typed_object{ dynamic_cast<TypedDataObject *>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }
    const std::unordered_map<std::string, std::unique_ptr<DataObjectBase>> & GetDataObjectMap(void) { return m_data_object_map; }

private:
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(const std::string & file_extension);
};
