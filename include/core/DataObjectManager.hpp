#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <mutex>

namespace rhbm_gem {

class DatabaseManager;
class FileProcessFactoryBase;
class FileProcessFactoryResolver;
class DataObjectBase;
class DataObjectVisitorBase;

class DataObjectManager
{
    std::unique_ptr<DatabaseManager> m_db_manager;
    std::shared_ptr<const FileProcessFactoryResolver> m_file_factory_resolver;
    std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> m_data_object_map;
    mutable std::mutex m_map_mutex; // protects m_data_object_map
    mutable std::mutex m_db_mutex;  // protects m_db_manager

public:
    DataObjectManager();
    explicit DataObjectManager(std::shared_ptr<const FileProcessFactoryResolver> file_factory_resolver);
    ~DataObjectManager();
    void ClearDataObjects();
    void SetDatabaseManager(const std::filesystem::path & dbname);
    void ProcessFile(const std::filesystem::path & filename, const std::string & key_tag);
    void ProduceFile(const std::filesystem::path & filename, const std::string & key_tag);
    bool HasDataObject(const std::string & key_tag) const;
    void LoadDataObject(const std::string & key_tag);
    void SaveDataObject(const std::string & key_tag, const std::string & renamed_key_tag="") const;
    void Accept(DataObjectVisitorBase * visitor, const std::vector<std::string> & key_tag_list={});
    std::shared_ptr<DataObjectBase> GetDataObject(const std::string & key_tag);
    std::shared_ptr<const DataObjectBase> GetDataObject(const std::string & key_tag) const;
    DatabaseManager * GetDatabaseManager() const;
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
    const std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> & GetDataObjectMap() const;

private:
    bool AddDataObject(const std::string & key_tag, std::shared_ptr<DataObjectBase> data_object);
    std::unique_ptr<FileProcessFactoryBase> CreateFileFactory(
        const std::filesystem::path & filename) const;
};

} // namespace rhbm_gem
