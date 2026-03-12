#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <mutex>
#include <functional>

namespace rhbm_gem {

class DataObjectBase;

class DataObjectManager
{
    struct Impl;
    std::unique_ptr<Impl> m_impl;
    std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> m_data_object_map;
    mutable std::mutex m_map_mutex; // protects m_data_object_map
    mutable std::mutex m_db_mutex;  // protects m_impl database state

public:
    struct IterateOptions
    {
        bool deterministic_order{ true };
    };

    DataObjectManager();
    ~DataObjectManager();
    DataObjectManager(const DataObjectManager &) = delete;
    DataObjectManager & operator=(const DataObjectManager &) = delete;
    void ClearDataObjects();
    void SetDatabaseManager(const std::filesystem::path & dbname);
    void ProcessFile(const std::filesystem::path & filename, const std::string & key_tag);
    void ProduceFile(const std::filesystem::path & filename, const std::string & key_tag);
    bool HasDataObject(const std::string & key_tag) const;
    void LoadDataObject(const std::string & key_tag);
    void SaveDataObject(const std::string & key_tag, const std::string & renamed_key_tag="") const;
    void ForEachDataObject(
        const std::function<void(DataObjectBase &)> & callback,
        const std::vector<std::string> & key_tag_list={});
    void ForEachDataObject(
        const std::function<void(DataObjectBase &)> & callback,
        const std::vector<std::string> & key_tag_list,
        const IterateOptions & options);
    void ForEachDataObject(
        const std::function<void(const DataObjectBase &)> & callback,
        const std::vector<std::string> & key_tag_list={}) const;
    void ForEachDataObject(
        const std::function<void(const DataObjectBase &)> & callback,
        const std::vector<std::string> & key_tag_list,
        const IterateOptions & options) const;
    std::shared_ptr<DataObjectBase> GetDataObject(const std::string & key_tag);
    std::shared_ptr<const DataObjectBase> GetDataObject(const std::string & key_tag) const;
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

private:
    bool AddDataObject(const std::string & key_tag, std::shared_ptr<DataObjectBase> data_object);
};

} // namespace rhbm_gem
