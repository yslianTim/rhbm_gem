#pragma once

#include <rhbm_gem/data/io/FileIO.hpp>

#include <memory>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <mutex>
#include <functional>
#include <type_traits>

namespace rhbm_gem {

class DataObjectBase;
class SQLitePersistence;

class DataObjectManager
{
    enum class StoreResult
    {
        Inserted,
        Replaced,
        RejectedNull
    };

    std::unique_ptr<SQLitePersistence> m_db_manager;
    std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> m_data_object_map;
    mutable std::mutex m_map_mutex; // protects m_data_object_map

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
    void OpenDatabase(const std::filesystem::path & dbname);
    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> ImportFileAs(
        const std::filesystem::path & filename,
        const std::string & key_tag)
    {
        auto data_object{ ReadFileAs<TypedDataObject>(filename) };
        data_object->SetKeyTag(key_tag);
        std::shared_ptr<TypedDataObject> shared_object{ std::move(data_object) };
        StoreDataObject(key_tag, shared_object);
        return shared_object;
    }
    void ExportToFile(
        const std::filesystem::path & filename,
        const std::string & key_tag) const;
    bool Contains(const std::string & key_tag) const;
    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> LoadFromDatabaseAs(const std::string & key_tag)
    {
        auto base_object{ LoadFromDatabase(key_tag) };
        auto typed_object{ std::dynamic_pointer_cast<TypedDataObject>(base_object) };
        if (!typed_object)
        {
            throw std::runtime_error("Invalid data type for " + key_tag);
        }
        StoreDataObject(key_tag, base_object);
        return typed_object;
    }
    void SaveToDatabase(
        const std::string & key_tag,
        const std::string & persisted_key = "") const;
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
    std::shared_ptr<TypedDataObject> Require(const std::string & key_tag)
    {
        auto base_object{ GetDataObject(key_tag) };
        auto typed_object{ std::dynamic_pointer_cast<TypedDataObject>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }
    template <typename TypedDataObject>
    std::shared_ptr<const TypedDataObject> Require(const std::string & key_tag) const
    {
        auto base_object{ GetDataObject(key_tag) };
        auto typed_object{ std::dynamic_pointer_cast<const TypedDataObject>(base_object) };
        if (!typed_object) throw std::runtime_error("Invalid data type for " + key_tag);
        return typed_object;
    }

private:
    template <typename TypedDataObject>
    static std::unique_ptr<TypedDataObject> ReadFileAs(const std::filesystem::path & filename)
    {
        if constexpr (std::is_same_v<TypedDataObject, ModelObject>)
        {
            return ReadModel(filename);
        }
        else if constexpr (std::is_same_v<TypedDataObject, MapObject>)
        {
            return ReadMap(filename);
        }
        else
        {
            static_assert(
                std::is_same_v<TypedDataObject, ModelObject>
                    || std::is_same_v<TypedDataObject, MapObject>,
                "ImportFileAs only supports ModelObject or MapObject.");
        }
    }

    std::shared_ptr<DataObjectBase> LoadFromDatabase(const std::string & key_tag);
    StoreResult StoreDataObject(
        const std::string & key_tag,
        std::shared_ptr<DataObjectBase> data_object);
};

} // namespace rhbm_gem
