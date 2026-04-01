#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

namespace rhbm_gem {

class MapObject;
class ModelObject;

namespace command_detail {

class CommandObjectCache
{
    using ObjectHandle = std::variant<std::shared_ptr<ModelObject>, std::shared_ptr<MapObject>>;

public:
    enum class ObjectKind
    {
        Model,
        Map
    };

    void Clear()
    {
        m_object_map.clear();
    }

    void PutModel(const std::string & key_tag, std::shared_ptr<ModelObject> data_object)
    {
        if (!data_object)
        {
            throw std::runtime_error(
                "CommandObjectCache::Put(): nullptr provided for key tag: " + key_tag);
        }
        m_object_map.insert_or_assign(key_tag, std::move(data_object));
    }

    void PutMap(const std::string & key_tag, std::shared_ptr<MapObject> data_object)
    {
        if (!data_object)
        {
            throw std::runtime_error(
                "CommandObjectCache::Put(): nullptr provided for key tag: " + key_tag);
        }
        m_object_map.insert_or_assign(key_tag, std::move(data_object));
    }

    ObjectKind GetKind(const std::string & key_tag) const
    {
        return std::visit(
            [](const auto & object_handle) -> ObjectKind {
                using HandleType = std::decay_t<decltype(object_handle)>;
                if constexpr (std::is_same_v<HandleType, std::shared_ptr<ModelObject>>)
                {
                    return ObjectKind::Model;
                }
                return ObjectKind::Map;
            },
            RequireEntry(key_tag));
    }

    std::shared_ptr<ModelObject> GetModel(const std::string & key_tag)
    {
        return GetTyped<ModelObject>(key_tag);
    }

    std::shared_ptr<MapObject> GetMap(const std::string & key_tag)
    {
        return GetTyped<MapObject>(key_tag);
    }

private:
    const ObjectHandle & RequireEntry(const std::string & key_tag) const
    {
        const auto iter{ m_object_map.find(key_tag) };
        if (iter == m_object_map.end())
        {
            throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
        }
        return iter->second;
    }

    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> GetTyped(const std::string & key_tag) const
    {
        const auto & object_handle{ RequireEntry(key_tag) };
        if (const auto typed_object{ std::get_if<std::shared_ptr<TypedDataObject>>(&object_handle) })
        {
            return *typed_object;
        }
        throw std::runtime_error("Invalid data type for " + key_tag);
    }

    std::unordered_map<std::string, ObjectHandle> m_object_map;
};

} // namespace command_detail

} // namespace rhbm_gem
