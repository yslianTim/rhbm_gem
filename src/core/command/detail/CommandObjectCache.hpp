#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace rhbm_gem {

class DataObjectBase;

namespace command_detail {

class CommandObjectCache
{
public:
    void Clear()
    {
        m_object_map.clear();
    }

    void Put(const std::string & key_tag, std::shared_ptr<DataObjectBase> data_object)
    {
        if (!data_object)
        {
            throw std::runtime_error(
                "CommandObjectCache::Put(): nullptr provided for key tag: " + key_tag);
        }
        m_object_map.insert_or_assign(key_tag, std::move(data_object));
    }

    std::shared_ptr<DataObjectBase> Get(const std::string & key_tag)
    {
        auto const_ptr{ std::as_const(*this).Get(key_tag) };
        return std::const_pointer_cast<DataObjectBase>(const_ptr);
    }

    std::shared_ptr<const DataObjectBase> Get(const std::string & key_tag) const
    {
        const auto iter{ m_object_map.find(key_tag) };
        if (iter == m_object_map.end())
        {
            throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
        }
        return iter->second;
    }

    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> Require(const std::string & key_tag)
    {
        auto base_object{ Get(key_tag) };
        auto typed_object{ std::dynamic_pointer_cast<TypedDataObject>(base_object) };
        if (!typed_object)
        {
            throw std::runtime_error("Invalid data type for " + key_tag);
        }
        return typed_object;
    }

private:
    std::unordered_map<std::string, std::shared_ptr<DataObjectBase>> m_object_map;
};

} // namespace command_detail

} // namespace rhbm_gem
