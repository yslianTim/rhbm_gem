#pragma once

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>

#include "DataObjectManager.hpp"

namespace rhbm_gem::command_data_access {

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> ProcessTypedFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    try
    {
        data_manager.ProcessFile(path, std::string(key_tag));
        return data_manager.GetTypedDataObject<TypedDataObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to process " + std::string(label) + " from '" + path.string()
            + "' as " + typeid(TypedDataObject).name() + ": " + ex.what());
    }
}

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> OptionalProcessTypedFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    if (path.empty())
    {
        return nullptr;
    }
    return ProcessTypedFile<TypedDataObject>(data_manager, path, key_tag, label);
}

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> LoadTypedObject(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label)
{
    try
    {
        data_manager.LoadDataObject(std::string(key_tag));
        return data_manager.GetTypedDataObject<TypedDataObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to load " + std::string(label) + " with key tag '"
            + std::string(key_tag) + "': " + ex.what());
    }
}

} // namespace rhbm_gem::command_data_access
