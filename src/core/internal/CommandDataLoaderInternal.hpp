#pragma once

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem::command_data_loader {

namespace detail {

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> ProcessTypedFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label,
    std::string_view object_type_name)
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
            + "' as " + std::string(object_type_name) + ": " + ex.what());
    }
}

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> LoadTypedObject(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label,
    std::string_view object_type_name)
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
            + std::string(key_tag) + "' as " + std::string(object_type_name)
            + ": " + ex.what());
    }
}

} // namespace detail

inline std::shared_ptr<ModelObject> ProcessModelFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    return detail::ProcessTypedFile<ModelObject>(
        data_manager,
        path,
        key_tag,
        label,
        "ModelObject");
}

inline std::shared_ptr<MapObject> ProcessMapFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    return detail::ProcessTypedFile<MapObject>(
        data_manager,
        path,
        key_tag,
        label,
        "MapObject");
}

inline std::shared_ptr<MapObject> OptionalProcessMapFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    if (path.empty())
    {
        return nullptr;
    }
    return ProcessMapFile(data_manager, path, key_tag, label);
}

inline std::shared_ptr<ModelObject> LoadModelObject(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label)
{
    return detail::LoadTypedObject<ModelObject>(
        data_manager,
        key_tag,
        label,
        "ModelObject");
}

} // namespace rhbm_gem::command_data_loader
