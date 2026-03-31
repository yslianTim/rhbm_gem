#include "internal/command/CommandDataLoader.hpp"

#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <memory>
#include <stdexcept>
#include <string>

namespace rhbm_gem {
namespace {

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

} // namespace

std::shared_ptr<ModelObject> LoadModelFromFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    return ProcessTypedFile<ModelObject>(
        data_manager,
        path,
        key_tag,
        label,
        "ModelObject");
}

std::shared_ptr<MapObject> LoadMapFromFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    return ProcessTypedFile<MapObject>(
        data_manager,
        path,
        key_tag,
        label,
        "MapObject");
}

std::shared_ptr<MapObject> MaybeLoadMapFromFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    if (path.empty())
    {
        return nullptr;
    }
    return LoadMapFromFile(data_manager, path, key_tag, label);
}

std::shared_ptr<ModelObject> LoadModelFromDatabase(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label)
{
    return LoadTypedObject<ModelObject>(
        data_manager,
        key_tag,
        label,
        "ModelObject");
}

} // namespace rhbm_gem
