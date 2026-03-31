#pragma once

#include <filesystem>
#include <memory>
#include <string_view>

namespace rhbm_gem {

class DataObjectManager;
class MapObject;
class ModelObject;

std::shared_ptr<ModelObject> LoadModelFromFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label);
std::shared_ptr<MapObject> LoadMapFromFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label);
std::shared_ptr<MapObject> MaybeLoadMapFromFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label);
std::shared_ptr<ModelObject> LoadModelFromDatabase(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label);

} // namespace rhbm_gem
