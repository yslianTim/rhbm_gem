#pragma once

#include <filesystem>
#include <memory>

namespace rhbm_gem {

class ModelObject;
class MapObject;

std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path & filename);
void WriteModel(
    const std::filesystem::path & filename,
    const ModelObject & model_object,
    int model_parameter = 0);

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path & filename);
void WriteMap(const std::filesystem::path & filename, const MapObject & map_object);

} // namespace rhbm_gem
