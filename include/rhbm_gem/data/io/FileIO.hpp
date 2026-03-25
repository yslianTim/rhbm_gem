#pragma once

#include <filesystem>
#include <memory>

namespace rhbm_gem {

class DataObjectBase;
class ModelObject;
class MapObject;

std::unique_ptr<DataObjectBase> ReadDataObject(const std::filesystem::path & filename);
void WriteDataObject(const std::filesystem::path & filename, const DataObjectBase & data_object);

std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path & filename);
void WriteModel(
    const std::filesystem::path & filename,
    const ModelObject & model_object,
    int model_parameter = 0);

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path & filename);
void WriteMap(const std::filesystem::path & filename, const MapObject & map_object);

} // namespace rhbm_gem
