#pragma once

#include <memory>

#include "FileFormatRegistry.hpp"

namespace rhbm_gem {

class ModelFileFormatBase;
class MapFileFormatBase;

std::unique_ptr<ModelFileFormatBase> CreateModelFileFormatBackend(ModelFormatBackend backend);
std::unique_ptr<MapFileFormatBase> CreateMapFileFormatBackend(MapFormatBackend backend);

} // namespace rhbm_gem
