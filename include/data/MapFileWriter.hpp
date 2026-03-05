#pragma once

#include <memory>
#include <string>
#include <array>
#include <ostream>
#include "MapFileFormatBase.hpp"

namespace rhbm_gem {

class MapObject;

class MapFileWriter
{
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;
    const MapObject * m_map_object;

public:
    explicit MapFileWriter(const std::string & filename, const MapObject * map_object);
    ~MapFileWriter() = default;
    void Write();

};

} // namespace rhbm_gem
