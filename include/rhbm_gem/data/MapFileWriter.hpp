#pragma once

#include <memory>
#include <string>
#include <array>
#include <ostream>

namespace rhbm_gem {

class MapObject;
class MapFileFormatBase;

class MapFileWriter
{
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;
    const MapObject * m_map_object;

public:
    explicit MapFileWriter(const std::string & filename, const MapObject * map_object);
    ~MapFileWriter();
    void Write();

};

} // namespace rhbm_gem
