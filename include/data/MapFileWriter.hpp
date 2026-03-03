#pragma once

#include <memory>
#include <string>
#include <array>
#include <ostream>
#include "FileWriterBase.hpp"
#include "MapFileFormatBase.hpp"

namespace rhbm_gem {

class MapObject;

class MapFileWriter : public FileWriterBase
{
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;
    const MapObject * m_map_object;

public:
    explicit MapFileWriter(const std::string & filename, const MapObject * map_object);
    ~MapFileWriter() = default;
    void Write(void) override;

private:
    void WriteHeader(std::ostream & stream);
    void WriteMapValueArray(std::ostream & stream);

};

} // namespace rhbm_gem
