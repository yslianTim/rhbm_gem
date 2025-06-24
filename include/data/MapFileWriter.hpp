#pragma once

#include <memory>
#include <string>
#include <array>
#include "FileWriterBase.hpp"

class MapFileFormatBase;
class MapObject;

class MapFileWriter : public FileWriterBase
{
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;
    const MapObject * m_map_object;

public:
    MapFileWriter(const std::string & filename, const MapObject * map_object);
    ~MapFileWriter();
    void Write(void) override;

private:
    void WriteHeader(void);
    void WriteMapValueArray(void);

};
