#pragma once

#include <memory>
#include <string>
#include <array>
#include "FileReaderBase.hpp"

class MapFileFormatBase;

class MapFileReader : public FileReaderBase
{
    bool m_successfully_read_file;
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;

public:
    MapFileReader(const std::string & filename);
    ~MapFileReader();
    void Read(void) override;
    bool IsSuccessfullyRead(void) override { return m_successfully_read_file; }
    std::unique_ptr<float[]> GetMapValueArray(void);
    std::array<int, 3> GetGridSizeArray(void);
    std::array<float, 3> GetGridSpacingArray(void);
    std::array<float, 3> GetOriginArray(void);

private:
    void ReadHeader(void);
    void ReadMapValueArray(void);

};