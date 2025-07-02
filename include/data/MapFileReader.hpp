#pragma once

#include <memory>
#include <string>
#include <array>
#include <fstream>
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
    bool IsSuccessfullyRead(void) const override { return m_successfully_read_file; }
    std::unique_ptr<float[]> GetMapValueArray(void);
    std::array<int, 3> GetGridSizeArray(void) const;
    std::array<float, 3> GetGridSpacingArray(void) const;
    std::array<float, 3> GetOriginArray(void) const;

private:
    bool ReadHeader(std::ifstream & stream);
    bool ReadMapValueArray(std::ifstream & stream);

};
