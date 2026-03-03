#pragma once

#include <memory>
#include <string>
#include <array>
#include <fstream>
#include "FileReaderBase.hpp"
#include "MapFileFormatBase.hpp"

namespace rhbm_gem {

class MapFileReader : public FileReaderBase
{
    bool m_successfully_read_file;
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;

public:
    explicit MapFileReader(const std::string & filename);
    ~MapFileReader() = default;
    void Read(void) override;
    bool IsSuccessfullyRead(void) const override { return m_successfully_read_file; }

    /**
     * @brief  Obtain the map value array read from file.
     *
     * Ownership of the returned pointer is transferred to the caller.
     * After a successful call the reader no longer stores the voxel data
     * and subsequent calls will return an empty pointer.
     */
    std::unique_ptr<float[]> GetMapValueArray(void);
    std::array<int, 3> GetGridSizeArray(void) const;
    std::array<float, 3> GetGridSpacingArray(void) const;
    std::array<float, 3> GetOriginArray(void) const;

private:
    bool ReadHeader(std::ifstream & stream);
    bool ReadMapValueArray(std::ifstream & stream);

};

} // namespace rhbm_gem
