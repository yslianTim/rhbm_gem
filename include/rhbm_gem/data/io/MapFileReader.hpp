#pragma once

#include <memory>
#include <string>
#include <array>

namespace rhbm_gem {

class MapFileFormatBase;
class FileFormatRegistry;

class MapFileReader
{
    std::string m_file_path;
    std::unique_ptr<MapFileFormatBase> m_file_format_helper;
    bool m_has_loaded_data;

public:
    explicit MapFileReader(const std::string & filename);
    MapFileReader(const std::string & filename, const FileFormatRegistry & file_format_registry);
    ~MapFileReader();
    void Read();

    /**
     * @brief  Obtain the map value array read from file.
     *
     * Ownership of the returned pointer is transferred to the caller.
     * After a successful call the reader no longer stores the voxel data
     * and subsequent calls will return an empty pointer.
     */
    std::unique_ptr<float[]> GetMapValueArray();
    std::array<int, 3> GetGridSizeArray() const;
    std::array<float, 3> GetGridSpacingArray() const;
    std::array<float, 3> GetOriginArray() const;

};

} // namespace rhbm_gem
