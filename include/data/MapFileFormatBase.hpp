#pragma once

#include <memory>
#include <string>
#include <istream>
#include <ostream>
#include <array>

namespace rhbm_gem {

class MapObject;

class MapFileFormatBase
{
public:
    virtual ~MapFileFormatBase() = default;
    virtual void Read(std::istream & stream, const std::string & source_name) = 0;
    virtual void Write(const MapObject & map_object, std::ostream & stream) = 0;

    /**
     * @brief  Retrieve the voxel data array owned by this object.
     * @return A unique pointer that takes ownership of the data.
     *
     * The caller becomes responsible for the returned array and the
     * internal pointer held by the implementation will be reset.  Calling
     * this method again after ownership has been transferred is allowed
     * but will typically return an empty pointer.
     */
    virtual std::unique_ptr<float[]> GetDataArray() = 0;
    virtual std::array<int, 3> GetGridSize() = 0;
    virtual std::array<float, 3> GetGridSpacing() = 0;
    virtual std::array<float, 3> GetOrigin() = 0;

};

} // namespace rhbm_gem
