#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <istream>
#include <ostream>
#include <array>

namespace rhbm_gem {

class MapFileFormatBase
{
public:
    virtual ~MapFileFormatBase() = default;
    virtual void InitHeader(void) = 0;
    virtual void LoadHeader(std::istream & stream) = 0;
    virtual void SaveHeader(std::ostream & stream) = 0;
    virtual void PrintHeader(void) const = 0;
    virtual void LoadDataArray(std::istream & stream) = 0;
    virtual void SaveDataArray(const float * data, size_t size, std::ostream & stream) = 0;

    /**
     * @brief  Retrieve the voxel data array owned by this object.
     * @return A unique pointer that takes ownership of the data.
     *
     * The caller becomes responsible for the returned array and the
     * internal pointer held by the implementation will be reset.  Calling
     * this method again after ownership has been transferred is allowed
     * but will typically return an empty pointer.
     */
    virtual std::unique_ptr<float[]> GetDataArray(void) = 0;
    virtual std::array<int, 3> GetGridSize(void) = 0;
    virtual std::array<float, 3> GetGridSpacing(void) = 0;
    virtual std::array<float, 3> GetOrigin(void) = 0;
    virtual void SetHeader(const std::array<int, 3> & grid_size,
                           const std::array<float, 3> & grid_spacing,
                           const std::array<float, 3> & origin) = 0;

};

} // namespace rhbm_gem
