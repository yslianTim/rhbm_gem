#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>

namespace rhbm_gem {

class MapObject
{
    std::string m_key_tag;
    size_t m_voxel_size;
    float m_map_value_mean, m_map_value_min, m_map_value_max, m_map_value_sd;
    std::array<int, 3> m_grid_size;
    std::array<float, 3> m_grid_spacing, m_origin, m_map_length;
    std::array<float, 3> m_overflow, m_underflow, m_upper_bound, m_lower_bound;
    std::unique_ptr<float[]> m_map_value_array;

public:
    MapObject();
    MapObject(const std::array<int, 3> & grid_size,
              const std::array<float, 3> & grid_spacing,
              const std::array<float, 3> & origin);
    MapObject(const std::array<int, 3> & grid_size,
              const std::array<float, 3> & grid_spacing,
              const std::array<float, 3> & origin,
              std::unique_ptr<float[]> map_value_array);
    ~MapObject();
    MapObject(const MapObject & other);
    void SetKeyTag(const std::string & label) { m_key_tag = label; }
    std::string GetKeyTag() const { return m_key_tag; }

    std::array<int, 3> GetGridSize() const { return m_grid_size; }
    std::array<float, 3> GetGridSpacing() const { return m_grid_spacing; }
    std::array<float, 3> GetOrigin() const { return m_origin; }
    std::size_t GetMapValueArraySize() const { return m_voxel_size; }
    const float * GetMapValueArray() const { return m_map_value_array.get(); }
    std::array<int, 3> GetIndexFromPosition(const std::array<float, 3> & position) const;
    std::array<int, 3> GetGridIndex(size_t global_index) const;
    std::array<float, 3> GetGridPosition(size_t global_index) const;
    float GetMapValue(size_t global_index) const;
    float GetMapValue(int index_x, int index_y, int index_z) const;
    float GetMapValueMean() const { return m_map_value_mean; }
    float GetMapValueMin() const { return m_map_value_min; }
    float GetMapValueMax() const { return m_map_value_max; }
    float GetMapValueSD() const { return m_map_value_sd; }
    void SetMapValueArray(std::unique_ptr<float[]> map_value_array);
    void MapValueArrayNormalization();

private:
    void RecomputeStatistics();
    void CheckIndex(int index_x, int index_y, int index_z) const;
    void CheckPosition(const std::array<float, 3> & position) const;
    size_t GetGlobalIndex(int index_x, int index_y, int index_z) const;
    void CalculateMapValueMean();
    void CalculateMapValueMin();
    void CalculateMapValueMax();
    void CalculateMapValueSD();
};

} // namespace rhbm_gem
