#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "DataObjectBase.hpp"

template <typename T> struct KDNode;

namespace rhbm_gem {

class GridNode;

class VoxelNode
{
    std::array<float, 3> m_position;
    float m_value;

public:
    VoxelNode(const std::array<float, 3> & position, float value) :
        m_position{ position }, m_value{ value } {}
    const std::array<float, 3> & GetPosition(void) const { return m_position; }
    float GetValue(void) const { return m_value; }
    void SetPosition(const std::array<float, 3> & position) { m_position = position; }
};

class MapObject : public DataObjectBase
{
    std::string m_key_tag;
    int m_thread_size;
    size_t m_voxel_size;
    float m_map_value_mean, m_map_value_min, m_map_value_max, m_map_value_sd;
    std::array<int, 3> m_grid_size;
    std::array<float, 3> m_grid_spacing, m_origin, m_map_length;
    std::array<float, 3> m_overflow, m_underflow, m_upper_bound, m_lower_bound;
    std::unique_ptr<float[]> m_map_value_array;
    std::unique_ptr<::KDNode<GridNode>> m_kd_tree_root;
    std::vector<GridNode> m_grid_node_list;

public:
    MapObject(void);
    MapObject(const std::array<int, 3> & grid_size,
              const std::array<float, 3> & grid_spacing,
              const std::array<float, 3> & origin);
    MapObject(const std::array<int, 3> & grid_size,
              const std::array<float, 3> & grid_spacing,
              const std::array<float, 3> & origin,
              std::unique_ptr<float[]> map_value_array);
    ~MapObject();
    MapObject(const MapObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
    void Update(void) override;
    void Accept(DataObjectVisitorBase * visitor) override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag(void) const override { return m_key_tag; }

    std::array<int, 3> GetGridSize(void) const { return m_grid_size; }
    std::array<float, 3> GetGridSpacing(void) const { return m_grid_spacing; }
    std::array<float, 3> GetOrigin(void) const { return m_origin; }
    std::size_t GetMapValueArraySize(void) const { return m_voxel_size; }
    const float * GetMapValueArray(void) const { return m_map_value_array.get(); }
    std::array<int, 3> GetIndexFromPosition(const std::array<float, 3> & position) const;
    std::array<int, 3> GetGridIndex(size_t global_index) const;
    std::array<float, 3> GetGridPosition(size_t global_index) const;
    size_t GetGlobalIndex(int index_x, int index_y, int index_z) const;
    float GetMapValue(size_t global_index) const;
    float GetMapValue(int index_x, int index_y, int index_z) const;
    float GetMapValueMean(void) const { return m_map_value_mean; }
    float GetMapValueMin(void) const { return m_map_value_min; }
    float GetMapValueMax(void) const { return m_map_value_max; }
    float GetMapValueSD(void) const { return m_map_value_sd; }
    ::KDNode<GridNode> * GetKDTreeRoot(void) const;
    void SetThreadSize(int value) { m_thread_size = value; }
    void SetMapValueArray(std::unique_ptr<float[]> map_value_array);
    void MapValueArrayNormalization(void);
    void BuildKDTreeRoot(void);

private:
    void CheckIndex(int index_x, int index_y, int index_z) const;
    void CheckPosition(const std::array<float, 3> & position) const;
    void CalculateMapValueMean(void);
    void CalculateMapValueMin(void);
    void CalculateMapValueMax(void);
    void CalculateMapValueSD(void);
    void BuildGridNodeList(void);
};

class GridNode
{
    size_t m_grid_index;
    const MapObject * m_map_object;

public:
    GridNode(size_t grid_index, const MapObject * map_object) :
        m_grid_index{ grid_index }, m_map_object{ map_object } {}
    size_t GetGridIndex(void) const { return m_grid_index; }
    std::array<float, 3> GetPosition(void) const
    {
        return m_map_object->GetGridPosition(m_grid_index);
    }
    float GetValue(void) const { return m_map_object->GetMapValue(m_grid_index); }
    void SetGridIndex(size_t grid_index) { m_grid_index = grid_index; }
};

} // namespace rhbm_gem
