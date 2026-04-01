#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

template <typename T> struct KDNode;

namespace rhbm_gem {

class MapObject;

class GridNode
{
    size_t m_grid_index;
    const MapObject * m_map_object;

public:
    GridNode(size_t grid_index, const MapObject * map_object) :
        m_grid_index{ grid_index }, m_map_object{ map_object } {}

    size_t GetGridIndex() const { return m_grid_index; }
    std::array<float, 3> GetPosition() const;
    float GetValue() const;
    void SetGridIndex(size_t grid_index) { m_grid_index = grid_index; }
};

class MapSpatialIndex
{
    const MapObject * m_map_object;
    std::unique_ptr<::KDNode<GridNode>> m_kd_tree_root;
    std::vector<GridNode> m_grid_node_list;

public:
    explicit MapSpatialIndex(const MapObject & map_object);
    ~MapSpatialIndex();

    void Invalidate();
    void Build(int thread_size);
    ::KDNode<GridNode> * GetRoot() const;

private:
    void BuildGridNodeList(int thread_size);
};

} // namespace rhbm_gem
