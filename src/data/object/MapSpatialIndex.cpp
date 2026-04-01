#include "data/object/MapSpatialIndex.hpp"

#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {

std::array<float, 3> GridNode::GetPosition() const
{
    return m_map_object->GetGridPosition(m_grid_index);
}

float GridNode::GetValue() const
{
    return m_map_object->GetMapValue(m_grid_index);
}

MapSpatialIndex::MapSpatialIndex(const MapObject & map_object) :
    m_map_object{ &map_object }
{
}

MapSpatialIndex::~MapSpatialIndex() = default;

void MapSpatialIndex::Invalidate()
{
    m_kd_tree_root.reset();
    m_grid_node_list.clear();
}

void MapSpatialIndex::Build(int thread_size)
{
    ScopeTimer timer("MapSpatialIndex::Build");
    if (m_kd_tree_root != nullptr)
    {
        return;
    }

    BuildGridNodeList(thread_size);
    if (m_grid_node_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No grids were found from the map.");
        return;
    }

    Logger::Log(
        LogLevel::Info,
        " /- Building KD-Tree from " +
            std::to_string(m_map_object->GetMapValueArraySize()) + " voxels...");
    m_kd_tree_root = KDTreeAlgorithm<GridNode>::BuildKDTree(
        m_grid_node_list,
        0,
        thread_size,
        true);
}

KDNode<GridNode> * MapSpatialIndex::GetRoot() const
{
    if (m_kd_tree_root == nullptr)
    {
        Logger::Log(
            LogLevel::Error,
            "MapSpatialIndex::GetRoot -> spatial index is not built yet. Call Build() first.");
        return nullptr;
    }
    return m_kd_tree_root.get();
}

void MapSpatialIndex::BuildGridNodeList(int thread_size)
{
    m_grid_node_list.clear();
    const auto voxel_size{ m_map_object->GetMapValueArraySize() };
    m_grid_node_list.reserve(voxel_size);

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        std::vector<GridNode> thread_local_list;
        thread_local_list.reserve(voxel_size / static_cast<size_t>(thread_size) + 1);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < voxel_size; i++)
        {
            thread_local_list.emplace_back(i, m_map_object);
        }

        #pragma omp critical
        {
            m_grid_node_list.insert(
                m_grid_node_list.end(),
                thread_local_list.begin(),
                thread_local_list.end());
        }
    }
#else
    for (size_t i = 0; i < voxel_size; i++)
    {
        m_grid_node_list.emplace_back(i, m_map_object);
    }
#endif
}

} // namespace rhbm_gem
