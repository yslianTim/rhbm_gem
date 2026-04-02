#include "data/detail/MapSpatialIndex.hpp"

#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {

namespace {

class GridNode
{
    size_t m_grid_index;
    const MapObject * m_map_object;

public:
    GridNode(size_t grid_index, const MapObject * map_object) :
        m_grid_index{ grid_index }, m_map_object{ map_object } {}

    size_t GetGridIndex() const { return m_grid_index; }
    std::array<float, 3> GetPosition() const
    {
        return m_map_object->GetGridPosition(m_grid_index);
    }
};

} // namespace

struct MapSpatialIndex::Impl
{
    std::unique_ptr<::KDNode<GridNode>> kd_tree_root;
    std::vector<GridNode> grid_node_list;
};

MapSpatialIndex::MapSpatialIndex(const MapObject & map_object, int thread_size) :
    m_map_object{ &map_object }, m_thread_size{ thread_size > 0 ? thread_size : 1 }
{
}

MapSpatialIndex::~MapSpatialIndex() = default;

void MapSpatialIndex::Build() const
{
    ScopeTimer timer("MapSpatialIndex::Build");
    if (m_impl != nullptr && m_impl->kd_tree_root != nullptr)
    {
        return;
    }

    if (m_impl == nullptr)
    {
        m_impl = std::make_unique<Impl>();
    }

    m_impl->grid_node_list.clear();
    const auto voxel_size{ m_map_object->GetMapValueArraySize() };
    m_impl->grid_node_list.reserve(voxel_size);

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_thread_size)
    {
        std::vector<GridNode> thread_local_list;
        thread_local_list.reserve(voxel_size / static_cast<size_t>(m_thread_size) + 1);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < voxel_size; i++)
        {
            thread_local_list.emplace_back(i, m_map_object);
        }

        #pragma omp critical
        {
            m_impl->grid_node_list.insert(
                m_impl->grid_node_list.end(),
                thread_local_list.begin(),
                thread_local_list.end());
        }
    }
#else
    for (size_t i = 0; i < voxel_size; i++)
    {
        m_impl->grid_node_list.emplace_back(i, m_map_object);
    }
#endif

    if (m_impl->grid_node_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No grids were found from the map.");
        return;
    }

    Logger::Log(
        LogLevel::Info,
        " /- Building KD-Tree from " +
            std::to_string(m_map_object->GetMapValueArraySize()) + " voxels...");
    m_impl->kd_tree_root = KDTreeAlgorithm<GridNode>::BuildKDTree(
        m_impl->grid_node_list,
        0,
        m_thread_size,
        true);
}

void MapSpatialIndex::CollectGridIndicesInRange(
    const std::array<float, 3> & center,
    float radius,
    std::vector<size_t> & grid_index_list) const
{
    grid_index_list.clear();
    Build();
    if (m_impl == nullptr || m_impl->kd_tree_root == nullptr)
    {
        Logger::Log(
            LogLevel::Error,
            "MapSpatialIndex::CollectGridIndicesInRange -> spatial index is not available.");
        return;
    }

    MapObject query_map_object({ 1, 1, 1 }, { 1.0f, 1.0f, 1.0f }, center);
    GridNode query_node(0, &query_map_object);
    std::vector<GridNode *> in_range_grid_node_list;
    KDTreeAlgorithm<GridNode>::RangeSearch(
        m_impl->kd_tree_root.get(),
        &query_node,
        radius,
        in_range_grid_node_list);

    grid_index_list.reserve(in_range_grid_node_list.size());
    for (const auto * grid_node : in_range_grid_node_list)
    {
        grid_index_list.emplace_back(grid_node->GetGridIndex());
    }
}

} // namespace rhbm_gem
