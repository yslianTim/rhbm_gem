#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace rhbm_gem {

class MapObject;

class MapSpatialIndex
{
    struct Impl;
    const MapObject * m_map_object;
    int m_thread_size;
    mutable std::unique_ptr<Impl> m_impl;

public:
    explicit MapSpatialIndex(const MapObject & map_object, int thread_size = 1);
    ~MapSpatialIndex();

    void Build() const;
    void CollectGridIndicesInRange(
        const std::array<float, 3> & center,
        float radius,
        std::vector<size_t> & grid_index_list) const;
};

} // namespace rhbm_gem
