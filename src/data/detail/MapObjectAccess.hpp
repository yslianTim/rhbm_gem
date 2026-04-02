#pragma once

#include <rhbm_gem/data/object/MapObject.hpp>

namespace rhbm_gem {

class MapSpatialIndex;

class MapObjectAccess
{
public:
    static MapSpatialIndex & SpatialIndex(MapObject & map_object)
    {
        return *map_object.m_spatial_index;
    }

    static const MapSpatialIndex & SpatialIndex(const MapObject & map_object)
    {
        return *map_object.m_spatial_index;
    }

    static void SetThreadSize(MapObject & map_object, int value)
    {
        map_object.m_thread_size = value;
    }

    static int ThreadSize(const MapObject & map_object)
    {
        return map_object.m_thread_size;
    }
};

} // namespace rhbm_gem
