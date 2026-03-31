#pragma once

#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <stdexcept>
#include <string>
#include <string_view>

namespace rhbm_gem::io_detail {

enum class TopLevelObjectKind
{
    Model,
    Map
};

inline const char * ResolveTopLevelObjectTypeName(const DataObjectBase & data_object) noexcept
{
    if (dynamic_cast<const ModelObject *>(&data_object) != nullptr) return "ModelObject";
    if (dynamic_cast<const MapObject *>(&data_object) != nullptr) return "MapObject";
    return "unsupported DataObject type";
}

inline TopLevelObjectKind ResolveTopLevelObjectKind(
    const DataObjectBase & data_object,
    std::string_view context)
{
    if (dynamic_cast<const ModelObject *>(&data_object) != nullptr)
    {
        return TopLevelObjectKind::Model;
    }
    if (dynamic_cast<const MapObject *>(&data_object) != nullptr)
    {
        return TopLevelObjectKind::Map;
    }
    throw std::runtime_error(
        std::string(context) + ": expected top-level DataObject but got "
        + ResolveTopLevelObjectTypeName(data_object) + ".");
}

inline std::string_view GetCatalogTypeName(TopLevelObjectKind kind) noexcept
{
    switch (kind)
    {
    case TopLevelObjectKind::Model:
        return "model";
    case TopLevelObjectKind::Map:
        return "map";
    }
    return "unknown";
}

inline TopLevelObjectKind ParseCatalogTypeName(std::string_view type_name)
{
    if (type_name == "model")
    {
        return TopLevelObjectKind::Model;
    }
    if (type_name == "map")
    {
        return TopLevelObjectKind::Map;
    }
    throw std::runtime_error("Unsupported data object type: " + std::string(type_name));
}

} // namespace rhbm_gem::io_detail
