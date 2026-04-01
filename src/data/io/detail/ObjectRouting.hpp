#pragma once

#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/DataObjectDispatch.hpp>

#include <stdexcept>
#include <string>
#include <string_view>

namespace rhbm_gem::io_internal {

enum class TopLevelObjectKind
{
    Model,
    Map
};

inline TopLevelObjectKind ResolveTopLevelObjectKind(
    const DataObjectBase & data_object,
    std::string_view context)
{
    if (AsModelObject(data_object) != nullptr)
    {
        return TopLevelObjectKind::Model;
    }
    if (AsMapObject(data_object) != nullptr)
    {
        return TopLevelObjectKind::Map;
    }
    throw std::runtime_error(
        std::string(context) + ": expected top-level DataObject but got unsupported DataObject type.");
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

} // namespace rhbm_gem::io_internal
