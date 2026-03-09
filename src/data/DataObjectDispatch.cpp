#include <rhbm_gem/data/DataObjectDispatch.hpp>

#include <stdexcept>
#include <string>
#include <string_view>

#include <rhbm_gem/data/AtomObject.hpp>
#include <rhbm_gem/data/BondObject.hpp>
#include <rhbm_gem/data/DataObjectBase.hpp>
#include <rhbm_gem/data/MapObject.hpp>
#include <rhbm_gem/data/ModelObject.hpp>

namespace rhbm_gem {
namespace {

const char * ResolveDataObjectTypeName(const DataObjectBase & data_object)
{
    if (dynamic_cast<const AtomObject *>(&data_object) != nullptr) return "AtomObject";
    if (dynamic_cast<const BondObject *>(&data_object) != nullptr) return "BondObject";
    if (dynamic_cast<const ModelObject *>(&data_object) != nullptr) return "ModelObject";
    if (dynamic_cast<const MapObject *>(&data_object) != nullptr) return "MapObject";
    return "unsupported DataObject type";
}

} // namespace

ModelObject * AsModelObject(DataObjectBase & data_object) noexcept
{
    return dynamic_cast<ModelObject *>(&data_object);
}

const ModelObject * AsModelObject(const DataObjectBase & data_object) noexcept
{
    return dynamic_cast<const ModelObject *>(&data_object);
}

MapObject * AsMapObject(DataObjectBase & data_object) noexcept
{
    return dynamic_cast<MapObject *>(&data_object);
}

const MapObject * AsMapObject(const DataObjectBase & data_object) noexcept
{
    return dynamic_cast<const MapObject *>(&data_object);
}

const ModelObject & ExpectModelObject(
    const DataObjectBase & data_object,
    std::string_view context)
{
    if (auto typed_object{ AsModelObject(data_object) })
    {
        return *typed_object;
    }
    throw std::runtime_error(
        std::string(context) + ": expected ModelObject but got "
        + ResolveDataObjectTypeName(data_object) + ".");
}

ModelObject & ExpectModelObject(
    DataObjectBase & data_object,
    std::string_view context)
{
    if (auto typed_object{ AsModelObject(data_object) })
    {
        return *typed_object;
    }
    throw std::runtime_error(
        std::string(context) + ": expected ModelObject but got "
        + ResolveDataObjectTypeName(data_object) + ".");
}

const MapObject & ExpectMapObject(
    const DataObjectBase & data_object,
    std::string_view context)
{
    if (auto typed_object{ AsMapObject(data_object) })
    {
        return *typed_object;
    }
    throw std::runtime_error(
        std::string(context) + ": expected MapObject but got "
        + ResolveDataObjectTypeName(data_object) + ".");
}

MapObject & ExpectMapObject(
    DataObjectBase & data_object,
    std::string_view context)
{
    if (auto typed_object{ AsMapObject(data_object) })
    {
        return *typed_object;
    }
    throw std::runtime_error(
        std::string(context) + ": expected MapObject but got "
        + ResolveDataObjectTypeName(data_object) + ".");
}

std::string GetCatalogTypeName(const DataObjectBase & data_object)
{
    if (dynamic_cast<const ModelObject *>(&data_object) != nullptr)
    {
        return "model";
    }
    if (dynamic_cast<const MapObject *>(&data_object) != nullptr)
    {
        return "map";
    }
    if (dynamic_cast<const AtomObject *>(&data_object) != nullptr)
    {
        throw std::runtime_error("GetCatalogTypeName(): AtomObject is not a top-level catalog type.");
    }
    if (dynamic_cast<const BondObject *>(&data_object) != nullptr)
    {
        throw std::runtime_error("GetCatalogTypeName(): BondObject is not a top-level catalog type.");
    }
    throw std::runtime_error("GetCatalogTypeName(): no catalog type resolved.");
}

} // namespace rhbm_gem
