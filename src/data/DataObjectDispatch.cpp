#include "DataObjectDispatch.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "DataObjectBase.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"

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

template <typename ExpectedType>
const ExpectedType & ExpectDataObjectType(
    const DataObjectBase & data_object,
    std::string_view context,
    const char * expected_type)
{
    if (auto typed_object{ dynamic_cast<const ExpectedType *>(&data_object) })
    {
        return *typed_object;
    }
    throw std::runtime_error(
        std::string(context) + ": expected " + expected_type + " but got "
        + ResolveDataObjectTypeName(data_object) + ".");
}

} // namespace

const ModelObject & ExpectModelObject(
    const DataObjectBase & data_object,
    std::string_view context)
{
    return ExpectDataObjectType<ModelObject>(data_object, context, "ModelObject");
}

ModelObject & ExpectModelObject(
    DataObjectBase & data_object,
    std::string_view context)
{
    const auto & const_ref{
        ExpectModelObject(static_cast<const DataObjectBase &>(data_object), context) };
    return const_cast<ModelObject &>(const_ref);
}

const MapObject & ExpectMapObject(
    const DataObjectBase & data_object,
    std::string_view context)
{
    return ExpectDataObjectType<MapObject>(data_object, context, "MapObject");
}

MapObject & ExpectMapObject(
    DataObjectBase & data_object,
    std::string_view context)
{
    const auto & const_ref{
        ExpectMapObject(static_cast<const DataObjectBase &>(data_object), context) };
    return const_cast<MapObject &>(const_ref);
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
