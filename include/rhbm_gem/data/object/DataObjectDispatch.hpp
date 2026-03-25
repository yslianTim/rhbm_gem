#pragma once

#include <string>
#include <string_view>

namespace rhbm_gem {

class DataObjectBase;
class ModelObject;
class MapObject;

ModelObject * AsModelObject(DataObjectBase & data_object) noexcept;
const ModelObject * AsModelObject(const DataObjectBase & data_object) noexcept;
MapObject * AsMapObject(DataObjectBase & data_object) noexcept;
const MapObject * AsMapObject(const DataObjectBase & data_object) noexcept;

const ModelObject & ExpectModelObject(
    const DataObjectBase & data_object,
    std::string_view context);
ModelObject & ExpectModelObject(
    DataObjectBase & data_object,
    std::string_view context);

const MapObject & ExpectMapObject(
    const DataObjectBase & data_object,
    std::string_view context);
MapObject & ExpectMapObject(
    DataObjectBase & data_object,
    std::string_view context);

std::string GetCatalogTypeName(const DataObjectBase & data_object);

} // namespace rhbm_gem
