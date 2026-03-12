#include "internal/io/sqlite/DataObjectDAOFactoryRegistry.hpp"
#include "internal/io/sqlite/MapObjectDAO.hpp"
#include "internal/io/sqlite/ModelObjectDAO.hpp"

#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

void RegisterDataObjectDaos(DataObjectDAOFactoryRegistry & registry)
{
    registry.RegisterFactory(
        typeid(ModelObject),
        "model",
        [](SQLiteWrapper * database)
        {
            return std::make_unique<ModelObjectDAO>(database);
        });
    registry.RegisterFactory(
        typeid(MapObject),
        "map",
        [](SQLiteWrapper * database)
        {
            return std::make_unique<MapObjectDAO>(database);
        });
}

} // namespace rhbm_gem
