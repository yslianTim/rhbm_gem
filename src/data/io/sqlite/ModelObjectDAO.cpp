#include "internal/io/sqlite/ModelObjectDAO.hpp"
#include "internal/io/sqlite/DataObjectDAOFactoryRegistry.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace
{
    rhbm_gem::DataObjectDAORegistrar<rhbm_gem::ModelObject, rhbm_gem::ModelObjectDAO>
        registrar_model_dao("model");
}
