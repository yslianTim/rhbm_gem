#include "internal/ModelObjectDAO.hpp"
#include "internal/DataObjectDAOFactoryRegistry.hpp"
#include <rhbm_gem/data/ModelObject.hpp>

namespace
{
    rhbm_gem::DataObjectDAORegistrar<rhbm_gem::ModelObject, rhbm_gem::ModelObjectDAO>
        registrar_model_dao("model");
}
