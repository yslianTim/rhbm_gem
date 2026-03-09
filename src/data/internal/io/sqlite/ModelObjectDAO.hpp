#pragma once

#include "ModelObjectDaoSqlite.hpp"

namespace rhbm_gem {

class ModelObjectDAO : public ModelObjectDAOSqlite
{
public:
    using ModelObjectDAOSqlite::ModelObjectDAOSqlite;
};

} // namespace rhbm_gem
