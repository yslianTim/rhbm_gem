#pragma once

#include "ModelObjectDaoSqlite.hpp"

namespace rhbm_gem {

class ModelObjectDAO : public ModelObjectDaoSqlite
{
public:
    using ModelObjectDaoSqlite::ModelObjectDaoSqlite;
};

} // namespace rhbm_gem
