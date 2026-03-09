#pragma once

#include <memory>
#include <string>

#include "DataObjectDAOBase.hpp"

namespace rhbm_gem {

class SQLiteWrapper;

class ModelObjectDaoSqlite : public DataObjectDAOBase
{
    SQLiteWrapper * m_database;

public:
    explicit ModelObjectDaoSqlite(SQLiteWrapper * db_manager);
    ~ModelObjectDaoSqlite();

    static void EnsureSchema(SQLiteWrapper & database);

    void Save(const DataObjectBase & obj, const std::string & key_tag) override;
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) override;
};

} // namespace rhbm_gem
