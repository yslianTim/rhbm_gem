#pragma once

#include <string>
#include <utility>

#include "internal/io/sqlite/SQLiteWrapper.hpp"

namespace rhbm_gem::model_io {

class SQLiteStatementBatch
{
    SQLiteWrapper & m_database;
    SQLiteWrapper::StatementGuard m_guard;

public:
    SQLiteStatementBatch(SQLiteWrapper & database, const std::string & sql) :
        m_database{ database },
        m_guard{ database }
    {
        m_database.Prepare(sql);
    }

    template <typename Binder>
    void Execute(Binder && binder)
    {
        std::forward<Binder>(binder)(m_database);
        m_database.StepOnce();
        m_database.Reset();
    }
};

} // namespace rhbm_gem::model_io
